#!/usr/bin/env python3
"""
Script: 03_curate_bioactivity.py
Description:
    Extracts and curates ALL compound bioactivity data from ChEMBL 37.
    Applies BROAD curation filters:
      - All organism targets (not restricted to human/rodent)
      - All standard activity types (IC50, Ki, EC50, Kd, etc.)
      - Includes mutant targets
      - Requires valid pChEMBL value OR numeric standard_value
      - Removes clear data quality flags (data_validity_comment exclusions)
    
    Enriches with:
      - Compound identifiers (ChEMBL ID, InChIKey, SMILES, canonical SMILES)
      - Compound properties (MW, LogP, HBD, HBA, PSA, RO5, QED, etc.)
      - Assay metadata (type, organism, confidence score, assay description)
      - Target metadata (ChEMBL target ID, name, type, taxonomy)
      - Document metadata (journal, year, DOI/PMID)
      - UniProt accession (via chembl_uniprot_mapping.txt)
      - Ligand efficiency metrics (LE, BEI, SEI, LLE)
    
    Outputs:
      data/curated/chembl_37_bioactivity_broad.tsv.gz   (main curated table)
      data/curated/chembl_37_bioactivity_stats.tsv       (summary statistics)

Usage: python scripts/03_curate_bioactivity.py
       python scripts/03_curate_bioactivity.py --chunksize 500000
"""

import sqlite3
import os
import sys
import gzip
import csv
import logging
import argparse
from datetime import datetime
from pathlib import Path
from collections import defaultdict

# --- Configuration -----------------------------------------------------------
CHEMBL_DB      = "data/raw/chembl_37.db"
UNIPROT_MAP    = "data/raw/chembl_uniprot_mapping.txt"
CURATED_DIR    = Path("data/curated")
LOG_DIR        = Path("logs")
LOGFILE        = LOG_DIR / f"curate_bioactivity_{datetime.now():%Y%m%d_%H%M%S}.log"

OUTPUT_FILE    = CURATED_DIR / "chembl_37_bioactivity_broad.tsv.gz"
STATS_FILE     = CURATED_DIR / "chembl_37_bioactivity_stats.tsv"

# --- Data quality exclusions (broad: remove only clear errors) ---------------
# These data_validity_comment values indicate erroneous/non-standard data
EXCLUDE_VALIDITY = {
    "Potential transcription error",
    "Outside typical range",
    "Non standard unit for type",
    "Author confirmed error",
    "Potential author error",
}

# Activity types to always flag (not exclude, but mark) for quality awareness
FLAGGED_TYPES = {"Ratio", "Fold", "log ratio"}

# --- Column order in output --------------------------------------------------
OUTPUT_COLUMNS = [
    # Activity identifiers & source tracking
    "activity_id",
    "assay_id",
    "doc_id",
    "record_id",
    "src_id",
    "src_short_name",

    # Compound identifiers
    "molregno",
    "chembl_id_compound",
    "compound_name",
    "compound_inchikey",
    "canonical_smiles",
    "max_phase",
    "molecule_type",
    "is_natural_product",
    "is_chemical_probe",

    # Compound properties
    "mw_freebase",
    "alogp",
    "hbd",
    "hba",
    "psa",
    "rtb",
    "ro3_pass",
    "num_ro5_violations",
    "full_mwt",
    "aromatic_rings",
    "heavy_atoms",
    "qed_weighted",
    "full_molformula",
    "np_likeness_score",

    # Original published measurements (Traceability - M03, M04)
    "original_type",
    "original_relation",
    "original_value",
    "original_units",
    "original_text_value",
    "original_upper_value",

    # Standardized bioactivity values (M01, M05)
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "standard_upper_value",
    "standard_flag",
    "standard_text_value",
    "pchembl_value",
    "activity_comment",
    "data_validity_comment",
    "potential_duplicate",
    "modality",
    "uo_units",
    "qudt_units",

    # BioAssay Ontology (BAO) annotations (M10)
    "bao_endpoint",
    "bao_format",

    # Ligand efficiency
    "bei",
    "le",
    "lle",
    "sei",

    # Assay metadata
    "assay_type",
    "assay_type_desc",
    "assay_test_type",
    "assay_category",
    "assay_description",
    "assay_organism",
    "assay_tax_id",
    "assay_strain",
    "assay_tissue",
    "assay_cell_type",
    "assay_subcellular_fraction",
    "confidence_score",
    "confidence_desc",
    "relationship_type",
    "relationship_desc",

    # Target metadata & variant tracking (C06)
    "tid",
    "variant_id",
    "chembl_id_target",
    "target_pref_name",
    "target_type",
    "target_organism",
    "target_tax_id",
    "species_group_flag",

    # Target annotations (UniProt)
    "uniprot_accession",
    "uniprot_name",
    "taxonomy_id",

    # Document metadata
    "doc_chembl_id",
    "journal",
    "year",
    "pubmed_id",
    "doi",
    "patent_id",
]


# --- Logging -----------------------------------------------------------------
LOG_DIR.mkdir(parents=True, exist_ok=True)
CURATED_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOGFILE),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)


def load_uniprot_mapping(filepath: str) -> dict:
    """
    Load chembl_uniprot_mapping.txt
    Format: UniProt_accession  UniProt_Name  ChEMBL_target_ID  taxonomy_id
    Returns dict: {chembl_target_id -> [(uniprot_acc, uniprot_name, tax_id), ...]}
    """
    mapping = defaultdict(list)
    if not os.path.exists(filepath):
        log.warning(f"UniProt mapping not found: {filepath}")
        return mapping
    with open(filepath) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                uniprot_acc, uniprot_name, chembl_tid, tax_id = parts[:4]
                mapping[chembl_tid].append((uniprot_acc, uniprot_name, tax_id))
    log.info(f"Loaded UniProt mapping: {len(mapping):,} ChEMBL targets with UniProt entries")
    return mapping


def main(chunksize: int = 200_000):
    log.info("=" * 70)
    log.info(" ChEMBL 37 Compound Bioactivity Curation (BROAD)")
    log.info("=" * 70)

    # Load UniProt mapping
    uniprot_map = load_uniprot_mapping(UNIPROT_MAP)

    # Connect
    if not os.path.exists(CHEMBL_DB):
        log.error(f"Database not found: {CHEMBL_DB}")
        log.error("Run scripts/01_download_chembl.sh first.")
        sys.exit(1)

    conn = sqlite3.connect(CHEMBL_DB)
    conn.row_factory = sqlite3.Row
    log.info(f"Connected to: {CHEMBL_DB}")

    # Performance pragma settings
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-1000000")  # ~1 GB cache
    conn.execute("PRAGMA mmap_size=30000000000")

    # -------------------------------------------------------------------------
    # MAIN QUERY: Broad bioactivity extraction
    # Uses CROSS JOIN to force sequential scan of activities table with
    # primary key lookups into dimension tables.
    # -------------------------------------------------------------------------
    MAIN_SQL = """
    SELECT
        -- Activity core & source tracking
        act.activity_id,
        act.assay_id,
        act.doc_id,
        act.record_id,
        act.src_id,
        src.src_short_name,

        -- Compound identifiers
        act.molregno,
        md.chembl_id                        AS chembl_id_compound,
        md.pref_name                        AS compound_name,
        cs.standard_inchi_key               AS compound_inchikey,
        cs.canonical_smiles,
        md.max_phase,
        md.molecule_type,
        md.natural_product                  AS is_natural_product,
        md.chemical_probe                   AS is_chemical_probe,

        -- Compound properties
        cp.mw_freebase,
        cp.alogp,
        cp.hbd,
        cp.hba,
        cp.psa,
        cp.rtb,
        cp.ro3_pass,
        cp.num_ro5_violations,
        cp.full_mwt,
        cp.aromatic_rings,
        cp.heavy_atoms,
        cp.qed_weighted,
        cp.full_molformula,
        cp.np_likeness_score,

        -- Original published measurements (Traceability - M03, M04)
        act.type                            AS original_type,
        act.relation                        AS original_relation,
        act.value                           AS original_value,
        act.units                           AS original_units,
        act.text_value                      AS original_text_value,
        act.upper_value                     AS original_upper_value,

        -- Standardized bioactivity values (M01, M05)
        act.standard_type,
        act.standard_relation,
        act.standard_value,
        act.standard_units,
        act.standard_upper_value,
        act.standard_flag,
        act.standard_text_value,
        act.pchembl_value,
        act.activity_comment,
        act.data_validity_comment,
        act.potential_duplicate,
        act.modality,
        act.uo_units,
        act.qudt_units,

        -- BioAssay Ontology (BAO) annotations (M10)
        act.bao_endpoint,
        asy.bao_format,

        -- Ligand efficiency
        le.bei,
        le.le,
        le.lle,
        le.sei,

        -- Assay metadata
        asy.assay_type,
        atl.assay_desc                      AS assay_type_desc,
        asy.assay_test_type,
        asy.assay_category,
        asy.description                     AS assay_description,
        asy.assay_organism,
        asy.assay_tax_id,
        asy.assay_strain,
        asy.assay_tissue,
        asy.assay_cell_type,
        asy.assay_subcellular_fraction,
        asy.confidence_score,
        csl.description                     AS confidence_desc,
        asy.relationship_type,
        rtl.relationship_desc,

        -- Target metadata & variant tracking (C06)
        asy.tid,
        asy.variant_id,
        td.chembl_id                        AS chembl_id_target,
        td.pref_name                        AS target_pref_name,
        td.target_type,
        td.organism                         AS target_organism,
        td.tax_id                           AS target_tax_id,
        td.species_group_flag,

        -- Document metadata
        doc.chembl_id                       AS doc_chembl_id,
        doc.journal,
        doc.year,
        doc.pubmed_id,
        doc.doi,
        doc.patent_id

    FROM activities act

    CROSS JOIN assays asy
        ON act.assay_id = asy.assay_id

    LEFT JOIN target_dictionary td
        ON asy.tid = td.tid

    LEFT JOIN molecule_dictionary md
        ON act.molregno = md.molregno

    LEFT JOIN compound_structures cs
        ON act.molregno = cs.molregno

    LEFT JOIN compound_properties cp
        ON act.molregno = cp.molregno

    LEFT JOIN ligand_eff le
        ON act.activity_id = le.activity_id

    LEFT JOIN docs doc
        ON act.doc_id = doc.doc_id

    LEFT JOIN source src
        ON act.src_id = src.src_id

    LEFT JOIN assay_type atl
        ON asy.assay_type = atl.assay_type

    LEFT JOIN confidence_score_lookup csl
        ON asy.confidence_score = csl.confidence_score

    LEFT JOIN relationship_type rtl
        ON asy.relationship_type = rtl.relationship_type

    WHERE act.molregno IS NOT NULL
      AND (
          act.standard_value IS NOT NULL
          OR act.pchembl_value IS NOT NULL
          OR act.activity_comment IS NOT NULL
          OR act.text_value IS NOT NULL
          OR act.standard_text_value IS NOT NULL
      )
    """

    # -------------------------------------------------------------------------
    # Count total for progress
    # -------------------------------------------------------------------------
    log.info("Counting total activities to process...")
    count_sql = """
        SELECT COUNT(*) FROM activities act
        WHERE (
            act.standard_value IS NOT NULL
            OR act.pchembl_value IS NOT NULL
            OR act.activity_comment IS NOT NULL
        )
    """
    total = conn.execute(count_sql).fetchone()[0]
    log.info(f"Total activities to process: {total:,}")

    # -------------------------------------------------------------------------
    # Stream and write in chunks
    # -------------------------------------------------------------------------
    log.info(f"Writing curated output: {OUTPUT_FILE}")
    log.info(f"Chunk size: {chunksize:,} rows")

    stats = defaultdict(int)
    stats["total_input"] = total

    with gzip.open(OUTPUT_FILE, "wt", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(
            fout,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()

        cursor = conn.execute(MAIN_SQL)
        columns = [desc[0] for desc in cursor.description]

        written = 0
        skipped_validity = 0
        chunk_num = 0

        while True:
            rows = cursor.fetchmany(chunksize)
            if not rows:
                break

            chunk_num += 1
            chunk_written = 0

            for row in rows:
                d = dict(zip(columns, row))

                # --- Curation filter: exclude clear data quality errors -------
                validity = d.get("data_validity_comment") or ""
                if validity in EXCLUDE_VALIDITY:
                    skipped_validity += 1
                    continue

                # --- Add UniProt mapping --------------------------------------
                chembl_tid = d.get("chembl_id_target") or ""
                uniprot_entries = uniprot_map.get(chembl_tid, [])
                if uniprot_entries:
                    # Use first entry; concat multiple accessions with ";"
                    d["uniprot_accession"] = ";".join(e[0] for e in uniprot_entries)
                    d["uniprot_name"]      = ";".join(e[1] for e in uniprot_entries)
                    d["taxonomy_id"]       = ";".join(e[2] for e in uniprot_entries)
                else:
                    d["uniprot_accession"] = ""
                    d["uniprot_name"]      = ""
                    d["taxonomy_id"]       = ""

                # --- Statistics trackers -------------------------------------
                if d.get("pchembl_value"):
                    stats["has_pchembl"] += 1
                if d.get("standard_type") in FLAGGED_TYPES:
                    stats["flagged_ratio_type"] += 1

                std_type = (d.get("standard_type") or "UNKNOWN").upper()
                stats[f"type_{std_type}"] += 1

                writer.writerow(d)
                chunk_written += 1

            written += chunk_written
            pct = 100 * written / total if total > 0 else 0
            log.info(
                f"  Chunk {chunk_num:4d}: +{chunk_written:>9,} rows | "
                f"Total: {written:>11,} / {total:>11,} ({pct:.1f}%)"
            )

    stats["total_written"] = written
    stats["skipped_quality_flag"] = skipped_validity

    conn.close()

    # -------------------------------------------------------------------------
    # Write statistics
    # -------------------------------------------------------------------------
    log.info(f"\nWriting statistics: {STATS_FILE}")
    with open(STATS_FILE, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in sorted(stats.items()):
            w.writerow([k, v])

    log.info("")
    log.info("=" * 70)
    log.info(" Bioactivity curation complete!")
    log.info(f"  Input activities  : {stats['total_input']:>12,}")
    log.info(f"  Skipped (quality) : {stats['skipped_quality_flag']:>12,}")
    log.info(f"  Written records   : {stats['total_written']:>12,}")
    log.info(f"  With pChEMBL value: {stats['has_pchembl']:>12,}")
    log.info(f"  Output file       : {OUTPUT_FILE}")
    log.info("=" * 70)
    log.info("Next step: python scripts/04_curate_target_interactions.py")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Curate ChEMBL 37 bioactivity data")
    parser.add_argument(
        "--chunksize", type=int, default=200_000,
        help="Number of rows to fetch per batch (default: 200000)"
    )
    args = parser.parse_args()
    main(chunksize=args.chunksize)
