#!/usr/bin/env python3
"""
Script: 04_curate_target_interactions.py
Description:
    Evidence-Tiered Compound-Target Interaction Engine for ChEMBL 37.
    
    Addresses Scientific & Pharmacological Audit Findings (C01-C12, M10-M22):
      1. Evidence Tiers: Separates direct single-protein binding from functional,
         complex, phenotypic, and ADMET assays.
      2. Endpoint Stratification: Never compares Ki, Kd, IC50, or EC50 against each other.
      3. Statistical Aggregation: Calculates median pChEMBL, min, max, exact count,
         censoring count (> and <), and independent doc/assay counts (NO selection bias).
      4. Target Variants: Preserves mutation/variant annotations (C06).
      5. Curated Drug Registry: Integrates mechanism of action, ATC, indications,
         and drug warnings directly onto interactions.
      6. Normalized Biology: Generates separate relational tables for target components,
         variants, and drug dictionary.

    Outputs:
      data/curated/chembl_37_compound_target_interactions_tiered.tsv.gz
      data/curated/chembl_37_target_components.tsv
      data/curated/chembl_37_target_variants.tsv
      data/curated/chembl_37_drug_dictionary.tsv
      data/curated/chembl_37_protein_classification.tsv
"""

import sqlite3
import os
import sys
import gzip
import csv
import logging
from datetime import datetime
from pathlib import Path
from collections import defaultdict
import statistics

# --- Configuration -----------------------------------------------------------
CHEMBL_DB            = "data/raw/chembl_37.db"
CURATED_DIR          = Path("data/curated")
LOG_DIR              = Path("logs")
LOGFILE              = LOG_DIR / f"curate_targets_{datetime.now():%Y%m%d_%H%M%S}.log"

TIERED_INTERACTIONS_FILE = CURATED_DIR / "chembl_37_compound_target_interactions_tiered.tsv.gz"
TARGET_COMPONENTS_FILE   = CURATED_DIR / "chembl_37_target_components.tsv"
TARGET_VARIANTS_FILE     = CURATED_DIR / "chembl_37_target_variants.tsv"
DRUG_DICT_FILE           = CURATED_DIR / "chembl_37_drug_dictionary.tsv"
PROTEIN_CLASS_FILE       = CURATED_DIR / "chembl_37_protein_classification.tsv"

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


def connect(db_path: str) -> sqlite3.Connection:
    if not os.path.exists(db_path):
        log.error(f"Database not found: {db_path}")
        sys.exit(1)
    # Open in URI mode for read-only safety (T08)
    conn = sqlite3.connect(f"file:{os.path.abspath(db_path)}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-1000000")
    conn.execute("PRAGMA mmap_size=30000000000")
    return conn


def export_target_components(conn: sqlite3.Connection):
    """
    Export normalized target-component mapping (M13, M15, M21).
    Preserves 1-to-1 association between each target, component, UniProt ID,
    homologue flag, and protein classification lineage.
    """
    log.info("Exporting normalized target components (M13, M15, M21)...")
    SQL = """
    SELECT
        td.tid,
        td.chembl_id                        AS chembl_id_target,
        td.pref_name                        AS target_pref_name,
        td.target_type,
        td.organism                         AS target_organism,
        td.tax_id                           AS target_tax_id,
        tc.component_id,
        tc.homologue,
        cs.accession                        AS uniprot_accession,
        cs.component_type,
        cs.description                      AS component_description,
        cs.db_source,
        pc.protein_class_id,
        pc.pref_name                        AS protein_class_name,
        pc.class_level                      AS protein_class_level
    FROM target_dictionary td
    INNER JOIN target_components tc ON td.tid = tc.tid
    LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
    LEFT JOIN component_class cc ON cs.component_id = cc.component_id
    LEFT JOIN protein_classification pc ON cc.protein_class_id = pc.protein_class_id
    ORDER BY td.tid, tc.component_id
    """
    rows = conn.execute(SQL).fetchall()
    log.info(f"  {len(rows):,} target-component records")
    tmp_file = TARGET_COMPONENTS_FILE.with_suffix(".tmp")
    with open(tmp_file, "w", newline="", encoding="utf-8") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows([dict(r) for r in rows])
    tmp_file.replace(TARGET_COMPONENTS_FILE)
    log.info(f"  Saved: {TARGET_COMPONENTS_FILE}")


def export_target_variants(conn: sqlite3.Connection):
    """
    Export target variants and resistance mutations (C06).
    """
    log.info("Exporting target variants and mutations (C06)...")
    SQL = """
    SELECT
        asy.tid,
        td.chembl_id                        AS chembl_id_target,
        td.pref_name                        AS target_pref_name,
        vs.variant_id,
        vs.mutation,
        vs.accession                        AS variant_accession,
        vs.isoform,
        vs.organism                         AS variant_organism,
        COUNT(DISTINCT asy.assay_id)        AS n_assays
    FROM assays asy
    INNER JOIN variant_sequences vs ON asy.variant_id = vs.variant_id
    INNER JOIN target_dictionary td ON asy.tid = td.tid
    GROUP BY asy.tid, vs.variant_id
    ORDER BY asy.tid, vs.variant_id
    """
    rows = conn.execute(SQL).fetchall()
    log.info(f"  {len(rows):,} target variant mappings")
    tmp_file = TARGET_VARIANTS_FILE.with_suffix(".tmp")
    with open(tmp_file, "w", newline="", encoding="utf-8") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows([dict(r) for r in rows])
    tmp_file.replace(TARGET_VARIANTS_FILE)
    log.info(f"  Saved: {TARGET_VARIANTS_FILE}")


def export_drug_dictionary(conn: sqlite3.Connection):
    """
    Export curated drug registry (M19, Section 4).
    Enriches with ATC codes, clinical indications, warnings, and mechanism of action.
    """
    log.info("Exporting curated drug registry (M19, Section 4)...")
    SQL = """
    SELECT
        md.molregno,
        md.chembl_id                        AS chembl_id_compound,
        md.pref_name                        AS drug_name,
        md.max_phase,
        md.first_approval,
        md.oral,
        md.parenteral,
        md.topical,
        md.black_box_warning,
        md.withdrawn_flag,
        cs.canonical_smiles,
        cs.standard_inchi_key               AS compound_inchikey,
        GROUP_CONCAT(DISTINCT atc.level5)   AS atc_codes,
        GROUP_CONCAT(DISTINCT di.mesh_heading) AS mesh_indications,
        GROUP_CONCAT(DISTINCT di.efo_term)  AS efo_indications,
        GROUP_CONCAT(DISTINCT dw.warning_type || ': ' || dw.warning_class) AS drug_warnings,
        GROUP_CONCAT(DISTINCT dm.mechanism_of_action || ' (' || dm.action_type || ')') AS drug_mechanisms
    FROM molecule_dictionary md
    LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
    LEFT JOIN molecule_atc_classification atc ON md.molregno = atc.molregno
    LEFT JOIN drug_indication di ON md.molregno = di.molregno
    LEFT JOIN drug_warning dw ON md.molregno = dw.molregno
    LEFT JOIN drug_mechanism dm ON md.molregno = dm.molregno
    WHERE md.max_phase >= 1 OR md.first_approval IS NOT NULL OR di.molregno IS NOT NULL
    GROUP BY md.molregno
    ORDER BY md.max_phase DESC, md.chembl_id
    """
    rows = conn.execute(SQL).fetchall()
    log.info(f"  {len(rows):,} curated clinical drug records")
    tmp_file = DRUG_DICT_FILE.with_suffix(".tmp")
    with open(tmp_file, "w", newline="", encoding="utf-8") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows([dict(r) for r in rows])
    tmp_file.replace(DRUG_DICT_FILE)
    log.info(f"  Saved: {DRUG_DICT_FILE}")


def export_protein_classification(conn: sqlite3.Connection):
    """Export hierarchical protein classification."""
    log.info("Exporting protein classification hierarchy...")
    SQL = """
    SELECT
        pc.protein_class_id,
        pc.parent_id,
        pc.pref_name,
        pc.short_name,
        pc.class_level,
        parent.pref_name AS parent_name
    FROM protein_classification pc
    LEFT JOIN protein_classification parent ON pc.parent_id = parent.protein_class_id
    ORDER BY pc.class_level, pc.protein_class_id
    """
    rows = conn.execute(SQL).fetchall()
    log.info(f"  {len(rows):,} protein class records")
    tmp_file = PROTEIN_CLASS_FILE.with_suffix(".tmp")
    with open(tmp_file, "w", newline="", encoding="utf-8") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows([dict(r) for r in rows])
    tmp_file.replace(PROTEIN_CLASS_FILE)
    log.info(f"  Saved: {PROTEIN_CLASS_FILE}")


def build_evidence_tiered_interactions(conn: sqlite3.Connection, chunksize: int = 250_000):
    """
    Builds the Evidence-Tiered Compound-Target Interaction Engine.
    
    Stratification: (molregno, tid, variant_id, evidence_tier, standard_type, standard_units)
    Evidence Tiers:
      - Tier A: Direct Single-Protein Biochemical Binding (Confidence >= 8, Assay B, Target 'SINGLE PROTEIN')
      - Tier B: Direct Single-Protein Functional Potency (Confidence >= 8, Assay F, Target 'SINGLE PROTEIN')
      - Tier C: Protein Complex / Family / Homologue Molecular Target
      - Tier D: Cellular / Tissue / Organism Phenotype (GI50, MIC, etc.)
      - Tier E: ADMET / Toxicity / Physicochemical
    
    Statistics:
      - n_activities, n_exact, n_censored_gt, n_censored_lt, n_inactive
      - median_pchembl, min_pchembl, max_pchembl
      - median_standard_value, standard_units
      - n_assays, n_docs
      - Integrated drug mechanism of action annotations
    """
    log.info("Building Evidence-Tiered Compound-Target Interaction Engine...")

    # Load target cache
    cur = conn.cursor()
    target_info = {}
    for r in cur.execute("""
        SELECT
            td.tid,
            td.chembl_id AS chembl_id_target,
            td.pref_name AS target_pref_name,
            td.target_type,
            td.organism AS target_organism,
            td.tax_id AS target_tax_id,
            GROUP_CONCAT(DISTINCT cs.accession) AS uniprot_accessions
        FROM target_dictionary td
        LEFT JOIN target_components tc ON td.tid = tc.tid
        LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
        GROUP BY td.tid
    """):
        target_info[r["tid"]] = dict(r)
    log.info(f"  Cached {len(target_info):,} target metadata records")

    # Load variant cache
    variant_info = {}
    for r in cur.execute("SELECT variant_id, mutation, accession FROM variant_sequences"):
        variant_info[r["variant_id"]] = dict(r)
    log.info(f"  Cached {len(variant_info):,} variant mutation records")

    # Load drug mechanism annotations indexed by (molregno, tid)
    drug_mechs = defaultdict(list)
    for r in cur.execute("""
        SELECT molregno, tid, mechanism_of_action, action_type, direct_interaction
        FROM drug_mechanism
        WHERE molregno IS NOT NULL AND tid IS NOT NULL
    """):
        drug_mechs[(r["molregno"], r["tid"])].append(
            f"{r['mechanism_of_action']} [{r['action_type'] or 'UNKNOWN'}]"
        )
    log.info(f"  Cached {len(drug_mechs):,} drug-target mechanism mappings")

    # Output schema
    OUTPUT_COLUMNS = [
        "molregno",
        "chembl_id_compound",
        "compound_name",
        "compound_inchikey",
        "canonical_smiles",
        "max_phase",
        "tid",
        "chembl_id_target",
        "target_pref_name",
        "target_type",
        "target_organism",
        "target_tax_id",
        "uniprot_accessions",
        "variant_id",
        "mutation",
        "evidence_tier",
        "evidence_tier_desc",
        "standard_type",
        "standard_units",
        "n_activities",
        "n_exact",
        "n_censored_gt",
        "n_censored_lt",
        "n_inactive_comments",
        "median_pchembl",
        "min_pchembl",
        "max_pchembl",
        "median_standard_value",
        "n_assays",
        "n_docs",
        "curated_drug_mechanisms",
    ]

    # Main extraction query using fk_act_molregno to stream rows pre-grouped by molregno (constant memory)
    SQL = """
    SELECT
        act.molregno,
        asy.tid,
        COALESCE(asy.variant_id, 0) AS variant_id,
        asy.assay_type,
        asy.confidence_score,
        td.target_type,
        act.standard_type,
        act.standard_relation,
        act.standard_value,
        act.standard_units,
        act.pchembl_value,
        act.activity_comment,
        act.assay_id,
        act.doc_id,
        md.chembl_id                        AS chembl_id_compound,
        md.pref_name                        AS compound_name,
        cs.standard_inchi_key               AS compound_inchikey,
        cs.canonical_smiles,
        md.max_phase
    FROM activities act INDEXED BY fk_act_molregno
    CROSS JOIN assays asy ON act.assay_id = asy.assay_id
    LEFT JOIN target_dictionary td ON asy.tid = td.tid
    LEFT JOIN molecule_dictionary md ON act.molregno = md.molregno
    LEFT JOIN compound_structures cs ON act.molregno = cs.molregno
    WHERE act.molregno IS NOT NULL
      AND asy.tid IS NOT NULL
      AND act.standard_type IS NOT NULL
      AND (
          act.standard_value IS NOT NULL
          OR act.pchembl_value IS NOT NULL
          OR act.activity_comment IS NOT NULL
      )
      AND (
          act.data_validity_comment IS NULL
          OR act.data_validity_comment NOT IN (
              'Potential transcription error',
              'Outside typical range',
              'Non standard unit for type',
              'Author confirmed error',
              'Potential author error'
          )
      )
    """

    tmp_file = CURATED_DIR / (TIERED_INTERACTIONS_FILE.name + ".tmp")
    log.info(f"Streaming and aggregating interactions to: {tmp_file}")

    cursor = conn.execute(SQL)
    written_count = 0
    row_batch = 0
    current_molregno = None

    with gzip.open(tmp_file, "wt", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        # Stream accumulator: group activities per compound
        accumulator = defaultdict(lambda: {
            "pchembls": [],
            "exact_vals": [],
            "n_exact": 0,
            "n_gt": 0,
            "n_lt": 0,
            "n_inactive": 0,
            "assays": set(),
            "docs": set(),
            "compound_meta": {},
        })

        def flush_accumulator():
            nonlocal written_count
            tier_desc_map = {
                "Tier A": "Direct Single-Protein Biochemical Binding",
                "Tier B": "Direct Single-Protein Functional Potency",
                "Tier C": "Protein Complex / Family / Homologue Molecular Target",
                "Tier D": "Cellular / Tissue / Organism Phenotype",
                "Tier E": "ADMET / Toxicity / Physicochemical / Unassigned",
            }
            for key, data in accumulator.items():
                molregno, tid, vid, tier, st, units = key
                cm = data["compound_meta"]
                tmeta = target_info.get(tid, {})
                vmeta = variant_info.get(vid, {})

                pchembls = data["pchembls"]
                exact_vals = data["exact_vals"]

                med_p = round(statistics.median(pchembls), 3) if pchembls else ""
                min_p = round(min(pchembls), 3) if pchembls else ""
                max_p = round(max(pchembls), 3) if pchembls else ""
                med_val = round(statistics.median(exact_vals), 4) if exact_vals else ""

                # Drug mechanism annotation
                mechs = "; ".join(drug_mechs.get((molregno, tid), []))

                out_row = {
                    "molregno": molregno,
                    "chembl_id_compound": cm.get("chembl_id_compound", ""),
                    "compound_name": cm.get("compound_name", ""),
                    "compound_inchikey": cm.get("compound_inchikey", ""),
                    "canonical_smiles": cm.get("canonical_smiles", ""),
                    "max_phase": cm.get("max_phase", ""),
                    "tid": tid,
                    "chembl_id_target": tmeta.get("chembl_id_target", ""),
                    "target_pref_name": tmeta.get("target_pref_name", ""),
                    "target_type": tmeta.get("target_type", ""),
                    "target_organism": tmeta.get("target_organism", ""),
                    "target_tax_id": tmeta.get("target_tax_id", ""),
                    "uniprot_accessions": tmeta.get("uniprot_accessions", ""),
                    "variant_id": vid if vid > 0 else "",
                    "mutation": vmeta.get("mutation", ""),
                    "evidence_tier": tier,
                    "evidence_tier_desc": tier_desc_map.get(tier, ""),
                    "standard_type": st,
                    "standard_units": units,
                    "n_activities": data["n_exact"] + data["n_gt"] + data["n_lt"] + data["n_inactive"],
                    "n_exact": data["n_exact"],
                    "n_censored_gt": data["n_gt"],
                    "n_censored_lt": data["n_lt"],
                    "n_inactive_comments": data["n_inactive"],
                    "median_pchembl": med_p,
                    "min_pchembl": min_p,
                    "max_pchembl": max_p,
                    "median_standard_value": med_val,
                    "n_assays": len(data["assays"]),
                    "n_docs": len(data["docs"]),
                    "curated_drug_mechanisms": mechs,
                }
                writer.writerow(out_row)
                written_count += 1
            accumulator.clear()

        while True:
            rows = cursor.fetchmany(chunksize)
            if not rows:
                break

            for r in rows:
                row_batch += 1
                m = r["molregno"]
                if m != current_molregno:
                    if current_molregno is not None:
                        flush_accumulator()
                    current_molregno = m

                # Evidence tier categorization (C02, C04, C05)
                atype = r["assay_type"]
                cscore = r["confidence_score"] or 0
                ttype = r["target_type"] or ""

                if atype == "B" and cscore >= 8 and ttype == "SINGLE PROTEIN":
                    tier = "Tier A"
                elif atype == "F" and cscore >= 8 and ttype == "SINGLE PROTEIN":
                    tier = "Tier B"
                elif ttype in ("PROTEIN COMPLEX", "PROTEIN FAMILY", "PROTEIN-PROTEIN INTERACTION", "CHIMERIC PROTEIN") or cscore in (4, 5, 6, 7):
                    tier = "Tier C"
                elif ttype in ("CELL-LINE", "ORGANISM", "TISSUE") or cscore <= 3:
                    tier = "Tier D"
                else:
                    tier = "Tier E"

                key = (
                    m,
                    r["tid"],
                    r["variant_id"],
                    tier,
                    r["standard_type"],
                    r["standard_units"] or "",
                )

                acc = accumulator[key]
                acc["assays"].add(r["assay_id"])
                if r["doc_id"]:
                    acc["docs"].add(r["doc_id"])

                if not acc["compound_meta"]:
                    acc["compound_meta"] = {
                        "chembl_id_compound": r["chembl_id_compound"],
                        "compound_name": r["compound_name"],
                        "compound_inchikey": r["compound_inchikey"],
                        "canonical_smiles": r["canonical_smiles"],
                        "max_phase": r["max_phase"],
                    }

                rel = r["standard_relation"]
                if rel in ("=", "~", "=="):
                    acc["n_exact"] += 1
                    if r["pchembl_value"] is not None:
                        acc["pchembls"].append(float(r["pchembl_value"]))
                    if r["standard_value"] is not None:
                        acc["exact_vals"].append(float(r["standard_value"]))
                elif rel in (">", ">="):
                    acc["n_gt"] += 1
                elif rel in ("<", "<="):
                    acc["n_lt"] += 1
                else:
                    acc["n_inactive"] += 1

            if row_batch % 1_000_000 == 0:
                log.info(f"  Ingested {row_batch:>10,} activities | Written tiered pairs: {written_count:>10,}")

        # Flush final compound's groups
        if accumulator:
            flush_accumulator()

    tmp_file.replace(TIERED_INTERACTIONS_FILE)
    log.info(f"  Total tiered interaction rows written: {written_count:,}")
    log.info(f"  Saved: {TIERED_INTERACTIONS_FILE}")
    return written_count


def main():
    log.info("=" * 70)
    log.info(" ChEMBL 37 Evidence-Tiered Interaction Curation")
    log.info("=" * 70)

    conn = connect(CHEMBL_DB)

    # 1. Normalized target-component registry (M13, M15, M21)
    export_target_components(conn)

    # 2. Target variants & mutations registry (C06)
    export_target_variants(conn)

    # 3. Curated drug & clinical registry (M19, Section 4)
    export_drug_dictionary(conn)

    # 4. Protein classification hierarchy (M14)
    export_protein_classification(conn)

    # 5. Evidence-Tiered Compound-Target Interaction Engine (C01-C12)
    n_interactions = build_evidence_tiered_interactions(conn)

    conn.close()

    log.info("")
    log.info("=" * 70)
    log.info(" Scientific target curation complete!")
    log.info(f"  Tiered Interactions   : {n_interactions:>10,}")
    log.info(f"  Interactions file     : {TIERED_INTERACTIONS_FILE}")
    log.info(f"  Target components     : {TARGET_COMPONENTS_FILE}")
    log.info(f"  Target variants       : {TARGET_VARIANTS_FILE}")
    log.info(f"  Drug dictionary       : {DRUG_DICT_FILE}")
    log.info(f"  Protein classes       : {PROTEIN_CLASS_FILE}")
    log.info("=" * 70)


if __name__ == "__main__":
    main()
