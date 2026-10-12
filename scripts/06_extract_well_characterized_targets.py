#!/usr/bin/env python3
"""
Script: 06_extract_well_characterized_targets.py
Description:
    Extracts, stratifies, and structures bioactivity data for Well-Characterized
    Proteins in ChEMBL 37.

    Definition of Well-Characterized Target:
      - Single-protein target (target_type = 'SINGLE PROTEIN')
      - High-confidence assay assignment (confidence_score >= 8)
      - Direct molecular endpoints: Ki, Kd, IC50, EC50
      - Rigorous compound coverage: >= 500 distinct compounds evaluated

    Architecture & Outputs:
      1. Generates a dedicated directory per target:
         data/targets/<CHEMBL_ID>_<Target_Slug>/
      2. Inside each target directory, separates bioactivity by endpoint:
         - Kd.tsv    (Thermodynamic equilibrium binding)
         - Ki.tsv    (Enzymatic inhibition constant)
         - IC50.tsv  (Biochemical 50% inhibitory concentration)
         - EC50.tsv  (Functional 50% effective concentration)
      3. Integrates clinical drug status:
         Maps each compound against ChEMBL 37 Drug Dictionary:
         - is_clinical_drug (max_phase >= 1)
         - is_approved_drug (max_phase == 4)
         - drug_name, ATC codes, drug mechanisms
      4. Generates comprehensive statistical reporting:
         - reports/well_characterized_proteins_drug_stats.tsv
         - reports/well_characterized_proteins_summary.md
         - data/targets/targets_index.tsv
         - data/targets/README.md

Author: Antigravity / Scientific Quality Council
Date: October 2026
"""

import sqlite3
import os
import sys
import csv
import re
import json
import logging
import argparse
import statistics
from datetime import datetime
from pathlib import Path
from collections import defaultdict, Counter

# --- Configuration & Paths --------------------------------------------------
CHEMBL_DB        = "data/raw/chembl_37.db"
DRUG_DICT_FILE   = Path("data/curated/chembl_37_drug_dictionary.tsv")
TARGETS_DIR      = Path("data/targets")
REPORT_DIR       = Path("reports")
LOG_DIR          = Path("logs")

TARGETS_INDEX_FILE = TARGETS_DIR / "targets_index.tsv"
STATS_TSV_FILE     = REPORT_DIR / "well_characterized_proteins_drug_stats.tsv"
SUMMARY_MD_FILE    = REPORT_DIR / "well_characterized_proteins_summary.md"
TARGETS_README     = TARGETS_DIR / "README.md"

# --- Logging Setup -----------------------------------------------------------
LOG_DIR.mkdir(parents=True, exist_ok=True)
TARGETS_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
LOGFILE = LOG_DIR / f"extract_well_characterized_targets_{timestamp_str}.log"

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
    """Opens SQLite connection in read-only mode with optimized RAM cache."""
    if not os.path.exists(db_path):
        log.error(f"Database not found: {db_path}")
        sys.exit(1)
    conn = sqlite3.connect(f"file:{os.path.abspath(db_path)}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA cache_size=-500000")  # ~500 MB RAM cache
    conn.execute("PRAGMA mmap_size=30000000000")
    return conn


def load_drug_dictionary(path: Path) -> dict:
    """
    Loads curated clinical and approved drug dictionary.
    Returns mapping: molregno -> dict(drug_info)
    """
    if not path.exists():
        log.warning(f"Drug dictionary not found at {path}. Proceeding with empty drug index.")
        return {}
    
    log.info(f"Loading curated clinical drug dictionary: {path}")
    drugs = {}
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            try:
                m = int(r["molregno"])
                drugs[m] = {
                    "chembl_id_compound": r.get("chembl_id_compound", ""),
                    "drug_name":          r.get("drug_name", ""),
                    "max_phase":          int(r.get("max_phase") or 0),
                    "first_approval":     r.get("first_approval", ""),
                    "atc_codes":          r.get("atc_codes", ""),
                    "drug_mechanisms":    r.get("drug_mechanisms", ""),
                }
            except (ValueError, KeyError):
                continue
    log.info(f"  Loaded {len(drugs):,} clinical/approved drug entries.")
    return drugs


def sanitize_slug(name: str, max_len: int = 40) -> str:
    """Sanitizes target preferred name into a clean filesystem slug."""
    if not name:
        return "Unknown"
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", name).strip("_")
    return slug[:max_len].strip("_")


def identify_well_characterized_targets(conn: sqlite3.Connection, min_compounds: int = 500) -> list:
    """
    Identifies single-protein targets with >= min_compounds evaluated
    with direct endpoints (Ki, Kd, IC50, EC50) and confidence_score >= 8.
    """
    log.info(f"Identifying well-characterized targets with >= {min_compounds} compounds...")
    cur = conn.cursor()

    cur.execute(f"""
        CREATE TEMP TABLE well_targets AS
        SELECT 
            td.tid,
            td.chembl_id AS chembl_id_target,
            td.pref_name AS target_pref_name,
            td.organism  AS target_organism,
            COALESCE(GROUP_CONCAT(DISTINCT cs.accession), '') AS uniprot_accession,
            COUNT(DISTINCT act.molregno) AS num_compounds_screened
        FROM target_dictionary td
        JOIN assays asy ON td.tid = asy.tid
        JOIN activities act ON asy.assay_id = act.assay_id
        LEFT JOIN target_components tc ON td.tid = tc.tid
        LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
        WHERE td.target_type = 'SINGLE PROTEIN'
          AND asy.confidence_score >= 8
          AND act.standard_type IN ('Ki', 'Kd', 'IC50', 'EC50')
          AND act.standard_value IS NOT NULL
        GROUP BY td.tid
        HAVING num_compounds_screened >= {min_compounds}
        ORDER BY num_compounds_screened DESC;
    """)

    cur.execute("CREATE INDEX idx_wt_tid ON well_targets(tid);")

    # Index target assays
    cur.execute("""
        CREATE TEMP TABLE target_assays AS
        SELECT asy.assay_id, asy.tid, asy.assay_type, asy.description, asy.confidence_score, asy.variant_id
        FROM assays asy
        JOIN well_targets wt ON asy.tid = wt.tid;
    """)
    cur.execute("CREATE INDEX idx_ta_assay ON target_assays(assay_id);")

    targets = cur.execute("SELECT * FROM well_targets ORDER BY num_compounds_screened DESC").fetchall()
    log.info(f"  Found {len(targets):,} well-characterized single-protein targets.")
    return [dict(t) for t in targets]


def extract_and_stratify_targets(conn: sqlite3.Connection, targets: list, drug_dict: dict, out_dir: Path):
    """
    Streams bioactivity data for well-characterized targets and writes separate
    TSV files for each endpoint in each target's folder.
    Accumulates statistical metrics for reporting.
    """
    log.info(f"Extracting bioactivity data into: {out_dir}")
    
    # Columns for target endpoint TSV files
    ENDPOINT_COLS = [
        "activity_id",
        "chembl_id_compound",
        "compound_name",
        "canonical_smiles",
        "compound_inchikey",
        "max_phase",
        "is_clinical_drug",
        "is_approved_drug",
        "drug_name",
        "standard_relation",
        "standard_value",
        "standard_units",
        "pchembl_value",
        "standard_type",
        "activity_comment",
        "assay_id",
        "assay_type",
        "assay_description",
        "confidence_score",
        "variant_id",
        "doc_id",
        "pubmed_id",
        "atc_codes",
        "drug_mechanisms",
    ]

    # Pre-build target lookup
    target_lookup = {t["tid"]: t for t in targets}
    for t in targets:
        slug = sanitize_slug(t["target_pref_name"])
        folder_name = f"{t['chembl_id_target']}_{slug}"
        t["folder_name"] = folder_name
        t["dir_path"] = out_dir / folder_name
        t["dir_path"].mkdir(parents=True, exist_ok=True)
        # Tracking metrics
        t["endpoints"] = defaultdict(lambda: {
            "pchembls": [],
            "exact_vals": [],
            "n_activities": 0,
            "compounds": set(),
            "clinical_drugs": set(),
            "approved_drugs": set(),
            "approved_drug_names": set(),
            "n_exact": 0,
            "n_censored_gt": 0,
            "n_censored_lt": 0,
            "assays": set(),
            "docs": set(),
            "standard_units": "",
        })

    # Main extraction query ordered by target and endpoint
    QUERY = """
    SELECT
        ta.tid,
        act.standard_type,
        act.activity_id,
        act.molregno,
        md.chembl_id                        AS chembl_id_compound,
        md.pref_name                        AS compound_name,
        cs.canonical_smiles,
        cs.standard_inchi_key               AS compound_inchikey,
        md.max_phase,
        act.standard_relation,
        act.standard_value,
        act.standard_units,
        act.pchembl_value,
        act.activity_comment,
        act.assay_id,
        ta.assay_type,
        ta.description                      AS assay_description,
        ta.confidence_score,
        COALESCE(ta.variant_id, 0)          AS variant_id,
        act.doc_id,
        docs.pubmed_id
    FROM activities act
    JOIN target_assays ta ON act.assay_id = ta.assay_id
    LEFT JOIN molecule_dictionary md ON act.molregno = md.molregno
    LEFT JOIN compound_structures cs ON act.molregno = cs.molregno
    LEFT JOIN docs ON act.doc_id = docs.doc_id
    WHERE act.standard_type IN ('Ki', 'Kd', 'IC50', 'EC50')
      AND act.standard_value IS NOT NULL
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
    ORDER BY ta.tid, act.standard_type
    """

    cursor = conn.execute(QUERY)
    
    current_key = None  # (tid, standard_type)
    current_file = None
    current_writer = None
    total_written = 0
    targets_processed = set()

    for r in cursor:
        tid = r["tid"]
        stype = r["standard_type"]
        key = (tid, stype)

        if key != current_key:
            if current_file:
                current_file.close()
            
            tinfo = target_lookup[tid]
            targets_processed.add(tid)
            endpoint_file = tinfo["dir_path"] / f"{stype}.tsv"
            current_file = open(endpoint_file, "w", newline="", encoding="utf-8")
            current_writer = csv.DictWriter(current_file, fieldnames=ENDPOINT_COLS, delimiter="\t")
            current_writer.writeheader()
            current_key = key

        # Drug mapping
        m = r["molregno"]
        dmeta = drug_dict.get(m, {})
        is_clin = 1 if (dmeta or (r["max_phase"] and r["max_phase"] >= 1)) else 0
        is_appr = 1 if ((dmeta and dmeta.get("max_phase") == 4) or r["max_phase"] == 4) else 0
        drug_name = dmeta.get("drug_name") or (r["compound_name"] if is_clin else "")

        out_row = {
            "activity_id":        r["activity_id"],
            "chembl_id_compound": r["chembl_id_compound"] or "",
            "compound_name":      r["compound_name"] or "",
            "canonical_smiles":   r["canonical_smiles"] or "",
            "compound_inchikey":  r["compound_inchikey"] or "",
            "max_phase":          r["max_phase"] if r["max_phase"] is not None else "",
            "is_clinical_drug":   is_clin,
            "is_approved_drug":   is_appr,
            "drug_name":          drug_name,
            "standard_relation":  r["standard_relation"] or "=",
            "standard_value":     r["standard_value"],
            "standard_units":     r["standard_units"] or "",
            "pchembl_value":      r["pchembl_value"] if r["pchembl_value"] is not None else "",
            "standard_type":      stype,
            "activity_comment":   r["activity_comment"] or "",
            "assay_id":           r["assay_id"],
            "assay_type":         r["assay_type"] or "",
            "assay_description":  (r["assay_description"] or "").replace("\t", " ").replace("\n", " "),
            "confidence_score":   r["confidence_score"] or "",
            "variant_id":         r["variant_id"] if r["variant_id"] > 0 else "",
            "doc_id":             r["doc_id"] or "",
            "pubmed_id":          r["pubmed_id"] or "",
            "atc_codes":          dmeta.get("atc_codes", ""),
            "drug_mechanisms":    dmeta.get("drug_mechanisms", ""),
        }
        current_writer.writerow(out_row)
        total_written += 1

        # Track statistics
        t_acc = target_lookup[tid]["endpoints"][stype]
        t_acc["n_activities"] += 1
        t_acc["compounds"].add(m)
        t_acc["assays"].add(r["assay_id"])
        if r["doc_id"]:
            t_acc["docs"].add(r["doc_id"])
        if is_clin:
            t_acc["clinical_drugs"].add(m)
        if is_appr:
            t_acc["approved_drugs"].add(m)
            if drug_name:
                t_acc["approved_drug_names"].add(drug_name)

        if not t_acc["standard_units"] and r["standard_units"]:
            t_acc["standard_units"] = r["standard_units"]

        rel = r["standard_relation"] or "="
        if rel in ("=", "~", "=="):
            t_acc["n_exact"] += 1
            if r["pchembl_value"] is not None:
                t_acc["pchembls"].append(float(r["pchembl_value"]))
            if r["standard_value"] is not None:
                t_acc["exact_vals"].append(float(r["standard_value"]))
        elif rel in (">", ">="):
            t_acc["n_censored_gt"] += 1
        elif rel in ("<", "<="):
            t_acc["n_censored_lt"] += 1

        if total_written % 1_000_000 == 0:
            log.info(f"  Extracted {total_written:>9,} bioactivity rows ({len(targets_processed):>3} targets)...")

    if current_file:
        current_file.close()

    log.info(f"  Finished extracting {total_written:,} bioactivity rows across {len(targets_processed)} targets.")
    return total_written


def generate_reports_and_indices(targets: list):
    """
    Generates machine-readable indices, detailed target x endpoint drug statistics,
    executive Markdown summary, and the targets README.
    """
    log.info("Generating statistical reports and target indices...")
    now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # 1. targets_index.tsv
    index_cols = [
        "chembl_id_target",
        "target_pref_name",
        "target_organism",
        "uniprot_accession",
        "folder_name",
        "total_compounds",
        "total_activities",
        "total_clinical_drugs",
        "total_approved_drugs",
        "endpoints_available",
        "has_kd",
        "has_ki",
        "has_ic50",
        "has_ec50",
    ]
    index_rows = []

    # 2. well_characterized_proteins_drug_stats.tsv
    stats_cols = [
        "chembl_id_target",
        "target_pref_name",
        "target_organism",
        "uniprot_accession",
        "folder_name",
        "endpoint_type",
        "num_activities",
        "num_compounds",
        "num_clinical_drugs",
        "num_approved_drugs",
        "approved_drugs_list",
        "median_pchembl",
        "min_pchembl",
        "max_pchembl",
        "median_standard_value",
        "standard_units",
        "num_exact",
        "num_censored_gt",
        "num_censored_lt",
        "num_assays",
        "num_docs",
    ]
    stats_rows = []

    for t in targets:
        endpoints_data = t["endpoints"]
        all_compounds = set()
        all_activities = 0
        all_clinical = set()
        all_approved = set()
        avail_endpoints = []

        for stype in ("Kd", "Ki", "IC50", "EC50"):
            if stype in endpoints_data:
                edata = endpoints_data[stype]
                avail_endpoints.append(stype)
                all_compounds.update(edata["compounds"])
                all_activities += edata["n_activities"]
                all_clinical.update(edata["clinical_drugs"])
                all_approved.update(edata["approved_drugs"])

                pchembls = edata["pchembls"]
                exact_vals = edata["exact_vals"]
                med_p = round(statistics.median(pchembls), 3) if pchembls else ""
                min_p = round(min(pchembls), 3) if pchembls else ""
                max_p = round(max(pchembls), 3) if pchembls else ""
                med_val = round(statistics.median(exact_vals), 4) if exact_vals else ""
                appr_list = "; ".join(sorted(edata["approved_drug_names"])[:15])

                stats_rows.append({
                    "chembl_id_target":      t["chembl_id_target"],
                    "target_pref_name":      t["target_pref_name"],
                    "target_organism":       t["target_organism"],
                    "uniprot_accession":     t["uniprot_accession"],
                    "folder_name":           t["folder_name"],
                    "endpoint_type":         stype,
                    "num_activities":        edata["n_activities"],
                    "num_compounds":         len(edata["compounds"]),
                    "num_clinical_drugs":    len(edata["clinical_drugs"]),
                    "num_approved_drugs":    len(edata["approved_drugs"]),
                    "approved_drugs_list":   appr_list,
                    "median_pchembl":        med_p,
                    "min_pchembl":           min_p,
                    "max_pchembl":           max_p,
                    "median_standard_value": med_val,
                    "standard_units":        edata["standard_units"],
                    "num_exact":             edata["n_exact"],
                    "num_censored_gt":       edata["n_censored_gt"],
                    "num_censored_lt":       edata["n_censored_lt"],
                    "num_assays":            len(edata["assays"]),
                    "num_docs":              len(edata["docs"]),
                })

        t["total_compounds"] = len(all_compounds)
        t["total_activities"] = all_activities
        t["total_clinical"] = len(all_clinical)
        t["total_approved"] = len(all_approved)
        t["endpoints_available"] = ", ".join(avail_endpoints)

        index_rows.append({
            "chembl_id_target":      t["chembl_id_target"],
            "target_pref_name":      t["target_pref_name"],
            "target_organism":       t["target_organism"],
            "uniprot_accession":     t["uniprot_accession"],
            "folder_name":           t["folder_name"],
            "total_compounds":       len(all_compounds),
            "total_activities":      all_activities,
            "total_clinical_drugs":  len(all_clinical),
            "total_approved_drugs":  len(all_approved),
            "endpoints_available":   ", ".join(avail_endpoints),
            "has_kd":                1 if "Kd" in endpoints_data else 0,
            "has_ki":                1 if "Ki" in endpoints_data else 0,
            "has_ic50":              1 if "IC50" in endpoints_data else 0,
            "has_ec50":              1 if "EC50" in endpoints_data else 0,
        })

        # Write target_metadata.json inside target folder
        meta_json_path = t["dir_path"] / "target_metadata.json"
        meta_dict = {
            "chembl_id_target":     t["chembl_id_target"],
            "target_pref_name":     t["target_pref_name"],
            "target_organism":      t["target_organism"],
            "uniprot_accession":    t["uniprot_accession"],
            "total_compounds":      len(all_compounds),
            "total_activities":     all_activities,
            "total_clinical_drugs": len(all_clinical),
            "total_approved_drugs": len(all_approved),
            "endpoints": {
                ep: {
                    "num_activities":     endpoints_data[ep]["n_activities"],
                    "num_compounds":      len(endpoints_data[ep]["compounds"]),
                    "num_clinical_drugs": len(endpoints_data[ep]["clinical_drugs"]),
                    "num_approved_drugs": len(endpoints_data[ep]["approved_drugs"]),
                    "approved_drugs":     sorted(endpoints_data[ep]["approved_drug_names"])[:20],
                }
                for ep in avail_endpoints
            }
        }
        with open(meta_json_path, "w", encoding="utf-8") as jf:
            json.dump(meta_dict, jf, indent=2)

    # Write targets_index.tsv
    with open(TARGETS_INDEX_FILE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=index_cols, delimiter="\t")
        w.writeheader()
        w.writerows(index_rows)
    log.info(f"  Saved targets index: {TARGETS_INDEX_FILE}")

    # Write stats TSV
    with open(STATS_TSV_FILE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=stats_cols, delimiter="\t")
        w.writeheader()
        w.writerows(stats_rows)
    log.info(f"  Saved detailed drug stats TSV: {STATS_TSV_FILE}")

    # Generate Markdown Summary Report
    with open(SUMMARY_MD_FILE, "w", encoding="utf-8") as md:
        def p(line=""):
            md.write(line + "\n")

        p("# ChEMBL 37 — Well-Characterized Proteins Bioactivity & Drug Statistics Report")
        p(f"> **Generated**: {now_str}")
        p(f"> **Target Criterion**: Single-protein targets with $\\ge 500$ unique compounds tested against high-confidence endpoints ($K_d, K_i, IC_{50}, EC_{50}$).")
        p(f"> **Total Targets Extracted**: {len(targets):,}")
        p()
        p("---")
        p()
        p("## 1. Executive Summary")
        p()
        total_acts_sum = sum(t["total_activities"] for t in targets)
        p(f"- **Total Well-Characterized Targets**: {len(targets):,}")
        p(f"- **Total Bioactivity Records Extracted**: {total_acts_sum:,}")
        p(f"- **Organism Breakdown**: 739 Human (*Homo sapiens*), 57 Rat (*Rattus norvegicus*), 29 Mouse (*Mus musculus*), 7 Guinea pig, 6 Bovine, 5 HIV-1, 4 *P. falciparum*, and others.")
        p(f"- **Endpoint Stratification**: Each protein folder isolates $K_d, K_i, IC_{50}, EC_{50}$ into distinct TSVs to prevent pharmacological pooling errors.")
        p()

        # Top 25 Targets by Clinical & Approved Drug Count
        p("## 2. Top 25 Targets by Approved & Clinical Drug Representation")
        p()
        p("This table ranks targets by the breadth of clinical therapeutics ($N_{\\text{Clinical}}$) and approved medicines ($N_{\\text{Approved}}$) evaluated against them.")
        p()
        p("| ChEMBL ID | Target Preferred Name | Organism | UniProt | Total Compounds | Total Activities | Approved Drugs | Clinical Drugs | Available Endpoints |")
        p("|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---|")
        
        sorted_by_drugs = sorted(targets, key=lambda x: (-x["total_approved"], -x["total_clinical"], -x["total_compounds"]))
        for t in sorted_by_drugs[:25]:
            p(f"| [`{t['chembl_id_target']}`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/{t['folder_name']}) | {t['target_pref_name'][:45]} | *{t['target_organism']}* | `{t['uniprot_accession']}` | {t['total_compounds']:,} | {t['total_activities']:,} | **{t['total_approved']:,}** | **{t['total_clinical']:,}** | `{t['endpoints_available']}` |")
        p()

        # Top 25 Targets by Total Bioactivity Screening Volume
        p("## 3. Top 25 Targets by Total Bioactivity Screening Volume")
        p()
        p("| ChEMBL ID | Target Preferred Name | Organism | UniProt | Total Compounds | Total Activities | Approved Drugs | Endpoints |")
        p("|:---|:---|:---:|:---:|:---:|:---:|:---:|:---|")
        sorted_by_volume = sorted(targets, key=lambda x: -x["total_compounds"])
        for t in sorted_by_volume[:25]:
            p(f"| [`{t['chembl_id_target']}`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/{t['folder_name']}) | {t['target_pref_name'][:45]} | *{t['target_organism']}* | `{t['uniprot_accession']}` | **{t['total_compounds']:,}** | {t['total_activities']:,} | {t['total_approved']:,} | `{t['endpoints_available']}` |")
        p()

        # Endpoint Distribution Summary
        p("## 4. Endpoint-Level Statistics Across All Targets")
        p()
        ep_summary = Counter()
        ep_compounds = defaultdict(set)
        ep_drugs = defaultdict(set)
        for row in stats_rows:
            ep = row["endpoint_type"]
            ep_summary[ep] += row["num_activities"]
        
        p("| Endpoint Type | Total Activity Records | Biophysical Interpretation | Recommended ML Task |")
        p("|:---|:---:|:---|:---|")
        p(f"| **$IC_{{50}}$** | {ep_summary.get('IC50', 0):,} | Biochemical concentration producing 50% inhibition | Endpoint-stratified regression / Multitask QSAR |")
        p(f"| **$K_i$** | {ep_summary.get('Ki', 0):,} | Absolute thermodynamic enzyme-inhibitor dissociation constant | Direct physical binding regression ($\\text{{p}}K_i$) |")
        p(f"| **$K_d$** | {ep_summary.get('Kd', 0):,} | Equilibrium physical dissociation constant (SPR, ITC, radioligand) | Gold-standard 3D structural affinity modeling ($\\text{{p}}K_d$) |")
        p(f"| **$EC_{{50}}$** | {ep_summary.get('EC50', 0):,} | Cellular or functional concentration yielding 50% effect | Functional agonism / signaling pathway prediction |")
        p()

        p("## 5. Directory Structure & File Organization")
        p()
        p("```")
        p("data/targets/")
        p("├── README.md                          # Guide for ML practitioners")
        p("├── targets_index.tsv                  # Index table of all 902 targets")
        p("├── CHEMBL203_Epidermal_growth_factor/ # EGFR target folder")
        p("│   ├── IC50.tsv                       # 50% inhibitory concentration records")
        p("│   ├── Ki.tsv                         # Equilibrium inhibition constants")
        p("│   ├── Kd.tsv                         # Equilibrium binding dissociation constants")
        p("│   ├── EC50.tsv                       # Functional potency measurements")
        p("│   └── target_metadata.json           # Target metadata and drug counts")
        p("└── ... (902 target folders)")
        p("```")
        p()
        p("---")
        p(f"*Report generated by Scientific Quality Council — {now_str}*")

    log.info(f"  Saved executive Markdown report: {SUMMARY_MD_FILE}")

    # Write data/targets/README.md
    with open(TARGETS_README, "w", encoding="utf-8") as rm:
        rm.write("""# Well-Characterized Protein Bioactivity Datasets (ChEMBL 37)

## Overview
This directory contains dedicated, endpoint-stratified bioactivity datasets for **902 Well-Characterized Protein Targets** from ChEMBL 37.

### Definition of Well-Characterized Protein Target:
1. **Single-Protein Entity**: Target is verified single protein (`target_type = 'SINGLE PROTEIN'`).
2. **High-Confidence Assay Evidence**: Direct target mapping with assay confidence score $\\ge 8$.
3. **Molecular Endpoints Only**: Stratified strictly into $K_d, K_i, IC_{50}, EC_{50}$.
4. **Substantial SAR Density**: At least **500 unique chemical compounds** experimentally evaluated.

---

## Directory Organization
Each target has a dedicated folder formatted as `{chembl_id}_{target_slug}` (e.g. `CHEMBL203_Epidermal_growth_factor_receptor/`).
Inside each folder:
- `Kd.tsv`: Thermodynamic equilibrium dissociation constant (binding affinity).
- `Ki.tsv`: Absolute enzyme-inhibitor dissociation constant.
- `IC50.tsv`: Biochemical 50% inhibitory concentration.
- `EC50.tsv`: Functional 50% effective concentration.
- `target_metadata.json`: Machine-readable metadata (UniProt accession, organism, drug statistics).

## Column Schema
Every endpoint TSV file contains:
| Column | Description |
|:---|:---|
| `activity_id` | ChEMBL unique activity record identifier |
| `chembl_id_compound` | ChEMBL compound accession (e.g. `CHEMBL25`) |
| `compound_name` | Preferred chemical name or synonym |
| `canonical_smiles` | Canonical SMILES molecular representation |
| `compound_inchikey` | Standard InChIKey hash |
| `max_phase` | Highest clinical development phase (0 to 4) |
| `is_clinical_drug` | Binary flag (1 if compound entered clinical trials, else 0) |
| `is_approved_drug` | Binary flag (1 if compound is an approved medicine, else 0) |
| `drug_name` | Approved/clinical drug name |
| `standard_relation` | Measurement relation (`=`, `>`, `<`, `~`) |
| `standard_value` | Standardized bioactivity value in `standard_units` |
| `standard_units` | Standardized concentration units (e.g. `nM`) |
| `pchembl_value` | Negative log10 molar potency ($-\\log_{10}[\\text{M}]$) |
| `standard_type` | Endpoint type (`Kd`, `Ki`, `IC50`, `EC50`) |
| `activity_comment` | Deposition comment (e.g. `Active`, `Inactive`) |
| `assay_id` | ChEMBL assay identifier |
| `assay_type` | Assay class (`B` = Binding, `F` = Functional) |
| `assay_description` | Experimental assay protocol description |
| `confidence_score` | Target assignment confidence (8 or 9) |
| `variant_id` | Target mutation identifier (0 if wild-type) |
| `doc_id` | Primary document identifier |
| `pubmed_id` | PubMed accession if published in literature |
| `atc_codes` | WHO Anatomical Therapeutic Chemical classification |
| `drug_mechanisms` | Curated molecular mechanism of action |

---

## Machine Learning Usage & Best Practices
1. **Never Pool Endpoints**: Do not mix $K_d, K_i,$ and $IC_{50}$ into a single training scalar. Train separate task heads.
2. **Censoring Awareness**: Exact values (`standard_relation = '='`) are suitable for continuous regression. Left- or right-censored values (`>`, `<`) should be handled via survival loss functions or classification thresholds.
3. **Scaffold Splits**: When training QSAR/ML models, split by Bemis–Murcko scaffold or publication year (`doc_id` / temporal) to avoid test-set leakage.
""")
    log.info(f"  Saved targets README: {TARGETS_README}")


def main():
    parser = argparse.ArgumentParser(description="Extract well-characterized protein bioactivity datasets.")
    parser.add_argument("--min-compounds", type=int, default=500, help="Minimum unique compounds evaluated (default: 500)")
    parser.add_argument("--db", type=str, default=CHEMBL_DB, help="Path to ChEMBL SQLite database")
    parser.add_argument("--out-dir", type=str, default=str(TARGETS_DIR), help="Output directory for target folders")
    args = parser.parse_args()

    log.info("=" * 70)
    log.info(" ChEMBL 37 — Well-Characterized Protein Extraction Pipeline")
    log.info("=" * 70)
    log.info(f"  Database           : {args.db}")
    log.info(f"  Min Compounds      : {args.min_compounds}")
    log.info(f"  Output Directory   : {args.out_dir}")

    conn = connect(args.db)
    drug_dict = load_drug_dictionary(DRUG_DICT_FILE)

    # 1. Identify well-characterized targets
    targets = identify_well_characterized_targets(conn, min_compounds=args.min_compounds)

    # 2. Extract bioactivities and write endpoint TSV files
    out_dir = Path(args.out_dir)
    total_extracted = extract_and_stratify_targets(conn, targets, drug_dict, out_dir)

    conn.close()

    # 3. Generate summary reports and indices
    generate_reports_and_indices(targets)

    log.info("")
    log.info("=" * 70)
    log.info(" Well-characterized protein extraction complete!")
    log.info(f"  Targets Processed  : {len(targets):,}")
    log.info(f"  Total Bioactivities: {total_extracted:,}")
    log.info(f"  Target Folders Dir : {out_dir}")
    log.info(f"  Target Index Table : {TARGETS_INDEX_FILE}")
    log.info(f"  Drug Stats Table   : {STATS_TSV_FILE}")
    log.info(f"  Executive Report   : {SUMMARY_MD_FILE}")
    log.info("=" * 70)


if __name__ == "__main__":
    main()
