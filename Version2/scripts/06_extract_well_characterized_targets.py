#!/usr/bin/env python3
"""
06_extract_well_characterized_targets.py — Target Partitioning & Endpoint Extraction (Version 2)
==================================================================================================
Audited for Scientific Quality Council standards.
Addresses Reviewer 2 findings:
  - C01: Enforces identical qualification and export predicate (confidence >= 8, valid records).
         Eliminates leakage of low-confidence or invalid records into exports.
  - C02: Preserves original measurement relations; never coerces missing relations to '='.
  - C03: Converts concentration units to a unified molar scale (nM) before computing medians;
         summarizes pChEMBL independently; separates incompatible non-molar dimensions.
  - H03: Retains standard_upper_value, data_validity_comment, potential_duplicate, doc_id, doi, year.
  - H04: Float-safe clinical phase parsing (handles '4.0', '1.0', '0.5', -1).
  - H07: Explicitly labels missing variants as 'Unspecified' rather than confirmed wild-type.
  - H08: Derives all report numbers and counts dynamically from execution without hardcoding.
"""

import os
import sys
import re
import csv
import json
import sqlite3
import argparse
import datetime
import statistics
from collections import defaultdict
from tqdm import tqdm

DB_PATH = "data/raw/chembl_37.db"
DEFAULT_OUT_DIR = "Version2/data/targets"
REPORTS_DIR = "Version2/reports"

ENDPOINTS = ["Kd", "Ki", "IC50", "EC50"]

def parse_phase(val):
    if val is None or val == "":
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None

def normalize_to_molar(val, units):
    if val is None or units is None:
        return None
    u = units.strip().lower()
    try:
        v = float(val)
        if v <= 0:
            return None
        if u in ("nm", "nanomolar", "nmol/l"):
            return v * 1e-9
        elif u in ("um", "micromolar", "umol/l", "µm"):
            return v * 1e-6
        elif u in ("mm", "millimolar", "mmol/l"):
            return v * 1e-3
        elif u in ("pm", "picomolar", "pmol/l"):
            return v * 1e-12
        elif u in ("m", "molar", "mol/l"):
            return v * 1.0
        return None  # Incompatible unit (e.g. ug/ml, %)
    except (ValueError, TypeError):
        return None

def sanitize_slug(name):
    clean = re.sub(r'[^A-Za-z0-9]+', '_', name or 'target').strip('_')
    return clean[:50]

def load_drug_dictionary(conn):
    print("Loading drug dictionary with float-safe phase parsing...")
    drugs = {}
    c = conn.cursor()
    c.execute("""
        SELECT molregno, chembl_id, pref_name, max_phase, therapeutic_flag
        FROM molecule_dictionary
        WHERE max_phase IS NOT NULL OR therapeutic_flag = 1;
    """)
    for molregno, chembl_id, pref_name, max_phase, th_flag in c.fetchall():
        phase_num = parse_phase(max_phase)
        is_approved = (phase_num is not None and phase_num >= 4.0)
        is_clinical = (phase_num is not None and phase_num >= 1.0) or (th_flag == 1)
        drugs[molregno] = {
            "chembl_id": chembl_id,
            "pref_name": pref_name or chembl_id,
            "max_phase": phase_num,
            "is_approved": is_approved,
            "is_clinical": is_clinical
        }
    print(f"Loaded {len(drugs):,} drug / clinical candidate entries.")
    return drugs

def main():
    parser = argparse.ArgumentParser(description="Extract well-characterized targets with rigorous qualification (V2)")
    parser.add_argument("--min-compounds", type=int, default=500, help="Minimum unique compounds required to qualify")
    parser.add_argument("--out-dir", type=str, default=DEFAULT_OUT_DIR, help="Output target directory")
    args = parser.parse_args()

    min_compounds = args.min_compounds
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(REPORTS_DIR, exist_ok=True)

    print(f"[{datetime.datetime.now()}] Starting Well-Characterized Target Extraction (V2)...")
    print(f"Configuration: min_compounds={min_compounds}, out_dir={out_dir}")

    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    c = conn.cursor()

    drug_dict = load_drug_dictionary(conn)

    # UNIFIED ELIGIBLE PREDICATE (used for both qualification and export)
    # 1. Target qualification query
    print("Step 1: Identifying qualified single-protein targets under strict validity predicate...")
    qual_query = f"""
    SELECT
        td.tid,
        td.chembl_id AS chembl_id_target,
        td.pref_name AS target_pref_name,
        td.organism AS target_organism,
        cs.accession AS uniprot_accession,
        pc.pref_name AS protein_class_name,
        count(DISTINCT act.molregno) AS num_compounds,
        count(*) AS num_activities
    FROM activities act
    JOIN assays a ON act.assay_id = a.assay_id
    JOIN target_dictionary td ON a.tid = td.tid
    LEFT JOIN target_components tc ON td.tid = tc.tid
    LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
    LEFT JOIN component_class cc ON tc.component_id = cc.component_id
    LEFT JOIN protein_classification pc ON cc.protein_class_id = pc.protein_class_id
    WHERE td.target_type = 'SINGLE PROTEIN'
      AND a.confidence_score >= 8
      AND act.standard_type IN ('Kd', 'Ki', 'IC50', 'EC50')
      AND act.standard_value IS NOT NULL
      AND (act.data_validity_comment IS NULL OR act.data_validity_comment NOT IN (
          'Potential transcription error', 'Outside typical range',
          'Non standard unit for type', 'Author confirmed error', 'Potential author error'
      ))
    GROUP BY td.tid
    HAVING num_compounds >= ?
    ORDER BY num_compounds DESC;
    """

    c.execute(qual_query, (min_compounds,))
    qualified_targets = c.fetchall()
    print(f"Total qualified targets identified: {len(qualified_targets):,}")

    target_stats_rows = []
    targets_index_rows = []
    total_exported_activities = 0

    # Step 2: Per-target extraction
    export_query = """
    SELECT
        act.activity_id,
        md.chembl_id AS chembl_id_compound,
        cs.canonical_smiles,
        cs.standard_inchi_key,
        act.standard_type,
        act.standard_relation,
        act.standard_value,
        act.standard_units,
        act.standard_upper_value,
        act.pchembl_value,
        act.activity_comment,
        act.data_validity_comment,
        act.potential_duplicate,
        a.variant_id,
        a.assay_id,
        a.chembl_id AS chembl_id_assay,
        a.assay_type,
        a.confidence_score,
        a.description AS assay_description,
        d.doc_id,
        d.year,
        d.doi,
        act.molregno
    FROM activities act
    JOIN molecule_dictionary md ON act.molregno = md.molregno
    LEFT JOIN compound_structures cs ON act.molregno = cs.molregno
    JOIN assays a ON act.assay_id = a.assay_id
    JOIN target_dictionary td ON a.tid = td.tid
    LEFT JOIN docs d ON act.doc_id = d.doc_id
    WHERE td.tid = ?
      AND td.target_type = 'SINGLE PROTEIN'
      AND a.confidence_score >= 8
      AND act.standard_type = ?
      AND act.standard_value IS NOT NULL
      AND (act.data_validity_comment IS NULL OR act.data_validity_comment NOT IN (
          'Potential transcription error', 'Outside typical range',
          'Non standard unit for type', 'Author confirmed error', 'Potential author error'
      ))
    ORDER BY act.molregno, act.activity_id;
    """

    tsv_header = [
        "activity_id", "chembl_id_compound", "canonical_smiles", "standard_inchi_key",
        "standard_type", "standard_relation", "standard_value", "standard_units",
        "standard_value_nM", "standard_upper_value", "pchembl_value",
        "activity_comment", "data_validity_comment", "potential_duplicate",
        "variant_annotation", "chembl_id_assay", "assay_type", "confidence_score",
        "assay_description", "doc_id", "year", "doi"
    ]

    for tid, target_chembl, target_name, organism, uniprot, prot_class, n_cmpds, n_acts in tqdm(qualified_targets, desc="Exporting Target Directories"):
        slug = sanitize_slug(target_name)
        folder_name = f"{target_chembl}_{slug}"
        target_dir = os.path.join(out_dir, folder_name)
        os.makedirs(target_dir, exist_ok=True)

        target_compounds_all = set()
        target_approved_drugs = set()
        target_clinical_drugs = set()
        endpoint_summary = {}

        for ep in ENDPOINTS:
            c.execute(export_query, (tid, ep))
            rows = c.fetchall()
            if not rows:
                continue

            ep_compounds = set()
            ep_approved = set()
            ep_clinical = set()
            ep_molar_nM_values = []
            ep_pchembl_values = []
            ep_exact_count = 0
            ep_approx_count = 0
            ep_upper_count = 0
            ep_lower_count = 0
            ep_unknown_count = 0

            ep_tsv_path = os.path.join(target_dir, f"{ep}.tsv")
            with open(ep_tsv_path, "w", encoding="utf-8") as ep_f:
                writer = csv.writer(ep_f, delimiter="\t", lineterminator="\n")
                writer.writerow(tsv_header)

                for r in rows:
                    (aid, cid, smi, inchikey, st_type, st_rel, st_val, st_units,
                     st_upper, pchembl, comment, val_comment, pot_dup, var_id,
                     assay_id, assay_chembl, atype, conf, desc, doc_id, year, doi, molregno) = r

                    target_compounds_all.add(molregno)
                    ep_compounds.add(molregno)

                    # Drug status
                    if molregno in drug_dict:
                        d_info = drug_dict[molregno]
                        if d_info["is_approved"]:
                            target_approved_drugs.add(d_info["pref_name"])
                            ep_approved.add(d_info["pref_name"])
                        if d_info["is_clinical"]:
                            target_clinical_drugs.add(d_info["pref_name"])
                            ep_clinical.add(d_info["pref_name"])

                    # Relation classification without coercing missing to '='
                    rel_clean = (st_rel or "").strip()
                    if rel_clean == "=":
                        ep_exact_count += 1
                    elif rel_clean == "~":
                        ep_approx_count += 1
                    elif rel_clean in ("<", "<="):
                        ep_upper_count += 1
                    elif rel_clean in (">", ">="):
                        ep_lower_count += 1
                    else:
                        ep_unknown_count += 1

                    # Unit conversion to normalized nM
                    molar_val = normalize_to_molar(st_val, st_units)
                    molar_nM_str = ""
                    if molar_val is not None:
                        val_nM = molar_val * 1e9
                        ep_molar_nM_values.append(val_nM)
                        molar_nM_str = f"{val_nM:.4f}"

                    if pchembl is not None:
                        try:
                            ep_pchembl_values.append(float(pchembl))
                        except (ValueError, TypeError):
                            pass

                    # Variant explicit labeling
                    var_label = f"Variant_{var_id}" if var_id else "Unspecified"

                    writer.writerow([
                        aid, cid, smi or "", inchikey or "",
                        st_type, rel_clean if rel_clean else "unspecified",
                        st_val, st_units or "", molar_nM_str, st_upper if st_upper is not None else "",
                        pchembl if pchembl is not None else "",
                        comment or "", val_comment or "", pot_dup or "",
                        var_label, assay_chembl, atype or "", conf,
                        desc or "", doc_id or "", year or "", doi or ""
                    ])

            total_exported_activities += len(rows)

            med_nM = f"{statistics.median(ep_molar_nM_values):.2f}" if ep_molar_nM_values else ""
            med_pchembl = f"{statistics.median(ep_pchembl_values):.2f}" if ep_pchembl_values else ""

            endpoint_summary[ep] = {
                "num_activities": len(rows),
                "num_compounds": len(ep_compounds),
                "num_approved_drugs": len(ep_approved),
                "num_clinical_drugs": len(ep_clinical),
                "median_value_nM": med_nM,
                "median_pchembl": med_pchembl,
                "exact_relations": ep_exact_count,
                "approx_relations": ep_approx_count,
                "upper_bounds": ep_upper_count,
                "lower_bounds": ep_lower_count
            }

            target_stats_rows.append([
                target_chembl, target_name or "", organism or "", uniprot or "",
                folder_name, ep, len(rows), len(ep_compounds), len(ep_clinical), len(ep_approved),
                ", ".join(sorted(ep_approved)) if ep_approved else "",
                med_pchembl, med_nM,
                ep_exact_count, ep_approx_count, ep_upper_count, ep_lower_count, ep_unknown_count
            ])

        # Write target_metadata.json
        meta_path = os.path.join(target_dir, "target_metadata.json")
        with open(meta_path, "w", encoding="utf-8") as jf:
            json.dump({
                "chembl_id_target": target_chembl,
                "target_pref_name": target_name,
                "target_organism": organism,
                "uniprot_accession": uniprot,
                "protein_classification": prot_class,
                "total_unique_compounds": len(target_compounds_all),
                "total_approved_drugs": len(target_approved_drugs),
                "approved_drugs_list": sorted(target_approved_drugs),
                "total_clinical_candidates": len(target_clinical_drugs),
                "endpoints": endpoint_summary
            }, jf, indent=2)

        targets_index_rows.append([
            target_chembl, target_name or "", organism or "", uniprot or "", prot_class or "",
            len(target_compounds_all), len(target_approved_drugs), len(target_clinical_drugs),
            folder_name
        ])

    conn.close()

    # Write target index
    index_tsv = os.path.join(out_dir, "targets_index.tsv")
    with open(index_tsv, "w", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["chembl_id_target", "target_pref_name", "target_organism", "uniprot_accession",
                         "protein_classification", "num_compounds", "num_approved_drugs",
                         "num_clinical_drugs", "folder_name"])
        writer.writerows(targets_index_rows)

    # Write stats report table
    stats_tsv = os.path.join(REPORTS_DIR, "well_characterized_proteins_drug_stats_v2.tsv")
    with open(stats_tsv, "w", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "chembl_id_target", "target_pref_name", "target_organism", "uniprot_accession",
            "folder_name", "endpoint_type", "num_activities", "num_compounds", "num_clinical_drugs",
            "num_approved_drugs", "approved_drugs_list", "median_pchembl", "median_value_nM",
            "num_exact", "num_approximate", "num_upper_bound", "num_lower_bound", "num_unknown_relation"
        ])
        writer.writerows(target_stats_rows)

    # Generate Markdown Summary dynamically
    summary_md = os.path.join(REPORTS_DIR, "well_characterized_proteins_summary_v2.md")
    with open(summary_md, "w", encoding="utf-8") as mf:
        mf.write("# ChEMBL 37 — Well-Characterized Proteins Profile (Version 2)\n")
        mf.write(f"> **Generated**: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        mf.write(f"> **Qualification Threshold**: $\\ge {min_compounds}$ unique tested compounds under strict single-protein validity predicate\n\n")
        mf.write(f"- **Total Qualified Single-Protein Targets**: **{len(qualified_targets):,}**\n")
        mf.write(f"- **Total Exported Bioactivity Determinations**: **{total_exported_activities:,}**\n")
        mf.write(f"- **Master Target Index**: [`{index_tsv}`]({index_tsv})\n")
        mf.write(f"- **Detailed Drug & Endpoint Table**: [`{stats_tsv}`]({stats_tsv})\n\n")
        mf.write("## Top 10 Data-Rich Targets\n\n")
        mf.write("| ChEMBL Target ID | Preferred Name | Organism | UniProt | Compounds Tested | Approved Drugs |\n")
        mf.write("|:---|:---|:---|:---:|:---:|:---:|\n")
        for r in targets_index_rows[:10]:
            mf.write(f"| `{r[0]}` | {r[1]} | *{r[2]}* | `{r[3]}` | {r[5]:,} | {r[6]} |\n")

    print(f"[{datetime.datetime.now()}] Extraction complete. {len(qualified_targets)} targets exported ({total_exported_activities:,} rows).")

if __name__ == "__main__":
    main()
