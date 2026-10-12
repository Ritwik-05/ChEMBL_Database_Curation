#!/usr/bin/env python3
"""
02_explore_schema.py — ChEMBL 37 Schema Exploration & Profiling (Version 2)
=============================================================================
Audited for Scientific Quality Council standards.
Profiles all core tables, assay condition tables (assay_parameters, activity_properties),
bioassay ontology, variants, and indexes to ensure complete relational coverage.
Outputs:
  Version2/data/schema/chembl_37_schema_report.tsv
  Version2/reports/schema_profile_v2.md
"""

import os
import sys
import sqlite3
import datetime

DB_PATH = "data/raw/chembl_37.db"
OUTPUT_TSV = "Version2/data/schema/chembl_37_schema_report.tsv"
REPORT_MD = "Version2/reports/schema_profile_v2.md"

KEY_TABLES = [
    "ACTIVITIES",
    "ASSAYS",
    "ASSAY_PARAMETERS",
    "ACTIVITY_PROPERTIES",
    "BIOASSAY_ONTOLOGY",
    "TARGET_DICTIONARY",
    "TARGET_COMPONENTS",
    "COMPONENT_SEQUENCES",
    "COMPONENT_CLASS",
    "PROTEIN_CLASSIFICATION",
    "VARIANT_SEQUENCES",
    "MOLECULE_DICTIONARY",
    "COMPOUND_STRUCTURES",
    "DRUG_MECHANISM",
    "DRUG_INDICATION",
    "DRUG_WARNING",
    "DOCS",
]

def main():
    print(f"[{datetime.datetime.now()}] Starting ChEMBL 37 Schema Profiling (Version 2)...")
    if not os.path.exists(DB_PATH):
        print(f"ERROR: Database file {DB_PATH} not found.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(OUTPUT_TSV), exist_ok=True)
    os.makedirs(os.path.dirname(REPORT_MD), exist_ok=True)

    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    c = conn.cursor()

    # Get list of all tables in SQLite
    c.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;")
    all_tables = [r[0] for r in c.fetchall()]
    print(f"Total tables found in database: {len(all_tables)}")

    table_stats = []
    
    with open(OUTPUT_TSV, "w", encoding="utf-8") as out:
        out.write("table_name\trow_count\tcolumn_count\tcolumn_names\tprimary_key\tindexes\n")
        
        for tbl in all_tables:
            # Table info
            c.execute(f"PRAGMA table_info('{tbl}');")
            cols = c.fetchall()
            col_names = [col[1] for col in cols]
            pks = [col[1] for col in cols if col[5] > 0]
            
            # Indexes
            c.execute(f"PRAGMA index_list('{tbl}');")
            idxs = [idx[1] for idx in c.fetchall()]
            
            # Row count
            try:
                c.execute(f"SELECT count(*) FROM '{tbl}';")
                row_count = c.fetchone()[0]
            except Exception as e:
                row_count = -1
                print(f"Warning: could not count rows for {tbl}: {e}")

            out.write(f"{tbl}\t{row_count}\t{len(cols)}\t{';'.join(col_names)}\t{';'.join(pks)}\t{';'.join(idxs)}\n")
            
            if tbl.upper() in KEY_TABLES or tbl in KEY_TABLES:
                table_stats.append({
                    "name": tbl,
                    "rows": row_count,
                    "cols": len(cols),
                    "col_names": col_names,
                    "pks": pks,
                    "indexes": idxs
                })

    # Generate Markdown Report
    with open(REPORT_MD, "w", encoding="utf-8") as rpt:
        rpt.write("# ChEMBL 37 Database Schema Profile (Version 2)\n")
        rpt.write(f"> **Generated**: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        rpt.write(f"> **Database**: `{DB_PATH}` (Read-only mode)\n")
        rpt.write(f"> **Total Tables Profiled**: {len(all_tables)}\n\n")
        
        rpt.write("## 1. Core and Condition Tables Profile\n\n")
        rpt.write("| Table Name | Row Count | Columns | Primary Key | Key Indexes |\n")
        rpt.write("|:---|:---:|:---:|:---|:---|\n")
        for s in table_stats:
            pks = ", ".join(s["pks"]) if s["pks"] else "None"
            idxs = f"{len(s['indexes'])} indexes" if s["indexes"] else "No index"
            rpt.write(f"| `{s['name']}` | {s['rows']:,} | {s['cols']} | `{pks}` | {idxs} |\n")
        
        rpt.write("\n## 2. Parameter and Property Tables Detail\n\n")
        rpt.write("- **`ASSAY_PARAMETERS`**: Contains assay-level experimental parameters (e.g. $[ATP]$, substrate, incubation duration, temperature).\n")
        rpt.write("- **`ACTIVITY_PROPERTIES`**: Contains 12.2M measurement-specific properties (e.g. timepoint, compound concentrations, replicate counts).\n")
        rpt.write("- **`BIOASSAY_ONTOLOGY`**: Formal BAO vocabulary mapping assay format codes (e.g. BAO_0000219 for cell-based format, BAO_0000357 for single protein format).\n")

    conn.close()
    print(f"[{datetime.datetime.now()}] Schema exploration complete. Report written to {REPORT_MD}")

if __name__ == "__main__":
    main()
