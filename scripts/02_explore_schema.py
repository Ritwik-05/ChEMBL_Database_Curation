#!/usr/bin/env python3
"""
Script: 02_explore_schema.py
Description: Connects to the ChEMBL 37 SQLite database, explores the schema,
             generates a comprehensive schema report with table sizes,
             column info, and key relationships relevant to bioactivity curation.
             Outputs: data/schema/chembl_37_schema_report.tsv
                      data/schema/chembl_37_table_sizes.tsv
Usage: python scripts/02_explore_schema.py
"""

import sqlite3
import os
import sys
import csv
import json
import logging
from datetime import datetime
from pathlib import Path

# --- Configuration -----------------------------------------------------------
CHEMBL_DB   = "data/raw/chembl_37.db"
SCHEMA_DIR  = Path("data/schema")
LOG_DIR     = Path("logs")
LOGFILE     = LOG_DIR / f"schema_explore_{datetime.now():%Y%m%d_%H%M%S}.log"

# Tables of primary interest for bioactivity + target interaction curation
KEY_TABLES = [
    "activities",
    "assays",
    "target_dictionary",
    "target_components",
    "component_sequences",
    "molecule_dictionary",
    "compound_structures",
    "compound_properties",
    "compound_records",
    "docs",
    "protein_classification",
    "protein_class_synonyms",
    "target_relations",
    "organism_class",
    "cell_dictionary",
    "tissue_dictionary",
    "action_type",
    "activity_properties",
    "activity_stds_lookup",
    "binding_sites",
    "site_components",
    "mechanism_refs",
    "drug_mechanism",
    "predicted_binding_domains",
    "ligand_eff",
    "metabolism",
]

# --- Logging -----------------------------------------------------------------
LOG_DIR.mkdir(parents=True, exist_ok=True)
SCHEMA_DIR.mkdir(parents=True, exist_ok=True)

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
        log.error("Run scripts/01_download_chembl.sh first.")
        sys.exit(1)
    log.info(f"Connecting to: {db_path}")
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def get_all_tables(conn: sqlite3.Connection) -> list[str]:
    cur = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
    )
    return [row[0] for row in cur.fetchall()]


def get_table_info(conn: sqlite3.Connection, table: str) -> list[dict]:
    cur = conn.execute(f"PRAGMA table_info('{table}')")
    return [dict(row) for row in cur.fetchall()]


def get_row_count(conn: sqlite3.Connection, table: str) -> int:
    try:
        cur = conn.execute(f"SELECT COUNT(*) FROM '{table}'")
        return cur.fetchone()[0]
    except Exception:
        return -1


def get_foreign_keys(conn: sqlite3.Connection, table: str) -> list[dict]:
    cur = conn.execute(f"PRAGMA foreign_key_list('{table}')")
    return [dict(row) for row in cur.fetchall()]


def sample_values(conn: sqlite3.Connection, table: str, column: str, n: int = 5) -> list:
    try:
        cur = conn.execute(
            f"SELECT DISTINCT \"{column}\" FROM \"{table}\" "
            f"WHERE \"{column}\" IS NOT NULL LIMIT {n}"
        )
        return [str(row[0]) for row in cur.fetchall()]
    except Exception:
        return []


def main():
    log.info("=" * 60)
    log.info(" ChEMBL 37 Schema Explorer")
    log.info("=" * 60)

    conn = connect(CHEMBL_DB)

    # 1. Get all tables
    all_tables = get_all_tables(conn)
    log.info(f"Total tables in DB: {len(all_tables)}")

    # 2. Table sizes
    log.info("\nCollecting table row counts...")
    table_sizes = []
    for t in all_tables:
        count = get_row_count(conn, t)
        is_key = t in KEY_TABLES
        table_sizes.append({"table": t, "row_count": count, "is_key_table": is_key})
        if count > 0:
            marker = " <-- KEY" if is_key else ""
            log.info(f"  {t:<45} {count:>12,} rows{marker}")

    sizes_file = SCHEMA_DIR / "chembl_37_table_sizes.tsv"
    with open(sizes_file, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["table", "row_count", "is_key_table"], delimiter="\t")
        w.writeheader()
        w.writerows(sorted(table_sizes, key=lambda x: x["row_count"], reverse=True))
    log.info(f"\nTable sizes saved to: {sizes_file}")

    # 3. Detailed column schema for key tables
    log.info("\nExtracting detailed schema for key tables...")
    schema_rows = []
    fk_rows = []

    for table in KEY_TABLES:
        if table not in all_tables:
            log.warning(f"  Table not found: {table}")
            continue

        cols = get_table_info(conn, table)
        fks  = get_foreign_keys(conn, table)
        fk_map = {fk["from"]: f"{fk['table']}.{fk['to']}" for fk in fks}

        for col in cols:
            samples = sample_values(conn, table, col["name"])
            schema_rows.append({
                "table":        table,
                "cid":          col["cid"],
                "column":       col["name"],
                "type":         col["type"],
                "notnull":      col["notnull"],
                "pk":           col["pk"],
                "fk_reference": fk_map.get(col["name"], ""),
                "sample_values": " | ".join(samples),
            })

        for fk in fks:
            fk_rows.append({
                "from_table":  table,
                "from_column": fk["from"],
                "to_table":    fk["table"],
                "to_column":   fk["to"],
            })

    schema_file = SCHEMA_DIR / "chembl_37_schema_report.tsv"
    with open(schema_file, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["table","cid","column","type","notnull","pk","fk_reference","sample_values"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(schema_rows)
    log.info(f"Schema report saved to: {schema_file}")

    fk_file = SCHEMA_DIR / "chembl_37_foreign_keys.tsv"
    with open(fk_file, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["from_table","from_column","to_table","to_column"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(fk_rows)
    log.info(f"Foreign keys saved to: {fk_file}")

    # 4. Activity type summary
    log.info("\nSummarizing activity types (standard_type)...")
    act_types = conn.execute("""
        SELECT standard_type, COUNT(*) as cnt
        FROM activities
        WHERE standard_type IS NOT NULL
        GROUP BY standard_type
        ORDER BY cnt DESC
        LIMIT 50
    """).fetchall()

    act_file = SCHEMA_DIR / "chembl_37_activity_type_summary.tsv"
    with open(act_file, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["standard_type", "count"])
        for row in act_types:
            w.writerow([row[0], row[1]])
    log.info(f"Top 50 activity types saved to: {act_file}")

    # 5. Target type summary
    log.info("\nSummarizing target types...")
    tgt_types = conn.execute("""
        SELECT target_type, COUNT(*) as cnt
        FROM target_dictionary
        GROUP BY target_type
        ORDER BY cnt DESC
    """).fetchall()

    tgt_file = SCHEMA_DIR / "chembl_37_target_type_summary.tsv"
    with open(tgt_file, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["target_type", "count"])
        for row in tgt_types:
            w.writerow([row[0], row[1]])
    log.info(f"Target types saved to: {tgt_file}")

    # 6. Organism summary
    log.info("\nSummarizing organisms in assays...")
    orgs = conn.execute("""
        SELECT assay_organism, COUNT(*) as cnt
        FROM assays
        WHERE assay_organism IS NOT NULL
        GROUP BY assay_organism
        ORDER BY cnt DESC
        LIMIT 30
    """).fetchall()

    org_file = SCHEMA_DIR / "chembl_37_assay_organism_summary.tsv"
    with open(org_file, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["assay_organism", "assay_count"])
        for row in orgs:
            w.writerow([row[0], row[1]])
    log.info(f"Top 30 organisms saved to: {org_file}")

    conn.close()

    log.info("")
    log.info("=" * 60)
    log.info(" Schema exploration complete!")
    log.info(f"  Reports: {SCHEMA_DIR}/")
    log.info("=" * 60)
    log.info("Next step: python scripts/03_curate_bioactivity.py")


if __name__ == "__main__":
    main()
