# ChEMBL Database Curation

A complete, reproducible pipeline to download, extract, and curate **ChEMBL 37** compound bioactivity and compound-target interaction data.

## 📊 Database Statistics (ChEMBL 37, May 2026)

| Metric | Count |
|--------|-------|
| Compounds | 2,921,148 |
| Activities | 24,527,044 |
| Assays | 1,970,438 |
| Targets | 18,552 |
| Documents | 101,100 |

## 🗂️ Project Structure

```
ChEMBL_Database_Curation/
├── run_pipeline.sh                  # Master pipeline runner
├── requirements.txt
├── scripts/
│   ├── 01_download_chembl.sh        # Download ChEMBL 37 SQLite + schema
│   ├── 02_explore_schema.py         # Schema exploration & reports
│   ├── 03_curate_bioactivity.py     # Lossless broad bioactivity curation
│   ├── 04_curate_target_interactions.py  # Evidence-tiered target interactions
│   ├── 05_generate_report.py        # Master curation report & manifest
│   ├── 06_extract_well_characterized_targets.py # Per-target endpoint extractor
│   └── 07_extract_assay_conditions.py # Assay format, conditions & protocol extractor
├── data/
│   ├── raw/                         # Downloaded ChEMBL files (gitignored)
│   ├── curated/                     # Curated master TSV files
│   │   └── assay_conditions/        # IC50 assay conditions & protocol index
│   ├── targets/                     # 902 Well-characterized protein directories
│   └── schema/                      # Schema documentation & reports
├── reports/                         # Scientific reports, statistics, manifests
└── logs/                            # Pipeline execution logs
```

## 🚀 Quick Start

### 1. Set up environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Run the full pipeline

```bash
bash run_pipeline.sh
```

This will:
1. Download ChEMBL 37 SQLite (~5.4 GB) + schema docs + UniProt mapping
2. Explore the schema and generate schema reports
3. Curate all bioactivity data → `data/curated/chembl_37_bioactivity_broad.tsv.gz`
4. Curate compound-target interactions → `data/curated/chembl_37_compound_target_interactions.tsv.gz`
5. Generate a summary Markdown report

### 3. Resume from a step (if interrupted)

```bash
bash run_pipeline.sh --from-step 3   # Resume from bioactivity curation
bash run_pipeline.sh --skip-download  # Skip download, DB already present
```

## 📁 Output Files

| File | Description |
|------|-------------|
| `data/curated/chembl_37_bioactivity_broad.tsv.gz` | **All bioactivity records** (23.2M) with compound, assay, target, document metadata |
| `data/curated/chembl_37_compound_target_interactions_tiered.tsv.gz` | **Evidence-tiered interactions** (16.9M) stratified by Tier A–E, endpoint, mutation |
| `data/curated/chembl_37_drug_dictionary.tsv` | Clinical & approved drug registry (16.9k) with indications, warnings, ATC, MoA |
| `data/curated/chembl_37_target_components.tsv` | Normalized 1-to-1 UniProt accessions & protein classifications |
| `data/curated/chembl_37_target_variants.tsv` | Target mutations and resistance variant mappings |
| `data/curated/chembl_37_protein_classification.tsv` | Full protein class hierarchy (Levels 1–4) |
| `data/curated/assay_conditions/ic50_assay_conditions.tsv.gz` | **3.44M IC50 determinations** enriched with assay format, [S]/[ATP], incubation time, and detection tech |
| `data/curated/assay_conditions/assay_protocols_conditions_index.tsv` | **293k Assay protocols index** with structured experimental conditions |
| `data/targets/` | **902 Target directories**, each containing separate `Kd.tsv`, `Ki.tsv`, `IC50.tsv`, `EC50.tsv` |
| `data/targets/targets_index.tsv` | Master index of all 902 well-characterized targets with compound & drug counts |
| `reports/assay_conditions_stats.tsv` | Statistical distributions of assay formats, detection methods, and incubation times |
| `reports/assay_conditions_summary.md` | Executive pharmacological summary on experimental IC50 confounders |
| `reports/well_characterized_proteins_drug_stats.tsv` | Detailed target × endpoint statistics with clinical & approved drug representation |
| `reports/well_characterized_proteins_summary.md` | Executive summary and ranking of top therapeutic targets |
| `reports/chembl_37_curation_summary.md` | Full curation report with evidence tiers and endpoint breakdowns |
| `reports/run_manifest.json` | Provenance and reproducibility manifest |

## 🔑 Key Columns (Bioactivity Table)

| Column | Description |
|--------|-------------|
| `chembl_id_compound` | ChEMBL compound ID (e.g. CHEMBL25) |
| `canonical_smiles` | Canonical SMILES structure |
| `compound_inchikey` | Standard InChIKey |
| `standard_type` | Activity type: IC50, Ki, EC50, Kd, MIC, ... |
| `standard_value` | Numeric value |
| `standard_units` | Units (nM, uM, ...) |
| `pchembl_value` | pChEMBL = -log10(IC50 in M), standardized across types |
| `chembl_id_target` | ChEMBL target ID |
| `target_pref_name` | Target preferred name |
| `target_type` | SINGLE PROTEIN / PROTEIN COMPLEX / ORGANISM / ... |
| `uniprot_accession` | UniProt accession(s) |
| `assay_organism` | Organism the assay was run in |
| `confidence_score` | Target assignment confidence (0–9; 9=directly assigned single protein) |
| `data_validity_comment` | Quality flag (filtered: transcription errors, out-of-range) |
| `modality` | e.g. "Targeted protein degradation" (new in ChEMBL 37) |
| `le` | Ligand efficiency |
| `bei` | Binding efficiency index |
| `lle` | Lipophilic ligand efficiency |

## ⚙️ Curation Filters (Broad Mode)

### Excluded
Records where `data_validity_comment` ∈:
- `Potential transcription error`
- `Outside typical range`
- `Non standard unit for type`
- `Author confirmed error`
- `Potential author error`

### Included (broad mode)
- ✅ All activity types (IC50, Ki, EC50, Kd, MIC, GI50, ...)
- ✅ All organism targets (human, rodent, pathogen, plant, ...)
- ✅ Mutant targets
- ✅ All assay types (Binding, Functional, ADME, Toxicity, ...)
- ✅ Qualitative records (activity_comment only, no numeric value)

## 📝 ChEMBL 37 Notable Changes
- New **ACTIVITIES.MODALITY** field for Targeted Protein Degradation annotation
- ~29,000 TPD bioactivity records newly annotated
- 169 compounds upgraded to MAX_PHASE=4 after manual curation
- ~400 tautomer duplicates merged
- Deprecated: `CURATION_LOOKUP`, `PROTEIN_CLASS_SYNONYMS` tables

## 📜 License

ChEMBL data is provided under the [Creative Commons Attribution-ShareAlike 3.0 Unported License](https://creativecommons.org/licenses/by-sa/3.0/).

Please cite: *Zdrazil B, et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity complementary capabilities. Nucleic Acids Research 51, D1373-D1380 (2023).*

## 🔗 Resources & Documentation
- [Comprehensive Scientific & Computational Methodology](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/METHODOLOGY.md) (`METHODOLOGY.md`)
- [ChEMBL FTP](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)
- [ChEMBL Interface Documentation](https://chembl.gitbook.io/chembl-interface-documentation/)
- [ChEMBL Schema Docs](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.html)
