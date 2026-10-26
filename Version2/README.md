# ChEMBL Database Curation — Version 2 (Audited & Re-engineered)

A fully audited, scientifically validated curation release of **ChEMBL 37** addressing all critical, major, and engineering findings raised in the Reviewer 2 Audit.

---

## 🌟 What is New in Version 2?

Version 2 was developed in a dedicated, isolated directory (`Version2/`) to preserve original Version 1 files while completely re-engineering the curation pipeline:

1. **Exact Target Qualification (Resolves C01)**:
   - Unified eligible SQL predicate applied across both target qualification and data export.
   - Identified exactly **880 true single-protein targets** ($\ge 500$ unique compounds evaluated against high-confidence $K_d, K_i, IC_{50}, EC_{50}$).
   - Low-confidence assays ($\text{confidence} < 8$) and flagged transcription/author errors are strictly barred from target folders.
2. **Censoring-Aware Measurement Relations (Resolves C02)**:
   - Preserves original relations (`=`, `~`, `<`, `<=`, `>`, `>=`).
   - Missing relations are labeled `unspecified`—never fabricated as exact `=`.
   - Activity status (Active, Inactive, Inconclusive) is decoupled from relation syntax.
3. **Molar Unit-Normalized Aggregation (Resolves C03)**:
   - All concentration metrics are normalized to a standard molar scale ($\text{nM}$) prior to calculating arithmetic medians.
   - Separates incompatible dimensions; maintains independent $p\text{Scale}$ summaries.
4. **Multi-Stage & Signed Experimental Extraction (Resolves C05, C06, C07, C09, C10)**:
   - Supports scientific exponential notation (`1e-3 mM ATP` $\rightarrow 1\ \mu\text{M}$ ATP).
   - Multi-stage durations (pre-incubation vs. reaction) parsed separately; seconds and days converted to normalized minutes before bucketing.
   - Optical negation handled ("no fluorescence; measured by LC-MS" maps to Mass Spectrometry).
   - Negative temperature preservation (`-20 °C` storage vs. assay reaction temperature).
   - Long-form parameter preservation without dictionary overwriting.
5. **Float-Safe Clinical Phase Parsing (Resolves H04)**:
   - Correctly parses `4.0`, `1.0`, `0.5`, `-1`, recovering over 2,300 previously dropped clinical candidates.
6. **Rejection Ledger & True Cryptographic Manifest (Resolves H01, H08, H09)**:
   - Generates an explicit `rejection_ledger.tsv` reconciling all 24.5M raw activities.
   - `run_manifest_v2.json` computes genuine SHA-256 checksums of all scripts and outputs.

---

## 🗂️ Version 2 Directory Structure

```
Version2/
├── README.md                                       # This file
├── METHODOLOGY_v2.md                               # Comprehensive scientific methodology
├── run_pipeline_v2.sh                              # Complete master execution script
├── scripts/                                        # Corrected execution scripts
│   ├── 02_explore_schema.py                        # Enhanced schema & parameter profiling
│   ├── 03_curate_bioactivity.py                    # Curation with explicit rejection accounting
│   ├── 04_curate_target_interactions.py           # Orthogonal evidence stratification
│   ├── 05_generate_report.py                       # Dynamic report & cryptographic manifest
│   ├── 06_extract_well_characterized_targets.py    # Strict single-protein target extractor
│   └── 07_extract_assay_conditions.py              # Multi-stage, negation-aware condition extractor
├── data/
│   ├── schema/
│   │   └── chembl_37_schema_report.tsv             # Column, index, and table profiles
│   ├── curated/
│   │   ├── chembl_37_bioactivity_curated.tsv.gz    # Clean curated activities with full provenance
│   │   ├── chembl_37_compound_target_associations_v2.tsv.gz # Stratified associations
│   │   ├── chembl_37_drug_dictionary.tsv           # Normalized drug & candidate dictionary
│   │   ├── chembl_37_target_components.tsv         # Normalized 1:1 UniProt mappings
│   │   ├── chembl_37_protein_classification.tsv    # Protein family hierarchy
│   │   ├── chembl_37_target_variants.tsv           # Target variant registry
│   │   └── assay_conditions/
│   │       ├── ic50_assay_conditions_v2.tsv.gz     # 3.44M IC50 records with validated conditions
│   │       └── assay_protocols_conditions_index_v2.tsv # Protocol condition index
│   └── targets/
│       ├── targets_index.tsv                       # Master index of 880 qualified targets
│       └── <CHEMBL_ID>_<Slug>/                     # 880 Target-specific folders
│           ├── Kd.tsv
│           ├── Ki.tsv
│           ├── IC50.tsv
│           ├── EC50.tsv
│           └── target_metadata.json
├── reports/
│   ├── changes_and_audit_report.md                 # Itemized disposition of Reviewer 2 audit
│   ├── bioactivity_rejection_ledger.tsv            # Raw activity attrition accounting table
│   ├── well_characterized_proteins_drug_stats_v2.tsv # Target × Endpoint drug statistics
│   ├── well_characterized_proteins_summary_v2.md   # Dynamic target profile report
│   ├── assay_conditions_stats_v2.tsv               # Dynamic assay condition distributions
│   ├── assay_conditions_summary_v2.md              # Dynamic assay condition report
│   ├── chembl_37_curation_summary_v2.md            # Master curation summary
│   ├── chembl_37_curation_summary_v2.tsv           # High-level pipeline metrics
│   └── run_manifest_v2.json                        # Cryptographic SHA-256 reproducibility manifest
└── logs/                                           # Pipeline execution logs
```

---

## 🚀 Quick Execution Guide

```bash
# Run the complete audited Version 2 pipeline
bash Version2/run_pipeline_v2.sh
```
