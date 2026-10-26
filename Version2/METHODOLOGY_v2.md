# Scientific & Computational Methodology (Version 2)
## Audited ChEMBL 37 Curation, Target-Specific Partitioning, and Assay Measurement Conditions

> **Version**: 2.0 (Audited & Re-engineered)  
> **Release Date**: October 2026  
> **Host Database**: ChEMBL 37 (EMBL-EBI SQLite, 29.0 GB raw)  
> **Repository Directory**: `Version2/`  


---

## 1. Executive Summary & Purpose of Version 2 Re-engineering

Version 1 established the initial database extraction architecture. However, an independent audit revealed critical scientific and numerical failure modes that compromised the reliability of the derived datasets for machine learning. 

Version 2 re-engineers the entire pipeline to address all reviewer critiques without overwriting Version 1 assets:
- **Unified Target Qualification & Export**: Eliminates the mismatch where targets qualified using invalid records and exported low-confidence assays. A single, strict predicate is enforced across both stages, identifying exactly **880 data-rich single-protein targets** ($\ge 500$ unique compounds).
- **Censoring & Relation Preservation**: Eliminates the fabrication of exact measurements (`=`). Distinguishes exact (`=`), approximate (`~`), upper bound (`<`, `<=`), lower bound (`>`, `>=`), and unknown relations.
- **Unit-Normalized Aggregation**: Concentrations are normalized to a consistent molar scale ($nM$) prior to arithmetic median calculations; incompatible units are never pooled; $p\text{Scale}$ summaries are maintained independently.
- **Multi-Stage & Signed Assay Condition Extraction**:
  - Handles scientific exponential notation (e.g., `1e-3 mM ATP` correctly parsed as $1\ \mu\text{M}$ rather than $3\ \text{mM}$).
  - Multi-stage durations (preincubation vs. reaction) are parsed separately, with seconds and days converted to normalized minutes before bucketing.
  - Optical negation (e.g., "no fluorescence detection; measured by LC-MS") correctly maps to mass spectrometry rather than false fluorescence intensity.
  - Signed temperature regex preserves negative storage conditions (e.g., `-20 °C`), distinguishing storage from reaction temperatures.
- **Relational Preservation of Parameters**: Structured parameters from `assay_parameters` are stored in long-form lists rather than overwritten dictionaries.
- **Cryptographic Provenance**: Every generated dataset, script, and environment parameter is verified using cryptographic SHA-256 digests in a comprehensive manifest (`run_manifest_v2.json`).

---

## 2. Itemized Resolution of  Findings

| Finding ID |  Critique | Version 1 Defect | Version 2 Resolution |
|:---|:---|:---|:---|
| **C01** | Target qualification vs. export eligibility mismatch | Qualification counted invalid records; export query omitted confidence filter. | Unified eligible-record SQL predicate applied across both qualification and export: `td.target_type = 'SINGLE PROTEIN' AND a.confidence_score >= 8 AND act.standard_value IS NOT NULL AND (validity is NULL or valid)`. 880 targets qualified. |
| **C02** | Relation syntax & inactivity handling | Missing relations coerced to `=`; `~` treated as exact; relations inferred as inactive. | Preserves original relations. Missing relations exported as `unspecified`. Activity status (Active, Inactive, Inconclusive) evaluated independently from relation syntax. |
| **C03** | Mixed-unit median aggregation | Raw standard values pooled and labeled with the first encountered unit. | Concentrations converted to common molar scale ($\text{nM}$) before computing medians; incompatible non-molar dimensions separated. |
| **C04** | Evidence Tier A overstatement | Tier A treated confidence 8 as direct binding and included inactive edges. | Deconstructed into 5 independent orthogonal axes: Assignment Confidence (9 direct vs 8 homologous), Format, Purpose, Relation, and Activity Status. |
| **C05** | Assay format over-inference & BAO errors | Cell lysate mapped to intact cells; BAO_0000221 misdocumented as biochemical binding. | Format mapped strictly from database ontology labels; cell lysate classified as `Biochemical (Cell-Free Lysate)`; BAO documentation corrected. |
| **C06** | Incubation duration & stage errors | 120 sec bucketed into 1–4 hours; days unsupported; multi-stage collapsed. | Numeric conversion to minutes for all units (sec, min, hr, day); pre-incubation and reaction stages parsed separately. |
| **C07** | Substrate/ATP concentration errors | `1e-3 mM ATP` parsed as `3 mM ATP`; ATP entity lost; $K_m$ mixed with substrate. | Scientific notation regex implemented; ATP identity preserved; $K_m$ extracted into a dedicated separate field. |
| **C08** | Structured parameter overwriting | Dict assignment overwrote duplicate parameter types; activity properties ignored. | Parameters stored in long-form list per assay; structured summaries record all parameters. |
| **C09** | Detection technology negation & lexicon | "No fluorescence; LC-MS" mapped to fluorescence; ADP-Glo missing. | Negation patterns handled; expanded lexicon including ADP-Glo, TR-FRET, FP, AlphaLISA, SPR, LC-MS. |
| **C10** | Temperature sign loss & storage confounding | `-20 °C` became `20 °C`; storage confused with assay temperature. | Signed regex matches negative temperatures; storage temperatures separated from reaction temperatures. |
| **H01** | "Lossless" claim inaccurate | Quality exclusions made without an attrition ledger. | Explicit `rejection_ledger.tsv` generated with exact row accounting for all raw activities. |
| **H03** | Critical modeling fields dropped | Upper bounds, validity comments, and duplicate flags stripped. | `standard_upper_value`, `data_validity_comment`, `potential_duplicate`, `doc_id`, `doi` retained in all exports. |
| **H04** | `int(max_phase)` silently dropped valid forms | Floats like `4.0`, `1.0`, `0.5` caused `ValueError` and were discarded. | Float-safe parser (`parse_phase`) correctly processes numeric strings without truncation. |
| **H08** | Hardcoded numbers in reports | Target counts (902) and percentages hardcoded in Markdown text. | All report metrics dynamically calculated from current execution outputs. |
| **H09** | Manifest lacked real SHA-256 hashes | Manifest wrote placeholders; succeeded even when input files were missing. | Genuine SHA-256 hashes computed; script fails with exit code 1 if any mandatory file is absent. |
| **H10** | SQL query lacked ORDER BY contract | Relied on index scan without explicit SQL contract. | Strict `ORDER BY act.molregno, act.activity_id` enforced in all queries. |

---

## 3. Directory Layout in Version 2

All Version 2 outputs are segregated into `Version2/`:

```
Version2/
├── run_pipeline_v2.sh                              # Master pipeline orchestrator
├── METHODOLOGY_v2.md                               # This document
├── scripts/
│   ├── 02_explore_schema.py                        # Enhanced schema & parameter profiling
│   ├── 03_curate_bioactivity.py                    # Curation with explicit rejection ledger
│   ├── 04_curate_target_interactions.py           # Orthogonal evidence stratification
│   ├── 05_generate_report.py                       # Manifest & dynamic SHA-256 reporting
│   ├── 06_extract_well_characterized_targets.py    # Strict single-protein target extractor
│   └── 07_extract_assay_conditions.py              # Multi-stage, negation-aware condition extractor
├── data/
│   ├── schema/
│   │   └── chembl_37_schema_report.tsv             # Detailed column & index inventory
│   ├── curated/
│   │   ├── chembl_37_bioactivity_curated.tsv.gz    # Clean bioactivity dataset with full provenance
│   │   ├── chembl_37_compound_target_associations_v2.tsv.gz # Stratified associations
│   │   ├── chembl_37_drug_dictionary.tsv           # Normalized drug & candidate registry
│   │   ├── chembl_37_target_components.tsv         # Normalized 1:1 UniProt accessions
│   │   ├── chembl_37_protein_classification.tsv    # Protein class hierarchy
│   │   ├── chembl_37_target_variants.tsv           # Target variant registry
│   │   └── assay_conditions/
│   │       ├── ic50_assay_conditions_v2.tsv.gz     # 3.44M IC50 records with validated conditions
│   │       └── assay_protocols_conditions_index_v2.tsv # Protocol condition index
│   └── targets/
│       ├── targets_index.tsv                       # Master index of 880 qualified targets
│       └── <CHEMBL_ID>_<Slug>/                     # 880 Target-specific directories
│           ├── Kd.tsv
│           ├── Ki.tsv
│           ├── IC50.tsv
│           ├── EC50.tsv
│           └── target_metadata.json
├── reports/
│   ├── bioactivity_rejection_ledger.tsv            # Raw activity attrition accounting table
│   ├── well_characterized_proteins_drug_stats_v2.tsv # Target × Endpoint drug statistics
│   ├── well_characterized_proteins_summary_v2.md   # Dynamic target profile report
│   ├── assay_conditions_stats_v2.tsv               # Dynamic assay condition distributions
│   ├── assay_conditions_summary_v2.md              # Dynamic assay condition report
│   ├── chembl_37_curation_summary_v2.md            # Master curation summary
│   ├── chembl_37_curation_summary_v2.tsv           # High-level pipeline metrics
│   └── run_manifest_v2.json                        # Cryptographic SHA-256 reproducibility manifest
└── logs/                                           # Execution logs
```

---

## 4. Machine Learning Guidance & Confounder Handling

1. **Assay Format Conditioning**:
   - Biochemical assays ($\sim 10\ \mu\text{M}$ ATP) and cell-based assays ($1\text{--}5\text{ mM}$ intracellular ATP) must never be pooled into an unconditioned target label.
   - For structural affinity prediction (e.g., docking, 3D binding models), filter by `assay_format_category == 'Biochemical (Single Protein)'` and `assignment_confidence_axis == 'DIRECT_SINGLE_PROTEIN (Conf 9)'`.
   - For phenotypic modeling, condition models on the biological format token.
2. **Censored & Approximate Loss Functions**:
   - Use survival or interval-censored regression (e.g., Tobit loss or thresholded margin ranking) for upper bound (`<`, `<=`) and lower bound (`>`, `>=`) measurements rather than treating them as exact numerical targets.
3. **Screening Out Optical Artifacts**:
   - When evaluating compounds active exclusively in `Fluorescence Intensity` assays, check for conjugated fluorophore substructures to prevent learning assay interference artifacts.

---

## 5. Execution Guide

To execute the complete Version 2 pipeline from the workspace root:

```bash
bash Version2/run_pipeline_v2.sh
```

Individual scripts can also be executed independently:
```bash
python3 Version2/scripts/02_explore_schema.py
python3 Version2/scripts/03_curate_bioactivity.py
python3 Version2/scripts/04_curate_target_interactions.py
python3 Version2/scripts/06_extract_well_characterized_targets.py
python3 Version2/scripts/07_extract_assay_conditions.py
python3 Version2/scripts/05_generate_report.py
```
