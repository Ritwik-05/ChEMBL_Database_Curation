# Scientific and Computational Methodology
## ChEMBL 37 Bioactivity Curation, Target-Specific Partitioning, and Assay Conditions Extraction

> **Version**: 1.0  
> **Release Date**: October 2026  
> **Host Database**: ChEMBL 37 (EMBL-EBI, SQLite Distribution, 29.0 GB raw)  
> **Repository**: `ChEMBL_Database_Curation`  
> **Compliance**: Scientific Quality Council Standards (`AGENTS.md`)

---

## 1. Executive Summary & Scientific Problem Statement

Quantitative bioactivity datasets in public repositories are frequently consumed by machine learning (ML) practitioners as tabular matrices of compound SMILES, protein IDs, and numerical activity values (such as $IC_{50}$ or $p\text{ChEMBL}$). However, treating raw bioactivity records as homogeneous physical binding measurements introduces profound systematic biases and false scientific inferences:

1. **Assay Format Confounding**: An $IC_{50}$ determined in an intact cell assay (e.g., cell viability or reporter gene assay) reflects cellular permeability, active efflux, and competition with millimolar intracellular ATP ($\sim 1\text{--}5\text{ mM}$). In contrast, a cell-free biochemical kinase assay typically uses $10\ \mu\text{M}$ ATP. The identical compound can exhibit an apparent $IC_{50}$ that differs by over 300-fold solely due to assay format.
2. **Target Assignment Ambiguity**: In screening databases, an assay associated with a target protein does not always represent a direct physical binding event against that purified protein. Assays may measure phenotypic pathway readouts, whole organism lethality, or responses in multi-protein complexes.
3. **Endpoint Incompatibility**: Thermodynamically rigorous equilibrium constants ($K_d, K_i$) must not be pooled directly with protocol-dependent functional metrics ($IC_{50}, EC_{50}$) or phenomenological thresholds ($MIC, GI_{50}$).
4. **Optical & Biophysical Interference**: Standard fluorescence intensity assays are susceptible to compound autofluorescence and optical quenching, causing non-specific false-positive readouts.

This curation framework provides a scientifically rigorous, fully traceable, and reproducible pipeline that:
- Curates **all 24.5 million activities** in ChEMBL 37 while purging transcription and reporting errors.
- Segregates data into **5 explicit evidence tiers** (Tiers A through E) separating direct single-protein binding from cellular, complex, and organism assays.
- Extracts and partitions data for **902 well-characterized single-protein targets** ($\ge 500$ tested compounds per target) into dedicated directories with separated endpoint tables ($K_d, K_i, IC_{50}, EC_{50}$).
- Parses and annotates **all 3.44 million $IC_{50}$ determinations** with experimental assay inputs: biological format, substrate/ATP concentration, incubation duration, and detection technology.
- Segregates datasets, logs, reports, and statistical tables into distinct, dedicated directories.

---

## 2. Architecture & File Structure

The workspace strictly separates executable code, raw data, curated datasets, statistical tables, markdown reports, and execution logs:

```
ChEMBL_Database_Curation/
├── README.md                                       # High-level overview & quick-start guide
├── METHODOLOGY.md                                  # Comprehensive scientific methodology (this document)
├── requirements.txt                                # Python dependencies (tqdm, pandas, etc.)
├── run_pipeline.sh                                 # Master pipeline execution script
├── scripts/                                        # Modular execution scripts
│   ├── 01_download_chembl.sh                       # Step 1: Download & verify ChEMBL 37 SQLite
│   ├── 02_explore_schema.py                        # Step 2: Database schema exploration & table mapping
│   ├── 03_curate_bioactivity.py                    # Step 3: Lossless broad bioactivity curation
│   ├── 04_curate_target_interactions.py           # Step 4: Evidence-tiered compound-target interaction curation
│   ├── 05_generate_report.py                       # Step 5: Master curation report & run manifest
│   ├── 06_extract_well_characterized_targets.py    # Step 6: Per-target folder & endpoint extraction
│   └── 07_extract_assay_conditions.py              # Step 7: IC50 experimental condition extraction
├── data/
│   ├── raw/                                        # Raw ChEMBL 37 database files (gitignored)
│   │   ├── chembl_37.db                            # SQLite database (29.0 GB)
│   │   ├── chembl_37_sqlite.tar.gz                 # Compressed archive
│   │   └── checksums.txt                           # Official MD5 checksums
│   ├── schema/                                     # Database schema definitions & column profiles
│   │   ├── chembl_37_schema_report.tsv             # Table row counts, column types, foreign keys
│   │   └── schema_documentation.html               # Official ChEMBL HTML schema guide
│   ├── curated/                                    # Master curated dataset files
│   │   ├── chembl_37_bioactivity_broad.tsv.gz      # Complete curated bioactivity table (23.2M rows)
│   │   ├── chembl_37_compound_target_interactions_tiered.tsv.gz # Tiered interactions (16.9M rows)
│   │   ├── chembl_37_drug_dictionary.tsv           # Approved & clinical drug registry (16.9k compounds)
│   │   ├── chembl_37_target_components.tsv         # Normalized 1:1 UniProt accessions
│   │   ├── chembl_37_target_variants.tsv           # Target mutations & resistance variants
│   │   ├── chembl_37_protein_classification.tsv    # Protein family hierarchy (Levels 1–4)
│   │   └── assay_conditions/                       # IC50 experimental conditions dataset
│   │       ├── ic50_assay_conditions.tsv.gz        # 3.44M IC50 records with experimental conditions
│   │       └── assay_protocols_conditions_index.tsv# 293k assay protocol profiles
│   └── targets/                                    # Target-specific data store
│       ├── targets_index.tsv                       # Master index of 902 well-characterized targets
│       ├── README.md                               # Target folder user guide
│       └── <CHEMBL_ID>_<Target_Slug>/              # 902 Individual protein target directories
│           ├── Kd.tsv                              # Dissociation constants
│           ├── Ki.tsv                              # Inhibition constants
│           ├── IC50.tsv                            # Half-maximal inhibitory concentrations
│           ├── EC50.tsv                            # Half-maximal effective concentrations
│           └── target_metadata.json                # Target profile & drug statistics
├── reports/                                        # Scientific audit reports & summary tables
│   ├── chembl_37_curation_summary.md               # Master pipeline curation report
│   ├── chembl_37_curation_summary.tsv              # High-level pipeline summary table
│   ├── well_characterized_proteins_summary.md      # Well-characterized target executive report
│   ├── well_characterized_proteins_drug_stats.tsv  # Target × Endpoint drug statistics table
│   ├── assay_conditions_summary.md                 # IC50 experimental condition report
│   ├── assay_conditions_stats.tsv                  # Detailed format, detection, and incubation stats
│   └── run_manifest.json                           # Pipeline run manifest & SHA256 checksums
└── logs/                                           # Verifiable execution logs
    ├── download_*.log                              # Download & checksum logs
    ├── schema_explore_*.log                        # Schema exploration logs
    ├── curate_bioactivity_*.log                    # Broad bioactivity extraction logs
    ├── curate_targets_*.log                        # Tiered target interaction logs
    ├── report_*.log                                # Report generation logs
    ├── extract_well_characterized_targets_*.log    # Target folder extraction logs
    └── extract_assay_conditions_*.log              # Assay condition extraction logs
```

---

## 3. Detailed Step-by-Step Methodology

### Step 1: Database Ingestion and Checksum Verification
- **Script**: [`scripts/01_download_chembl.sh`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/scripts/01_download_chembl.sh)
- **Data Source**: Official EMBL-EBI FTP server (`https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/`)
- **Procedure**:
  1. Downloads `chembl_37_sqlite.tar.gz`, `schema_documentation.html`, and `chembl_uniprot_mapping.txt`.
  2. Verifies the downloaded file against official MD5 checksums:
     ```bash
     md5sum -c checksums.txt
     ```
  3. Extracts `chembl_37.db` (29.0 GB uncompressed SQLite database) into `data/raw/`.
  4. Opens SQLite in read-only immutable mode (`file:chembl_37.db?mode=ro`) to prevent accidental modifications to raw data.
- **Log**: [`logs/download_20261009_160217.log`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/logs/download_20261009_160217.log)

---

### Step 2: Relational Schema & Ontology Exploration
- **Script**: [`scripts/02_explore_schema.py`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/scripts/02_explore_schema.py)
- **Objective**: Establish table schemas, relationships, foreign keys, row counts, and indexing strategies.
- **Key Relationships Mapped**:
  - `ACTIVITIES`: Primary table (24,527,044 rows). Foreign keys: `molregno` (compound), `assay_id` (assay), `doc_id` (document).
  - `ASSAYS`: Assay definitions (1,970,438 rows). Foreign key: `tid` (target). Contains `bao_format` (BioAssay Ontology format) and `description`.
  - `ASSAY_PARAMETERS`: Specific experimental conditions (460,883 rows) linking `assay_id` to `standard_type` (e.g., `[ATP]`, `TIME`, `SUBSTRATE`).
  - `TARGET_DICTIONARY`: Target definitions (18,552 rows). Contains `target_type` (SINGLE PROTEIN, PROTEIN COMPLEX, etc.).
  - `TARGET_COMPONENTS`: Bridge table mapping `tid` to `component_id` (UniProt accessions).
  - `COMPONENT_CLASS`: Hierarchical protein classification (Levels 1–4).
  - `MOLECULE_DICTIONARY`: Compound definitions (2,921,148 rows) with `max_phase`, `chembl_id`.
  - `COMPOUND_STRUCTURES`: Chemical structures with `canonical_smiles`, `standard_inchi_key`.
- **Outputs**: [`data/schema/chembl_37_schema_report.tsv`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/schema/chembl_37_schema_report.tsv)
- **Log**: [`logs/schema_explore_20261009_180409.log`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/logs/schema_explore_20261009_180409.log)

---

### Step 3: Lossless Broad Bioactivity Curation
- **Script**: [`scripts/03_curate_bioactivity.py`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/scripts/03_curate_bioactivity.py)
- **Design Philosophy**: Lossless broad curation. Preserves all endpoints, all organisms, mutants, and qualitative measurements, while strictly filtering out verified data-entry and transcription errors.
- **Filtering Rules**:
  - Excluded records where `activities.data_validity_comment` matches:
    - `Potential transcription error`
    - `Outside typical range`
    - `Non standard unit for type`
    - `Author confirmed error`
    - `Potential author error`
  - Retained: 23,246,904 valid bioactivity determinations (filtered out 1,280,140 flagged or corrupted rows).
- **Joined Attributes**:
  - Compound identifiers: `chembl_id_compound`, `canonical_smiles`, `compound_inchikey`.
  - Activity details: `standard_type`, `standard_relation`, `standard_value`, `standard_units`, `pchembl_value`, `activity_comment`, `modality` (Targeted Protein Degradation).
  - Assay details: `chembl_id_assay`, `assay_type`, `assay_organism`, `confidence_score`.
  - Target details: `chembl_id_target`, `target_pref_name`, `target_type`, `target_organism`.
  - Document details: `journal`, `year`, `doi`, `pubmed_id`.
  - Ligand efficiency metrics: `le` (Ligand Efficiency), `bei` (Binding Efficiency Index), `lle` (Lipophilic Ligand Efficiency).
- **Performance Architecture**: Chunked streaming from SQLite (100,000 rows per batch) with gzip-compressed tab-separated streaming to maintain bounded memory footprint (<500 MB RAM).
- **Outputs**: [`data/curated/chembl_37_bioactivity_broad.tsv.gz`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/curated/chembl_37_bioactivity_broad.tsv.gz) (1.2 GB compressed)
- **Log**: [`logs/curate_bioactivity_20261009_183026.log`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/logs/curate_bioactivity_20261009_183026.log)

---

### Step 4: Evidence-Tiered Compound-Target Interaction Curation
- **Script**: [`scripts/04_curate_target_interactions.py`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/scripts/04_curate_target_interactions.py)
- **Scientific Rationale (Addressing Reviewer 2 Audit)**:
  An assay associated with a target does not necessarily demonstrate direct target engagement. To prevent downstream false-positive pharmacology networks, records are stratified into five mutually exclusive confidence tiers:
  
  | Evidence Tier | Target Type | Assay Type | Assignment Confidence | Biological Interpretation |
  |:---:|:---|:---|:---:|:---|
  | **Tier A** | `SINGLE PROTEIN` | `B` (Binding) | High ($\ge 7$) | **High-Confidence Direct Binding**: Direct molecular interaction against a defined single protein. |
  | **Tier B** | `PROTEIN COMPLEX` | `B` (Binding) | Moderate–High | **Protein Complex Engagement**: Interaction with a multi-protein assembly. |
  | **Tier C** | `SINGLE PROTEIN` (Mutant) | Any | Any | **Mutant / Variant Interaction**: Explicit amino acid substitution or resistance mutation. |
  | **Tier D** | Any | `F` (Functional) | Any | **Functional / Cellular Response**: Pathway activation, reporter expression, or second messenger readout. |
  | **Tier E** | `ORGANISM`, `TISSUE`, `CELL-LINE` | Any | Low–Any | **Phenotypic / Whole-Organism Response**: Tissue contraction, organism lethality, or phenotypic screens. |

- **Auxiliary Curation Tables Generated**:
  1. `data/curated/chembl_37_target_components.tsv`: Normalized 1:1 mapping between `tid`, `chembl_id_target`, and UniProt accession. Resolves multi-protein complexes into constituent subunits with individual UniProt IDs.
  2. `data/curated/chembl_37_drug_dictionary.tsv`: Complete registry of 16,929 investigational and approved drugs with `max_phase`, `first_approval`, `oral`/`parenteral` flags, black box warnings, and WHO ATC codes.
  3. `data/curated/chembl_37_target_variants.tsv`: Complete registry of target mutations, isoforms, and resistance variants.
  4. `data/curated/chembl_37_protein_classification.tsv`: ChEMBL protein family hierarchy across 4 levels (e.g., *Enzyme* $\rightarrow$ *Kinase* $\rightarrow$ *Protein Kinase* $\rightarrow$ *TK group*).
- **Outputs**: [`data/curated/chembl_37_compound_target_interactions_tiered.tsv.gz`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/curated/chembl_37_compound_target_interactions_tiered.tsv.gz) (16,903,245 rows)
- **Log**: [`logs/curate_targets_20261012_105703.log`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/logs/curate_targets_20261012_105703.log)

---

### Step 5: Master Curation Reporting & Run Manifest
- **Script**: [`scripts/05_generate_report.py`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/scripts/05_generate_report.py)
- **Objective**: Synthesize statistical distributions, compute SHA256 checksums, and generate the formal reproducibility run manifest.
- **Outputs**:
  - Markdown Report: [`reports/chembl_37_curation_summary.md`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/reports/chembl_37_curation_summary.md)
  - Summary TSV: [`reports/chembl_37_curation_summary.tsv`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/reports/chembl_37_curation_summary.tsv)
  - Run Manifest: [`reports/run_manifest.json`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/reports/run_manifest.json)
- **Log**: [`logs/report_20261012_110641.log`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/logs/report_20261012_110641.log)

---

### Step 6: Well-Characterized Proteins Extraction & Directory Partitioning
- **Script**: [`scripts/06_extract_well_characterized_targets.py`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/scripts/06_extract_well_characterized_targets.py)
- **Definition of Well-Characterized Protein**:
  - Target type must be `SINGLE PROTEIN`.
  - Target must have at least **500 unique tested chemical compounds** with high-confidence quantitative affinity/potency measurements ($K_d, K_i, IC_{50}, EC_{50}$).
  - Total qualified targets identified in ChEMBL 37: **902 single-protein targets**.
- **Execution Architecture**:
  - Avoided SQLite disk-sorting bottlenecks by streaming directly through `chembl_37_bioactivity_broad.tsv.gz` using Python's `csv.DictReader` with an in-memory lookup index.
  - Streaming throughput: ~120,000 rows/second with low memory consumption (<650 MB RAM).
- **Per-Target Folder Structure**:
  For each of the 902 targets, a dedicated directory is generated at `data/targets/<CHEMBL_ID>_<Target_Slug>/`:
  - `Kd.tsv`: Dissociation constants ($K_d$).
  - `Ki.tsv`: Inhibition constants ($K_i$).
  - `IC50.tsv`: Half-maximal inhibitory concentrations ($IC_{50}$).
  - `EC50.tsv`: Half-maximal effective concentrations ($EC_{50}$).
  - `target_metadata.json`: JSON document recording target accession, organism, protein class, total determinations, approved drug count, clinical drug count, and names of approved drugs.
- **Statistical Table Generated**:
  [`reports/well_characterized_proteins_drug_stats.tsv`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/reports/well_characterized_proteins_drug_stats.tsv)
  - 2,959 target × endpoint rows detailing:
    - Target ID, preferred name, UniProt accession, organism, protein family.
    - Endpoint type ($K_d, K_i, IC_{50}, EC_{50}$).
    - Total bioactivity records ($N_{\text{records}}$).
    - Unique compounds tested ($N_{\text{compounds}}$).
    - Approved drug count ($N_{\text{approved\_drugs}}$).
    - Clinical candidate count ($N_{\text{clinical\_drugs}}$).
    - Comma-delimited list of approved drug generic names.
- **Outputs**:
  - 902 target directories under `data/targets/` (total 2,496,289 bioactivity determinations).
  - Target Index: [`data/targets/targets_index.tsv`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/targets_index.tsv)
  - Target Guide: [`data/targets/README.md`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/README.md)
  - Executive Report: [`reports/well_characterized_proteins_summary.md`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/reports/well_characterized_proteins_summary.md)
- **Log**: [`logs/extract_well_characterized_targets_20261012_114642.log`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/logs/extract_well_characterized_targets_20261012_114642.log)

---

### Step 7: Assay Inputs & Experimental Measurement Conditions Extraction
- **Script**: [`scripts/07_extract_assay_conditions.py`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/scripts/07_extract_assay_conditions.py)
- **Scientific Objective**:
  $IC_{50}$ values are conditional measurements. To enable condition-aware machine learning, this step extracts experimental parameters for **all 3,435,223 $IC_{50}$ determinations** across **293,458 distinct assays** in ChEMBL 37.
- **Information Extraction Pipeline**:
  1. **BioAssay Ontology (BAO) Format**: Mapped from `assays.bao_format` into standardized biological formats:
     - `CELL-BASED (INTACT WHOLE CELL)`: BAO_0000219
     - `BIOCHEMICAL (PURIFIED PROTEIN)`: BAO_0000357
     - `BIOCHEMICAL BINDING`: BAO_0000221
     - `BIOCHEMICAL (CELL MEMBRANE PREP)`: BAO_0000224
     - `ORGANISM-BASED (WHOLE ORGANISM)`: BAO_0000218
     - `FUNCTIONAL CELL-BASED`: BAO_0000010
     - `TISSUE-BASED (EX VIVO)`: BAO_0000222
     - `BIOCHEMICAL (CELL-FREE / MICROSOMAL)`: BAO_0000054
  2. **Structured Parameters (`assay_parameters` table)**:
     - 460,883 parameter rows queried for explicit substrate concentrations ($[S]$), cofactor concentrations ($[ATP]$), radioligand concentrations, and incubation durations.
  3. **NLP Free-Text Description Parsing**:
     - Fast regex pattern matching applied across `assays.description` to extract:
       - **Reaction Duration**: Matched against patterns such as `\b(\d+(?:\.\d+)?)\s*(min|mins|minute|minutes|h|hr|hrs|hour|hours|day|days|s|sec|seconds)\b`. Bucketed into: `Not Stated`, `<= 15 min`, `16 - 60 min`, `1 - 4 hours`, `4 - 24 hours`, `> 24 hours`.
       - **Detection Methodology**: Categorized into 12 distinct biophysical technologies:
         - `TR-FRET / HTRF`
         - `Fluorescence Polarization (FP)`
         - `Fluorescence Intensity`
         - `Luminescence / Luciferase`
         - `AlphaScreen / AlphaLISA`
         - `Radioligand Binding / SPA`
         - `Colorimetric / MTT / ELISA`
         - `Mass Spectrometry`
         - `Electrophysiology / Patch Clamp`
         - `Biochemical Binding Assay`
         - `Functional Assay`
         - `Unspecified Assay Readout`
       - **Substrate / Cofactor Detection**: Pattern matching for ATP concentrations (e.g., `10 uM ATP`, `1 mM ATP`), substrate peptides, or competitor radioligands.
- **Output Datasets & Reports**:
  - Enriched Bioactivity Dataset: [`data/curated/assay_conditions/ic50_assay_conditions.tsv.gz`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/curated/assay_conditions/ic50_assay_conditions.tsv.gz) (126 MB compressed, 3,435,223 rows).
  - Assay Protocols Index: [`data/curated/assay_conditions/assay_protocols_conditions_index.tsv`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/curated/assay_conditions/assay_protocols_conditions_index.tsv) (36 MB, 293,458 assay protocols).
  - Statistical Table: [`reports/assay_conditions_stats.tsv`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/reports/assay_conditions_stats.tsv).
  - Scientific Report: [`reports/assay_conditions_summary.md`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/reports/assay_conditions_summary.md).
- **Log**: [`logs/extract_assay_conditions_20261012_141443.log`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/logs/extract_assay_conditions_20261012_141443.log)

---

## 4. Key Scientific Findings & Guidance for Machine Learning

### 1. The $IC_{50}$ Format Confounder
- **Discovery**: 39.9% of all $IC_{50}$ measurements in ChEMBL are intact whole-cell assays, whereas 49.6% are cell-free biochemical assays (purified protein, membrane preps, or cell-free lysates).
- **Impact on ML**: Training a universal compound-protein prediction model on pooled $IC_{50}$ values without format conditioning forces the model to learn contradictory physical relationships. A cell-based $IC_{50}$ reflects cell permeability, active transport, and competition against 1–5 mM intracellular ATP. A cell-free $IC_{50}$ reflects direct binding under artificial buffer conditions (often 10 $\mu\text{M}$ ATP).
- **Recommendation**:
  - When training structural affinity models (e.g., docking, 3D binding affinity, $K_d/K_i$ predictors), use **Tier A** biochemical binding records or filter $IC_{50}$ by `assay_format_category == 'Biochemical (Purified Protein)'`.
  - For phenotypic or whole-cell potency models, condition your network on the format token:  
    $$\text{Affinity} = f(\text{SMILES}, \text{Protein Sequence}, \text{Assay Format}, [ATP], \text{Incubation Time})$$

### 2. Incubation Time and Irreversible Inhibitors
- **Discovery**: Over 145,000 $IC_{50}$ determinations involve incubation times exceeding 24 hours (primarily cytotoxicity/antiproliferation screens), while biochemical assays typically run for 15 to 60 minutes.
- **Impact on ML**: For covalent or slow-off-rate inhibitors, apparent $IC_{50}$ decreases as incubation time increases because inhibition accumulates over time ($IC_{50}^{\text{app}} \propto 1 / t$). Pooling 1-hour and 72-hour assays distorts potency rankings.

### 3. Optical Artifacts in Single-Wavelength Readouts
- **Discovery**: 10.4% (357,540) of $IC_{50}$ determinations utilize standard single-wavelength fluorescence intensity.
- **Impact on ML**: Compounds containing conjugated aromatic systems, polyphenols, or rhodanine cores frequently exhibit intrinsic fluorescence or quench assay fluorophores. The explicit `detection_method` column enables filtering out putative false positives that appear active solely in fluorescence intensity formats.

---

## 5. Execution Guide & Reproducibility

### Environment Setup
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Complete Pipeline Execution
```bash
# Execute complete pipeline (Steps 1 to 7)
bash run_pipeline.sh
```

### Running Individual Steps
```bash
# Step 1: Download ChEMBL 37 SQLite + schema
bash scripts/01_download_chembl.sh

# Step 2: Explore database schema
python3 scripts/02_explore_schema.py

# Step 3: Curate broad bioactivity
python3 scripts/03_curate_bioactivity.py

# Step 4: Curate tiered target interactions
python3 scripts/04_curate_target_interactions.py

# Step 5: Generate master report and manifest
python3 scripts/05_generate_report.py

# Step 6: Extract well-characterized protein directories and endpoint tables
python3 scripts/06_extract_well_characterized_targets.py

# Step 7: Extract assay formats, conditions, and detection technologies
python3 scripts/07_extract_assay_conditions.py
```

---

## 6. Verification and Audit

All datasets and code changes were evaluated against the Scientific Quality Council release criteria:
1. **Input Integrity**: Validated against official EMBL-EBI MD5 checksums. Zero rows lost to unhandled exceptions.
2. **Scientific Validity**: Preserved physical distinction between thermodynamic parameters ($K_d, K_i$) and functional/format-dependent readouts ($IC_{50}, EC_{50}$).
3. **Statistical Validity**: Verified row-count conservation across transformations; no silent imputation or duplicate fabrication.
4. **Code / Numerical Correctness**: Memory-bounded streaming execution with streaming gzip generation.
5. **Traceability**: All steps produce timestamps, row counts, parameters, and outputs in dedicated log files under `logs/`.
