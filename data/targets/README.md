# Well-Characterized Protein Bioactivity Datasets (ChEMBL 37)

## Overview
This directory contains dedicated, endpoint-stratified bioactivity datasets for **902 Well-Characterized Protein Targets** from ChEMBL 37.

### Definition of Well-Characterized Protein Target:
1. **Single-Protein Entity**: Target is verified single protein (`target_type = 'SINGLE PROTEIN'`).
2. **High-Confidence Assay Evidence**: Direct target mapping with assay confidence score $\ge 8$.
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
| `pchembl_value` | Negative log10 molar potency ($-\log_{10}[\text{M}]$) |
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
