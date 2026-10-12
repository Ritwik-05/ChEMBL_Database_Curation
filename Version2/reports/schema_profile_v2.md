# ChEMBL 37 Database Schema Profile (Version 2)
> **Generated**: 2026-10-12 15:07:49
> **Database**: `data/raw/chembl_37.db` (Read-only mode)
> **Total Tables Profiled**: 72

## 1. Core and Condition Tables Profile

| Table Name | Row Count | Columns | Primary Key | Key Indexes |
|:---|:---:|:---:|:---|:---|
| `activities` | 24,527,044 | 29 | `activity_id` | 20 indexes |
| `activity_properties` | 12,213,211 | 14 | `ap_id` | 8 indexes |
| `assay_parameters` | 460,048 | 13 | `assay_param_id` | 11 indexes |
| `assays` | 1,970,438 | 24 | `assay_id` | 8 indexes |
| `bioassay_ontology` | 311 | 2 | `bao_id` | 1 indexes |
| `component_class` | 13,180 | 3 | `comp_class_id` | 2 indexes |
| `component_sequences` | 12,986 | 10 | `component_id` | 2 indexes |
| `compound_structures` | 2,897,819 | 5 | `molregno` | 7 indexes |
| `docs` | 101,100 | 19 | `doc_id` | 8 indexes |
| `drug_indication` | 60,055 | 8 | `drugind_id` | 3 indexes |
| `drug_mechanism` | 7,561 | 14 | `mec_id` | 1 indexes |
| `drug_warning` | 2,304 | 11 | `warning_id` | 1 indexes |
| `molecule_dictionary` | 2,921,148 | 28 | `molregno` | 6 indexes |
| `protein_classification` | 905 | 7 | `protein_class_id` | 2 indexes |
| `target_components` | 17,284 | 4 | `targcomp_id` | 2 indexes |
| `target_dictionary` | 18,552 | 7 | `tid` | 7 indexes |
| `variant_sequences` | 2,836 | 8 | `variant_id` | 2 indexes |

## 2. Parameter and Property Tables Detail

- **`ASSAY_PARAMETERS`**: Contains assay-level experimental parameters (e.g. $[ATP]$, substrate, incubation duration, temperature).
- **`ACTIVITY_PROPERTIES`**: Contains 12.2M measurement-specific properties (e.g. timepoint, compound concentrations, replicate counts).
- **`BIOASSAY_ONTOLOGY`**: Formal BAO vocabulary mapping assay format codes (e.g. BAO_0000219 for cell-based format, BAO_0000357 for single protein format).
