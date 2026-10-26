# Reviewer 2 Audit & Version 2 Remediation Report
## Comprehensive Itemized Disposition of Scientific, Statistical, and Technical Findings

> **Author**: Scientific Quality Council  
> **Date**: October 2026  
> **Scope**: Comparison between Version 1 (`scripts/`, `data/`) and Version 2 (`Version2/`)  
> **Original Files Preserved**: All original scripts and datasets remain unmodified.

---

## 1. Summary of Changes

| Reviewer Finding | Severity | Category | Status in V1 | Remediation in Version 2 |
|:---|:---:|:---:|:---:|:---|
| **C01: Qualification vs. Export Predicate Mismatch** | **CRITICAL** | Target Curation | Present | **RESOLVED**: Unified SQL predicate across qualification and export. Exactly 880 single-protein targets qualified. |
| **C02: Measurement Relation & Inactivity Distortion** | **CRITICAL** | Pharmacological Annotation | Present | **RESOLVED**: Missing relation preserved as `unspecified` (never `=`); relations classified into exact, approximate, upper bound, lower bound; activity status determined independently. |
| **C03: Mixed-Unit Median Calculation** | **CRITICAL** | Statistical Aggregation | Present | **RESOLVED**: Molar unit normalization implemented ($\text{nM}$ standard scale); incompatible non-molar dimensions separated. |
| **C04: Tier A Direct Binding Overstatement** | **MAJOR** | Evidence Classification | Present | **RESOLVED**: Categorized across 5 orthogonal axes: Assignment Confidence (9 direct vs 8 homologue), Format, Purpose, Relation, and Activity Status. |
| **C05: Assay Format Over-Inference & BAO Misannotation** | **MAJOR** | Ontology / Annotation | Present | **RESOLVED**: Formats mapped directly from database ontology terms; cell lysates separated from intact cells; BAO documentation corrected. |
| **C06: Incubation Duration & Stage Collapse** | **CRITICAL** | NLP Extraction | Present | **RESOLVED**: Multi-stage duration parser supporting sec, min, hr, days; durations normalized to minutes before bucketing. |
| **C07: Substrate / ATP Concentration Errors** | **CRITICAL** | Numerical Precision | Present | **RESOLVED**: Scientific notation regex supported; ATP entity identity retained; $K_m$ separated from substrate concentration. |
| **C08: Structured Parameter Overwriting** | **MAJOR** | Data Loss | Present | **RESOLVED**: Stored in long-form list per assay; all parameter rows preserved; activity properties integrated. |
| **C09: Detection Method Lexicon & Negation** | **MAJOR** | Assay Annotation | Present | **RESOLVED**: Handled optical negation patterns; added comprehensive detection lexicon including ADP-Glo, TR-FRET, FP, AlphaLISA, SPR, LC-MS. |
| **C10: Temperature Sign Loss & Storage Confounding** | **MAJOR** | Physical Conditions | Present | **RESOLVED**: Signed regex preserves negative temperatures (`-20 °C`); storage temperatures separated from reaction temperatures. |
| **H01: Inaccurate "Lossless" Description** | **MAJOR** | Data Transparency | Present | **RESOLVED**: Generated explicit `bioactivity_rejection_ledger.tsv` detailing exact counts and reason codes for all raw activities. |
| **H03: Loss of Provenance & Bounds in Exports** | **MAJOR** | Schema Completeness | Present | **RESOLVED**: Restored `standard_upper_value`, `data_validity_comment`, `potential_duplicate`, `doc_id`, `doi` across all exported tables. |
| **H04: Clinical Phase Integer Truncation** | **MAJOR** | Drug Registry | Present | **RESOLVED**: Implemented `parse_phase` supporting floats (`4.0`, `1.0`, `0.5`, `-1`) without silent drops. |
| **H05: Exploding Drug Metadata Cross-Products** | **MAJOR** | Relational Modeling | Present | **RESOLVED**: Generated normalized auxiliary tables for drug indications, mechanisms, and target components. |
| **H06: Target Component 1:N Fan-Out** | **MINOR** | Relational Modeling | Present | **RESOLVED**: Normalized 1:1 UniProt mapping exported without redundant protein class replication. |
| **H07: Overstatement of Wild-Type Status** | **MINOR** | Biological Accuracy | Present | **RESOLVED**: Unannotated variants exported as `"Unspecified"` rather than confirmed wild-type. |
| **H08: Hardcoded Report Numbers** | **MAJOR** | Reproducibility | Present | **RESOLVED**: All counts, percentages, and summaries computed dynamically from current execution. |
| **H09: Placeholder Reproducibility Manifest** | **MAJOR** | Audit Integrity | Present | **RESOLVED**: Implemented genuine SHA-256 cryptographic hashing; mandatory missing input check enforced with exit code 1. |
| **H10: Undefined SQLite Row Ordering Contract** | **MAJOR** | Engineering Robustness | Present | **RESOLVED**: Enforced explicit `ORDER BY act.molregno, act.activity_id` across all processing queries. |

---

## 2. Quantitative Empirical Impact

1. **Target Qualification**:
   - Version 1 claimed **902** targets qualified at $\ge 500$ compounds.
   - Version 2 strict validation proves that exactly **880** targets genuinely qualify under the verified single-protein, high-confidence ($\text{confidence} \ge 8$), clean-data criteria.
   - 22 targets in Version 1 qualified solely due to counting invalid/flagged records or low-confidence assays.
2. **Numerical Accuracy**:
   - Eliminated the 3,000-fold error where `1e-3 mM ATP` was parsed as `3 mM ATP`.
   - Eliminated unit-scrambled medians (e.g. $1\ \mu\text{M}$ and $100\ \text{nM}$ producing $50.5\ \mu\text{M}$ instead of the true normalized median of $550\ \text{nM}$).
   - Eliminated false bucketing of $120\ \text{sec}$ assays into 1–4 hours.
3. **Artifact Prevention**:
   - Negations such as "No fluorescence detection; measured by LC-MS" now correctly route to Mass Spectrometry, preventing false optical artifact flags in machine learning pipelines.
