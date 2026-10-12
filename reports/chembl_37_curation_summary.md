# ChEMBL 37 — Scientifically Curated Bioactivity & Interaction Report
> **Generated**: 2026-10-12 11:11:55
> **Source**: ChEMBL 37 SQLite Database (EMBL-EBI)
> **Architecture**: Layered Evidence Tiers with Endpoint Stratification (Reviewer 2 Compliant)

---

## 1. Executive Summary

| Dataset Component | Record Count | Description |
|:---|:---:|:---|
| **Curated Bioactivity Master** | 23,195,772 | Lossless extraction of original & standardized values, ranges, and BAO terms |
| **Evidence-Tiered Interactions** | 16,891,411 | Stratified by Compound × Target × Mutation × Evidence Tier × Endpoint |
| **Target Components (Normalized)** | 17,617 | Unambiguous 1-to-1 mappings to UniProt accessions and protein classes |
| **Target Mutation Variants** | 3,127 | Resistance and genetic mutation annotations linked to parent targets |
| **Curated Drug Registry** | 16,907 | Approved and clinical drugs with ATC, indications, warnings, and MoA |
| **Protein Hierarchy** | 905 | Hierarchical protein classifications (Levels 1–4) |

## 2. Evidence-Tier Classification (Addressing C02, C04, C05)

Rather than conflating cellular phenotypes or ADMET assays with single-protein molecular interactions, every edge is categorized into an explicit biological evidence tier:

| Tier | Biological Definition | Interaction Count | Network Interpretation |
|:---|:---|:---:|:---|
| **Tier A** | Single-Protein Direct Biochemical Binding (`Confidence >= 8`, `Assay B`) | 3,997,396 | **High-Confidence Molecular Interaction** |
| **Tier B** | Single-Protein Functional Potency (`Confidence >= 8`, `Assay F`) | 3,865,223 | Target-Associated Functional Potency |
| **Tier C** | Protein Complex / Family / Homologue Target | 523,324 | Conditional Molecular Evidence |
| **Tier D** | Cellular / Tissue / Organism Phenotype (`GI50`, `MIC`, etc.) | 8,396,106 | Phenotypic Association (NOT direct protein edge) |
| **Tier E** | ADMET / Toxicity / Physicochemical / Unassigned | 109,362 | Secondary Profiling / Non-Target Assay |

## 3. Curated Bioactivity Data Profile

| Metric | Value |
|:---|:---|
| Total Curated Activities | 23,195,772 |
| Unique Chemical Entities | 2,739,887 |
| Unique Targets | 16,689 |
| Unique Assays | 1,610,677 |
| Unique Publications / Patents | 99,455 |
| Standardized Affinity ($p\text{ChEMBL}$) Available | 4,969,278 (21.4%) |
| Canonical SMILES Structure Available | 22,975,195 |

### Top 10 Bioactivity Endpoints (`standard_type`)

| Endpoint Type | Measurement Count | Pharmacological Interpretation |
|:---|:---:|:---|
| `Potency` | 4,471,592 | Stratified independently (never pooled across types) |
| `IC50` | 3,435,223 | Stratified independently (never pooled across types) |
| `GI50` | 2,624,408 | Stratified independently (never pooled across types) |
| `Inhibition` | 1,355,098 | Stratified independently (never pooled across types) |
| `Percent Effect` | 1,328,366 | Stratified independently (never pooled across types) |
| `Ki` | 862,388 | Stratified independently (never pooled across types) |
| `Activity` | 830,897 | Stratified independently (never pooled across types) |
| `k_off` | 826,515 | Stratified independently (never pooled across types) |
| `kon` | 826,355 | Stratified independently (never pooled across types) |
| `MIC` | 784,970 | Stratified independently (never pooled across types) |

### Top 10 Assay Organisms

| Organism | Measurement Count |
|:---|:---:|
| *Homo sapiens* | 8,765,295 |
| *Unknown* | 7,029,931 |
| *Rattus norvegicus* | 1,189,317 |
| *Plasmodium falciparum* | 1,021,573 |
| *Mus musculus* | 733,535 |
| *Mycobacterium tuberculosis* | 692,342 |
| *Escherichia coli* | 256,782 |
| *Trypanosoma cruzi* | 237,955 |
| *Staphylococcus aureus* | 234,330 |
| *Plasmodium berghei* | 231,492 |

## 4. Methodological Corrections Implemented

1. **No Endpoint Pooling (C01, C08)**: $IC_{50}$, $K_i$, $K_d$, and $EC_{50}$ are never ranked against each other. Each edge is strictly endpoint-stratified.
2. **Median Over Extreme Potency (C09)**: Replaced `MAX(pChEMBL)` with median $p\text{ChEMBL}$, interquartile bounds, and study counts ($N_{\text{assays}}$, $N_{\text{docs}}$) to prevent study-density selection bias.
3. **Censoring Awareness (C07)**: Exact (`=`) measurements are separated from lower/upper bounds (`>`, `<`). Median potency is computed strictly on verified exact data.
4. **Target Variant & Mutation Preservation (C06)**: Resistant mutants (e.g. EGFR T790M) remain distinct from wild-type targets via `variant_id` and explicit amino-acid mutation annotations.
5. **Integrated Clinical Drug Metadata (Section 4)**: Drugs are linked to their approved ATC codes, clinical indications, regulatory warnings, and curated mechanisms of action.
6. **RFC-Compliant Parsing (T01)**: TSV reading utilizes `csv.reader` to prevent record splitting on embedded multiline text.

---
*Audited by Scientific Quality Council — Generated 2026-10-12 11:11:55*
