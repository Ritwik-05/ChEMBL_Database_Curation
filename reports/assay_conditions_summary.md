# ChEMBL 37 — IC50 Assay Experimental Conditions & Measurement Profile
> **Generated**: 2026-10-12 14:26:59
> **Population**: Complete ChEMBL 37 $IC_{50}$ Bioactivity Records (3,435,223 determinations across 293,458 assays)
> **Dataset Path**: [`data/curated/assay_conditions/ic50_assay_conditions.tsv.gz`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/curated/assay_conditions/ic50_assay_conditions.tsv.gz)

---

## 1. Executive Summary

An $IC_{50}$ value is not an intrinsic thermodynamic constant; it is an experimental readout governed by assay format, substrate/cofactor concentration, incubation time, and detection methodology. This report describes **how $IC_{50}$ was measured** across all 3.44 million determinations in ChEMBL 37.

| Overall Metric | Count | Percentage | Pharmacological Interpretation |
|:---|:---:|:---:|:---|
| **Total $IC_{50}$ Determinations** | 3,435,223 | 100.0% | Complete curated population in ChEMBL 37 |
| **Distinct Assay Protocols** | 293,458 | — | Independent experimental protocol definitions |
| **Standardized Affinity ($p\text{ChEMBL}$)** | 2,489,470 | 72.5% | Values with $-\log_{10}[\text{M}]$ available |
| **Explicit Detection Technology Identified** | 1,595,085 | 46.4% | Optical / biophysical readout technology mapped |
| **Incubation Time Documented** | 1,021,187 | 29.7% | Explicit reaction duration recorded |
| **Substrate / Cofactor Documented** | 214,086 | 6.2% | Explicit $[S], [ATP],$ or radioligand captured |

## 2. Assay Format Distribution (The Primary Confounder)

A critical pharmacological discovery: **over 40% of all $IC_{50}$ values in ChEMBL are measured in whole cells**, reflecting cellular phenotypic uptake and pathway response rather than direct physical enzyme inhibition.

| Biological Format Category | Determination Count | % of All $IC_{50}$ | Biophysical & Pharmacological Interpretation |
|:---|:---:|:---:|:---|
| **Cell-based (Intact Whole Cell)** | 1,372,146 | 39.9% | Cellular whole-cell inhibition |
| **Biochemical (Purified Protein)** | 1,121,740 | 32.7% | Direct target inhibition in cell-free system |
| **Biochemical Binding** | 461,500 | 13.4% | Direct target inhibition in cell-free system |
| **Organism-based (Whole Organism)** | 156,628 | 4.6% | In vivo / ex vivo whole-organism response |
| **Functional Cell-based** | 149,473 | 4.4% | Cellular whole-cell inhibition |
| **Biochemical (Cell Membrane Prep)** | 92,858 | 2.7% | Direct target inhibition in cell-free system |
| **Tissue-based (Ex Vivo Organ / Tissue)** | 44,400 | 1.3% | In vivo / ex vivo whole-organism response |
| **Biochemical (Cell-free / Microsomal)** | 28,332 | 0.8% | Direct target inhibition in cell-free system |
| **Unclassified In Vitro** | 8,146 | 0.2% | In vivo / ex vivo whole-organism response |

## 3. Detection Technology & Instrumentation Readout

Detection technology governs vulnerability to optical interference (e.g., compound autofluorescence, turbidity, quenching):

| Readout Methodology | Measurement Count | % Share | Confounders & Characteristics |
|:---|:---:|:---:|:---|
| **Unspecified Assay Readout** | 1,042,651 | 30.4% | Standard metabolic or optical readout |
| **Functional Assay** | 689,028 | 20.1% | Standard metabolic or optical readout |
| **Colorimetric / MTT / ELISA** | 512,972 | 14.9% | Standard metabolic or optical readout |
| **Fluorescence Intensity** | 357,540 | 10.4% | Susceptible to compound autofluorescence |
| **TR-FRET / HTRF** | 263,131 | 7.7% | Homogeneous proximity readout (TR-FRET/HTRF) |
| **Radioligand Binding / SPA** | 160,097 | 4.7% | Gold-standard radioactive competition |
| **Luminescence / Luciferase** | 150,742 | 4.4% | Standard metabolic or optical readout |
| **Biochemical Binding Assay** | 108,459 | 3.2% | Standard metabolic or optical readout |
| **Fluorescence Polarization (FP)** | 52,350 | 1.5% | Susceptible to compound autofluorescence |
| **AlphaScreen / AlphaLISA** | 43,434 | 1.3% | Standard metabolic or optical readout |
| **Mass Spectrometry** | 26,457 | 0.8% | Standard metabolic or optical readout |
| **Electrophysiology / Patch Clamp** | 25,660 | 0.7% | Standard metabolic or optical readout |

## 4. Reaction Incubation Time Distribution

Incubation duration is critical for identifying tight-binding, slow off-rate ($k_{\text{off}}$), or irreversible covalent inhibitors (where apparent $IC_{50}$ decreases with time):

| Incubation Duration Bucket | Count | % Share |
|:---|:---:|:---:|
| `Not Stated` | 2,414,036 | 70.3% |
| `16 - 60 min` | 443,341 | 12.9% |
| `<= 15 min` | 244,233 | 7.1% |
| `> 24 hours` | 145,554 | 4.2% |
| `1 - 4 hours` | 138,676 | 4.0% |
| `4 - 24 hours` | 49,383 | 1.4% |

## 5. Machine Learning Guidance for Condition-Conditioned Modeling

1. **Do not pool Cell-based $IC_{50}$ with Biochemical $IC_{50}$**: In whole-cell assays, intracellular ATP is millimolar ($\sim 1	ext{--}5\,	ext{mM}$), whereas cell-free kinase assays typically use $10\,\mu	ext{M}$ ATP. The same kinase inhibitor can show a $300\times$ shift in $IC_{50}$ purely due to ATP competition.
2. **Condition-Aware Modeling**: Condition your neural network on the extracted feature vector: `Model(SMILES, FASTA, Format, [ATP], Time) -> IC50`.
3. **Filter by Optical Technology**: When screening for potential false positives, flag hits derived solely from single-wavelength Fluorescence Intensity assays that contain fluorophore-like substructures.

---
*Report generated by Scientific Quality Council — 2026-10-12 14:26:59*
