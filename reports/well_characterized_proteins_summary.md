# ChEMBL 37 — Well-Characterized Proteins Bioactivity & Drug Statistics Report
> **Generated**: 2026-10-12 11:51:23
> **Target Criterion**: Single-protein targets with $\ge 500$ unique compounds tested against high-confidence endpoints ($K_d, K_i, IC_50, EC_50$).
> **Total Targets Extracted**: 902

---

## 1. Executive Summary

- **Total Well-Characterized Targets**: 902
- **Total Bioactivity Records Extracted**: 2,496,289
- **Organism Breakdown**: 739 Human (*Homo sapiens*), 57 Rat (*Rattus norvegicus*), 29 Mouse (*Mus musculus*), 7 Guinea pig, 6 Bovine, 5 HIV-1, 4 *P. falciparum*, and others.
- **Endpoint Stratification**: Each protein folder isolates $K_d, K_i, IC_50, EC_50$ into distinct TSVs to prevent pharmacological pooling errors.

## 2. Top 25 Targets by Approved & Clinical Drug Representation

This table ranks targets by the breadth of clinical therapeutics ($N_{\text{Clinical}}$) and approved medicines ($N_{\text{Approved}}$) evaluated against them.

| ChEMBL ID | Target Preferred Name | Organism | UniProt | Total Compounds | Total Activities | Approved Drugs | Clinical Drugs | Available Endpoints |
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| [`CHEMBL6020`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL6020_Bile_salt_export_pump) | Bile salt export pump | *Homo sapiens* | `O95342` | 440 | 634 | **232** | **296** | `Ki, IC50` |
| [`CHEMBL240`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL240_Voltage_gated_inwardly_rectifying_potass) | Voltage-gated inwardly rectifying potassium c | *Homo sapiens* | `Q12809` | 15,994 | 19,741 | **227** | **478** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL340`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL340_Cytochrome_P450_3A4) | Cytochrome P450 3A4 | *Homo sapiens* | `P08684` | 9,953 | 12,434 | **192** | **466** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL289`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL289_Cytochrome_P450_2D6) | Cytochrome P450 2D6 | *Homo sapiens* | `P10635` | 6,835 | 7,636 | **170** | **403** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL3397`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL3397_Cytochrome_P450_2C9) | Cytochrome P450 2C9 | *Homo sapiens* | `P11712` | 6,035 | 6,556 | **148** | **373** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL222`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL222_Sodium_dependent_noradrenaline_transport) | Sodium-dependent noradrenaline transporter | *Homo sapiens* | `P23975` | 4,130 | 5,507 | **131** | **186** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL1867`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL1867_Alpha_2A_adrenergic_receptor) | Alpha-2A adrenergic receptor | *Homo sapiens* | `P08913` | 1,255 | 1,750 | **130** | **176** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL228`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL228_Sodium_dependent_serotonin_transporter) | Sodium-dependent serotonin transporter | *Homo sapiens* | `P31645` | 6,027 | 7,956 | **128** | **177** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL3622`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL3622_Cytochrome_P450_2C19) | Cytochrome P450 2C19 | *Homo sapiens* | `P33261` | 3,790 | 4,103 | **125** | **303** | `Ki, IC50, EC50` |
| [`CHEMBL224`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL224_5_hydroxytryptamine_receptor_2A) | 5-hydroxytryptamine receptor 2A | *Homo sapiens* | `P28223` | 6,949 | 9,818 | **124** | **197** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL234`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL234_D_3_dopamine_receptor) | D(3) dopamine receptor | *Homo sapiens* | `P35462` | 6,666 | 10,119 | **122** | **190** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL225`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL225_5_hydroxytryptamine_receptor_2C) | 5-hydroxytryptamine receptor 2C | *Homo sapiens* | `P28335` | 4,674 | 6,377 | **120** | **181** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL203`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL203_Epidermal_growth_factor_receptor) | Epidermal growth factor receptor | *Homo sapiens* | `P00533` | 14,536 | 28,200 | **118** | **351** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL3356`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL3356_Cytochrome_P450_1A2) | Cytochrome P450 1A2 | *Homo sapiens* | `P05177` | 4,548 | 4,876 | **117** | **300** | `Ki, IC50, EC50` |
| [`CHEMBL217`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL217_D_2_dopamine_receptor) | D(2) dopamine receptor | *Homo sapiens* | `P14416` | 11,028 | 17,317 | **115** | **208** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL1833`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL1833_5_hydroxytryptamine_receptor_2B) | 5-hydroxytryptamine receptor 2B | *Homo sapiens* | `P41595` | 2,535 | 3,542 | **114** | **168** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL1841`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL1841_Tyrosine_protein_kinase_Fyn) | Tyrosine-protein kinase Fyn | *Homo sapiens* | `P06241` | 2,240 | 3,560 | **110** | **302** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL1942`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL1942_Alpha_2B_adrenergic_receptor) | Alpha-2B adrenergic receptor | *Homo sapiens* | `P18089` | 744 | 1,107 | **109** | **149** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL238`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL238_Sodium_dependent_dopamine_transporter) | Sodium-dependent dopamine transporter | *Homo sapiens* | `Q01959` | 4,144 | 5,709 | **107** | **155** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL216`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL216_Muscarinic_acetylcholine_receptor_M1) | Muscarinic acetylcholine receptor M1 | *Homo sapiens* | `P11229` | 3,457 | 4,416 | **107** | **137** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL287`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL287_Sigma_non_opioid_intracellular_receptor) | Sigma non-opioid intracellular receptor 1 | *Homo sapiens* | `Q99720` | 3,666 | 4,825 | **103** | **160** | `Kd, Ki, IC50` |
| [`CHEMBL211`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL211_Muscarinic_acetylcholine_receptor_M2) | Muscarinic acetylcholine receptor M2 | *Homo sapiens* | `P08172` | 2,693 | 3,507 | **97** | **126** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL4302`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL4302_ATP_dependent_translocase_ABCB1) | ATP-dependent translocase ABCB1 | *Homo sapiens* | `P08183` | 2,769 | 4,137 | **96** | **133** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL1916`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL1916_Alpha_2C_adrenergic_receptor) | Alpha-2C adrenergic receptor | *Homo sapiens* | `P18825` | 860 | 1,358 | **96** | **132** | `Kd, Ki, IC50, EC50` |
| [`CHEMBL245`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL245_Muscarinic_acetylcholine_receptor_M3) | Muscarinic acetylcholine receptor M3 | *Homo sapiens* | `P20309` | 3,213 | 4,246 | **95** | **120** | `Kd, Ki, IC50, EC50` |

## 3. Top 25 Targets by Total Bioactivity Screening Volume

| ChEMBL ID | Target Preferred Name | Organism | UniProt | Total Compounds | Total Activities | Approved Drugs | Endpoints |
|:---|:---|:---:|:---:|:---:|:---:|:---:|:---|
| [`CHEMBL240`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL240_Voltage_gated_inwardly_rectifying_potass) | Voltage-gated inwardly rectifying potassium c | *Homo sapiens* | `Q12809` | **15,994** | 19,741 | 227 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL2971`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL2971_Tyrosine_protein_kinase_JAK2) | Tyrosine-protein kinase JAK2 | *Homo sapiens* | `O60674` | **15,439** | 23,478 | 62 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL203`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL203_Epidermal_growth_factor_receptor) | Epidermal growth factor receptor | *Homo sapiens* | `P00533` | **14,536** | 28,200 | 118 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL5251`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL5251_Tyrosine_protein_kinase_BTK) | Tyrosine-protein kinase BTK | *Homo sapiens* | `Q06187` | **13,414** | 23,709 | 75 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL279`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL279_Vascular_endothelial_growth_factor_recep) | Vascular endothelial growth factor receptor 2 | *Homo sapiens* | `P35968` | **13,321** | 18,557 | 57 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL2835`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL2835_Tyrosine_protein_kinase_JAK1) | Tyrosine-protein kinase JAK1 | *Homo sapiens* | `P23458` | **12,924** | 19,741 | 71 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL1163125`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL1163125_Bromodomain_containing_protein_4) | Bromodomain-containing protein 4 | *Homo sapiens* | `O60885` | **12,514** | 23,776 | 12 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL1741186`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL1741186_Nuclear_receptor_ROR_gamma) | Nuclear receptor ROR-gamma | *Homo sapiens* | `P51449` | **12,129** | 18,922 | 3 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL217`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL217_D_2_dopamine_receptor) | D(2) dopamine receptor | *Homo sapiens* | `P14416` | **11,028** | 17,317 | 115 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL3130`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL3130_Phosphatidylinositol_4_5_bisphosphate_3) | Phosphatidylinositol 4,5-bisphosphate 3-kinas | *Homo sapiens* | `O00329` | **10,849** | 15,592 | 35 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL4822`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL4822_Beta_secretase_1) | Beta-secretase 1 | *Homo sapiens* | `P56817` | **10,483** | 16,257 | 5 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL3778`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL3778_Interleukin_1_receptor_associated_kinase) | Interleukin-1 receptor-associated kinase 4 | *Homo sapiens* | `Q9NWZ3` | **10,187** | 13,210 | 65 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL2148`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL2148_Tyrosine_protein_kinase_JAK3) | Tyrosine-protein kinase JAK3 | *Homo sapiens* | `P52333` | **10,143** | 14,697 | 46 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL253`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL253_Cannabinoid_receptor_2) | Cannabinoid receptor 2 | *Homo sapiens* | `P34972` | **10,007** | 13,661 | 9 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL205`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL205_Carbonic_anhydrase_2) | Carbonic anhydrase 2 | *Homo sapiens* | `P00918` | **9,974** | 14,326 | 67 | `Kd, Ki, IC50` |
| [`CHEMBL340`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL340_Cytochrome_P450_3A4) | Cytochrome P450 3A4 | *Homo sapiens* | `P08684` | **9,953** | 12,434 | 192 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL325`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL325_Histone_deacetylase_1) | Histone deacetylase 1 | *Homo sapiens* | `Q13547` | **9,771** | 12,364 | 20 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL3553`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL3553_Non_receptor_tyrosine_protein_kinase_TYK) | Non-receptor tyrosine-protein kinase TYK2 | *Homo sapiens* | `P29597` | **9,615** | 16,666 | 73 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL233`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL233_Mu_type_opioid_receptor) | Mu-type opioid receptor | *Homo sapiens* | `P35372` | **9,601** | 13,977 | 90 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL4005`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL4005_Phosphatidylinositol_4_5_bisphosphate_3) | Phosphatidylinositol 4,5-bisphosphate 3-kinas | *Homo sapiens* | `P42336` | **9,481** | 12,399 | 34 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL218`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL218_Cannabinoid_receptor_1) | Cannabinoid receptor 1 | *Homo sapiens* | `P21554` | **9,370** | 12,928 | 14 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL261`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL261_Carbonic_anhydrase_1) | Carbonic anhydrase 1 | *Homo sapiens* | `P00915` | **8,896** | 12,236 | 56 | `Kd, Ki, IC50` |
| [`CHEMBL4792`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL4792_Orexin_receptor_type_2) | Orexin receptor type 2 | *Homo sapiens* | `O43614` | **8,719** | 13,289 | 4 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL251`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL251_Adenosine_receptor_A2a) | Adenosine receptor A2a | *Homo sapiens* | `P29274` | **8,395** | 10,673 | 25 | `Kd, Ki, IC50, EC50` |
| [`CHEMBL4296`](file:///home/cbb/My_GIT_Projects/ChEMBL_Database_Curation/data/targets/CHEMBL4296_Sodium_channel_protein_type_9_subunit_al) | Sodium channel protein type 9 subunit alpha | *Homo sapiens* | `Q15858` | **8,379** | 10,609 | 13 | `Kd, Ki, IC50, EC50` |

## 4. Endpoint-Level Statistics Across All Targets

| Endpoint Type | Total Activity Records | Biophysical Interpretation | Recommended ML Task |
|:---|:---:|:---|:---|
| **$IC_{50}$** | 1,589,907 | Biochemical concentration producing 50% inhibition | Endpoint-stratified regression / Multitask QSAR |
| **$K_i$** | 597,050 | Absolute thermodynamic enzyme-inhibitor dissociation constant | Direct physical binding regression ($\text{p}K_i$) |
| **$K_d$** | 108,185 | Equilibrium physical dissociation constant (SPR, ITC, radioligand) | Gold-standard 3D structural affinity modeling ($\text{p}K_d$) |
| **$EC_{50}$** | 201,147 | Cellular or functional concentration yielding 50% effect | Functional agonism / signaling pathway prediction |

## 5. Directory Structure & File Organization

```
data/targets/
├── README.md                          # Guide for ML practitioners
├── targets_index.tsv                  # Index table of all 902 targets
├── CHEMBL203_Epidermal_growth_factor/ # EGFR target folder
│   ├── IC50.tsv                       # 50% inhibitory concentration records
│   ├── Ki.tsv                         # Equilibrium inhibition constants
│   ├── Kd.tsv                         # Equilibrium binding dissociation constants
│   ├── EC50.tsv                       # Functional potency measurements
│   └── target_metadata.json           # Target metadata and drug counts
└── ... (902 target folders)
```

---
*Report generated by Scientific Quality Council — 2026-10-12 11:51:23*
