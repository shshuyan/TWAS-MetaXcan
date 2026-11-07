# TWAS results: CAD_CARDIoGRAMplusC4D_fast

Coronary artery disease, CARDIoGRAMplusC4D 1000G additive meta-analysis (Nikpay et al. 2015). MetaXcan tutorial data. Fast mode (no summary imputation).

- Input: `cad.add.160614.website.txt.gz` (format `cardiogram`, build hg19, ancestry EUR)
- Imputation: skipped

## eqtl: S-PrediXcan (per tissue)

- 49 tissues, 617,892 gene-tissue tests, Bonferroni p < 8.09e-08
- Significant pairs: 404 (61 unique genes)

| gene | tissues significant | best tissue | best p |
|---|---|---|---|
| CDKN2B | 23 | Artery_Aorta | 3.13e-93 |
| PHACTR1 | 3 | Artery_Tibial | 1.41e-42 |
| LPA | 1 | Testis | 4.41e-37 |
| CDKN2A | 6 | Brain_Cerebellum | 3.92e-27 |
| PSRC1 | 16 | Whole_Blood | 4.63e-23 |
| WDR12 | 14 | Brain_Nucleus_accumbens_basal_ganglia | 7.78e-19 |
| NBEAL1 | 19 | Artery_Tibial | 4.38e-18 |
| ICA1L | 35 | Kidney_Cortex | 1.44e-17 |
| AP000318.2 | 1 | Artery_Aorta | 7.58e-16 |
| MRPS6 | 1 | Artery_Aorta | 1.27e-15 |
| SLC5A3 | 2 | Artery_Aorta | 1.27e-15 |
| MORF4L1 | 6 | Artery_Aorta | 1.84e-15 |
| ADAMTS7 | 8 | Brain_Cerebellum | 3.68e-15 |
| LDLR | 3 | Minor_Salivary_Gland | 2.88e-14 |
| LIPA | 16 | Whole_Blood | 5.24e-14 |
| LINC00310 | 1 | Artery_Aorta | 5.63e-14 |
| TCF21 | 18 | Small_Intestine_Terminal_Ileum | 5.71e-14 |
| PLPP3 | 3 | Brain_Cerebellum | 9.83e-14 |
| AIDA | 10 | Brain_Hypothalamus | 7.47e-12 |
| MIA3 | 13 | Colon_Sigmoid | 9.60e-12 |
| CARF | 38 | Esophagus_Gastroesophageal_Junction | 1.19e-11 |
| KCNE2 | 2 | Artery_Aorta | 1.47e-11 |
| CELSR2 | 5 | Skin_Sun_Exposed_Lower_leg | 3.35e-11 |
| FAM177B | 7 | Small_Intestine_Terminal_Ileum | 6.78e-11 |
| ABO | 2 | Brain_Spinal_cord_cervical_c-1 | 7.14e-11 |
| ATP2B1 | 1 | Artery_Tibial | 8.53e-11 |
| GGCX | 36 | Cells_EBV-transformed_lymphocytes | 2.18e-10 |
| SMARCA4 | 2 | Esophagus_Gastroesophageal_Junction | 4.19e-10 |
| KIAA1462 | 5 | Heart_Left_Ventricle | 4.38e-10 |
| ZC3HC1 | 2 | Heart_Atrial_Appendage | 7.46e-10 |

## eqtl: S-MultiXcan (joint across tissues)

- 22,014 genes tested, Bonferroni p < 2.27e-06; significant: 70

| gene | p | tissues used | best single-tissue p |
|---|---|---|---|
| CDKN2B | 1.24e-39 | 29 | 3.13e-93 |
| PHACTR1 | 2.01e-36 | 35 | 1.41e-42 |
| LPA | 8.02e-36 | 9 | 4.41e-37 |
| CDKN2A | 1.20e-31 | 48 | 3.92e-27 |
| PSRC1 | 1.72e-16 | 47 | 4.63e-23 |
| LDLR | 4.46e-16 | 3 | 2.88e-14 |
| MORF4L1 | 9.46e-16 | 12 | 1.84e-15 |
| ICA1L | 1.03e-15 | 46 | 1.44e-17 |
| AP000318.2 | 1.82e-15 | 5 | 7.58e-16 |
| WDR12 | 7.51e-15 | 34 | 7.78e-19 |
| TGFB1 | 3.18e-14 | 34 | 7.52e-09 |
| EDNRA | 3.23e-14 | 30 | 3.13e-09 |
| NBEAL1 | 3.84e-14 | 33 | 4.38e-18 |
| TCF21 | 4.35e-14 | 27 | 5.71e-14 |
| SLC5A3 | 2.01e-13 | 19 | 1.27e-15 |
| ADAMTS7 | 2.45e-13 | 33 | 3.68e-15 |
| LINC00310 | 6.34e-13 | 39 | 5.63e-14 |
| MRPS6 | 1.61e-12 | 13 | 1.27e-15 |
| GUCY1A3 | 1.82e-11 | 37 | 7.24e-09 |
| PLPP3 | 4.65e-11 | 24 | 9.83e-14 |
| KCNE2 | 5.26e-11 | 43 | 1.47e-11 |
| LIPA | 6.33e-10 | 40 | 5.24e-14 |
| ABCG8 | 7.83e-10 | 18 | 8.41e-07 |
| AS3MT | 8.91e-10 | 49 | 6.42e-07 |
| CCDC97 | 2.55e-09 | 48 | 3.60e-08 |
| APOE | 2.75e-09 | 42 | 9.44e-08 |
| TOMM40 | 2.92e-09 | 17 | 6.83e-08 |
| MIA3 | 3.33e-09 | 48 | 9.60e-12 |
| GGCX | 3.42e-09 | 49 | 2.18e-10 |
| CARF | 3.75e-09 | 46 | 1.19e-11 |

