# Figure and table index

Generated for the working v0.6 manuscript package.

## Main figures

| Figure | Asset stem | Primary content | Required format |
|---:|---|---|---|
| 1 | `FIGURE_1_TOP_BENCHMARK_AUC` | Top benchmark AUCs across matched-500 and GSE174188 CD4 | PNG + PDF |
| 2 | `FIGURE_2_P9_TRANSFER_AUC` | Source-only CD4 transfer summary | PNG + PDF |
| 3 | `FIGURE_3_COVARIATE_SENSITIVITY` | Covariate-only and residualized AUCs | PNG + PDF |
| 4 | `FIGURE_4_REPEATED_CV_STABILITY` | Repeated-CV stability | PNG + PDF |
| 5 | `FIGURE_5_PERMUTATION_NULLS` | Formal permutation nulls | PNG + PDF |
| 6 | `FIGURE_6_CELL_COUNT_SENSITIVITY` | Cell-count sensitivity | PNG + PDF |
| 7 | `FIGURE_7_COVARIATE_DELTA` | Residualized methods versus covariates-only | PNG + PDF |
| 8 | `FIGURE_8_KME_DIAGNOSTICS` | KME bandwidth and prediction-correlation diagnostics | PNG + PDF |
| 9 | `FIGURE_9_THREE_COHORT_METHOD_HEATMAP` | Three-cohort method AUC heatmap | PNG + PDF |
| 10 | `FIGURE_10_ISG_SIGNATURE_VALIDATION` | ISG donor-score validation | PNG + PDF |
| 11 | `FIGURE_11_TOP_SLE_HIGH_GENES` | Top SLE-high donor-level genes | PNG + PDF |
| 12 | `FIGURE_12_TOP_SINGLE_GENE_AUC` | Top single-gene donor discriminators | PNG + PDF |
| 13 | `FIGURE_13_COVARIATE_DECOMPOSITION` | GSE174188 covariate group decomposition | PNG + PDF |
| 14 | `FIGURE_14_PAIRED_AUC_DELTA` | Paired bootstrap ΔAUC with Holm correction | PNG + PDF |
| 15 | `FIGURE_15_MATCHED500_SCGPT_CONTEXT` | Matched-500 scGPT mean context | PNG + PDF |
| 16 | `FIGURE_16_P9_STRUCTURED_TRANSFER` | Structured source-only transfer extension | PNG + PDF |

## Main tables embedded in manuscript

| Table | Content |
|---:|---|
| 1 | Top benchmark results |
| 2 | Primary method summary |
| 3 | P9 source-only transfer summary |
| 4 | P8.4 cell-count sensitivity summary |
| 5 | Covariate sensitivity |
| 6 | KME identity diagnostics |

## Supplementary tables

| Table | File | Content |
|---|---|---|
| S1 | `SUPP_TABLE_S1_P5_MATCHED500_METHOD_SUMMARY.csv` | P5 matched-500 method summary |
| S2 | `SUPP_TABLE_S2_P6_GSE174188_METHOD_SUMMARY.csv` | P6 GSE174188 method summary |
| S3 | `SUPP_TABLE_S3_P6_ALL_COHORTS_METHOD_SUMMARY.csv` | P6 three-cohort method summary |
| S4 | `SUPP_TABLE_S4_P8_3_REPEATED_CV.csv` | P8.3 repeated-CV summary |
| S5 | `SUPP_TABLE_S5_P8_4_CELL_COUNT_FULL.csv` | P8.4 cell-count full results |
| S6 | `SUPP_TABLE_S6_P8_6_COVARIATE_SUMMARY.csv` | P8.6 covariate sensitivity summary |
| S7 | `SUPP_TABLE_S7_P8_6_SCORE_DELTA.csv` | P8.6 score-delta bootstrap |
| S8 | `SUPP_TABLE_S8_P8_7_KME_DIAGNOSTICS.csv` | P8.7 KME diagnostics |
| S9 | `SUPP_TABLE_S9_P9_TRANSFER_SUMMARY.csv` | P9 transfer summary |
| S10 | `SUPP_TABLE_S10_P9_TRANSFER_PREDICTIONS.csv` | P9 transfer donor predictions |
| S11 | `SUPP_TABLE_S11_P8_2_PERMUTATION_SUMMARY.csv` | P8.2 permutation summary |
| S12 | `SUPP_MECH_GSE174188_CD4_donor_level_gene_discrimination.csv` | Donor-level gene discrimination |
| S13 | `SUPP_MECH_TOP100_SLE_HIGH_GENES.csv`; `SUPP_MECH_TOP100_CONTROL_HIGH_GENES.csv`; `SUPP_MECH_TOP100_DISCRIMINATIVE_GENES.csv` | Top gene lists |
| S14 | `SUPP_MECH_GSE174188_CD4_donor_ISG_scores.csv` | Donor-level ISG scores |
| S15 | `SUPP_TABLE_S15_COVARIATE_DECOMPOSITION.csv` | Covariate decomposition and encoding sensitivity |
| S16 | `SUPP_TABLE_S16_PAIRED_AUC_DELTA_HOLM.csv` | Paired ΔAUC with Holm correction |
| S17 | `SUPP_TABLE_S17_MATCHED500_SCGPT_CONTEXT.csv` | Matched-500 scGPT mean context |
| S18 | `SUPP_TABLE_S18_P9_STRUCTURED_TRANSFER.csv` | P9 structured transfer extension |

## Figure legend requirements

Before journal submission, every figure legend should explicitly state:

- donor count and cohort;
- error-bar or interval definition;
- whether a value is fixed-fold, repeated-CV, bootstrap or transfer;
- expansion of method abbreviations at first use;
- whether the analysis is formal, sensitivity, or hypothesis-generating.

## Manifest

Figure checksums are in `figures/MANIFEST_SHA256.tsv`. Supplementary table checksums are in `supplementary_tables/MANIFEST_SHA256.tsv`.
