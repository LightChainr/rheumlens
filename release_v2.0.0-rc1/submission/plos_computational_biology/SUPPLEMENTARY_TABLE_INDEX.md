# Supplementary table index

The stable workbook `Supplementary_Tables_v2.0.0-rc1.xlsx` contains the following
worksheets. Every analytical sheet is a direct, formatted copy or compact derivation
of a versioned TSV in `results/` or `inputs/`.

| Sheet | Contents | Primary source |
|---|---|---|
| README | Workbook scope, notation and provenance | Release metadata |
| S1_Cohorts | Cohort sizes, class counts, analysis compartment and primary cell cap | Donor labels and methods lock |
| S2_CoreSummary | Locked repeated-CV and residualisation summary | `results/locked_validity_audit/summary.tsv` |
| S3_RepeatMetrics | Repeat-level AUC, average precision, Brier score and ECE | `results/locked_validity_audit/repeat_metrics.tsv` |
| S4_BatchPredict | Representation-to-batch predictability summary | `results/locked_validity_audit/representation_batch_predictability_summary.tsv` |
| S5_Strata | Headline restriction and matched-control results | `results/HEADLINE_stratification_table.tsv`, `results/pure_wave_stratum.tsv` |
| S6_GSE135Matched | Repeat-level mixed-batch restriction analysis | `results/locked_validity_audit/gse135779_mixed_batch_matched.tsv` |
| S7_Transfer | Strict source-only cross-cohort metrics | `results/strict_source_only_transfer/strict_source_only_transfer_metrics.tsv` |
| S8_PairedDeLong | Paired target-donor AUC comparisons | `results/strict_source_only_transfer/strict_source_only_transfer_paired_delong.tsv` |
| S9_Harmonization | Source-only feature matching, scaling and PCA provenance | `results/strict_source_only_transfer/strict_source_only_feature_harmonization.tsv` |
| S10_LearnedPooling | Independent-target learned-pooling metrics | `results/learned_pooling/learned_pooling_metrics.tsv` |
| S11_LearnedTests | Paired learned-pooling comparisons | `results/learned_pooling/learned_pooling_paired_tests.tsv` |
| S12_GSE135Metadata | Restored GSE135779 donor metadata | `results/gse135779_metadata/gse135779_donor_metadata_restored.tsv` |
| S13_CovariateAudit | Repeated covariate and residualisation summary | `results/covariate_audit/summary.tsv` |
| S14_CovariateRepeats | Repeat-level covariate audit metrics | `results/covariate_audit/repeat_level_metrics.tsv` |

Intervals labelled as repeat percentiles describe split sensitivity. Transfer
intervals are percentile intervals from 5,000 case/control-stratified target-donor
bootstrap resamples and are conditional on one fitted source model.
