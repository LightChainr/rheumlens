# Supplementary table index

The stable workbook `Supplementary_Tables_v2.0.0-rc3.xlsx` contains the following
worksheets. Every analytical sheet is a direct, formatted copy or compact derivation
of a versioned TSV in `results/`, `inputs/` or `identifiability_extension/results/`.

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
| S15_SimSummary | Simulation summary across design-label dependence and signal settings | `identifiability_extension/results/simulation/simulation_summary.tsv` |
| S16_DesignNulls | Observed and complete-pipeline null results for GSE135779 design models | `identifiability_extension/results/gse135779_null/design_null_summary.tsv` |
| S17_DesignNullDraws | Full label-permutation and design-column-shuffle null draws | `identifiability_extension/results/gse135779_null/design_null_distributions.tsv.gz` |
| S18_LabelInfo | Design-adjusted label-information fraction in GSE135779 | `identifiability_extension/results/gse135779_null/design_information_fraction.tsv` |
| S19_CD4Composition | Repeated-CV disease performance from GSE174188 CD4-subtype composition | `identifiability_extension/results/composition/composition_summary.tsv` |
| S20_CD4Restriction | Wave-restricted and size-matched CD4-composition analyses | `identifiability_extension/results/composition/composition_restriction_metrics.tsv` |
| S21_CD4WavePredict | Predictability of study wave from CD4-composition features | `identifiability_extension/results/composition/composition_wave_predictability.tsv` |
| S22_DesignOverlap | Discrete support and design-label overlap diagnostics in GSE135779 | `identifiability_extension/results/gse135779_null/discrete_overlap_summary.tsv` |
| S23_RestrictionTests | Sign-flip tests for residualisation and matched-restriction contrasts | `identifiability_extension/results/gse135779_null/restriction_sign_flip_tests.tsv` |
| S24_AdjustSummary | Summary of linear, nonlinear, batch-location and overlap-weighted adjustment analyses | `results/plos_robustness/residualisation_sensitivity_summary.tsv` |
| S25_AdjustRepeats | Repeat-level metrics for adjustment sensitivity analyses | `results/plos_robustness/residualisation_sensitivity_metrics.tsv` |
| S26_NonlinearDesign | Nonlinear and linear design-only disease prediction | `results/plos_robustness/nonlinear_design_only_summary.tsv` |
| S27_BatchNegative | Internal batch-exposure negative control | `results/plos_robustness/batch_exposure_negative_control.tsv` |
| S28_AttenuationCI | Bootstrap sensitivity intervals for residualisation-versus-restriction attenuation differences | `results/plos_robustness/attenuation_difference_bootstrap.tsv` |
| S29_PositiveControl | Synthetic positive-control results | `results/plos_robustness/synthetic_positive_control.tsv` |
| S30_LeakageSensitivity | Confound-leakage simulation summary | `results/plos_robustness/confound_leakage_summary.tsv` |

Intervals labelled as repeat percentiles describe split sensitivity. Transfer
intervals are percentile intervals from 5,000 case/control-stratified target-donor
bootstrap resamples and are conditional on one fitted source model.
