# RheumLens: a donor-level audit benchmark for single-cell foundation-model disease representations in systemic lupus erythematosus

**Running title:** Donor-level audit of single-cell foundation models in SLE  
**Manuscript version:** working v0.6, assembled 2026-06-23  
**Authors:** [AUTHOR NAMES TO BE ADDED]  
**Affiliations:** [AFFILIATIONS TO BE ADDED]  
**Corresponding author:** [NAME, EMAIL TO BE ADDED]

## Abstract

Single-cell foundation models are increasingly used as reusable feature extractors for disease studies. Their donor-level value remains difficult to quantify because many evaluations mix cells from the same donor across train/test splits, omit strong expression baselines, or report high discrimination without testing covariate sensitivity. We built RheumLens, a donor-level benchmark and audit workflow for systemic lupus erythematosus (SLE) single-cell RNA-seq cohorts. The workflow uses donor-disjoint folds, out-of-fold prediction, source-only transfer, paired resampling, representation provenance checks, cell-count sensitivity, covariate sensitivity, and kernel diagnostics.

Across GSE135779 matched-500, GSE174188 CD4, and GSE285773 CD4 analyses, disease signal was strong and reproducible. In GSE135779, raw pseudobulk reached AUC 0.972, focus_lite@scGPT reached 0.967, and Geneformer mean reached 0.926 after full provenance repair. In GSE174188 CD4, raw pseudobulk reached AUC 0.988, donor expression PCA 0.986, donor mean HVG 0.984, scGPT mean 0.978, and kme_multiscale@scGPT 0.978. Source-only transfer between GSE285773 and GSE174188 showed cross-study signal: GSE285773→GSE174188 AUCs were 0.730 for scGPT mean, 0.790 for donor expression PCA, and 0.806 for donor mean HVG; GSE174188→GSE285773 AUCs were 0.881, 0.931, and 0.900, respectively. Covariate sensitivity was substantial: measured covariates alone achieved AUC 0.846 in GSE174188 CD4, while residualized scGPT mean, donor expression PCA, and donor mean HVG reached AUCs of 0.868, 0.853, and 0.870.

RheumLens shows that single-cell foundation embeddings recover reproducible SLE-associated donor signal, while expression baselines remain highly competitive. The strongest manuscript finding is the audit itself: high disease AUCs require donor-level validation, expression baselines, provenance controls, and covariate sensitivity before biological or model-superiority claims are made.

**Keywords:** single-cell RNA-seq; foundation model; scGPT; Geneformer; systemic lupus erythematosus; donor-level classification; benchmark; covariate sensitivity

## 1. Introduction

Single-cell foundation models convert transcriptomes into compact representations intended to transfer across tissues, cell types, and tasks. Models such as scGPT and Geneformer are now practical tools for cell-level embedding extraction, and their use in disease datasets is expanding. The disease setting creates a specific statistical problem: thousands of cells are measured, but the independent unit is usually the donor. Evaluations that split cells instead of donors can report inflated performance. Evaluations without conventional expression baselines can also mistake recoverable expression signal for foundation-model-specific value.

Systemic lupus erythematosus provides a demanding test case. SLE has a strong, reproducible interferon-associated transcriptional program and multiple public single-cell cohorts with case-control donors. These properties make it possible to test whether frozen embeddings recover disease signal, whether that signal transfers across studies, and whether more complex donor-level methods add information beyond straightforward expression summaries.

RheumLens was developed as a donor-level benchmark and audit framework for this setting. The project evaluates scGPT, Geneformer, expression pseudobulk, donor expression PCA, HVG baselines, KME-style distribution embeddings, FOCUS-style methods, attention/MIL-style models, and covariate-adjusted sensitivity analyses. The final analysis emphasizes out-of-fold donor predictions, fixed folds, source-only transfer, provenance checks, and explicit claim boundaries.

## 2. Results

### 2.1 RheumLens assembled a donor-level benchmark spanning matched cells, large adult CD4 data, and pediatric CD4 transfer

The final local analysis package contains the core processed data, embeddings, fold files, out-of-fold predictions, bootstrap summaries, and audit reports needed for manuscript assembly. GSE135779 matched-500 provides a balanced, controlled setting with 44 donors and 500 cells per donor. GSE174188 CD4 provides the largest adult CD4 case-control analysis with 261 donors. GSE285773 provides an independent pediatric CD4 cohort with 26 donors. All formal analyses use donor-level labels and donor-disjoint evaluation.

Figure 1 summarizes the leading benchmark results. The matched-500 and GSE174188 CD4 analyses both show high donor-level disease discrimination across multiple method families.

![Figure 1. Top benchmark AUCs.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_1_TOP_BENCHMARK_AUC.png)

**Table 1. Top benchmark results.**

| analysis              | method_id                      |   n_donors |   roc_auc |   pr_auc |   brier |
|:----------------------|:-------------------------------|-----------:|----------:|---------:|--------:|
| GSE135779 matched-500 | raw_pseudobulk                 |         44 |     0.972 |    0.991 |   0.192 |
| GSE135779 matched-500 | focus_lite@scgpt               |         44 |     0.967 |    0.99  |   0.088 |
| GSE135779 matched-500 | attention_mil@geneformer       |         44 |     0.964 |    0.99  |   0.089 |
| GSE135779 matched-500 | red@geneformer                 |         44 |     0.964 |    0.989 |   0.106 |
| GSE135779 matched-500 | tail_fractions@geneformer      |         44 |     0.961 |    0.987 |   0.102 |
| GSE135779 matched-500 | focus_lite@geneformer          |         44 |     0.961 |    0.987 |   0.104 |
| GSE135779 matched-500 | red@scgpt                      |         44 |     0.956 |    0.987 |   0.092 |
| GSE135779 matched-500 | prototype_histogram@geneformer |         44 |     0.942 |    0.981 |   0.104 |
| GSE135779 matched-500 | moments_mean_var@geneformer    |         44 |     0.942 |    0.981 |   0.135 |
| GSE135779 matched-500 | moments_mean_var@scgpt         |         44 |     0.926 |    0.977 |   0.143 |
| GSE135779 matched-500 | geneformer_mean                |         44 |     0.926 |    0.974 |   0.157 |
| GSE135779 matched-500 | quantiles@geneformer           |         44 |     0.917 |    0.974 |   0.161 |
| GSE174188 CD4         | raw_pseudobulk                 |        261 |     0.988 |    0.993 |   0.054 |
| GSE174188 CD4         | donor_expression_pca           |        261 |     0.986 |    0.992 |   0.054 |
| GSE174188 CD4         | donor_mean_hvg                 |        261 |     0.984 |    0.991 |   0.055 |
| GSE174188 CD4         | scgpt_mean                     |        261 |     0.978 |    0.988 |   0.062 |
| GSE174188 CD4         | kme_multiscale@scgpt           |        261 |     0.978 |    0.986 |   0.056 |
| GSE174188 CD4         | attention_mil                  |        261 |     0.97  |    0.984 |   0.064 |
| GSE174188 CD4         | cckme_u_weighted               |        261 |     0.97  |    0.979 |   0.065 |
| GSE174188 CD4         | cckme_u@scgpt                  |        261 |     0.97  |    0.978 |   0.064 |
| GSE174188 CD4         | cckme_u                        |        261 |     0.97  |    0.978 |   0.064 |
| GSE174188 CD4         | topk_mil                       |        261 |     0.969 |    0.982 |   0.066 |
| GSE174188 CD4         | moments_mean_var@scgpt         |        261 |     0.967 |    0.965 |   0.066 |
| GSE174188 CD4         | tail_fractions@scgpt           |        261 |     0.962 |    0.968 |   0.064 |

### 2.2 Foundation embeddings recover SLE-associated donor signal, while expression baselines remain highly competitive

In GSE135779 matched-500, raw pseudobulk achieved AUC 0.972. scGPT-derived methods were also strong: focus_lite@scGPT reached 0.967 and scGPT mean reached 0.846. After Geneformer provenance repair, Geneformer mean reached AUC 0.926 and focus_lite@geneformer reached 0.961. This corrected the earlier failure mode caused by an unknown-provenance Geneformer NPZ and incompatible extraction environment.

In GSE174188 CD4, raw pseudobulk achieved AUC 0.988, donor expression PCA 0.986, donor mean HVG 0.984, scGPT mean 0.978, and kme_multiscale@scGPT 0.978. These results show robust disease signal in frozen embeddings and direct expression summaries. Expression baselines were among the strongest methods in the final analysis.

**Table 2. Primary method anchors.**

| analysis              | method_id                |   n_donors |   roc_auc |   pr_auc |   brier |
|:----------------------|:-------------------------|-----------:|----------:|---------:|--------:|
| GSE135779 matched-500 | raw_pseudobulk           |         44 |     0.972 |    0.991 |   0.192 |
| GSE135779 matched-500 | focus_lite@scgpt         |         44 |     0.967 |    0.99  |   0.088 |
| GSE135779 matched-500 | attention_mil@geneformer |         44 |     0.964 |    0.99  |   0.089 |
| GSE135779 matched-500 | geneformer_mean          |         44 |     0.926 |    0.974 |   0.157 |
| GSE135779 matched-500 | donor_expression_pca     |         44 |     0.915 |    0.972 |   0.112 |
| GSE135779 matched-500 | donor_mean_hvg           |         44 |     0.893 |    0.955 |   0.181 |
| GSE135779 matched-500 | kme_multiscale@scgpt     |         44 |     0.873 |    0.958 |   0.198 |
| GSE135779 matched-500 | scgpt_mean               |         44 |     0.846 |    0.951 |   0.164 |
| GSE174188 CD4         | raw_pseudobulk           |        261 |     0.988 |    0.993 |   0.054 |
| GSE174188 CD4         | donor_expression_pca     |        261 |     0.986 |    0.992 |   0.054 |
| GSE174188 CD4         | donor_mean_hvg           |        261 |     0.984 |    0.991 |   0.055 |
| GSE174188 CD4         | scgpt_mean               |        261 |     0.978 |    0.988 |   0.062 |
| GSE174188 CD4         | kme_multiscale@scgpt     |        261 |     0.978 |    0.986 |   0.056 |
| GSE174188 CD4         | attention_mil            |        261 |     0.97  |    0.984 |   0.064 |
| GSE174188 CD4         | focus_lite@scgpt         |        261 |     0.958 |    0.974 |   0.08  |

### 2.3 Source-only CD4 transfer supports cross-study signal

We repaired the local GSE174188 CD4 expression matrix metadata by recovering true gene symbols from the server-side h5ad file. After this repair, GSE174188 and GSE285773 shared 18,916 expression features, enabling expression-space transfer as well as scGPT transfer.

In source-only P9 transfer, training on GSE285773 and testing on GSE174188 CD4 produced AUC 0.730 for scGPT mean, 0.790 for donor expression PCA, and 0.806 for donor mean HVG. The reverse direction, training on GSE174188 CD4 and testing on GSE285773, produced AUC 0.881 for scGPT mean, 0.931 for donor expression PCA, and 0.900 for donor mean HVG. The reverse target has only 26 donors, so intervals are wide, but the signal is clear.

![Figure 2. P9 source-only CD4 transfer.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_2_P9_TRANSFER_AUC.png)

**Table 3. Source-only transfer results.**

| direction                  | method_id            |   n_donors |   n_cases |   n_controls |   auc |   auc_ci_low |   auc_ci_high |   pr_auc |   pr_auc_ci_low |   pr_auc_ci_high |   brier |
|:---------------------------|:---------------------|-----------:|----------:|-------------:|------:|-------------:|--------------:|---------:|----------------:|-----------------:|--------:|
| GSE285773_to_GSE174188_CD4 | scgpt_mean           |        261 |       162 |           99 | 0.73  |        0.67  |         0.791 |    0.845 |           0.808 |            0.88  |   0.443 |
| GSE285773_to_GSE174188_CD4 | donor_expression_pca |        261 |       162 |           99 | 0.79  |        0.737 |         0.844 |    0.896 |           0.868 |            0.924 |   0.252 |
| GSE285773_to_GSE174188_CD4 | donor_mean_hvg       |        261 |       162 |           99 | 0.806 |        0.755 |         0.858 |    0.893 |           0.862 |            0.922 |   0.313 |
| GSE174188_CD4_to_GSE285773 | scgpt_mean           |         26 |        16 |           10 | 0.881 |        0.75  |         1     |    0.935 |           0.856 |            0.99  |   0.615 |
| GSE174188_CD4_to_GSE285773 | donor_expression_pca |         26 |        16 |           10 | 0.931 |        0.833 |         1     |    0.969 |           0.917 |            1     |   0.183 |
| GSE174188_CD4_to_GSE285773 | donor_mean_hvg       |         26 |        16 |           10 | 0.9   |        0.78  |         1     |    0.953 |           0.887 |            1     |   0.476 |

### 2.4 Repeated cross-validation and permutation testing confirmed that the primary signals are stable

Repeated cross-validation on GSE174188 CD4 produced tightly concentrated AUC distributions for the three primary methods. Across 30 repeats, mean AUC was 0.977 for scGPT mean, 0.986 for donor expression PCA, and 0.975 for kme_multiscale@scGPT. These values agreed with the fixed-fold P8.2 reference runs.

Formal donor-label permutation tests were available as full null distributions for donor expression PCA and kme_multiscale@scGPT. Both observed AUCs were far outside the 1,000-replicate null distributions, with empirical P=0.000999. This establishes that the primary GSE174188 CD4 signals are not artifacts of the fixed fold assignment.

![Figure 4. Repeated CV stability.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_4_REPEATED_CV_STABILITY.png)

![Figure 5. Formal permutation nulls.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_5_PERMUTATION_NULLS.png)

### 2.5 Cell-count sensitivity confirms that the primary signal is stable across donor sampling levels

P8.4 tested common-support sampling across cell-count levels. The analysis used nested sampling manifests and a corrected V5 runner that invoked the registered method objects rather than collapsing all methods into mean pooling. At common-support level 50, the 10-repeat mean AUC was 0.953 for scGPT mean, 0.970 for donor expression PCA, and 0.948 for kme_multiscale@scGPT. Across levels 25, 50, 100, 200, and 500, disease discrimination remained high. The result supports stability of the donor-level signal under cell-count perturbation.

![Figure 6. Cell-count sensitivity.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_6_CELL_COUNT_SENSITIVITY.png)

**Table 4. Cell-count sensitivity summary.**

| estimand        |   level | method               |   n_repeats |   mean_auc |   sd_auc |   min_auc |   max_auc |   mean_pr_auc |
|:----------------|--------:|:---------------------|------------:|-----------:|---------:|----------:|----------:|--------------:|
| available_donor |      25 | donor_expression_pca |          10 |      0.965 |    0.002 |     0.962 |     0.969 |         0.979 |
| available_donor |      25 | kme_multiscale@scgpt |          10 |      0.937 |    0.008 |     0.924 |     0.947 |         0.962 |
| available_donor |      25 | scgpt_mean           |          10 |      0.944 |    0.008 |     0.933 |     0.958 |         0.967 |
| available_donor |      50 | donor_expression_pca |          10 |      0.974 |    0.003 |     0.969 |     0.978 |         0.984 |
| available_donor |      50 | kme_multiscale@scgpt |          10 |      0.95  |    0.01  |     0.93  |     0.965 |         0.969 |
| available_donor |      50 | scgpt_mean           |          10 |      0.958 |    0.008 |     0.946 |     0.972 |         0.975 |
| available_donor |     100 | donor_expression_pca |          10 |      0.978 |    0.004 |     0.971 |     0.983 |         0.987 |
| available_donor |     100 | kme_multiscale@scgpt |          10 |      0.961 |    0.01  |     0.944 |     0.977 |         0.976 |
| available_donor |     100 | scgpt_mean           |          10 |      0.966 |    0.006 |     0.958 |     0.977 |         0.98  |
| available_donor |     200 | donor_expression_pca |          10 |      0.983 |    0.004 |     0.978 |     0.988 |         0.99  |
| available_donor |     200 | kme_multiscale@scgpt |          10 |      0.967 |    0.007 |     0.959 |     0.98  |         0.979 |
| available_donor |     200 | scgpt_mean           |          10 |      0.97  |    0.004 |     0.964 |     0.977 |         0.982 |
| available_donor |     500 | donor_expression_pca |          10 |      0.983 |    0.003 |     0.979 |     0.986 |         0.989 |
| available_donor |     500 | kme_multiscale@scgpt |          10 |      0.972 |    0.004 |     0.965 |     0.979 |         0.982 |
| available_donor |     500 | scgpt_mean           |          10 |      0.974 |    0.002 |     0.971 |     0.976 |         0.984 |
| common_support  |      25 | donor_expression_pca |          10 |      0.96  |    0.003 |     0.956 |     0.965 |         0.974 |
| common_support  |      25 | kme_multiscale@scgpt |          10 |      0.936 |    0.01  |     0.922 |     0.949 |         0.957 |
| common_support  |      25 | scgpt_mean           |          10 |      0.937 |    0.01  |     0.926 |     0.955 |         0.96  |
| common_support  |      50 | donor_expression_pca |          10 |      0.97  |    0.005 |     0.961 |     0.978 |         0.98  |
| common_support  |      50 | kme_multiscale@scgpt |          10 |      0.948 |    0.011 |     0.927 |     0.964 |         0.966 |

### 2.6 Covariate sensitivity materially changes interpretation of high AUCs

The covariate audit is central to the final manuscript. In GSE174188 CD4, measured donor/sample covariates alone achieved AUC 0.846. After residualization against the frozen covariate set, scGPT mean reached AUC 0.868, donor expression PCA 0.853, and donor mean HVG 0.870. Residualized methods were not cleanly separated from the covariates-only baseline by the paired bootstrap analysis. This result does not eliminate biological disease signal; it shows that measured donor/sample variables capture a substantial component of the apparent discrimination and must be reported with the benchmark.

![Figure 3. Covariate sensitivity.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_3_COVARIATE_SENSITIVITY.png)

Residualized method scores were only modestly above the covariates-only baseline. The AUC increment was +0.022 for scGPT mean, +0.007 for donor expression PCA, and +0.023 for donor mean HVG; the corresponding bootstrap intervals included zero. This result motivates presenting covariate sensitivity as a primary audit endpoint rather than a supplementary afterthought.

![Figure 7. Residualized methods versus covariates-only.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_7_COVARIATE_DELTA.png)

**Table 5. Covariate sensitivity.**

| method_id                   |   n_donors |   n_sle |   n_control |   auc |   pr_auc |   brier |   feature_dim |   q025 |   q975 |
|:----------------------------|-----------:|--------:|------------:|------:|---------:|--------:|--------------:|-------:|-------:|
| covariates_only             |        261 |     162 |          99 | 0.846 |    0.905 |   0.162 |            13 |  0.798 |  0.892 |
| scgpt_mean_covres           |        261 |     162 |          99 | 0.868 |    0.896 |   0.148 |            13 |  0.819 |  0.911 |
| donor_expression_pca_covres |        261 |     162 |          99 | 0.853 |    0.87  |   0.156 |            13 |  0.801 |  0.899 |
| donor_mean_hvg_covres       |        261 |     162 |          99 | 0.87  |    0.896 |   0.169 |            13 |  0.822 |  0.913 |

### 2.7 KME diagnostics resolved identity collapse and established a bounded interpretation

P8.7 confirmed that KME predictions are distinct from scGPT mean predictions after the V5 runner fix. Across cell-count levels, KME had zero exact prediction identity with scGPT mean, stable bandwidth values around 7.55–7.65, and effective rank approximately 2.6–2.8. The low effective rank indicates a concentrated kernel spectrum. The diagnostic result supports use of KME as a distinct distribution-level method in the benchmark and sets the interpretation boundary for mechanism or superiority claims.

![Figure 8. KME diagnostics.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_8_KME_DIAGNOSTICS.png)

**Table 6. KME identity diagnostics.**

|   level |   max_abs_diff |   pearson_r |   spearman_r |   n_exact_equal |
|--------:|---------------:|------------:|-------------:|----------------:|
|      25 |          0.661 |       0.942 |        0.956 |               0 |
|      50 |          0.739 |       0.949 |        0.967 |               0 |
|     100 |          0.574 |       0.968 |        0.976 |               0 |
|     200 |          0.614 |       0.963 |        0.974 |               0 |
|     500 |          0.69  |       0.957 |        0.976 |               0 |

### 2.8 Three-cohort method heatmap summarizes the final benchmark landscape

The three-cohort method heatmap provides a compact view of method behavior across GSE135779, GSE174188 CD4, and GSE285773. Expression methods dominate the upper range in the larger CD4 analyses, while foundation-embedding methods remain competitive and provide a reusable representation layer for distribution-level and attention-style aggregators. This pattern supports the manuscript's main empirical message: disease signal is strong and reproducible, and the audit framework is required to determine when complex representations add value beyond expression.

![Figure 9. Three-cohort method heatmap.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_9_THREE_COHORT_METHOD_HEATMAP.png)

### 2.9 Interferon and high-discrimination genes anchor the benchmark in known SLE biology

To add a direct biological anchor to the benchmark, we computed donor-level CD4 pseudobulk expression in GSE174188 after repairing gene symbols from the source h5ad file. A 20-gene interferon-stimulated gene (ISG) score was higher in SLE donors than controls and reached donor-level AUC 0.906, Cohen d=1.37 and Welch P=3.6×10^-29. All 20 predefined ISG genes were present in the expression matrix.

The top SLE-high donor-level genes included EPSTI1, ISG15, IFI44, IFI44L, IFI6, XAF1, IRF7 and USP18, consistent with a strong type I interferon program. Additional high-ranking genes included PSME2, LY6E, TNFRSF4 and LGALS3BP. This result supports the clinical interpretation that much of the benchmarked signal tracks a canonical lupus inflammatory axis rather than an opaque model-specific latent factor.

![Figure 10. ISG signature validation.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_10_ISG_SIGNATURE_VALIDATION.png)

![Figure 11. Top SLE-high donor-level genes.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_11_TOP_SLE_HIGH_GENES.png)

![Figure 12. Top single-gene donor discriminators.](/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/FIGURE_12_TOP_SINGLE_GENE_AUC.png)

The mechanism analysis also clarifies the clinical meaning of the covariate audit. The adult GSE174188 CD4 dataset contains strong SLE/control transcriptional separation, but donor/sample covariates alone also classify well. These two findings are compatible: clinically meaningful disease biology and cohort-structured covariates can both contribute to high AUC. RheumLens therefore treats mechanism validation and covariate sensitivity as paired requirements for interpreting single-cell disease classifiers.

### 2.10 FOCUS and surrogate-null audits were handled by trigger logic

FOCUS results exist in P5 and P6. The earlier broad superiority phrasing was withdrawn during audit. P8.8 was therefore treated as a sensitivity-control trigger audit rather than a formal superiority test. P8.5 surrogate-null computation was skipped because KME did not show a positive incremental advantage over scGPT mean in the audited analyses. These skip decisions are recorded in the failed/skipped registry and preserve a clear boundary between completed formal evidence and exploratory method code.

## 3. Discussion

RheumLens establishes a practical donor-level audit standard for single-cell foundation-model disease analysis. The benchmark confirms that frozen scGPT and repaired Geneformer embeddings carry SLE-associated donor signal. It also shows that expression baselines are very strong. Raw pseudobulk, donor expression PCA, and donor mean HVG frequently match or exceed embedding-based methods in AUC.

The source-only transfer results strengthen the manuscript. They show that the donor-level signal is not restricted to within-cohort cross-validation. Transfer between pediatric and adult CD4 datasets remains informative, with expression summaries strongest in the current formal run and scGPT mean retaining clear cross-study signal.

The covariate analysis should anchor the interpretation. High AUCs in GSE174188 CD4 are partly recoverable from measured donor/sample covariates. This finding raises the technical standard for future claims. Disease-prediction performance in single-cell cohorts should be reported alongside covariate-only baselines and residualized sensitivity analyses.

The mechanism supplement gives the benchmark a biological axis. The strongest donor-level CD4 signal in GSE174188 is an interferon-centered lupus program, with EPSTI1, ISG15, IFI44, IFI44L, IFI6, XAF1, IRF7 and USP18 among the highest SLE-associated genes. This is clinically plausible because interferon activation is a recurrent SLE feature and is relevant to disease stratification and therapeutic interpretation. The pediatric-to-adult transfer results suggest that part of this transcriptional program is stable across age groups, while the covariate audit shows that clinical use would require explicit control for donor, demographic and sample-level structure.

The Geneformer repair and KME runner audit are also important. Geneformer initially appeared defective because an unknown-provenance embedding had poor disease performance. Full extraction repair restored expected performance. KME initially collapsed to mean-pooling behavior because the runner ignored registered method aggregators. The corrected V5 pipeline resolved this. These episodes support one of the paper's main engineering claims: model-comparison papers need provenance gates and implementation audits, not only output metrics.

The main limitation is that formal claims remain retrospective and cohort-specific. GSE285773 has only 26 donors, and GSE174188 covariates strongly affect interpretation. The benchmark evaluates frozen or fixed representations rather than full task-specific fine-tuning. Clinical deployment, longitudinal prediction, treatment response, and causal biology require additional datasets and prospective designs.

## 4. Methods

### 4.1 Data sources and processed objects

The formal snapshot contains processed assets for GSE135779 matched-500, GSE174188 CD4, and GSE285773. GSE135779 matched-500 includes scGPT and Geneformer embeddings and expression matrices for 44 donors. GSE174188 CD4 includes 380,275 cells from 261 donors. GSE285773 includes 262,908 cells from 26 donors. GSE174188 feature names were repaired from the source h5ad `var.feature_name` field before P9 expression transfer.

### 4.2 Donor-level evaluation

All formal classification metrics are donor-level. Cell-level representations are aggregated to donor-level features by the registered method. Cross-validation uses donor-disjoint folds. Out-of-fold predictions are pooled across held-out donors before AUC, PR-AUC, and Brier score calculation. Source-only transfer fits all preprocessing and classifiers on the source cohort and applies the fitted pipeline to the target cohort.

### 4.3 Method families

The benchmark includes raw pseudobulk, donor expression PCA, donor mean HVG, scGPT mean, Geneformer mean, KME/multiscale kernel mean embeddings, FOCUS-style methods, RED/tail/moment summaries, attention/MIL methods, and additional learned donor-set models. The manuscript tables emphasize methods with completed formal audit status.

### 4.4 Uncertainty estimates

Within-cohort sensitivity analyses use paired donor bootstrap where available. P8.2 formal permutation used 1,000 donor-label permutations for primary GSE174188 CD4 methods and empirical P-values with one-count correction. P8.3 repeated cross-validation used 30 stratified donor-level repeats to quantify fold-resampling stability. P9 source-only transfer reports AUC intervals using a Hanley-McNeil large-sample approximation and PR-AUC intervals using stratified donor bootstrap. Bootstrap distributions, permutation distributions, repeated-CV summaries and out-of-fold predictions are archived with SHA256 manifests.

### 4.5 Covariate sensitivity

P8.6 used a frozen covariate set derived from authoritative donor metadata. The final V3.1 analysis includes covariates-only prediction and residualized versions of scGPT mean, donor expression PCA, and donor mean HVG. The sensitivity analysis is interpreted as an audit of robustness to measured covariates.

### 4.6 Mechanism-oriented gene analysis

For GSE174188 CD4, log-normalized cell-level expression was aggregated to donor-level pseudobulk means after gene-symbol repair from the source h5ad `var.feature_name` field. Donor-level SLE-control gene differences were summarized by mean difference, Welch t statistic, Benjamini-Hochberg FDR, Cohen d and single-gene donor AUC. Duplicate gene symbols were collapsed by retaining the most discriminative feature. A predefined 20-gene ISG score was computed as the donor-level mean of IFI44L, IFI44, IFI27, ISG15, MX1, MX2, OAS1, OAS2, OAS3, OASL, IFIT1, IFIT2, IFIT3, RSAD2, SIGLEC1, STAT1, IRF7, USP18, XAF1 and EPSTI1.

### 4.7 Reproducibility

The local final run directory is `/Volumes/Mac Data/Research/RheumLens_20260622/mac_final_results/run_20260623T045526Z`. The manuscript asset directory contains all plot-ready tables, figures, and SHA256 manifests. Server-derived GSE174188 feature-name repair is documented in `metadata/GSE174188_feature_names_from_h5ad.tsv`.

## 5. Data and code availability

Source datasets are publicly available from GEO or associated repositories. Processed embeddings, fold files, out-of-fold predictions, bootstrap summaries, audit reports, and plot-ready tables will be deposited in the project repository or archival storage before submission. Large cell-level embeddings should be distributed according to source-data terms and repository size limits.

## 6. Current manuscript figure legends

**Figure 1. Top benchmark AUCs across matched-500 and GSE174188 CD4 analyses.** Horizontal bars show leading method ROC-AUC values. The figure highlights high donor-level disease discrimination across expression and foundation-embedding methods.

**Figure 2. Source-only CD4 cross-study transfer.** Models were trained on one cohort and evaluated on the other without target-label fitting. Error bars show AUC intervals.

**Figure 3. Covariate sensitivity in GSE174188 CD4.** Covariates-only and residualized method AUCs are shown with bootstrap intervals. Measured covariates account for a substantial component of donor-level discrimination.

**Figure 4. Repeated donor-level cross-validation stability.** Bars show mean AUC across 30 repeated donor-level cross-validation runs; error bars show the standard deviation. Black points show the corresponding fixed-fold P8.2 reference values.

**Figure 5. Formal donor-label permutation null distributions.** Histograms show 1,000 donor-label permutation AUCs for donor expression PCA and kme_multiscale@scGPT. Red lines show observed AUCs. Dashed black lines show the 97.5th percentile of each null distribution.

**Figure 6. Cell-count sensitivity.** Common-support donor sampling was evaluated at 25, 50, 100, 200 and 500 cells per donor. Points show mean AUC across 10 repeats; error bars show 2.5th–97.5th percentile intervals.

**Figure 7. Residualized methods versus covariates-only.** Bars show AUC increment over the covariates-only baseline for residualized scGPT mean, donor expression PCA and donor mean HVG. Error bars show bootstrap intervals.

**Figure 8. KME diagnostics.** Left: bandwidth stability across cell-count levels. Right: Pearson and Spearman correlations between KME and scGPT mean donor predictions after identity-collapse repair.

**Figure 9. Three-cohort method heatmap.** Heatmap summarizing method AUCs across completed cohorts and modalities in the final P6/P7 integrated benchmark table.

**Figure 10. ISG signature validation.** Donor-level 20-gene interferon-stimulated gene score in GSE174188 CD4, stratified by control and SLE labels. The score achieved AUC 0.906, Cohen d=1.37 and Welch P=3.6×10^-29.

**Figure 11. Top SLE-high donor-level CD4 genes.** Bars show Cohen d for the top SLE-high genes from GSE174188 CD4 donor-level pseudobulk expression. Interferon-associated genes are highlighted.

**Figure 12. Top single-gene donor discriminators.** Bars show single-gene donor-level discriminative AUC for the top genes in GSE174188 CD4.

## 7. Supplementary material plan

- Supplementary Table S1: full method summary for P5 matched-500.
- Supplementary Table S2: full method summary for P6 GSE174188 CD4.
- Supplementary Table S3: integrated P6 all-cohort method summary.
- Supplementary Table S4: P8.3 repeated-CV summary.
- Supplementary Table S5: P8.4 cell-count sensitivity full table.
- Supplementary Table S6: P8.6 covariate summary.
- Supplementary Table S7: P8.6 residualized-method increment over covariates-only.
- Supplementary Table S8: P8.7 KME diagnostics.
- Supplementary Table S9: P9 source-only transfer summary.
- Supplementary Table S10: P9 source-only transfer donor-level predictions.
- Supplementary Table S11: P8.2 formal permutation summary.
- Supplementary Table S12: GSE174188 CD4 donor-level gene discrimination table.
- Supplementary Table S13: top SLE-high and control-high donor-level genes.
- Supplementary Table S14: donor-level ISG scores and ISG summary.
- Supplementary Methods: Geneformer provenance repair, KME V5 runner audit, FOCUS trigger audit, and skipped-analysis registry.

## 8. References to add before submission

- Baechler, E.C. et al. (2003). Interferon-inducible gene expression signature in peripheral blood cells of patients with severe lupus. *Proceedings of the National Academy of Sciences USA*, 100, 2610–2615.
- Bennett, L. et al. (2003). Interferon and granulopoiesis signatures in systemic lupus erythematosus blood. *Journal of Experimental Medicine*, 197, 711–723.
- Collins, G.S. et al. (2024). TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ*.
- Cui, H. et al. (2024). scGPT: toward building a foundation model for single-cell multi-omics using generative AI. *Nature Methods*, 21, 1470–1480.
- Hanley, J.A. and McNeil, B.J. (1982). The meaning and use of the area under a receiver operating characteristic curve. *Radiology*, 143, 29–36.
- Nehar-Belaid, D. et al. (2020). Mapping systemic lupus erythematosus heterogeneity at the single-cell level. *Nature Immunology*, 21, 1094–1106.
- Perez, R.K. et al. (2022). Single-cell RNA-seq reveals cell type-specific molecular and genetic associations to lupus. *Science*, 376, eabf1970.
- Squair, J.W. et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12, 5692.
- Theodoris, C.V. et al. (2023). Transfer learning enables predictions in network biology. *Nature*, 618, 616–624.
- GSE285773 GEO record: Single cell RNA profiling of blood CD4+ T cells identifies distinct helper and dysfunctional regulatory clusters in children with SLE [CD4+ T cells]. Full publication details to be verified before submission.
