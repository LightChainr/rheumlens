# RheumLens manuscript extension analyses

Date: 2026-06-23

## 1. GSE174188 covariate decomposition

The formal-compatible covariate model achieved OOF AUC 0.848, reproducing accepted P8.6 V3.1. In this formal-compatible encoding, `processing_cohort` is numeric-coded by the accepted provider implementation and `sex` is one-hot categorical. A sensitivity encoding that one-hot encodes `processing_cohort` achieved AUC 0.930; this is reported as sensitivity only and is not substituted for the accepted P8.6 result.

Group-only models showed:

- technical depth/QC covariates: AUC 0.776;
- processing covariates: AUC 0.685;
- demographic covariates: AUC 0.653.

This supports the reviewer-facing point that the strong covariates-only signal is not a black box; it can be decomposed into technical/QC, processing, and demographic components. The result remains a sensitivity analysis and is not causal.

Top single-covariate AUCs:

| variable                   |      auc |   pr_auc |    brier |
|:---------------------------|---------:|---------:|---------:|
| cells_per_donor            | 0.689113 | 0.774567 | 0.224596 |
| age_from_development_stage | 0.639481 | 0.692074 | 0.240889 |
| mean_pct_mito              | 0.614042 | 0.753492 | 0.244491 |
| n_suspensions              | 0.568213 | 0.711666 | 0.231261 |
| mean_umi_per_cell          | 0.531425 | 0.615345 | 0.249677 |
| sex                        | 0.520358 | 0.639942 | 0.245819 |
| n_samples                  | 0.520264 | 0.656932 | 0.243121 |
| n_libraries                | 0.481419 | 0.601204 | 0.247415 |
| processing_cohort          | 0.459409 | 0.581697 | 0.250583 |
| mean_genes_per_cell        | 0.44644  | 0.554891 | 0.251915 |

Cross-cohort covariate extension status:

| cohort               |   n_donors |   n_cases |   n_controls | available_covariates                                                | missing_material_covariates                             | formal_adjustment_status   |
|:---------------------|-----------:|----------:|-------------:|:--------------------------------------------------------------------|:--------------------------------------------------------|:---------------------------|
| GSE135779_matched500 |         44 |        33 |           11 | age,age_group                                                       | batch,ancestry,treatment,site,chemistry,processing_date | limited_metadata_only      |
| GSE285773            |         26 |        16 |           10 | cells_per_donor,mean_umi_per_cell,mean_genes_per_cell,mean_pct_mito | batch,ancestry,treatment,site,chemistry,processing_date | limited_metadata_only      |

## 2. Paired method comparison

Paired donor-bootstrap ΔAUC tests were run within each cohort against scGPT mean where both methods had donor-matched OOF predictions. Holm correction was applied across successful comparisons.

Significant Holm-adjusted comparisons:

| cohort               | method_a         | method_b   |    auc_a |    auc_b |   delta_auc |     ci_low |     ci_high |   bootstrap_p_two_sided |    holm_p |
|:---------------------|:-----------------|:-----------|---------:|---------:|------------:|-----------:|------------:|------------------------:|----------:|
| GSE135779_matched500 | raw_pseudobulk   | scgpt_mean | 0.972452 | 0.84573  |   0.126722  |  0.0385675 |  0.23416    |              0.00179982 | 0.0323968 |
| GSE135779_matched500 | focus_lite@scgpt | scgpt_mean | 0.966942 | 0.84573  |   0.121212  |  0.0358127 |  0.225895   |              0.0019998  | 0.0323968 |
| GSE174188_CD4        | focus_lite@scgpt | scgpt_mean | 0.957912 | 0.978364 |  -0.0204514 | -0.0366629 | -0.00679636 |              0.00179982 | 0.0323968 |

Full comparison table:

| cohort               | method_a             | method_b   | status   |   n_donors |    auc_a |    auc_b |    delta_auc |      ci_low |     ci_high |   bootstrap_p_two_sided |   mean_delta |   sd_delta |    holm_p |
|:---------------------|:---------------------|:-----------|:---------|-----------:|---------:|---------:|-------------:|------------:|------------:|------------------------:|-------------:|-----------:|----------:|
| GSE135779_matched500 | raw_pseudobulk       | scgpt_mean | success  |         44 | 0.972452 | 0.84573  |  0.126722    |  0.0385675  |  0.23416    |              0.00179982 |  0.126425    | 0.0496713  | 0.0323968 |
| GSE135779_matched500 | donor_expression_pca | scgpt_mean | success  |         44 | 0.914601 | 0.84573  |  0.0688705   |  0.0110193  |  0.146006   |              0.0133987  |  0.0689906   | 0.0345187  | 0.187581  |
| GSE135779_matched500 | donor_mean_hvg       | scgpt_mean | success  |         44 | 0.892562 | 0.84573  |  0.046832    | -0.0385675  |  0.140496   |              0.262774   |  0.047143    | 0.0443843  | 1         |
| GSE135779_matched500 | geneformer_mean      | scgpt_mean | success  |         44 | 0.92562  | 0.84573  |  0.0798898   |  0.0110193  |  0.170799   |              0.0185981  |  0.0798165   | 0.0406102  | 0.241776  |
| GSE135779_matched500 | focus_lite@scgpt     | scgpt_mean | success  |         44 | 0.966942 | 0.84573  |  0.121212    |  0.0358127  |  0.225895   |              0.0019998  |  0.120847    | 0.0485406  | 0.0323968 |
| GSE135779_matched500 | red@scgpt            | scgpt_mean | success  |         44 | 0.955923 | 0.84573  |  0.110193    |  0.030303   |  0.209366   |              0.00419958 |  0.10988     | 0.0461297  | 0.0629937 |
| GSE135779_matched500 | tail_fractions@scgpt | scgpt_mean | success  |         44 | 0.909091 | 0.84573  |  0.0633609   | -0.0220386  |  0.157025   |              0.168983   |  0.0633372   | 0.0465983  | 1         |
| GSE174188_CD4        | raw_pseudobulk       | scgpt_mean | success  |        261 | 0.987904 | 0.978364 |  0.00953984  | -0.00237093 |  0.0235082  |              0.118988   |  0.00955907  | 0.00657262 | 1         |
| GSE174188_CD4        | donor_expression_pca | scgpt_mean | success  |        261 | 0.985659 | 0.978364 |  0.00729517  | -0.00486345 |  0.0203267  |              0.232177   |  0.00727171  | 0.0062619  | 1         |
| GSE174188_CD4        | donor_mean_hvg       | scgpt_mean | success  |        261 | 0.983976 | 0.978364 |  0.00561167  | -0.00717047 |  0.0188926  |              0.383962   |  0.00556351  | 0.00654879 | 1         |
| GSE174188_CD4        | kme_multiscale@scgpt | scgpt_mean | success  |        261 | 0.977865 | 0.978364 | -0.000498815 | -0.00904103 |  0.00748223 |              0.934707   | -0.000472515 | 0.00420979 | 1         |
| GSE174188_CD4        | focus_lite@scgpt     | scgpt_mean | success  |        261 | 0.957912 | 0.978364 | -0.0204514   | -0.0366629  | -0.00679636 |              0.00179982 | -0.0204237   | 0.00762192 | 0.0323968 |
| GSE174188_CD4        | attention_mil        | scgpt_mean | success  |        261 | 0.969946 | 0.978364 | -0.00841751  | -0.0230702  |  0.00405287 |              0.19678    | -0.00833675  | 0.00682709 | 1         |
| GSE174188_CD4        | cckme_u_weighted     | scgpt_mean | success  |        261 | 0.969759 | 0.978364 | -0.00860456  | -0.0232573  |  0.00299289 |              0.172583   | -0.00854825  | 0.00668939 | 1         |
| GSE285773            | raw_pseudobulk       | scgpt_mean | success  |         26 | 0.95625  | 0.925    |  0.03125     | -0.0875     |  0.15625    |              0.630137   |  0.0306688   | 0.0592389  | 1         |
| GSE285773            | donor_expression_pca | scgpt_mean | success  |         26 | 0.96875  | 0.925    |  0.04375     | -0.025      |  0.1375     |              0.30477    |  0.0430481   | 0.0416347  | 1         |
| GSE285773            | donor_mean_hvg       | scgpt_mean | success  |         26 | 0.9625   | 0.925    |  0.0375      | -0.05       |  0.1375     |              0.435556   |  0.03692     | 0.0459995  | 1         |
| GSE285773            | kme_multiscale@scgpt | scgpt_mean | success  |         26 | 0.8625   | 0.925    | -0.0625      | -0.16875    |  0.0125     |              0.172183   | -0.0622562   | 0.0462293  | 1         |

## 3. GSE135779 matched-500 scGPT mean underperformance

In GSE135779 matched-500, scGPT mean AUC was 0.846. The top method in the ranking was raw_pseudobulk with AUC 0.972.

The pattern supports a focused discussion that mean pooling can underuse donor-distribution structure in this cohort: distributional/structured summaries such as focus, RED, tail fractions, and moments often outperformed simple scGPT mean. This should be framed as cohort-specific evidence, not a universal mean-pooling failure claim.

Top GSE135779 methods:

| method_id                      |      auc |   delta_vs_scgpt_mean | method_family                      |
|:-------------------------------|---------:|----------------------:|:-----------------------------------|
| raw_pseudobulk                 | 0.972452 |             0.126722  | expression_or_isg                  |
| focus_lite@scgpt               | 0.966942 |             0.121212  | structured_scgpt_or_distributional |
| attention_mil@geneformer       | 0.964187 |             0.118457  | geneformer_or_structured           |
| red@geneformer                 | 0.964187 |             0.118457  | geneformer_or_structured           |
| tail_fractions@geneformer      | 0.961433 |             0.115702  | geneformer_or_structured           |
| focus_lite@geneformer          | 0.961433 |             0.115702  | geneformer_or_structured           |
| red@scgpt                      | 0.955923 |             0.110193  | structured_scgpt_or_distributional |
| prototype_histogram@geneformer | 0.942149 |             0.0964187 | geneformer_or_structured           |
| moments_mean_var@geneformer    | 0.942149 |             0.0964187 | geneformer_or_structured           |
| moments_mean_var@scgpt         | 0.92562  |             0.0798898 | structured_scgpt_or_distributional |
| geneformer_mean                | 0.92562  |             0.0798898 | geneformer_or_structured           |
| quantiles@geneformer           | 0.917355 |             0.0716253 | geneformer_or_structured           |
| cckme_u@scgpt                  | 0.917355 |             0.0716253 | structured_scgpt_or_distributional |
| donor_expression_pca           | 0.914601 |             0.0688705 | expression_or_isg                  |
| quantiles@scgpt                | 0.909091 |             0.0633609 | structured_scgpt_or_distributional |

## Output files

- `tables/EXT_P8_6_GSE174188_covariate_decomposition.csv`
- `tables/EXT_P8_6_GSE174188_single_covariate_auc.csv`
- `tables/EXT_P8_6_cross_cohort_covariate_availability.csv`
- `tables/EXT_paired_auc_delta_holm.csv`
- `tables/EXT_GSE135779_method_ranking_scgpt_mean_underperformance.csv`
- `figures/EXT_FIG_covariate_group_auc.png`
- `figures/EXT_FIG_paired_delta_vs_scgpt_mean.png`
- `figures/EXT_FIG_GSE135779_scgpt_mean_context.png`
