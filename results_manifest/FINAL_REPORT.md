# RheumLens final integration report

## Formal stage state

- P4 Geneformer provenance: accepted.
- P5 matched-500: accepted (27/27).
- P6 three-cohort benchmarks: accepted with modality-specific boundaries.
- P8.3/P8.4/P8.6/P8.7: accepted with documented caveats.
- P8.5: skipped because no positive KME increment triggered the null audit.
- P8.8: skipped after FOCUS superiority claim retraction; exploratory only.
- P9: source-only CD4 cross-study transfer completed by this run.

## P9 results

| direction                  | method_id            |      auc |   pr_auc |    brier |   calibration_intercept |   calibration_slope |   n_donors |   n_cases |   n_controls | auc_ci_method   |   auc_ci_low |   auc_ci_high |    auc_se | pr_auc_ci_method           |   pr_auc_ci_low |   pr_auc_ci_high |
|:---------------------------|:---------------------|---------:|---------:|---------:|------------------------:|--------------------:|-----------:|----------:|-------------:|:----------------|-------------:|--------------:|----------:|:---------------------------|----------------:|-----------------:|
| GSE285773_to_GSE174188_CD4 | scgpt_mean           | 0.730453 | 0.844889 | 0.443241 |                 1.80664 |            0.1741   |        261 |       162 |           99 | Hanley-McNeil   |     0.670368 |      0.790537 | 0.0306558 | stratified_donor_bootstrap |        0.808418 |         0.880371 |
| GSE285773_to_GSE174188_CD4 | donor_expression_pca | 0.790248 | 0.896082 | 0.252029 |                 1.97206 |            1.46075  |        261 |       162 |           99 | Hanley-McNeil   |     0.736718 |      0.843778 | 0.0273117 | stratified_donor_bootstrap |        0.867603 |         0.923666 |
| GSE285773_to_GSE174188_CD4 | donor_mean_hvg       | 0.806397 | 0.892524 | 0.313122 |                 1.71113 |            0.403405 |        261 |       162 |           99 | Hanley-McNeil   |     0.75492  |      0.857874 | 0.0262642 | stratified_donor_bootstrap |        0.861522 |         0.921864 |
| GSE174188_CD4_to_GSE285773 | scgpt_mean           | 0.88125  | 0.935056 | 0.615318 |               218.291   |           15.7873   |         26 |        16 |           10 | Hanley-McNeil   |     0.750462 |      1        | 0.0667299 | stratified_donor_bootstrap |        0.856068 |         0.990189 |
| GSE174188_CD4_to_GSE285773 | donor_expression_pca | 0.93125  | 0.96875  | 0.183072 |                -1.77845 |            0.687568 |         26 |        16 |           10 | Hanley-McNeil   |     0.832806 |      1        | 0.0502276 | stratified_donor_bootstrap |        0.916667 |         1        |
| GSE174188_CD4_to_GSE285773 | donor_mean_hvg       | 0.9      | 0.952679 | 0.476457 |                18.0714  |            1.37628  |         26 |        16 |           10 | Hanley-McNeil   |     0.780337 |      1        | 0.0610538 | stratified_donor_bootstrap |        0.887031 |         1        |

## P9 method availability

All configured P9 transfer methods completed successfully.

## Global claim boundary

Performance is retrospective and cohort-specific. Covariate sensitivity materially attenuated AUC; none of the results establishes causal biology or clinical utility.
