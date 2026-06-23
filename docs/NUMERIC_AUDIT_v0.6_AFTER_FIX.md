# RheumLens manuscript v0.6 numeric audit after fix

- Manuscript: `/Users/lc/Documents/RheumLens/manuscript/RheumLens_manuscript_working_v0.6.md`
- Checks: 24
- PASS: 24
- CHECK: 0

| group            | claim                                           | expected   | status   | evidence                                                                  |
|:-----------------|:------------------------------------------------|:-----------|:---------|:--------------------------------------------------------------------------|
| P5 matched-500   | raw_pseudobulk                                  | 0.972      | PASS     | source auc=0.972452                                                       |
| P5 matched-500   | focus_lite@scgpt                                | 0.967      | PASS     | source auc=0.966942                                                       |
| P5 matched-500   | geneformer_mean                                 | 0.926      | PASS     | source auc=0.925620                                                       |
| P6 GSE174188     | raw_pseudobulk                                  | 0.988      | PASS     | source auc=0.987904                                                       |
| P6 GSE174188     | donor_expression_pca                            | 0.986      | PASS     | source auc=0.985659                                                       |
| P6 GSE174188     | donor_mean_hvg                                  | 0.984      | PASS     | source auc=0.983976                                                       |
| P6 GSE174188     | scgpt_mean                                      | 0.978      | PASS     | source auc=0.978364                                                       |
| P6 GSE174188     | kme_multiscale@scgpt                            | 0.978      | PASS     | source auc=0.977865                                                       |
| P8.4 cell-count  | common_support level 50 scgpt_mean              | 0.953      | PASS     | formal 10-repeat mean_auc=0.952637; manuscript expected 0.953             |
| P8.4 cell-count  | common_support level 50 donor_expression_pca    | 0.970      | PASS     | formal 10-repeat mean_auc=0.969614; manuscript expected 0.970             |
| P8.4 cell-count  | common_support level 50 kme_multiscale@scgpt    | 0.948      | PASS     | formal 10-repeat mean_auc=0.947976; manuscript expected 0.948             |
| P8.6 covariates  | covariates_only                                 | 0.846      | PASS     | source auc=0.846490                                                       |
| P8.6 covariates  | scgpt_mean_covres                               | 0.868      | PASS     | source auc=0.868250                                                       |
| P8.6 covariates  | donor_expression_pca_covres                     | 0.853      | PASS     | source auc=0.853348                                                       |
| P8.6 covariates  | donor_mean_hvg_covres                           | 0.870      | PASS     | source auc=0.869809                                                       |
| P9 transfer      | GSE285773_to_GSE174188_CD4 scgpt_mean           | 0.730      | PASS     | source auc=0.730453                                                       |
| P9 transfer      | GSE285773_to_GSE174188_CD4 donor_expression_pca | 0.790      | PASS     | source auc=0.790248                                                       |
| P9 transfer      | GSE285773_to_GSE174188_CD4 donor_mean_hvg       | 0.806      | PASS     | source auc=0.806397                                                       |
| P9 transfer      | GSE174188_CD4_to_GSE285773 scgpt_mean           | 0.881      | PASS     | source auc=0.881250                                                       |
| P9 transfer      | GSE174188_CD4_to_GSE285773 donor_expression_pca | 0.931      | PASS     | source auc=0.931250                                                       |
| P9 transfer      | GSE174188_CD4_to_GSE285773 donor_mean_hvg       | 0.900      | PASS     | source auc=0.900000                                                       |
| Mechanism        | ISG donor score AUC                             | 0.906      | PASS     | source mechanism report                                                   |
| Mechanism        | Cohen                                           | 1.372      | PASS     | source mechanism report                                                   |
| Regression guard | old P8.4 smoke values absent from P8.4 prose    | absent     | PASS     | prevents single-repeat smoke values from replacing formal 10-repeat means |
