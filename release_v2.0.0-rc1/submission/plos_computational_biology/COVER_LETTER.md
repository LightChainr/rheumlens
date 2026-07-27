# Draft cover letter

Dear Editors,

We submit the Research Article **“When design predicts disease: residualisation
failure and generalisation gaps in donor-level single-cell classifiers”** for
consideration in *PLOS Computational Biology*.

Patient-level single-cell studies increasingly report high donor-level
cross-validated discrimination after compressing thousands of cells into one
representation. The unresolved question is whether that separation reflects
disease, study design, or a mixture that will fail outside the originating cohort.
We address this question by reconstructing donor-level design metadata in two widely
reused systemic lupus erythematosus cohorts and placing frozen Geneformer,
HVG-pseudobulk and PCA-pseudobulk representations under one repeated,
donor-isolated pipeline.

The result is a methodological contradiction with broad consequences. Complete
measured design predicts disease at AUC 0.953 and 0.952 in the two cohorts, while
the molecular representations recover recorded batch identity at AUC up to 0.9997.
Yet disease discrimination largely survives restriction to overlapping designs with
size- and label-matched controls. Conventional fold-contained residualisation gives
the opposite answer, removing up to 0.430 AUC and reducing all representations to
chance in one cohort. Independent cross-batch and cross-cohort transfer, rather than
residualisation, exposes the reproducible generalisation gap. Three learned pooling
architectures do not repair it.

These results establish that design predictability, residualisation, restriction and
external transfer are not interchangeable confounding diagnostics. We convert the
finding into four executable checks for patient-level single-cell studies: report
design-only disease AUC, representation-to-design predictability,
design-restricted performance with matched donor controls, and source-only external
performance. The contribution is therefore not another model leaderboard. It is an
empirically supported validity standard for computational studies whose statistical
unit is a patient but whose measurements are individual cells.

The work is suited to *PLOS Computational Biology* because it combines a timely
machine-learning problem with a general lesson about computational inference from
high-dimensional biological data. All analyses use public data. Code, restored
metadata, donor-level predictions, source tables, figures, software details and
checksums are being assembled as a versioned open release.

The manuscript is original, is not under consideration elsewhere, and has been
approved by all authors. The authors declare no competing interests.

Sincerely,

Hongyu Ying, Dandan Yun and Dan Liu  
Department of Rheumatology and Immunology  
Shanghai Pudong Hospital, Fudan University Pudong Medical Center  
Shanghai 201399, China

