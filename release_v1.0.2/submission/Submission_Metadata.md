# Submission Metadata

## Destination

- Journal: Applied Sciences
- Article type: Article
- Section: Computing and Artificial Intelligence
- Special Issue: Research on Computational Biology and Bioinformatics
- Special Issue URL: https://www.mdpi.com/journal/applsci/special_issues/SSY60TBE8W

## Title

Cross-Cohort Benchmarking of Frozen Geneformer and Expression Pseudobulk for Donor-Level SLE Classification

## Authors

1. Hongyu Ying - Department of Rheumatology and Immunology, Shanghai Pudong Hospital, Fudan University Pudong Medical Center, Shanghai 201399, China
2. Dandan Yun - same affiliation
3. Dan Liu - same affiliation; corresponding author; danliu600@126.com

## Abstract

Patient-level classifiers derived from single-cell RNA sequencing must preserve disease information across cohorts, yet most representation comparisons remain within-dataset. We benchmarked mean-pooled frozen Geneformer embeddings against highly variable gene (HVG) and principal-component pseudobulk for systemic lupus erythematosus classification in three cohorts comprising 44, 261, and 26 donors. Evaluation used 20 repeated donor-level five-fold cross-validations and reciprocal external transfer between the two CD4-positive cohorts, with scaling, feature selection, dimensionality reduction, and regularization fitted exclusively in training donors. Internal ordering varied by cohort after fold-contained scaling and regularization selection, although all representations retained strong discrimination. When the 261-donor cohort served as target, HVG and PCA pseudobulk achieved AUCs of 0.919 and 0.926, compared with 0.884 for Geneformer. Source-only fusion and seven fixed distribution-preserving pooling alternatives produced no stable bidirectional improvement. Source-size curves indicated that the pseudobulk advantage widened as training donors accumulated, while cell-budget analyses showed little gain beyond 500-1000 cells per donor. Interferon-program residualization reduced Geneformer discrimination, linking predictive signal to structured SLE expression. External scores measured discrimination rather than calibrated clinical risk. Under frozen, fixed-pooling conditions, source-fitted pseudobulk provided the stronger transferable donor representation and a practical benchmark for future learned patient models.

Abstract word count: 196

## Keywords

single-cell RNA sequencing; systemic lupus erythematosus; Geneformer; pseudobulk; patient-level prediction; cross-cohort transfer; foundation models; bioinformatics benchmarking

## Featured Application

This study provides an externally evaluated workflow for selecting donor-level representations in clinical single-cell transcriptomics and a reproducible benchmark for testing patient-level foundation-model embeddings under cohort shift.

## Suggested Short Title

Cross-Cohort SLE Donor Representations

## Data and Code

- Stable Zenodo concept DOI: https://doi.org/10.5281/zenodo.20813922
- GitHub: https://github.com/LightChainr/rheumlens
