# Strict Source-only Transfer Recalculation

## Design

All learned preprocessing, feature selection, scaling, PCA, and classifier regularization selection was fitted on source donors only. Target cohorts contributed only already-normalized donor vectors and unlabeled feature rows.

## Results

- SLE_GSE174188_CD4_to_SLE_GSE285773_CD4 / Frozen Geneformer: AUC 0.900 (0.762-0.994); Brier 0.615.
- SLE_GSE174188_CD4_to_SLE_GSE285773_CD4 / Source-HVG pseudobulk: AUC 0.944 (0.844-1.000); Brier 0.544.
- SLE_GSE174188_CD4_to_SLE_GSE285773_CD4 / Source-PCA pseudobulk: AUC 0.981 (0.925-1.000); Brier 0.180.
- SLE_GSE285773_CD4_to_SLE_GSE174188_CD4 / Frozen Geneformer: AUC 0.884 (0.843-0.920); Brier 0.237.
- SLE_GSE285773_CD4_to_SLE_GSE174188_CD4 / Source-HVG pseudobulk: AUC 0.919 (0.884-0.946); Brier 0.155.
- SLE_GSE285773_CD4_to_SLE_GSE174188_CD4 / Source-PCA pseudobulk: AUC 0.926 (0.895-0.953); Brier 0.159.

## Paired comparisons

- SLE_GSE174188_CD4_to_SLE_GSE285773_CD4 / source_hvg_pseudobulk: Geneformer minus baseline AUC -0.044; paired DeLong BH p=0.4914.
- SLE_GSE174188_CD4_to_SLE_GSE285773_CD4 / source_pca_pseudobulk: Geneformer minus baseline AUC -0.081; paired DeLong BH p=0.1061.
- SLE_GSE285773_CD4_to_SLE_GSE174188_CD4 / source_hvg_pseudobulk: Geneformer minus baseline AUC -0.035; paired DeLong BH p=0.01513.
- SLE_GSE285773_CD4_to_SLE_GSE174188_CD4 / source_pca_pseudobulk: Geneformer minus baseline AUC -0.042; paired DeLong BH p=0.0006027.

## Boundary

This recomputation does not infer single-patient prospective deployment performance beyond the evaluated two-cohort transfer setting.
