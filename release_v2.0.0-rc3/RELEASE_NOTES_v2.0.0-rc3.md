# v2.0.0-rc3 release notes

This release candidate converts the RC2 design-identifiability study into a
PLOS Computational Biology-ready review object.

## Added

- Nonlinear complete-design prediction using random forest and gradient boosting.
- Fold-contained random-forest residualisation for Geneformer, HVG pseudobulk and
  PCA pseudobulk.
- Batch-location adjustment and overlap-weighted sensitivity estimands.
- Internal batch-exposure negative control.
- Paired split-sensitivity intervals comparing residualisation attenuation with
  matched design restriction.
- Synthetic positive-control scenario with five prespecified pass checks.
- Gaussian confound-leakage sensitivity experiment.
- Figure 6 rebuilt as a six-panel robustness figure.
- Supplementary Tables S24-S30.
- Executable robustness validator with 16 checks.

## Revised

- Abstract, Introduction, Results, Discussion, Methods and Data Availability.
- Author Summary and cover letter for PLOS Computational Biology.
- Confound-control and single-cell benchmark references.
- Reproducibility documentation, metadata and package validator.

## Release identity

GitHub tag: `v2.0.0-rc3`. Stable Zenodo concept DOI:
`10.5281/zenodo.20813922`. Zenodo assigns a version-specific DOI from the GitHub
release.
