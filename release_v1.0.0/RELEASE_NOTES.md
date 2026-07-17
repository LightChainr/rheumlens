# Version 1.0.0

This release consolidates the manuscript, traceable donor-level evidence, and editable publication figures for the cross-cohort SLE donor-representation benchmark.

## Scientific additions

- Reciprocal source-only transfer of frozen Geneformer, HVG pseudobulk, and PCA pseudobulk.
- Shared source-internal fold sensitivity and paired target-donor comparisons.
- Source-size learning curves, matched 500/1000-cell analyses, and regularization paths.
- Fixed distribution-preserving pooling and source-only early-fusion stress tests.
- Repeated interferon conditioning with expression-matched module controls.
- Donor-mean embedding geometry, spectral rank, and cell-budget geometry preservation.

## Public research assets

- Complete manuscript snapshot in Markdown, DOCX, and PDF.
- Eleven editable, multi-panel SVG figures with PNG review copies.
- Thirty-seven figure-source records with an explicit provenance index.
- Analysis and figure-building scripts plus release-level validation reports.
- SHA256 manifests for the versioned evidence snapshot and repository.

## Software corrections

- Restored the previously omitted `rheumlens.data` package.
- Added pickle-free dense and sparse NPZ persistence.
- Added raw-count and embedding validation.
- Added source-ordered shared-feature alignment for cross-cohort expression transfer.
- Added regression tests for sparse persistence and cross-cohort feature ordering.

## Verification

- Python 3.12 clean-environment install and test suite: 13 passed.
- Synthetic end-to-end smoke workflow: passed as part of the test suite.
- Editable SVG/XML, PNG preview, Markdown image-link, and figure-source validation: passed.
- Historical tracked manifests and release checksums: passed.

The stable concept DOI is [10.5281/zenodo.20813922](https://doi.org/10.5281/zenodo.20813922). Zenodo will assign a version-specific DOI to this GitHub release.
