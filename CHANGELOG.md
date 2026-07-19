# Changelog

## 1.0.2 - 2026-07-19

### Added

- A 20-repeat donor-stratified internal benchmark with training-fold feature scaling and nested regularization selection.
- Supplementary Table S19 with available case-control cohort characteristics, standardized differences, missingness, and donor-level QC.
- Supplementary Table S20 with corrected internal sensitivity results.
- Exact Geneformer extraction details and a self-contained Applied Sciences submission snapshot.

### Changed

- Reframed internal AUCs as preprocessing-sensitive and retained strict reciprocal external transfer as the primary representation-selection evidence.
- Distinguished external discrimination from calibrated clinical risk using prevalence-null Brier references and explicit scope language.
- Synchronized package metadata and runtime `__version__` at 1.0.2.

### Fixed

- Corrected the primary internal comparison, which previously used unscaled feature blocks with fixed `C=1`.
- Restored sequential main-table numbering, complete Supplementary Tables S1-S20, and clean repository-link rendering.
- Archived as version DOI [10.5281/zenodo.21436893](https://doi.org/10.5281/zenodo.21436893).

## 1.0.1 - 2026-07-17

### Fixed

- Stored the release manuscript DOCX as a regular Git object so Zenodo's automatically generated source archive contains the complete document instead of a Git LFS pointer.
- Retained all analyses, figures, source tables, checksums, and scientific conclusions unchanged from version 1.0.0.
- Archived as version DOI [10.5281/zenodo.21412436](https://doi.org/10.5281/zenodo.21412436).

## 1.0.0 - 2026-07-17

### Added

- Final integrated manuscript snapshot in Markdown, DOCX, and PDF.
- Eleven editable multi-panel SVG figures and review PNGs.
- Exact figure-source tables and a figure-to-source provenance index.
- Shared source-internal fold sensitivity for reciprocal external transfer.
- Donor-mean embedding geometry and cap500-cap1000 geometry preservation analyses.
- Direct cap500 pooling stress test, fairness sensitivities, and regularization-path analysis.
- Release-level evidence disposition, document visual QA, and checksum validation.

### Changed

- Reframed the public project around cross-cohort patient-representation selection rather than an audit-branded manuscript identity.
- Updated the Zenodo title, description, keywords, release date, and citation metadata.
- Updated the Python package version to 1.0.0 while retaining the existing package name for compatibility.
- Replaced pending DOI language with the stable Zenodo concept DOI.

### Scope

The release supports claims about the evaluated SLE cohorts, frozen Geneformer representations, expression pseudobulk, and fixed donor-level aggregation. It does not claim clinical deployment readiness or general superiority over fine-tuned or learned patient models.
