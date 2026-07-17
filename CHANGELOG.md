# Changelog

## 1.0.1 - 2026-07-17

### Fixed

- Stored the release manuscript DOCX as a regular Git object so Zenodo's automatically generated source archive contains the complete document instead of a Git LFS pointer.
- Retained all analyses, figures, source tables, checksums, and scientific conclusions unchanged from version 1.0.0.

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
