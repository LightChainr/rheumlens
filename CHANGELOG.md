# Changelog

## 3.0.0-rc1 - 2026-09-08

### Added

- A nine-comparison, five-dataset public screen (809 donors) spanning lupus, COVID-19, influenza, cytomegalovirus infection, and sepsis-versus-COVID-19.
- Cohort-specific permutation nulls for design-only AUC and cross-fitted residual label variance.
- Free and collection-preserving complete-pipeline permutations, including an explicit undefined state when no collection stratum contains both labels.
- Seven simulation regimes across seven design-diagnosis association levels and two outcome types (19,600 simulated cohorts).
- A prespecified CMV negative control, variable-group localization, and decision-tree walkthroughs tied to the same collection strata used by the restricted permutation.

### Changed

- Removed the identifiability-boundary framing and recast the reported quantities as scoped association and validation checks.
- Reframed the five checks as complementary diagnostics rather than a required validity standard.
- Separated collection variables from demographics/case mix and sample-quality summaries throughout the manuscript.
- Expanded the empirical scope beyond SLE while keeping the detailed lupus representation/use, residualisation/restriction, and source-only transfer analyses.
- Updated repository, citation, Zenodo and claim-boundary metadata to the PLOS Computational Biology resubmission framing.

### Submission snapshot

- Prepared the 2026-09-08 PLOS Computational Biology submission snapshot under `ploscb_submission_2026-09-08/`.
- Generated a 72-page, double-spaced, continuously line-numbered submission PDF with eight main figures embedded.

## 2.0.0-rc3 - 2026-07-27

### Added

- A quantitative design-adjusted label-information scale, 11,200 simulations and
  an exact observational-equivalence construction.
- Complete-pipeline nulls, nonlinear design-only prediction and fold-contained
  random-forest residualisation.
- An internal batch-exposure negative control, attenuation-difference intervals,
  a synthetic positive control and confound-leakage sensitivity experiment.
- Eight main figures, Supplementary Tables S1-S30 and a complete PLOS
  Computational Biology review object.
- A five-check validity standard spanning estimability, full-pipeline nulls,
  design exposure, overlap restriction and source-only external transfer.

### Changed

- Reframed the manuscript from a representation comparison to a general
  identifiability problem in patient-level single-cell classification.
- Rebuilt Figure 6 as a six-panel robustness analysis.
- Updated the abstract, Author Summary, cover letter, references, reproducibility
  documentation and release validator.

### Validated

- Four theory tests, 18 identifiability-extension checks, 16 PLOS robustness
  checks and the 203-file release manifest passed.
- Archived as version DOI
  [10.5281/zenodo.21618824](https://doi.org/10.5281/zenodo.21618824).

## 2.0.0-rc1 - 2026-07-27

### Added

- A two-cohort design-validity reconstruction with restored GSE135779 batch,
  collection-year, demographic and sequencing metadata.
- Repeated donor-level design-only prediction, representation-to-design recovery,
  fold-contained residualisation, matched design restriction and source-only
  transfer.
- Six new main figures, locked donor-level predictions, source tables, methods
  provenance, verified references and a package-level SHA256 manifest.
- A four-check validity standard for patient-level single-cell classification.

### Changed

- Reframed the scientific object from a Geneformer-versus-pseudobulk leaderboard to
  attribution and transportability under design-label entanglement.
- Restricted learned-pooling claims to independent-target evidence.
- Updated release-candidate creator metadata to Hongyu Ying, Dandan Yun and Dan Liu.

### Preserved

- Version 1.0.2 and DOI
  [10.5281/zenodo.21436893](https://doi.org/10.5281/zenodo.21436893) remain unchanged.
- The `rheumlens` repository and Python package names remain for backward
  compatibility and are not used as manuscript terminology.

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
