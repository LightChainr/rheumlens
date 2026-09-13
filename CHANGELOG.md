# Changelog

## 3.0.0-rc2 - 2026-09-13

### Fixed

- Every fitted pipeline now runs through `scripts/cohorts/pipeline_core.py` (shipped as `ploscb_submission_2026-09-08/02_SOURCE/analysis_scripts/pipeline_core.py`). Imputation, one-hot encoding, the gene filter, standardisation, PCA and ridge residualisation are fitted on training donors only in the screen, the walkthrough and the calibration simulation.
- The metadata matrix's median imputation and one-hot levels were fitted on the whole cohort; they are now fitted per training fold (`run_design_screen.py --refresh-design` recomputes only the metadata columns). Metadata-only AUCs moved by at most 0.012; expression columns are unchanged and the 8/9 and 4/7 counts hold.
- Calibration arm A no longer fits PCA on all donors; arm B's unadjusted and residualised arms are the same call with one argument added. Both re-run (300 replicates per strength).
- The degeneracy gate has one criterion (no label-mixed collection stratum); a metadata-only AUC of 1 is no longer a trigger. Route A no longer recommends adjustment; restriction feasibility is a fixed rule.
- Table S8 (hyperparameters) is generated from the `PipelineSpec`/`ForestSpec` objects that ran.
- Bibliography corrections (references 19, 30, 31, 34), all checked against Crossref.

### Added

- `run_incremental.py`: AUC of metadata plus expression over metadata alone, five seeds (Table S5).
- `run_calibration.py --shard/--merge-from` for splitting the simulation across machines.
- Point-by-point response to the previous decision (`01_UPLOAD/Response_to_Previous_Review.md`).

### Changed

- Title: "Recorded metadata predicts the phenotype label across public single-cell cohorts and changes what internal validation can establish". Terminology unified (phenotype label, collection-stratified permutation, fixed-hyperparameter pipeline).
- Contribution reframed as an empirical characterisation of benchmark validity: four questions (presence, source, validation, transport) and three validation regimes (low recorded association, associated but conditionally evaluable, structurally non-overlapping).
- Results reordered so the nine-comparison screen (Section 3) precedes the simulations that interpret it (Section 4); Figures 2 and 3 swapped; supporting figures and tables renumbered into citation order (Figures S1-S7, Tables S1-S11).
- Subtractive language edit of the main text (about 27% shorter), removing meta-commentary and revision history; no number changed.
- Added a retrospective CELLxGENE Discover eligibility audit (`audit/cellxgene_eligibility_audit.py`, Table S9): 35 datasets in 15 collections meet the screen's dataset-level requirements, including all five analysed collections. Performed after the analysis and not used for selection.

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
- Generated a double-spaced, continuously line-numbered submission PDF with the figures embedded. Rebuilt at 60 pages with all 15 figures after the fold-contained re-run; it is produced from the manuscript markdown by `build_pdf.py` rather than from a separate LaTeX copy of the body.
- Committed the research object itself under `ploscb_submission_2026-09-08/`: manuscript source, cover letter, PDF, all 15 figures, Supporting Tables S1-S9, analysis and figure scripts, result tables, the five per-seed screen outputs, cohort registry, donor-level interface tables and QA records.
- Added `REPRODUCE.md`, `SHA256SUMS` and `verify.sh`, which re-runs the 24 number checks and the structural checks against the released files alone.
- Moved the screen's top-variance gene filter and PCA inside the cross-fitting fold and regenerated all 45 (comparison x seed) runs. Every design-only quantity is unchanged bit for bit; the diagnosis AUC and the two permutation p-values change, and the frozen and tuned classifiers now agree to within 0.066 AUC instead of 0.518.

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
