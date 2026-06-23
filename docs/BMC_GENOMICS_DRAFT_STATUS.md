# BMC Genomics draft status

Date: 2026-06-23

## Deliverables

The BMC-oriented manuscript package is stored in:

`manuscript_bmc/`

Key files:

- `RheumLens_BMC_Genomics_v0.1.md` (internal filename; manuscript text marked v0.2)
- `RheumLens_BMC_Genomics_v0.1.docx` (internal filename; header/text marked v0.2)
- `BMC_SUBMISSION_CHECKLIST.md`
- `BMC_ASSET_BUILD_REPORT.md`
- `MANIFEST_SHA256.tsv`
- `tables_bmc/`
- `figures_bmc/`
- `formal_key_analyses/`
- `additional_files/`

## Positioning

The BMC version reframes the project from a foundation-model benchmark to a single-cell transcriptomic disease-representation study:

> Donor-level analysis of public SLE single-cell transcriptomic cohorts shows that frozen foundation-model embeddings recover reproducible disease signal, while expression summaries remain highly competitive and high AUCs are materially shaped by donor/sample covariates and interferon biology.

## Main BMC figures

1. Cohort overview.
2. Primary donor-level representation anchors.
3. Adult-pediatric CD4 source-only transfer.
4. GSE174188 CD4 interferon signature.
5. Model disease scores versus ISG score.
6. Formal GSE174188 covariate sensitivity.
7. Formal cross-cohort ISG scoring.
8. Formal GSE174188 incremental models.
9. Leave-processing-cohort-out robustness.
10. Covariate-matched subset sensitivity.
11. Calibration and Brier-score summary.

Additional Figure S1 reports the one-hot covariate+ISG+representation-score sensitivity model and is not used as the formal covariate estimate.

## Main BMC tables

- `BMC_TABLE_1_COHORT_CHARACTERISTICS.csv`
- `BMC_TABLE_PRIMARY_METHOD_SUMMARY.csv`
- `BMC_TABLE_STRUCTURED_TRANSFER_SUMMARY.csv`
- `BMC_TABLE_ISG_CROSS_COHORT_SUMMARY.csv`
- `BMC_TABLE_ISG_MODEL_SCORE_CORRELATION.csv`
- `BMC_TABLE_FORMAL_COVARIATE_SENSITIVITY.csv`
- `BMC_TABLE_COVARIATE_ISG_INCREMENTAL_MODELS.csv` as sensitivity only
- `formal_key_analyses/tables/FORMAL_CROSS_COHORT_ISG_SUMMARY.csv`
- `formal_key_analyses/tables/FORMAL_GSE174188_INCREMENTAL_MODEL_SUMMARY.csv`
- `formal_key_analyses/tables/FORMAL_GSE174188_LEAVE_PROCESSING_COHORT_OUT_SUMMARY.csv`
- `formal_key_analyses/tables/FORMAL_GSE174188_PROPENSITY_MATCHED_SUBSET_SUMMARY.csv`
- `formal_key_analyses/tables/FORMAL_CALIBRATION_SUMMARY.csv`

## Verification

- BMC asset build script ran successfully.
- BMC DOCX generated successfully.
- DOCX render QA completed for v0.2: 16 rendered pages, 11 embedded figures, 1 embedded table.
- Current BMC package manifest: `manuscript_bmc/MANIFEST_SHA256.tsv`.

## Remaining blockers before BMC submission

- Fill authors, affiliations and corresponding author.
- Add final ethics, funding, author contribution, competing-interest and acknowledgement text.
- Mint Zenodo/figshare DOI values and replace placeholder data availability language.
- Finalize references in BMC style with DOI/PMID checks.
- Optional reviewer-response item: add graphical calibration curves for source-only transfer; Brier/calibration summary is already included.
