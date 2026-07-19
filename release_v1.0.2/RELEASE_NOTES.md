# Version 1.0.2

Version 1.0.2 is the submission-stage scientific and reproducibility update for the cross-cohort SLE donor-representation benchmark.

## Scientific corrections

- Replaced the primary unscaled fixed-`C=1` internal comparison with a 20-repeat donor-stratified analysis using training-fold feature scaling and nested regularization selection.
- Established that all three retained representations discriminate strongly within cohorts, while their internal ordering is preprocessing-sensitive.
- Preserved reciprocal source-only external transfer as the primary representation-selection evidence; the central pseudobulk transfer conclusion is unchanged.
- Distinguished external AUC from calibrated clinical risk and reported prevalence-null Brier references.

## New public assets

- Applied Sciences manuscript in Markdown, DOCX, and PDF.
- Supplementary Figures S1-S35 and indexed Supplementary Tables S1-S20.
- Cohort-characteristic and missingness summary in Supplementary Table S19.
- Corrected internal sensitivity results in Supplementary Table S20.
- Exact Geneformer extraction, tokenization, truncation, hidden-layer, dimensionality, and fold-scaling details.
- Eight final main figures, editable SVG sources, analysis scripts, source metrics, and QA records.

## Verification

- Final journal upload bundle: all automated package checks passed.
- Main manuscript and cover letter: rendered and visually reviewed.
- DOI resolution for the retained high-risk references: verified.
- Release and repository SHA256 manifests: generated and verified.
- Python package version metadata: synchronized at 1.0.2.

The stable Zenodo concept DOI is [10.5281/zenodo.20813922](https://doi.org/10.5281/zenodo.20813922). A version-specific DOI can be documented after Zenodo mints the v1.0.2 record.

