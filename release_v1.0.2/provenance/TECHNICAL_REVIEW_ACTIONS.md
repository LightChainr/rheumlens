# Technical Review Actions

Date: 18 July 2026

## Overall assessment

The external critique is substantially correct. The manuscript was scientifically close to submission, but calibration language, cohort reporting, implementation detail, version control, and several editorial defects required correction before upload.

## Corrections completed

1. Clinical-risk language was removed or narrowed. External AUC is now described as donor ranking/discrimination, not calibrated individual risk. The manuscript reports prevalence-null Brier scores and states that prospective validation, target-population recalibration, and decision-curve analysis were not performed.
2. Supplementary Table S19 now reports available case-control cohort characteristics, standardized differences, missingness, and donor-level QC. Unavailable clinical fields are identified rather than imputed.
3. Geneformer extraction is specified at the level needed to reproduce the retained representation: checkpoint revision, tokenizer and mapping assets, count normalization and gene ranking, CLS/SEP token handling, 4,094-gene truncation, final hidden-layer first-token extraction, and the 1,152-dimensional output.
4. The internal benchmark was rerun with fold-contained feature scaling and inner-CV regularization selection. The archived unscaled fixed-C=1 analysis is retained only as an implementation sensitivity.
5. Supplementary Tables S1-S20 are provided as one indexed workbook. Tables S19-S20 contain the new cohort-characteristic and scaled internal-sensitivity results.
6. Main-table numbering is now sequential (Tables 1-5). Literal `X` artifacts after repository links were removed. Random seeds are explicitly described as integer identifiers rather than dates.
7. The final DOCX and PDF were rendered and visually checked. The upload bundle passed automated checks for required sections, figure counts and dimensions, placeholders, archive integrity, file size, DOI/repository references, and the exact five-file upload set.

## Additional methodological issue identified locally

The previous primary internal comparison used fixed `C=1` without standardizing every final feature block. Because Geneformer and expression-derived features have different scales, that setup could alter regularization strength across representations. The corrected analysis uses training-fold scaling and nested selection of `C` from 10^-4 to 10^4.

Corrected mean ROC-AUC across 20 repeated donor-stratified five-fold analyses:

| Cohort | Geneformer | HVG pseudobulk | PCA pseudobulk |
|---|---:|---:|---:|
| GSE135779 | 0.904 | 0.949 | 0.946 |
| GSE174188 CD4 | 0.981 | 0.978 | 0.966 |
| GSE285773 CD4 | 0.872 | 0.872 | 0.890 |

The revised interpretation is that all three representations discriminate well internally, while their within-cohort ordering is preprocessing-sensitive. Strict external transfer therefore remains the more stable basis for representation selection. The reciprocal external-transfer results and the central pseudobulk conclusion were not changed by this internal sensitivity correction.

## Version-freeze status

The public v1.0.1 GitHub/Zenodo release remains an internally consistent archived prior release. Version 1.0.2 adds the fold-scaled analysis, Supplementary Tables S19-S20, revised manuscript, exact extraction details, and corrected package metadata. The release snapshot uses the stable Zenodo concept DOI so it remains valid before and after Zenodo mints the version-specific v1.0.2 DOI. The version DOI should be documented in a follow-up metadata commit after minting.
