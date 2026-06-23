# Data and code availability plan

This document is the public-facing availability checklist for the RheumLens manuscript package.

## Code

- Public repository: <https://github.com/LightChainr/rheumlens>
- Archived software DOI: pending Zenodo release.
- Reproducible environment files: `environment.yml`, `Dockerfile`, `Singularity.def`.
- Minimal verification command:

```bash
conda env create -f environment.yml
conda activate rheumlens
bash scripts/reproduce_minimal.sh
bash scripts/verify_manifests.sh
```

## Source data

Raw data are obtained from public source datasets or associated repositories:

- GSE135779
- GSE174188
- GSE285773

## Processed data release

The processed-data archive should be DOI-tagged before submission and should include the assets listed in `docs/ARCHIVAL_DATA_MANIFEST_TEMPLATE.tsv`. At minimum, it should include:

1. donor-level fold definitions;
2. out-of-fold prediction tables;
3. bootstrap, permutation and repeated-CV distributions;
4. covariate sensitivity outputs;
5. source-only transfer predictions;
6. figure-ready and supplement-ready tables;
7. SHA256 manifests for every deposited object.

Large cell-level matrices and foundation-model embeddings are handled separately from the Git repository. They may be deposited as large archival objects or regenerated from documented recipes, depending on repository limits and source-data redistribution terms.

## Metadata repair

The GSE174188 feature-name repair table is derived from the source h5ad `var.feature_name` field. It should be deposited as `metadata/GSE174188_feature_names_from_h5ad.tsv` or regenerated from the source h5ad during preprocessing.

## Pre-submission checklist

- Replace `TBD` DOI placeholders in `README.md`, `CITATION.cff` and `.zenodo.json`.
- Fill `docs/ARCHIVAL_DATA_MANIFEST_TEMPLATE.tsv` with final SHA256 and size values.
- Verify repository manifests with `bash scripts/verify_manifests.sh`.
- Confirm that all manuscript figure legends identify sample size, uncertainty definition and abbreviation expansion.
