# Data and code availability

This document describes the public assets for the cross-cohort SLE donor-representation benchmark.

## Code

- Public repository: <https://github.com/LightChainr/rheumlens>
- Stable Zenodo concept DOI: <https://doi.org/10.5281/zenodo.20813922>.
- Version 1.0.2 DOI: <https://doi.org/10.5281/zenodo.21436893>.
- Submission-stage research snapshot: `release_v1.0.2/`.
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

The versioned release includes or indexes:

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

## Release verification

- Verify repository manifests with `bash scripts/verify_manifests.sh`.
- Verify the Python package and tests with `pytest -q`.
- Cite version DOI `10.5281/zenodo.21436893` when referring to the submission-stage snapshot.
