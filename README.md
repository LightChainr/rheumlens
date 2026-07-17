# Cross-cohort SLE donor-representation benchmark

> Frozen Geneformer and expression pseudobulk for patient-level classification under cohort shift

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20813922.svg)](https://doi.org/10.5281/zenodo.20813922)

This repository provides code, donor-level results, editable figures, and manuscript assets for a cross-cohort benchmark of patient representations derived from single-cell RNA sequencing. The public repository name and Python package name remain `rheumlens` for backward compatibility; they are not used as a manuscript concept or method name.

## Scientific question

Single-cell RNA sequencing measures cells, whereas clinical classification is performed at the patient level. We test which donor representation best preserves systemic lupus erythematosus (SLE) discrimination when cohort, population, and technical context change.

The primary comparison is between:

- mean-pooled frozen Geneformer cell embeddings;
- source-selected highly variable gene pseudobulk;
- source-fitted PCA pseudobulk.

All supervised analyses use donors as the independent unit. Reciprocal external transfer uses source-only feature selection, scaling, dimensionality reduction, and regularization selection.

## Main findings

- All three representation families retain strong within-cohort donor discrimination.
- Pseudobulk provides the stronger external-transfer representation in the evaluated SLE cohorts.
- In the 261-donor independent target, Geneformer, HVG pseudobulk, and PCA pseudobulk achieved AUCs of 0.884, 0.919, and 0.926.
- A shared source-internal fold sensitivity preserved the larger-target pseudobulk advantage.
- Source-only early fusion and seven fixed distribution-preserving pooling alternatives produced no stable bidirectional improvement.
- The pseudobulk advantage widened as source donors accumulated, while per-donor performance plateaued between 500 and 1000 cells in the evaluated Geneformer configurations.
- Geneformer discrimination was associated with an interferon-related expression axis, but matched structured expression modules produced similar attenuation.
- Donor-mean Geneformer embeddings occupied a low-dimensional subspace relative to sampled cell embeddings; this is a descriptive geometric result, not a causal explanation of transfer performance.

## Cohorts

| Dataset | Donors | Analysis population |
|---|---:|---|
| GSE135779 | 44 | PBMC-derived donor benchmark |
| GSE174188 | 261 | CD4-positive T cells |
| GSE285773 | 26 | CD4-positive T cells |

## Version 1.0.1 archive

Version 1.0.1 is a packaging-only correction that stores the manuscript DOCX directly in the repository archive rather than as a Git LFS pointer. The research snapshot remains under `release_v1.0.0/` and contains the manuscript and all newly consolidated public assets:

```text
release_v1.0.0/
├── manuscript/           # Markdown, DOCX, and visually verified PDF
├── figures/svg/          # 11 editable multi-panel SVG figures
├── figures/previews/     # PNG review copies
├── source_data/          # 37 traceable figure-source records
├── scripts/              # Figure, transfer, pooling, and validation scripts
├── provenance/           # Source mapping, evidence boundaries, and QA reports
└── MANIFEST_SHA256.txt   # Release-level checksums
```

The broader repository retains the original benchmarking package, tests, historical supplementary tables, and earlier fixed-split analyses for provenance.

## Reproduce and validate

```bash
conda env create -f environment.yml
conda activate rheumlens
pip install -e ".[dev,io]"

pytest -q
bash scripts/reproduce_minimal.sh
bash scripts/verify_manifests.sh
```

The release-specific plotting and validation scripts are under `release_v1.0.0/scripts/`. Several scripts use project-relative archived inputs recorded in `release_v1.0.0/provenance/figure_source_index.tsv`; the exact source tables are copied under `release_v1.0.0/source_data/`.

## Scope

The results apply to the evaluated public SLE cohorts, frozen Geneformer embeddings, and fixed donor-level aggregation. They do not establish clinical diagnostic readiness, causal disease mechanisms, formal cell-count sufficiency, or general inferiority of fine-tuned or learned patient representations. See [CLAIM_BOUNDARY.md](CLAIM_BOUNDARY.md) and the release evidence disposition for the complete scope.

## Data and code availability

- Repository: <https://github.com/LightChainr/rheumlens>
- Stable concept DOI: <https://doi.org/10.5281/zenodo.20813922>
- Current version DOI (`v1.0.1`): <https://doi.org/10.5281/zenodo.21412436>
- Superseded version DOI (`v1.0.0`): <https://doi.org/10.5281/zenodo.21412278>
- Previous version DOI (`v0.1.1`): <https://doi.org/10.5281/zenodo.20813923>
- Raw datasets: GEO GSE135779, GSE174188, and GSE285773

Large raw matrices and cell-level embedding archives are not redistributed in Git. Their provenance, hashes, and regeneration context are retained in the project records. Donor-level and plot-ready outputs needed for the public manuscript figures are included in the versioned research archive.

## Citation

Use the versioned Zenodo citation generated for the release or the stable concept DOI above. Machine-readable citation metadata are provided in [CITATION.cff](CITATION.cff).

## License

MIT. See [LICENSE](LICENSE).

## Affiliation

Shanghai Pudong Hospital / 上海市浦东医院
