# RheumLens

> Donor-level audit benchmark for single-cell foundation-model disease representations in systemic lupus erythematosus

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/TBD.svg)](https://doi.org/TBD)

RheumLens is a donor-level benchmark and audit workflow for evaluating single-cell foundation-model disease representations. It uses donor-disjoint folds, out-of-fold prediction, source-only transfer, paired resampling, representation provenance checks, cell-count sensitivity, covariate sensitivity, and kernel diagnostics to produce auditable classification results for systemic lupus erythematosus (SLE) single-cell RNA-seq cohorts.

## Key findings

- **Expression baselines are highly competitive.** Raw pseudobulk, donor expression PCA, and donor mean HVG frequently match or exceed foundation-embedding methods in donor-level AUC.
- **Foundation embeddings carry reproducible SLE signal.** Frozen scGPT and repaired Geneformer embeddings recover donor-level disease discrimination across matched-500, adult CD4, and pediatric CD4 cohorts.
- **Covariate sensitivity matters.** In GSE174188 CD4, measured donor/sample covariates alone achieve AUC 0.846. Residualized methods are only modestly above this baseline — covariate sensitivity should be a primary audit endpoint.
- **Cross-study transfer is informative.** Source-only transfer between pediatric (GSE285773) and adult (GSE174188) CD4 datasets shows clear cross-cohort signal.
- **Known biology anchors the benchmark.** Donor-level gene analysis recovers a canonical type I interferon SLE axis with ISG score AUC 0.906 (Cohen d=1.37).

## Datasets

| Dataset | Donors | Cell type | Description |
|---|---|---|---|
| GSE135779 matched-500 | 44 | PBMC subsets | Balanced, controlled setting |
| GSE174188 CD4 | 261 | CD4⁺ T cells | Adult SLE case-control |
| GSE285773 | 26 | CD4⁺ T cells | Pediatric SLE case-control |

## Installation

Recommended reproducible setup:

```bash
conda env create -f environment.yml
conda activate rheumlens
pip install -e ".[dev,io]"
```

Container recipes are provided as `Dockerfile` and `Singularity.def` for CPU-side manuscript verification and table/figure audits. GPU-specific embedding extraction may require a CUDA/PyTorch image matched to the local driver.

Minimal editable install:

```bash
pip install -e .
```

For optional dependencies:

```bash
# Full analysis stack (including GPU support)
pip install -e ".[full]"

# IO only (AnnData, Parquet)
pip install -e ".[io]"

# Development
pip install -e ".[dev]"
```

## Quick start

```bash
# Run tests
pytest -q

# Run a minimal reproduction check
bash scripts/reproduce_minimal.sh

# Verify checked-in SHA256 manifests
bash scripts/verify_manifests.sh

# Run P9/P10 integration (requires processed data)
python scripts/mac_p9_p10.py --help
```

## Project structure

```
rheumlens/
├── src/rheumlens/           # Core Python package
│   ├── aggregators/         # Donor-level aggregation methods (KME, RED, FOCUS, moments, etc.)
│   ├── bag_models/          # MIL/attention/donor-set models
│   ├── data/                # Data loading, validation, operations
│   ├── estimators/          # Kernel and linear estimators
│   ├── evaluation/          # Engine, metrics, bootstrap, permutation, diagnostics
│   ├── postprocessors/      # Covariate residualization
│   ├── providers/           # Expression and embedding providers
│   └── utils/               # Array ops, manifest tools
├── tests/                   # Unit and smoke tests
├── scripts/                 # Reproduction and integration scripts
├── configs/                 # Project configuration YAML
├── manuscript/              # Working manuscript (v0.6)
├── figures/                 # Plot-ready manuscript figures (PNG+PDF)
├── supplementary_tables/    # Supplementary CSV tables (S1–S18)
├── docs/                    # Audit documents and assembly status
├── results_manifest/        # Data asset manifests with SHA256
├── pyproject.toml           # Package metadata and dependencies
├── CLAIM_BOUNDARY.md        # Allowed and disallowed interpretation claims
└── LICENSE                  # MIT license
```

## Reproducibility

The formal manuscript claims are supported by:

- Donor-level fixed folds with out-of-fold predictions
- 1,000 donor-label permutation tests (empirical P=0.000999)
- 30 repeated donor-level cross-validation runs
- Cell-count sensitivity across 25–500 cells per donor
- Covariate sensitivity with residualized methods
- SHA256-manifested data assets and audit reports

Large processed matrices and cell-level embeddings are **not** included in this repository. They are tracked by SHA256 in `results_manifest/RHEUMLENS_KEY_ASSETS.csv` and should be distributed through an archival data release.

## Data and code availability

- Code repository: <https://github.com/LightChainr/rheumlens>
- Archived code DOI: pending Zenodo release.
- Processed-data DOI: pending Zenodo/figshare release.
- Raw source datasets: GEO GSE135779, GSE174188 and GSE285773.

The archival data release should include donor-level fold files, out-of-fold predictions, repeated-CV summaries, bootstrap and permutation distributions, covariate sensitivity outputs, source-only transfer predictions, plot-ready tables and SHA256 manifests. Large cell-level matrices and foundation-model embeddings are distributed as separate archival objects or regenerated from documented recipes, subject to source-data redistribution terms. The GSE174188 feature-name repair table is derived from the source h5ad `var.feature_name` field and is included as a documented metadata asset when available.

## Claim boundary

See [CLAIM_BOUNDARY.md](CLAIM_BOUNDARY.md) for the full list of allowed and disallowed interpretation claims. In brief:

- **Allowed**: reproducible disease signal, competitive expression baselines, cross-study transfer signal, covariate sensitivity findings, ISG biology anchors.
- **Disallowed**: clinical deployment readiness, causal mechanism discovery, foundation-model superiority claims, generalization to untested platforms.

## Citation

If you use RheumLens in your research, please cite:

```
RheumLens: a donor-level audit benchmark for single-cell foundation-model disease representations in systemic lupus erythematosus. Code DOI and manuscript citation pending.
```

## License

MIT — see [LICENSE](LICENSE) for details.

## Affiliation

Shanghai Pudong Hospital / 上海市浦东医院
