# Reproduce the release candidate

Run commands from the extracted `v2.0.0-rc3` candidate directory.

## Validate the package

```bash
python3 scripts/validate_rc_package.py .
```

This verifies required files, creator metadata, GSE135779 donor count, forbidden
stale claims, eight figure triplets, the RC3 identifiability and robustness
outputs, PLOS submission artifacts and every generated SHA256 entry. Transient build files
and Python bytecode caches are excluded from the release manifest.

## Rerun the identifiability extension

```bash
cd identifiability_extension
python3 -m pytest -q
python3 scripts/validate_upgrade.py
cd ..
```

The simulations and complete-pipeline null distributions are already packaged.
Their manifests record the commands, seeds and source inputs used for the locked
run.

## Restore GSE135779 metadata

```bash
python3 scripts/design_validity/gse135779_metadata_audit.py
```

The script uses only `inputs/public_metadata/` and
`inputs/donor_level/SLE_GSE135779/donor_labels.tsv`.

## Rerun the locked two-cohort audit

```bash
N_JOBS=6 python3 scripts/design_validity/locked_validity_audit.py
python3 scripts/design_validity/augment_gse135779_design_residualisation.py
python3 scripts/design_validity/validate_design_validity_release.py
```

The audit uses the packaged donor-level matrices and design metadata. Increase
`N_JOBS` only when sufficient local CPU capacity is available.

## Rerun strict source-only transfer

```bash
python3 scripts/strict_source_only_transfer/recompute_strict_source_only_transfer.py
python3 scripts/strict_source_only_transfer/plot_strict_source_only_transfer.py
```

Every learned transformation is fitted on source donors. The target contributes no
labels or fitted statistics.

## Rerun the PLOS robustness extension

```bash
python3 scripts/plos_robustness/run_plos_robustness.py
python3 scripts/plos_robustness/validate_plos_robustness.py
python3 scripts/plos_robustness/make_plos_robustness_figure.py
```

This block reruns nonlinear residualisation, nonlinear design-only prediction, the
batch-exposure negative control, attenuation-difference resampling, the synthetic
positive control and the confound-leakage sensitivity experiment. The validator
checks 16 prespecified structural and numerical conditions.

## Rebuild locked figures

```bash
python3 scripts/design_validity/make_locked_validity_figures.py
```

Reproduced figures are written under `figures/reproduced_locked/` so the archived
main figures are never overwritten.

## Optional inputs

The following environment variables override portable defaults:

- `RHEUMLENS_RELEASE_ROOT`
- `RHEUMLENS_INPUT_ROOT`
- `RHEUMLENS_RESULTS_ROOT`
- `RHEUMLENS_FIGURES_ROOT`
- `RHEUMLENS_GSE174188_H5AD`
- `RHEUMLENS_LEARNED_POOLING_INPUT_ROOT`
- `RHEUMLENS_CELL_EMBEDDING_ROOT`
- `RHEUMLENS_TRANSFER_PREDICTIONS`

The raw GSE174188 H5AD is required only to rebuild the donor covariate table from
cell-level metadata. Cell-level cap-500 embeddings are required only to retrain the
exploratory learned-pooling networks. Their locked downstream results are already
included.
