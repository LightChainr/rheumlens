# Design validity in donor-level single-cell classifiers

This directory is the local `v2.0.0-rc1` release candidate for:

> **When design predicts disease: residualisation failure and generalisation gaps
> in donor-level single-cell classifiers**

It preserves the public `v1.0.2` benchmark and starts a new scientific release line.
This candidate has **not** been tagged, pushed to GitHub or deposited in Zenodo.

## Scientific contribution

Across two independently design-auditable systemic lupus erythematosus cohorts, the
release evaluates whether high donor-level discrimination is attributable and
transportable. It combines:

1. disease prediction from recorded design variables;
2. recovery of batch identity from donor representations;
3. fold-contained residualisation;
4. design-restricted analysis with size- and label-matched controls;
5. cross-batch and strict source-only cross-cohort transfer; and
6. independent-target evaluation of learned set pooling.

The release supports a four-check validity standard: report design-only disease AUC,
representation-to-design predictability, design-restricted performance with matched
controls, and source-only external performance.

## Contents

- `manuscript/manuscript.md`: integrated scientific manuscript.
- `figures/main/`: six main figures in SVG, PDF and PNG.
- `results/locked_validity_audit/`: repeated donor-level validity audit and OOF
  predictions.
- `results/gse135779_metadata/`: restored donor-level batch, year, demographic and
  sequencing metadata with provenance.
- `results/strict_source_only_transfer/`: external transfer predictions, metrics,
  feature harmonisation and paired tests.
- `results/covariate_audit/`: repeated covariate/residualisation supporting results.
- `results/learned_pooling/`: learned and deterministic pooling predictions and
  paired tests.
- `scripts/`: analysis code grouped by scientific block.
- `docs/`: methods lock, technical audit, verified references, novelty audit,
  release gate and 2026 journal strategy.
- `submission/plos_computational_biology/`: journal-specific working documents;
  these are not part of the journal-neutral scientific claim.

## Locked headline values

- Complete measured design predicts disease at AUC 0.953 in GSE174188 and 0.952
  in GSE135779.
- Donor representations recover recorded batch identity at AUC up to 0.9997.
- Restriction to observed design overlap produces small matched changes, whereas
  residualisation can remove 0.230–0.430 AUC and reduce performance to chance.
- Cross-batch and cross-cohort transfer expose the reproducible generalisation gap.
- Learned pooling does not repair the independent-target gap.

Exact values and validation hashes are recorded in
`results/locked_validity_audit/release_validation.json`.

## Reproducibility boundary

The release candidate contains compact donor-level outputs, metadata, predictions,
source tables and analysis scripts. Public raw sequencing matrices remain at their
GEO accessions. Large intermediate cell embeddings are not duplicated here. Exact
Geneformer checkpoint, tokenizer, extraction, sampling, software and random-seed
details are locked in `docs/METHODS_LOCK_20260726.md`.

Some scripts retain the original project-relative paths used to generate the locked
outputs. They are included for provenance. Before public `v2.0.0`, paths must be
parameterised or wrapped by one repository-level runner, and the package must pass
an independent clean-directory execution check.

## Release status

- Scientific reconstruction: **passed**
- Manuscript numeric validation: **passed**
- Main figures: **locked**
- Target-journal route: **PLOS Computational Biology → GigaScience →
  Bioinformatics Advances → Applied Sciences special issue**
- Public GitHub/Zenodo release: **not yet published**

