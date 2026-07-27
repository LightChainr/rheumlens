# Design validity in donor-level single-cell classifiers

This directory is the local `v2.0.0-rc1` release candidate for:

> **When design predicts disease: residualisation failure and generalisation gaps
> in donor-level single-cell classifiers**

It preserves the public `v1.0.2` benchmark and starts a new scientific release line.
This candidate is prepared for the `release/v2.0.0-rc1` prerelease branch. It has
not been tagged as a final release or deposited as a new Zenodo version.

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
- `inputs/donor_level/`: donor labels, pseudobulk matrices and primary frozen
  Geneformer donor embeddings for all three cohorts.
- `inputs/public_metadata/` and `inputs/design_metadata/`: compact inputs required
  for the two-cohort design audit.
- `scripts/`: analysis code grouped by scientific block.
- `docs/`: methods lock, technical audit, verified references, novelty audit,
  release gate, clean-copy reproduction record and 2026 journal strategy.
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

The release candidate contains the donor-level inputs required to rerun the locked
validity audit and strict source-only transfer. Public raw sequencing matrices remain
at their GEO accessions. Cell-level embeddings used to train the exploratory learned
pooling networks are not duplicated in Git; their locked predictions and paired
tests are included, and the optional cell-input location is configurable through
`RHEUMLENS_LEARNED_POOLING_INPUT_ROOT`.

All release-facing paths default to the extracted package and can be overridden
through the environment variables documented in `REPRODUCE.md`. Exact Geneformer
checkpoint, tokenizer, extraction, sampling, software and random-seed details are
locked in `docs/METHODS_LOCK_20260726.md`.

## Release status

- Scientific reconstruction: **passed**
- Manuscript numeric validation: **passed**
- Main figures: **locked**
- Target-journal route: **PLOS Computational Biology → GigaScience →
  Bioinformatics Advances → Applied Sciences special issue**
- Final GitHub tag and Zenodo version: **not yet published**
- Clean-copy donor-level reproduction: **passed**
- PLOS initial-submission PDF and stable supplementary workbook: **passed**

See `REPRODUCE.md` for validation and analysis commands.
