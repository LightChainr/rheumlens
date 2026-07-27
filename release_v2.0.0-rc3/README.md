# Identifiability in patient-level single-cell classifiers

This directory is the versioned scientific `v2.0.0-rc3` research object for:

> **When study design predicts disease: an identifiability limit for patient-level
> single-cell classifiers**

It preserves the validated RC2 scientific object and adds nonlinear adjustment,
internal negative-control, attenuation-uncertainty, positive-control and
confound-leakage sensitivity analyses requested during PLOS-specific review. The
immutable GitHub tag identifies the code release; the stable Zenodo concept DOI is
`10.5281/zenodo.20813922`.

## Scientific contribution

Across two independently design-auditable systemic lupus erythematosus cohorts, the
candidate evaluates whether high donor-level discrimination is estimable,
attributable and transportable. It combines:

1. a design-adjusted label-information fraction and an exact non-identifiability
   construction;
2. 11,200 donor-level simulations;
3. disease prediction from recorded design variables and complete-pipeline nulls;
4. recovery of batch identity from donor representations;
5. fold-contained residualisation;
6. design-restricted analysis with size- and label-matched controls;
7. a CD4-subtype composition replication;
8. cross-batch and strict source-only cross-cohort transfer; and
9. independent-target evaluation of learned set pooling;
10. nonlinear design-only and fold-contained residualisation analyses;
11. a batch-exposure negative control; and
12. synthetic positive-control and confound-leakage sensitivity experiments.

The candidate supports a five-check validity standard: quantify design-label
estimability, test the complete pipeline under null data, report
representation-to-design predictability, pair residualisation with matched overlap
restriction, and require source-only external performance.

## Contents

- `manuscript/manuscript.md`: integrated scientific manuscript.
- `figures/main/`: eight main figures in SVG, PDF and PNG.
- `identifiability_extension/`: new simulation, null, composition, theory, tests,
  manifests and two figure groups.
- `results/locked_validity_audit/`: repeated donor-level validity audit and OOF
  predictions.
- `results/gse135779_metadata/`: restored donor-level batch, year, demographic and
  sequencing metadata with provenance.
- `results/strict_source_only_transfer/`: external transfer predictions, metrics,
  feature harmonisation and paired tests.
- `results/covariate_audit/`: repeated covariate/residualisation supporting results.
- `results/learned_pooling/`: learned and deterministic pooling predictions and
  paired tests.
- `results/plos_robustness/`: nonlinear adjustment, internal negative control,
  attenuation intervals, positive control and leakage-sensitivity outputs.
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
- The GSE135779 complete-design AUC exceeds 1,000 full-pipeline label permutations
  (`p=0.000999`) and leaves 20.3% design-adjusted label information.
- CD4 composition predicts disease at AUC 0.899-0.915 and falls to 0.689-0.709
  after wave residualisation, while matched wave restriction has a much smaller
  effect.
- Across simulation, design-adjusted label information falls to zero at perfect
  design-label collinearity.
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
- Identifiability extension checks: **18/18 passed**
- PLOS robustness checks: **16/16 passed**
- Manuscript numeric validation: **passed**
- Eight main figure groups: **validated and numbered**
- GitHub release object: **v2.0.0-rc3**
- Target-journal route: **PLOS Computational Biology → GigaScience →
  Bioinformatics Advances → Applied Sciences special issue**
- Zenodo integration: **version DOI assigned from the GitHub release**
- Clean-copy donor-level reproduction: **passed**
- PLOS Author Summary: **198 words; within the 150-200-word requirement**
- PLOS cover letter: **554 words**
- Stable supplementary workbook: **S1-S30; passed**
- RC3 manuscript rendering and complete package hash refresh: **passed**

See `REPRODUCE.md` for validation and analysis commands.
