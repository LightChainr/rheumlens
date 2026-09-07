# rheumlens — design validity for patient-level single-cell classifiers

Do high donor-level AUCs in single-cell disease classifiers measure biology, or
do they measure how the cohort was collected?

This repository holds the analysis behind that question: the locked v2 study
(three SLE cohorts) and the v3 multi-disease extension now in progress.

## Repository layout

```
scripts/                 analysis pipeline
  design_validity/       locked v2 validity audit
  covariate_audit/       matched covariate analysis
  learned_pooling/       pooling comparison + source-only transfer protocol
  plos_robustness/       robustness and sensitivity analyses
  strict_source_only_transfer/
  cohorts/               v3 multi-disease extension  <-- new
inputs/                  donor-level tables, design metadata, public metadata
results/                 locked outputs
figures/                 manuscript figures
identifiability_extension/  simulation study
cohorts/registry.yaml    v3 cohort definitions
docs/                    methods notes and review context
manuscript/, submission/ v2 manuscript and the PLOS submission package
```

Large matrices are not tracked. `inputs/donor_level/` holds donor-level
summaries only; cohort `.h5ad` files are downloaded into `data/h5ad/` and are
git-ignored.

## Status

**v2 (locked).** Three SLE cohorts. Submitted to PLOS Computational Biology as
`PCOMPBIOL-D-26-01856` and rejected on 2026-09-03 after four reviews. The
analysis, results and figures are preserved; the manuscript's central claim is
not. See `docs/REVIEW_RESPONSE_CONTEXT.md`.

**v3 (in progress).** Two changes:

1. *A design-confounding spectrum instead of one disease.* Five cohorts spanning
   four diseases, chosen so that expected design confounding runs from weak
   (a standardised CMV atlas) through medium (SLE) to strong (two COVID-19
   cohorts), plus COMBAT, where COVID-19, influenza and controls share one
   acquisition process — so the disease contrast can change while acquisition
   is held approximately constant.
2. *A design-preserving null.* Reviewer #4 showed that standard label
   permutation destroys the design–label association and therefore tests
   leakage rather than design exploitation. The screen now permutes labels
   within design strata as well.

To run v3, see [`TASK_BRIEF.md`](TASK_BRIEF.md).

## Reproducing v2

`REPRODUCE.md` documents the locked pipeline; `SHA256SUMS` and
`RELEASE_MANIFEST.json` cover the archived inputs and results. Environments are
pinned in `environment-analysis-lock.yml` (CPU analysis) and
`environment-geneformer-lock.yml` (embedding generation, GPU).

## History

This repository previously carried four full release snapshots
(`release_v1.0.0`, `release_v1.0.2`, `release_v2.0.0-rc1`, `release_v2.0.0-rc3`)
plus the v1-era working tree, totalling ~650 MB. The v2.0.0-rc3 content is now
the repository root and the rest was removed on the `restructure/multi-disease-v3`
branch. **Nothing was rewritten out of git history** — every deleted path is
recoverable from its tag (`v1.0.0`, `v1.0.2`, `v2.0.0-rc3`, …), and the
published Zenodo archives are unaffected.

## Citation

See `CITATION.cff`. Claim scope is stated in `CLAIM_BOUNDARY.md`.
