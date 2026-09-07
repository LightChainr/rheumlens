# Task brief: multi-disease design screen

**Audience:** the compute agent running this repository on Huawei Cloud.
**Prepared:** 2026-09-07
**Context:** v2 of this work was rejected by PLOS Computational Biology
(`PCOMPBIOL-D-26-01856`). Two of the four reviewers, plus the academic editor,
rejected the manuscript's central claim because *all* empirical evidence came
from one disease (SLE) and because the proposed "required checks" were never
validated beyond it. This task supplies the missing evidence.

Read `docs/REVIEW_RESPONSE_CONTEXT.md` before starting if you want the full
reviewer reasoning. The short version is below.

---

## What you are producing

One table. It is the scientific deliverable:

| cohort | disease | expected | I_D | design AUC | disease AUC | p_standard | p_design_preserving |
|---|---|---|---|---|---|---|---|
| CMV_HIHA | CMV | **weak** | high | ~0.5 | ? | ? | ? |
| SLE_GSE174188 | SLE | medium | ? | ? | ? | ? | ? |
| COMBAT_COVID | COVID-19 | controlled | ? | ? | ? | ? | ? |
| COMBAT_INFLUENZA | influenza | controlled | ? | ? | ? | ? | ? |
| COMBAT_CROSS | COVID vs flu | controlled | ? | ? | ? | ? | ? |
| COVID_STEPHENSON | COVID-19 | **strong** | low | high | ? | ? | ? |
| COVID_REN | COVID-19 | **strong** | low | high | ? | ? | ? |

The "expected" column is a prediction registered **before** the run. Report what
you actually get, including disagreements with the prediction. A cohort that
contradicts the expectation is a finding, not a failure — and `CMV_HIHA`
contradicting it would be the single most important result of the run, because a
framework that flags every cohort as confounded is not falsifiable.

---

## Why the design-preserving null matters

PLOS reviewer #4 identified a real defect in the v2 method, not a presentation
problem:

> a standard complete-pipeline label permutation destroys the observed
> design–label association, so it primarily tests leakage, overfitting, and
> pipeline validity rather than whether the classifier can exploit acquisition
> structure.

`run_design_screen.py` therefore computes **two** null distributions:

* `p_standard` — permute the label freely (the v2 test).
* `p_design_preserving` — permute the label **within design strata**, so the
  design–label association survives while disease-specific molecular signal is
  destroyed.

Interpretation:

| p_standard | p_design_preserving | meaning |
|---|---|---|
| small | small | discrimination exceeds what design alone explains |
| small | **large** | apparent discrimination is attributable to design |
| large | — | pipeline problem; investigate before reporting anything |

Do not drop the second null. It is the reason this run is worth doing.

---

## Pipeline

```bash
# 0. environment
conda env create -f environment-analysis-lock.yml && conda activate rheumlens-analysis
pip install anndata pyyaml            # if absent from the lock file

# 1. download (38.2 GB, resumable — safe to re-run after an interruption)
python scripts/cohorts/fetch_cohort.py --all --data-dir /data/h5ad

# 2. build donor-level tables (one per contrast)
for C in COVID_REN COVID_STEPHENSON CMV_HIHA; do
  python scripts/cohorts/build_donor_level.py --cohort $C --data-dir /data/h5ad
done
for X in COMBAT_COVID COMBAT_INFLUENZA COMBAT_CROSS; do
  python scripts/cohorts/build_donor_level.py --cohort COMBAT --contrast $X --data-dir /data/h5ad
done

# 3. gate — do not proceed past a FAIL
python scripts/cohorts/validate_cohort.py --all

# 4. screen
python scripts/cohorts/run_design_screen.py --all --n-perm 1000
```

Outputs land in `results/multi_cohort_design_screen/`:
`design_screen.tsv`, `confounding_spectrum.tsv`, `screen_config.json`.

---

## Resources

| stage | needs | rough cost |
|---|---|---|
| download | network, 40 GB disk | ~1 h |
| build_donor_level | 16-32 GB RAM, CPU only | ~30-60 min per cohort |
| design screen | many CPU cores | ~2-6 h for 1000 permutations x 7 contrasts |

**No GPU is required.** Geneformer embeddings are *not* needed for this task —
the screen runs on pseudobulk. Adding embeddings later would only extend Check 3
(representation-level design prediction); it does not gate any result here.

Peak memory is set by the pseudobulk accumulator (`n_donor x n_gene` float64,
~50 MB) plus one 100k-cell chunk. If a cohort still exhausts memory, lower
`CHUNK` in `build_donor_level.py`.

---

## Things that will probably bite you

1. **`COVID_REN` contains lung and saliva cells.** The tissue filter in
   `cohorts/registry.yaml` (`tissue_filter: blood`) handles this, but verify the
   surviving cell count looks like blood only.
2. **`COVID_REN` mixes two 10x chemistries** (3' v3 and 5' v2). `assay` is
   therefore a genuine design variable and is added to the design block
   automatically. Do not "correct" this away — it is exactly the phenomenon
   under study.
3. **`COVID_STEPHENSON` has a third disease label**
   (`respiratory system disorder`) that is neither case nor control; those
   donors are excluded by the registry's `exclude_labels`.
4. **Raw counts location.** CELLxGENE usually puts normalised values in `.X` and
   raw counts in `.raw.X`, but not always. `pick_raw_layer()` sniffs for
   non-negative integers and prints which source it chose — **check that line in
   the log for every cohort.** A silent switch to normalised values would
   invalidate every QC covariate.
5. **Age parsing.** `development_stage` arrives as ontology labels
   (`'25-year-old stage'`, `'fifth decade stage'`, `'90 year-old and over stage'`).
   `parse_age()` handles these three shapes and returns NaN otherwise; the build
   log reports how many donors got an age. If that count is low for a cohort,
   say so — the demographic block is weaker there.
6. **`COMBAT_CROSS` has no `normal` donors by construction** (COVID vs
   influenza). That is intended.

---

## Reporting back

For each cohort report: donors kept vs registry count, cases/controls, cells
used, which raw-count source was chosen, batch-like columns detected, age
recovery count, and the screen row. Flag anything where the observed donor count
differs from the registry by more than ~10% — that means the filters removed
more than expected and the cause needs to be understood before the number is
used in a manuscript.

Do not tune thresholds, drop cohorts, or reorder the spectrum to make it look
cleaner. The spectrum's credibility rests on it being pre-registered here and
reported as observed.
