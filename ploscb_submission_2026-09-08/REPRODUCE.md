# Reproduce the 2026-09-08 PLOS Computational Biology submission

This directory is the research object for

> **Recorded study metadata predicts the diagnosis in eight of nine public
> single-cell case-control comparisons**

It is self-contained apart from the public expression matrices, which stay at
their original accessions (GEO and CZ CELLxGENE Discover) rather than being
mirrored here.

## What is here

| Path | Contents |
|---|---|
| `01_UPLOAD/` | The submitted manuscript, cover letter, PDF, 8 main and 7 supporting figures (PNG + SVG), and Supporting Tables S1–S9 |
| `02_SOURCE/analysis_scripts/` | The screen, calibration, simulation and walkthrough, plus the staging, table-building and checking tools |
| `02_SOURCE/cohort_scripts/` | Cohort fetch, donor-level build and validation |
| `02_SOURCE/result_tables/` | Every result table behind Figures 2–8, including the 19,600-cohort simulation output |
| `02_SOURCE/per_seed/` | The five per-seed screen outputs, one directory per seed, rather than only their summary |
| `02_SOURCE/figure_scripts/` | The R sources for all 15 figures |
| `verify.sh` | Re-runs both manuscript checks against the files in this directory |
| `03_QA/` | The checklist, the package SHA256 manifest, and the point-by-point review-response records |
| `cohorts/registry.yaml` | Dataset identifiers and donor counts, verified against the CELLxGENE Discover curation API on 2026-09-07 |
| `inputs/design_metadata/` | The donor-level covariate tables the screen actually fits |
| `inputs/donor_level/*/` | Donor labels and build summaries |

## Two tiers of reproduction

The expensive analyses need a 16-vCPU machine and hours; everything downstream
of their outputs rebuilds on a laptop in about a minute. The result tables in
`02_SOURCE/result_tables/` and the per-seed outputs in `02_SOURCE/per_seed/` are
the boundary between the two, so everything downstream can be checked without
re-running the analyses.

### Tier 1 — re-run the checks against the files published here

```bash
bash verify.sh
```

This assembles the layout the checkers expect from this directory in a temporary
tree and runs both of them, so the scripts published here stay byte-identical to
the ones used in the analysis:

- `check_manuscript_numbers.py` re-derives 24 numbers from the per-seed screen
  outputs, the simulation output and the walkthrough, and requires each one
  verbatim in the manuscript;
- `check_manuscript_structure.py` runs ~30 structural checks: length, citation
  order, figure and table coverage, retired vocabulary in the rendered SVGs and
  in the released tables, machine-specific paths, salted `hash()` calls, and the
  permutation counts stated in the text against `N_PERM` in the code.

Expected output is `24/24 assertions found verbatim in the manuscript` and
`all structural checks passed`. A number that drifted between the analysis and
the manuscript fails a check rather than surviving as a sentence nobody re-read.

To rebuild the tables and figures as well, rather than only check them, run
`02_SOURCE/analysis_scripts/build_all.sh` from the analysis workspace. It expects
that workspace's directory layout, not this package's, and it needs R.

### Tier 2 — re-run the analyses themselves

First rebuild the donor-level inputs, which downloads the public matrices:

```bash
python3 02_SOURCE/cohort_scripts/fetch_cohort.py --all
for c in SLE_GSE174188 COVID_REN COVID_STEPHENSON COMBAT CMV_HIHA; do
  python3 02_SOURCE/cohort_scripts/build_donor_level.py --cohort "$c"
done
python3 02_SOURCE/cohort_scripts/validate_cohort.py --all
```

Then:

```bash
# the nine-comparison screen, once per seed
for s in 20260907 20260908 20260909 20260910 20260911; do
  python3 02_SOURCE/analysis_scripts/run_design_screen.py --all --seed "$s"
done

python3 02_SOURCE/analysis_scripts/run_extended_simulation.py  # 19,600 cohorts
python3 02_SOURCE/analysis_scripts/run_calibration.py          # type-I calibration
python3 02_SOURCE/analysis_scripts/run_mediator_arm.py         # mediator arm
python3 02_SOURCE/analysis_scripts/run_walkthrough.py          # decision-tree traces
```

Set `RHEUMLENS_JOBS` to the number of worker processes. On a 16-vCPU container
with a 14.5-core quota, 14 workers with `OMP_NUM_THREADS=1` is what these runs
were done on.

## Determinism

Every seed is derived with `hashlib.blake2b`, not Python's built-in `hash()`,
which is salted per process for strings and would make string-keyed draws
irreproducible across runs. `check_manuscript_structure.py` parses each script
with `ast` and fails if a real `hash()` call reappears.

Permutation counts are fixed in the code: the complete-pipeline permutations use
1,000 draws (p floor 1/1001 ≈ 0.0010), and the V_D and design-only AUC nulls use
200 draws (p floor 1/201 ≈ 0.0050).

## A note on column names

The analysis layer and the manuscript use different vocabulary.
`02_SOURCE/analysis_scripts/si_names.py` holds the single rename map and applies
it when a supplementary table is exported, so `I_D` in the scripts is
`V_D_insample` in Table S8, and `p_standard` is `p_free`. The map also lists the
retired names, and exporting a table whose header still carries one fails the build.

## Environment

Python 3.11 with numpy, pandas, scikit-learn, scipy, pyarrow and joblib; R 4.6
with ggplot2, patchwork, svglite and ragg. Figure TIFFs additionally need
`rsvg-convert` and ImageMagick (`02_SOURCE/analysis_scripts/build_tiffs.py`).

## Integrity

`SHA256SUMS` covers every file in this directory. Verify with:

```bash
shasum -a 256 -c SHA256SUMS
```
