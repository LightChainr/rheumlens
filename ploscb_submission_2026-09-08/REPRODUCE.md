# Reproduce the PLOS Computational Biology submission

This directory is the research object for

> **Recorded metadata predicts the phenotype label in eight of nine public
> single-cell cohort comparisons**

It is self-contained apart from the public expression matrices, which stay at
their original accessions (GEO and CZ CELLxGENE Discover) rather than being
mirrored here.

## What is here

| Path | Contents |
|---|---|
| `01_UPLOAD/` | The submitted manuscript, cover letter, response to the previous review, PDF, 8 main and 7 supporting figures (PNG + SVG, plus 300 dpi TIFF), and Supporting Tables S1–S10 |
| `02_SOURCE/analysis_scripts/` | The screen, calibration, simulation and walkthrough, plus the staging, table-building and checking tools |
| `02_SOURCE/cohort_scripts/` | Cohort fetch, donor-level build and validation |
| `02_SOURCE/result_tables/` | Every result table behind Figures 2–8, including the 19,600-cohort simulation output |
| `02_SOURCE/per_seed/` | The five per-seed screen outputs, one directory per seed, rather than only their summary |
| `02_SOURCE/per_seed_incremental/` | The five per-seed outputs of `run_incremental.py` behind Table S5 and the "+ over metadata" column of Table 1 |
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

`verify.sh` checks the released files against the manuscript. It is not a full
rebuild: it does not redraw the figures or regenerate the supporting tables. That is
`02_SOURCE/analysis_scripts/build_all.sh`, which needs R and the analysis workspace's
directory layout rather than this package's, because the scripts staged here keep
their workspace-relative paths. Reproducing the figures from scratch therefore means
reconstructing that workspace, not running one command in this directory.

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
# the nine-comparison screen, once per seed, each into its own directory.
# --seed on its own now implies results/screen/seed_<seed>/; name --out anyway,
# because without it older copies of the script wrote every seed to the same
# design_screen.tsv and four of the five runs were silently overwritten.
# --workers parallelises the permutation stage only. Every permuted label is
# drawn first, in the order the single-process version drew it, so the output is
# byte-identical to --workers 1; only the wall time changes.
for s in 20260907 20260908 20260909 20260910 20260911; do
  OMP_NUM_THREADS=1 python3 02_SOURCE/analysis_scripts/run_design_screen.py \
      --all --seed "$s" --workers 14 --out "results/screen/seed_$s"
done

# A (cohort, seed) job is independent of every other, so the sweep also splits
# across machines. Run the 45 jobs wherever, one directory each, then join them
# without re-running anything - --merge-from applies the same spectrum and
# degeneracy code a single-process run applies, so a merged sweep cannot drift
# from an unsplit one:
#   python3 run_design_screen.py --merge-from out/seed_$s/*/design_screen.tsv \
#       --seed "$s" --out "results/screen/seed_$s"

# The increment of expression over the recorded metadata (Table S5), same seeds:
for s in 20260907 20260908 20260909 20260910 20260911; do
  python3 02_SOURCE/analysis_scripts/run_incremental.py --all --seed "$s" \
      --out "results/incremental/seed_$s"
done

# A metadata-only refresh: recompute only the metadata columns of an existing
# per-seed output and carry the expression columns across unchanged.
#   python3 run_design_screen.py --all --seed "$s" --refresh-design \
#       --out "results/screen/seed_$s"

python3 02_SOURCE/analysis_scripts/run_extended_simulation.py  # 19,600 cohorts
# Calibration arms A and B. Cells split across machines with --shard i/n and are
# joined with --merge-from, which aggregates exactly as an unsplit run does:
#   python3 run_calibration.py --shard 0/4 ...   (one per machine)
#   python3 run_calibration.py --merge-from raw_0 raw_1 raw_2 raw_3
python3 02_SOURCE/analysis_scripts/run_calibration.py
python3 02_SOURCE/analysis_scripts/run_mediator_arm.py         # mediator arm
python3 02_SOURCE/analysis_scripts/run_walkthrough.py          # decision-tree traces
```

Set `RHEUMLENS_JOBS` to the number of worker processes. On a 16-vCPU container
with a 14.5-core quota, 14 workers with `OMP_NUM_THREADS=1` is what these runs
were done on.

## One implementation for every fitted pipeline

`02_SOURCE/analysis_scripts/pipeline_core.py` holds every classifier and every
preprocessing step the paper fits: the screen, the increment analysis, the
decision-tree walkthrough and the calibration simulation all import it. Imputation,
one-hot encoding, the top-variance gene filter, standardisation, PCA and ridge
residualisation are fitted on the training donors of each fold. An unadjusted and a
residualised arm are the same call with one argument added. Each pipeline's settings
live in a `PipelineSpec`/`ForestSpec` object, and Table S8 is generated from those
objects by `build_hyperparameter_table.py`, which fails if a spec is added without
being printed. Scripts import the module by name, so keep it in the same directory.

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
`V_D_insample` in Table S9, and `p_standard` is `p_free`. The map also lists the
retired names, and exporting a table whose header still carries one fails the build.

## Environment

The screen runs in the pinned analysis environment, `environment-analysis-lock.yml`
(conda-forge: Python 3.13.7, numpy 2.4.6, scipy 1.18.0, scikit-learn 1.9.0,
pandas 3.0.3, pyarrow 21.0.0, joblib 1.5.3):

```bash
micromamba create -y -p ./env -f environment-analysis-lock.yml
```

That environment is not merely declared, it is checked. Before the released
screen was regenerated, the *previous* version of `run_design_screen.py` was run
inside it on a linux-aarch64 container and reproduced the previously released
`design_screen.tsv` byte for byte. Every difference between the old and the new
result tables is therefore attributable to the code change described in
Section 9.6 of the manuscript, not to the machine it ran on.

Downstream rebuilding needs R 4.6 with ggplot2, patchwork, svglite and ragg.
Figure TIFFs additionally need `rsvg-convert` and ImageMagick
(`02_SOURCE/analysis_scripts/build_tiffs.py`).

## Integrity

`SHA256SUMS` covers every file in this directory. Verify with:

```bash
shasum -a 256 -c SHA256SUMS
```
