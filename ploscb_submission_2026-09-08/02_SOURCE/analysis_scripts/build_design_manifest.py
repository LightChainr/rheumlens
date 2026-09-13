#!/usr/bin/env python3
"""Generate the design-variable manifest from the analysis code itself.

Table S2 in the previous package was maintained by hand and listed every column
of the covariate file, classified by a naming heuristic. That inventory was wrong
in two ways that matter: it listed `case_control` under "sample quality", and it
put the processing-wave fractions there too. Neither ever entered a model matrix
- `design_matrix()` in run_design_screen.py takes a fixed whitelist - but a table
that disagrees with the code cannot be used to rule leakage out.

This script imports the real `design_matrix` and `design_strata` and reports, per
cohort and per block, the raw columns consumed and the expanded columns produced.
Anything present in the covariate file but absent from every block is listed
explicitly, with the reason.
"""
from __future__ import annotations
import importlib.util, json, sys
import os
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(os.environ.get(
    "RHEUMLENS_REPO",
    Path(__file__).resolve().parents[2] / "20_repo_restructure_20260907"))
OUT  = Path(__file__).resolve().parents[1] / "supplementary"
OUT.mkdir(exist_ok=True)

spec = importlib.util.spec_from_file_location(
    "screen", REPO / "scripts" / "cohorts" / "run_design_screen.py")
screen = importlib.util.module_from_spec(spec)
sys.modules["screen"] = screen
spec.loader.exec_module(screen)

# The whitelists, taken from the module the screen imports rather than scraped out
# of a function body with a regex. The regex version broke the moment the matrix
# construction moved into pipeline_core, which is the right failure but a slow way
# to learn it; these are module constants and there is nothing to parse.
core = screen.core
QC   = list(core.QC_COLS)
DNUM = list(core.DEMO_NUM)
DCAT = list(core.DEMO_CAT)

LABEL_COLS = {"y_true", "case_control", "disease", "diagnosis", "label"}

rows, unused_rows, checks = [], [], []
cohorts = sorted(p.name for p in (REPO/"inputs"/"donor_level").iterdir() if p.is_dir())

for coh in cohorts:
    cov_path = REPO/"inputs"/"design_metadata"/f"{coh.lower()}_donor_covariates.tsv"
    if not cov_path.exists():
        continue
    cov = pd.read_csv(cov_path, sep="\t", dtype={"donor_id": str}).set_index("donor_id")
    blocks = screen.design_matrix(cov)

    qc   = [c for c in QC   if c in cov.columns]
    dnum = [c for c in DNUM if c in cov.columns]
    dcat = [c for c in DCAT if c in cov.columns]
    bcat = [c for c in cov.columns if c.startswith("batch__")]
    if "assay" in cov.columns and cov["assay"].nunique() > 1:
        bcat = bcat + ["assay"]
    raw = {"qc": qc, "demographic": dnum+dcat, "batch": bcat,
           "all": qc+dnum+dcat+bcat}

    for blk, M in blocks.items():
        for c in raw[blk]:
            n_exp = (1 if c in qc+dnum
                     else max(cov[c].astype(str).nunique()-1, 0))
            rows.append(dict(cohort=coh, block=blk, raw_column=c,
                             kind=("numeric" if c in qc+dnum else "categorical"),
                             n_levels=(np.nan if c in qc+dnum
                                       else int(cov[c].astype(str).nunique())),
                             n_expanded_columns=int(n_exp)))
        # the manifest's expanded count must equal the matrix the code builds
        want = sum(1 if c in qc+dnum else max(cov[c].astype(str).nunique()-1,0)
                   for c in raw[blk])
        checks.append(dict(cohort=coh, block=blk, manifest_columns=int(want),
                           model_matrix_columns=int(M.shape[1]),
                           agrees=bool(want == M.shape[1])))

    used = set(raw["all"])
    for c in cov.columns:
        if c in used:
            continue
        if c in LABEL_COLS:
            why = "diagnosis label or a direct proxy - never a predictor"
        elif c.startswith("processing_fraction"):
            why = "per-donor share of a processing wave; the wave itself enters as batch__processing_wave, so the fractions would duplicate it"
        elif c in ("cells_per_donor", "mean_umi_per_cell", "mean_genes_per_cell"):
            why = f"superseded by log_{c}, which is on the whitelist; the raw and logged forms would be redundant"
        elif c == "mean_cell_pct_mito":
            why = "aggregate_pct_mito is the whitelisted mitochondrial summary; the per-cell mean would be redundant"
        elif c in ("n_libraries", "n_samples", "n_suspensions", "n_processing_cohorts"):
            why = "per-donor counts of library preparations; not on the prespecified whitelist"
        elif c in ("assay", "suspension_type"):
            why = "constant in this cohort, so it carries no information and expands to zero columns"
        elif c == "disease_label":
            why = "human-readable form of the diagnosis label - never a predictor"
        else:
            why = "not on the prespecified whitelist"
        unused_rows.append(dict(cohort=coh, column=c, in_any_design_block="no",
                                reason=why))

man = pd.DataFrame(rows).sort_values(["cohort","block","raw_column"])
unu = pd.DataFrame(unused_rows).sort_values(["cohort","column"])
chk = pd.DataFrame(checks).sort_values(["cohort","block"])

man.to_csv(OUT/"Table_S2_design_manifest.tsv", sep="\t", index=False)
unu.to_csv(OUT/"Table_S2b_columns_not_used.tsv", sep="\t", index=False)
chk.to_csv(OUT/"Table_S2c_manifest_agreement_check.tsv", sep="\t", index=False)

leak = man[man.raw_column.isin(LABEL_COLS)]
print(f"cohorts: {chk.cohort.nunique()}   block rows: {len(chk)}")
print(f"manifest vs model matrix: {int(chk.agrees.sum())}/{len(chk)} agree")
if not chk.agrees.all():
    print(chk[~chk.agrees].to_string(index=False))
print(f"label columns found in any design block: {len(leak)}"
      + ("" if len(leak)==0 else "\n"+leak.to_string(index=False)))
