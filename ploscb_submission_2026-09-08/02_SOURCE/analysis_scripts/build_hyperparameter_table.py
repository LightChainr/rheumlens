#!/usr/bin/env python3
"""Write Table S10 from the pipeline objects that actually ran.

The hand-maintained version of this table described one configuration as though it
covered every fitted model in the paper. Pre-submission review checked it against the
code and found six disagreements - inner-fold count, forest size and depth, `max_iter`,
repeat budget, and where imputation and encoding were fitted - because the settings it
listed belonged to the lupus deep dive and the multi-cohort screen runs different ones.

So the table is generated. Every row below either comes from a `PipelineSpec` /
`ForestSpec` in `scripts/cohorts/pipeline_core.py`, which is the object the runners
import, or sits in EXTRA with the file and line that fixes it. A parameter cannot drift
from its description without this script changing what it prints.

    python3 tools/build_hyperparameter_table.py [--out supplementary/Table_S10_hyperparameters.tsv]
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPO = Path(os.environ.get("RHEUMLENS_REPO",
                           ROOT.parent / "20_repo_restructure_20260907"))
CORE_DIR = REPO / "scripts" / "cohorts"

# Settings that are not part of a classifier spec: the residualisers, the feature
# construction and the permutation budgets. Each names where it is fixed, so a
# reader can check the row against the line rather than against a claim.
EXTRA = [
    ("ridge residualiser", "used by the walkthrough, the lupus deep dive and "
     "calibration arm B", "alpha", "1.0",
     "fixed (Section 9.13)"),
    ("ridge residualiser", "used by the walkthrough, the lupus deep dive and "
     "calibration arm B", "adjustment matrix scaling",
     "standardised on the training donors of each fold",
     "fixed (Section 9.13)"),
    ("random forest residualiser", "nonlinear residualisation arm, GSE135779",
     "n_estimators / max_depth / min_samples_leaf", "96 / 4 / 3",
     "fixed (Section 9.13)"),
    ("histogram gradient boosting", "nonlinear metadata-only arm, lupus deep dive",
     "max_iter / learning_rate / max_leaf_nodes / min_samples_leaf",
     "150 / 0.05 / 7 / 5", "fixed (Section 9.5)"),
    ("pseudobulk construction", "all expression analyses", "HVG count",
     "4,000, ranked by training-fold variance",
     "fixed (Section 9.5)"),
    ("metadata matrix", "all metadata-only analyses",
     "numeric imputation / categorical encoding",
     "training-fold median / one-hot from training-fold levels, drop-first, "
     "unseen held-out levels encoded as the reference level",
     "fixed (Section 9.5)"),
    ("permutation budget", "metadata-only AUC and V_D nulls", "n_perm",
     "200 (resolution 1/201)", "fixed (Section 9.4)"),
    ("permutation budget", "free and collection-stratified tests", "n_perm",
     "1,000 (resolution 1/1001)", "fixed (Section 9.6)"),
    ("split seeds", "whole screen", "master seeds",
     "20260907-20260911, five runs; every quoted range is the min-max over them",
     "fixed (Section 9.8)"),
]

# Purpose text for each pipeline, keyed to the section and display numbering of the
# submitted manuscript (the descriptions carried by the specs predate the final order).
PURPOSE = {
    "screen_frozen": "expression AUC and both permutation nulls (Sections 3.2-3.4, Table 1)",
    "screen_metadata_frozen": "metadata-only AUC and its 200-draw permutation null (Sections 3.2-3.3, Table 1)",
    "screen_tuned": "descriptive expression AUC reported beside the fixed-hyperparameter statistic (S3 Table)",
    "screen_metadata_tuned": "descriptive metadata-only AUC, tuned logistic (S3 Table)",
    "screen_metadata_forest": "nonlinear metadata-only arm, checking that a linear model does not miss a "
                              "low-dimensional metadata-label relationship (S3 Table)",
    "walkthrough": "decision-tree walkthrough on CMV, Ren and influenza (Section 6.1, S7 Table); 20 repeats "
                   "because a single split is too noisy to compare adjusted and unadjusted arms",
    "simulation": "calibration arms (Sections 4.4-4.5, S5 Table)",
    "extended_simulation": "19,600-cohort regime sweep (Sections 4.1-4.3, Fig 3, S2 Table); a check on "
                           "estimator behaviour, not a mimic of the empirical pipeline",
    "sle_deep_dive": "GSE174188 residualisation, restriction and cross-cohort arms (Sections 5.3-5.5)",
    "sle_deep_dive_forest": "nonlinear residualiser and nonlinear metadata arm in GSE174188 (Section 5.3)",
}

ORDER = ["screen_frozen", "screen_metadata_frozen", "screen_tuned",
         "screen_metadata_tuned", "screen_metadata_forest", "walkthrough",
         "simulation", "extended_simulation", "sle_deep_dive",
         "sle_deep_dive_forest"]


def load_core():
    sys.path.insert(0, str(CORE_DIR))
    spec = importlib.util.spec_from_file_location(
        "pipeline_core", CORE_DIR / "pipeline_core.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["pipeline_core"] = m
    spec.loader.exec_module(m)
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "supplementary" /
                                         "Table_S10_hyperparameters.tsv"))
    args = ap.parse_args()

    core = load_core()
    by_name = {s.name: s for s in core.ALL_SPECS}
    missing = [n for n in ORDER if n not in by_name]
    if missing:
        sys.exit(f"pipeline_core no longer defines: {missing}")
    unlisted = [s.name for s in core.ALL_SPECS if s.name not in ORDER]
    if unlisted:
        sys.exit("pipeline_core defines specs this table does not print, so the "
                 f"table would be incomplete: {unlisted}")

    rows = []
    for name in ORDER:
        d = by_name[name].describe()
        purpose = PURPOSE[name]
        d.pop("used for")
        pipeline = d.pop("pipeline")
        for setting, value in d.items():
            rows.append({"pipeline": pipeline, "used for": purpose,
                         "setting": setting, "value": value,
                         "fixed or selected":
                             "selected inside each outer training fold"
                             if setting == "regularisation" and "inner" in value
                             else "fixed"})
    for pipeline, purpose, setting, value, fixed in EXTRA:
        rows.append({"pipeline": pipeline, "used for": purpose,
                     "setting": setting, "value": value,
                     "fixed or selected": fixed})

    df = pd.DataFrame(rows)
    out = Path(args.out)
    header = (
        "# Table S10. Every fitted pipeline in the paper, and its settings.\n"
        "# The four analyses are separate pipelines and are not interchangeable: the\n"
        "# multi-cohort screen, the decision-tree walkthrough and the calibration\n"
        "# simulation share one classifier and differ only in repeat budget, while the\n"
        "# GSE174188/GSE135779 deep dive predates them and runs a wider search.\n")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write(header)
        df.to_csv(fh, sep="\t", index=False)
    print(f"wrote {out}  ({len(df)} rows, {df.pipeline.nunique()} pipelines)")


if __name__ == "__main__":
    main()
