#!/usr/bin/env python3
"""Repeated complete fold-contained CV sensitivity for RheumLens."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import io
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from fold_contained_benchmark import COHORTS, run_oof


ROOT = Path(__file__).resolve().parents[1]
FEATURES = ROOT / "data" / "donor_features"
OUT = ROOT / "results" / "evidence_package" / "repeated_cv"


def evaluate(cohort: str, seeds: np.ndarray) -> pd.DataFrame:
    z = np.load(FEATURES / f"{cohort}_donor_features.npz", allow_pickle=True)
    y = z["y"].astype(int)
    expression = z["expression"].astype(np.float64)
    scgpt = z["scgpt"].astype(np.float64)
    rows = []
    for i, seed in enumerate(seeds, 1):
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=int(seed))
        folds = list(skf.split(np.arange(len(y)), y))
        with redirect_stdout(io.StringIO()):
            predictions, _ = run_oof(expression, scgpt, y, folds)
        auc = {name: roc_auc_score(y, p) for name, p in predictions.items()}
        rows.append({
            "cohort": cohort, "repeat": i, "seed": int(seed),
            "scgpt_auc": auc["scgpt"],
            "expression_pca_auc": auc["expression_pca"],
            "donor_mean_hvg_auc": auc["hvg_pseudobulk"],
            "pca_minus_scgpt": auc["expression_pca"] - auc["scgpt"],
            "hvg_minus_scgpt": auc["hvg_pseudobulk"] - auc["scgpt"],
        })
        print(f"{cohort} repeat {i}/{len(seeds)}", flush=True)
    return pd.DataFrame(rows)


def summarize(values: pd.DataFrame) -> pd.DataFrame:
    metrics = ["scgpt_auc", "expression_pca_auc", "donor_mean_hvg_auc",
               "pca_minus_scgpt", "hvg_minus_scgpt"]
    rows = []
    for cohort, group in values.groupby("cohort"):
        for metric in metrics:
            x = group[metric]
            rows.append({
                "cohort": cohort, "metric": metric, "n_repeats": len(x),
                "mean": x.mean(), "sd": x.std(ddof=1), "median": x.median(),
                "q025": x.quantile(0.025), "q25": x.quantile(0.25),
                "q75": x.quantile(0.75), "q975": x.quantile(0.975),
                "min": x.min(), "max": x.max(),
                "fraction_gt_zero": float((x > 0).mean()) if "minus" in metric else np.nan,
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", action="append", choices=COHORTS)
    parser.add_argument("--repeats", type=int, default=30)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    seeds = np.random.SeedSequence(20260619).generate_state(args.repeats)
    selected = args.cohort or list(COHORTS)
    for cohort in selected:
        evaluate(cohort, seeds).to_csv(OUT / f"{cohort}_repeated_cv.csv", index=False)
    paths = [OUT / f"{c}_repeated_cv.csv" for c in COHORTS]
    if all(p.exists() for p in paths):
        all_values = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
        all_values.to_csv(OUT / "all_cohorts_repeated_cv.csv", index=False)
        summarize(all_values).to_csv(OUT / "all_cohorts_repeated_cv_summary.csv", index=False)


if __name__ == "__main__":
    main()
