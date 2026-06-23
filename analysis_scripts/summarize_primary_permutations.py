#!/usr/bin/env python3
"""Summarize existing primary-model permutation null files."""

from pathlib import Path

import pandas as pd
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "evidence_package"
rows = []
for cohort in ("GSE135779", "GSE285773", "GSE174188"):
    null = pd.read_csv(OUT / f"{cohort}_primary_model_permutation_null.csv.gz")
    pred = pd.read_csv(ROOT / "results" / "fold_contained" / f"{cohort}_oof_predictions.csv")
    for method, column in (
        ("scgpt", "scgpt"),
        ("expression_pca", "expression_pca"),
        ("donor_mean_hvg", "hvg_pseudobulk"),
    ):
        observed = roc_auc_score(pred["label"], pred[column])
        exceedances = int((null[method] >= observed).sum())
        rows.append({
            "cohort": cohort, "method": method, "observed_auc": observed,
            "n_perm": len(null), "exceedances": exceedances,
            "empirical_p": (exceedances + 1) / (len(null) + 1),
            "null_mean": null[method].mean(), "null_sd": null[method].std(ddof=1),
        })
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "all_primary_models_permutation_tests.csv", index=False)
print(summary.to_string(index=False))
