#!/usr/bin/env python3
"""Classifier-C and donor-aggregation sensitivity for all RheumLens cohorts."""

from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import fold_contained_benchmark as fb


ROOT = Path(__file__).resolve().parents[1]
FEATURES = ROOT / "data/donor_features"
OUT = ROOT / "results/evidence_package/model_sensitivity"
CS = [0.01, 0.1, 1.0, 10.0, 100.0]


def oof_logistic(x, y, folds, c=1.0):
    p = np.full(len(y), np.nan)
    for tr, te in folds:
        m = make_pipeline(StandardScaler(), LogisticRegression(
            C=c, penalty="l2", solver="liblinear", class_weight="balanced",
            max_iter=10000, random_state=42)).fit(x[tr], y[tr])
        p[te] = m.predict_proba(x[te])[:, 1]
    return p


def paired_bootstrap(y, p_new, p_ref, n=10000):
    rng = np.random.default_rng(20260619)
    case, ctrl = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    auc, delta = np.empty(n), np.empty(n)
    for i in range(n):
        idx = np.r_[rng.choice(case, len(case), True), rng.choice(ctrl, len(ctrl), True)]
        auc[i] = roc_auc_score(y[idx], p_new[idx])
        delta[i] = auc[i] - roc_auc_score(y[idx], p_ref[idx])
    return np.quantile(auc, [.025, .975]), np.quantile(delta, [.025, .975])


def median_embeddings(cohort, donors):
    cfg = fb.COHORTS[cohort]
    if cohort == "GSE174188":
        a = ad.read_h5ad(cfg["embedding_h5ad"], backed="r")
        matrix = a.X
    else:
        a = ad.read_h5ad(cfg["h5ad"], backed="r")
        matrix = a.obsm[cfg["embedding"].split(":", 1)[1]]
    cell_donor = a.obs["donor_id"].astype(str).to_numpy()
    output = np.empty((len(donors), matrix.shape[1]), dtype=np.float32)
    for i, donor in enumerate(donors):
        idx = np.flatnonzero(cell_donor == str(donor))
        output[i] = np.median(np.asarray(matrix[idx, :], dtype=np.float32), axis=0)
    a.file.close()
    return output


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    reg_rows, agg_rows = [], []
    old_classifier = fb.classifier
    for cohort, cfg in fb.COHORTS.items():
        z = np.load(FEATURES / f"{cohort}_donor_features.npz", allow_pickle=True)
        donors, y = z["donors"].astype(str), z["y"].astype(int)
        expression = z["expression"].astype(np.float64)
        scgpt_mean = z["scgpt"].astype(np.float64)
        folds = fb.folds_from_json(Path(cfg["folds"]), donors, y)
        for c in CS:
            fb.classifier = lambda c=c: LogisticRegression(
                C=c, penalty="l2", solver="liblinear", class_weight="balanced",
                max_iter=10000, random_state=42)
            pred, _ = fb.run_oof(expression, scgpt_mean, y, folds)
            for method, p in pred.items():
                reg_rows.append({"cohort": cohort, "C": c, "method": method,
                                 "auc": roc_auc_score(y, p)})
        fb.classifier = old_classifier

        median = median_embeddings(cohort, donors)
        p_median = oof_logistic(median, y, folds)
        ref = pd.read_csv(ROOT / f"results/fold_contained/{cohort}_oof_predictions.csv")
        ref = ref.set_index("donor_id").reindex(donors)
        p_ref = ref["scgpt"].to_numpy()
        ci, dci = paired_bootstrap(y, p_median, p_ref)
        agg_rows.append({"cohort": cohort, "aggregation": "exact_cell_median",
                         "auc": roc_auc_score(y, p_median), "ci_low": ci[0], "ci_high": ci[1],
                         "median_minus_mean": roc_auc_score(y, p_median)-roc_auc_score(y, p_ref),
                         "delta_ci_low": dci[0], "delta_ci_high": dci[1]})
        pd.DataFrame({"donor_id": donors, "label": y,
                      "median_scgpt_probability": p_median,
                      "mean_scgpt_probability": p_ref}).to_csv(
            OUT / f"{cohort}_aggregation_oof.csv", index=False)
        print(cohort, agg_rows[-1], flush=True)
    fb.classifier = old_classifier
    pd.DataFrame(reg_rows).to_csv(OUT / "regularization_C_scan.csv", index=False)
    pd.DataFrame(agg_rows).to_csv(OUT / "aggregation_sensitivity_summary.csv", index=False)
    (OUT / "protocol.json").write_text(json.dumps({
        "C_values": CS, "folds": "archived fixed donor folds",
        "aggregation": "exact per-dimension median across all retained cells",
        "bootstrap": "10000 outcome-stratified paired donor resamples",
        "seed": 20260619}, indent=2))


if __name__ == "__main__":
    main()
