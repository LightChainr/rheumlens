#!/usr/bin/env python3
"""Fold-contained 15-gene ISG residual sensitivity in every RheumLens cohort."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from fold_contained_benchmark import (
    COHORTS, classifier, folds_from_json, ifn_score,
    residual_permutation_test,
)


ROOT = Path(__file__).resolve().parents[1]
FEATURES = ROOT / "data" / "donor_features"
FROZEN = ROOT / "results" / "fold_contained"
OUT = ROOT / "results" / "evidence_package" / "ifn_residual_all_cohorts"


def paired_bootstrap(y: np.ndarray, unadjusted: np.ndarray, residual: np.ndarray,
                     n_boot: int, seed: int = 20260619):
    rng = np.random.default_rng(seed)
    cases, controls = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    a = np.empty(n_boot)
    b = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.concatenate([
            rng.choice(cases, len(cases), replace=True),
            rng.choice(controls, len(controls), replace=True),
        ])
        a[i] = roc_auc_score(y[idx], unadjusted[idx])
        b[i] = roc_auc_score(y[idx], residual[idx])
    rows = []
    for name, observed, values in (
        ("scgpt", roc_auc_score(y, unadjusted), a),
        ("scgpt_ifn_residual", roc_auc_score(y, residual), b),
        ("scgpt_minus_residual", roc_auc_score(y, unadjusted) - roc_auc_score(y, residual), a - b),
    ):
        rows.append({"comparison": name, "estimate": observed,
                     "ci_low": np.quantile(values, 0.025),
                     "ci_high": np.quantile(values, 0.975),
                     "bootstrap_replicates": n_boot})
    return pd.DataFrame(rows), {"scgpt": a, "scgpt_ifn_residual": b,
                                "scgpt_minus_residual": a - b}


def run(cohort: str, n_perm: int = 1000, n_boot: int = 10000) -> pd.DataFrame:
    z = np.load(FEATURES / f"{cohort}_donor_features.npz", allow_pickle=True)
    donors = z["donors"].astype(str)
    y = z["y"].astype(int)
    expression = z["expression"].astype(np.float64)
    scgpt = z["scgpt"].astype(np.float64)
    genes = z["genes"].astype(str)
    score, available = ifn_score(expression, genes)
    folds = folds_from_json(Path(COHORTS[cohort]["folds"]), donors, y)

    probability = np.full(len(y), np.nan)
    for tr, te in folds:
        resid = LinearRegression().fit(score[tr, None], scgpt[tr])
        xtr = scgpt[tr] - resid.predict(score[tr, None])
        xte = scgpt[te] - resid.predict(score[te, None])
        model = make_pipeline(StandardScaler(), classifier()).fit(xtr, y[tr])
        probability[te] = model.predict_proba(xte)[:, 1]
    observed = float(roc_auc_score(y, probability))
    p_perm, null = residual_permutation_test(
        scgpt, score, y, folds, observed_auc=observed, n_perm=n_perm, seed=20260619
    )

    frozen = pd.read_csv(FROZEN / f"{cohort}_oof_predictions.csv")
    unadjusted = frozen["scgpt"].to_numpy()
    ci, boot = paired_bootstrap(y, unadjusted, probability, n_boot=n_boot)
    ci.insert(0, "cohort", cohort)
    ci.to_csv(OUT / f"{cohort}_ifn_residual_bootstrap_ci.csv", index=False)
    np.savez_compressed(OUT / f"{cohort}_ifn_residual_bootstrap_distributions.npz", **boot)
    np.save(OUT / f"{cohort}_ifn_residual_permutation_null.npy", null)
    pd.DataFrame({
        "cohort": cohort, "donor_id": donors, "label": y,
        "isg_score": score, "scgpt_probability": unadjusted,
        "scgpt_ifn_residual_probability": probability,
    }).to_csv(OUT / f"{cohort}_ifn_residual_oof_predictions.csv", index=False)

    r, r_p = stats.pearsonr(score, unadjusted)
    case, control = score[y == 1], score[y == 0]
    welch = stats.ttest_ind(case, control, equal_var=False)
    pooled = np.sqrt(((len(case) - 1) * case.var(ddof=1) + (len(control) - 1) * control.var(ddof=1)) /
                     (len(case) + len(control) - 2))
    d = (case.mean() - control.mean()) / pooled
    row = pd.DataFrame([{
        "cohort": cohort,
        "n_donors": len(y),
        "n_isg": len(available),
        "isg_genes": ";".join(available),
        "scgpt_auc": roc_auc_score(y, unadjusted),
        "residual_auc": observed,
        "auc_attenuation": roc_auc_score(y, unadjusted) - observed,
        "permutations": n_perm,
        "permutation_p": p_perm,
        "isg_vs_scgpt_r": r,
        "isg_vs_scgpt_p": r_p,
        "isg_sle_vs_hc_cohen_d": d,
        "isg_sle_vs_hc_welch_p": welch.pvalue,
    }])
    row.to_csv(OUT / f"{cohort}_ifn_residual_summary.csv", index=False)
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", action="append", choices=COHORTS)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--bootstrap", type=int, default=10000)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    selected = args.cohort or list(COHORTS)
    rows = [run(c, args.permutations, args.bootstrap) for c in selected]
    existing = [OUT / f"{c}_ifn_residual_summary.csv" for c in COHORTS]
    if all(p.exists() for p in existing):
        pd.concat([pd.read_csv(p) for p in existing], ignore_index=True).to_csv(
            OUT / "all_cohorts_ifn_residual_summary.csv", index=False
        )
    print(pd.concat(rows, ignore_index=True).to_string(index=False))


if __name__ == "__main__":
    main()
