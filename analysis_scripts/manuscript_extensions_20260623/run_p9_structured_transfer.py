#!/usr/bin/env python3
"""P9 extension: source-only transfer for non-neural structured methods.

This extends the Mac-side P9 table beyond mean/PCA/HVG while keeping a strict
source-only fit boundary. Target labels are used only for metrics after scores
are produced.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


ROOT = Path("/Users/lc/Documents/RheumLens")
SNAP = Path("/Volumes/Mac Data/Research/RheumLens_20260622/server_snapshot/rheumlens")
CODE = ROOT / "release/rheumlens_open_repo_v0"
OUT = ROOT / "manuscript/extension_results_20260623/P9_structured_transfer"
SEED = 20260623

METHODS = [
    "scgpt_mean",
    "kme_multiscale@scgpt",
    "moments_mean_var@scgpt",
    "quantiles@scgpt",
    "tail_fractions@scgpt",
    "red@scgpt",
    "focus_lite@scgpt",
    "donor_expression_pca",
    "donor_mean_hvg",
]
PAIRS = [("GSE285773", "GSE174188_CD4"), ("GSE174188_CD4", "GSE285773")]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def dataset_paths(cohort: str) -> dict[str, Path]:
    if cohort == "GSE285773":
        return {
            "scgpt": SNAP / "embeddings/scgpt/GSE285773_v1/GSE285773_scgpt.npz",
            "lognorm": SNAP / "data/processed/GSE285773/lognorm.npz",
        }
    if cohort == "GSE174188_CD4":
        return {
            "scgpt": SNAP / "embeddings/scgpt/GSE174188_CD4_v1/GSE174188_CD4_scgpt.npz",
            "lognorm": SNAP / "data/processed/GSE174188_CD4/lognorm.npz",
        }
    raise KeyError(cohort)


def metrics(frame: pd.DataFrame) -> dict[str, float]:
    y = frame["y_true"].to_numpy(int)
    score = frame["score"].to_numpy(float)
    p = np.clip(score, 1e-6, 1 - 1e-6)
    out = {
        "auc": float(roc_auc_score(y, score)),
        "pr_auc": float(average_precision_score(y, score)),
        "brier": float(brier_score_loss(y, p)),
        "n_donors": int(len(y)),
        "n_cases": int(y.sum()),
        "n_controls": int((1 - y).sum()),
    }
    try:
        logit = np.log(p / (1 - p)).reshape(-1, 1)
        cal = LogisticRegression(C=1e6, solver="lbfgs").fit(logit, y)
        out["calibration_intercept"] = float(cal.intercept_[0])
        out["calibration_slope"] = float(cal.coef_[0, 0])
    except Exception:
        out["calibration_intercept"] = np.nan
        out["calibration_slope"] = np.nan
    return out


def bootstrap_ci(frame: pd.DataFrame, reps: int = 5000) -> dict[str, float]:
    rng = np.random.default_rng(SEED)
    y = frame["y_true"].to_numpy(int)
    score = frame["score"].to_numpy(float)
    case = np.flatnonzero(y == 1)
    control = np.flatnonzero(y == 0)
    aucs = np.empty(reps)
    for i in range(reps):
        idx = np.r_[rng.choice(case, len(case), True), rng.choice(control, len(control), True)]
        aucs[i] = roc_auc_score(y[idx], score[idx])
    return {
        "auc_ci_low": float(np.quantile(aucs, 0.025)),
        "auc_ci_high": float(np.quantile(aucs, 0.975)),
        "auc_bootstrap_reps": reps,
    }


def query_bank() -> dict:
    # Use the default registry ISG query if the locked query bank is unavailable.
    qb = SNAP / "configs/query_bank.locked.yaml"
    if qb.exists():
        # The registry only needs a dict with queries; avoid adding a PyYAML
        # dependency by using defaults for now unless a future version needs the
        # full locked bank.
        return {}
    return {}


def main() -> None:
    sys.path.insert(0, str(CODE / "src"))
    from rheumlens.data.io import load_npz_dataset
    from rheumlens.evaluation.transfer import run_source_target
    from rheumlens.registry import build_method

    OUT.mkdir(parents=True, exist_ok=True)
    cache = {}
    rows = []
    preds = []
    failures = []
    defaults = {
        "n_hvg": 2000,
        "donor_pca_dim": 25,
        "cell_pca_dim": 50,
        "cell_pca_max_fit_cells": 200_000,
        "kme_rff_dim": 256,
        "kme_max_bandwidth_points": 4000,
        "kme_max_diagnostic_points": 512,
        "focus_topk_fraction": 0.05,
    }
    for source, target in PAIRS:
        for method_id in METHODS:
            reg = build_method(method_id, defaults=defaults, query_bank=query_bank(), seed=SEED)
            key = reg.data_key
            try:
                for cohort, data_key in [(source, key), (target, key)]:
                    ck = (cohort, data_key)
                    if ck not in cache:
                        cache[ck] = load_npz_dataset(dataset_paths(cohort)[data_key])
                expression_source = expression_target = None
                if reg.expression_key is not None:
                    for cohort in (source, target):
                        ck = (cohort, reg.expression_key)
                        if ck not in cache:
                            cache[ck] = load_npz_dataset(dataset_paths(cohort)[reg.expression_key])
                    expression_source = cache[(source, reg.expression_key)]
                    expression_target = cache[(target, reg.expression_key)]
                frame = run_source_target(
                    cache[(source, key)],
                    cache[(target, key)],
                    reg.method,
                    expression_source=expression_source,
                    expression_target=expression_target,
                    random_state=SEED,
                )
                frame["direction"] = f"{source}_to_{target}"
                frame["status"] = "SUCCESS"
                preds.append(frame)
                row = {
                    "direction": f"{source}_to_{target}",
                    "method_id": method_id,
                    "data_key": key,
                    "expression_key": reg.expression_key or "",
                    "status": "SUCCESS",
                    **metrics(frame),
                    **bootstrap_ci(frame),
                }
                rows.append(row)
                print(f"OK {row['direction']} {method_id} AUC={row['auc']:.3f}", flush=True)
            except Exception as exc:
                failures.append(
                    {
                        "direction": f"{source}_to_{target}",
                        "method_id": method_id,
                        "data_key": key,
                        "expression_key": reg.expression_key or "",
                        "status": "FAILED",
                        "error": repr(exc),
                    }
                )
                print(f"FAIL {source}->{target} {method_id}: {exc!r}", flush=True)
    pred_df = pd.concat(preds, ignore_index=True) if preds else pd.DataFrame()
    pred_df.to_csv(OUT / "EXT_P9_structured_transfer_predictions.csv", index=False)
    summary = pd.DataFrame(rows).sort_values(["direction", "auc"], ascending=[True, False])
    summary.to_csv(OUT / "EXT_P9_structured_transfer_summary.csv", index=False)
    pd.DataFrame(failures).to_csv(OUT / "EXT_P9_structured_transfer_failures.csv", index=False)
    report = "# P9 structured source-only transfer extension\n\n"
    report += "Non-neural structured methods were fit on source donors only and scored on target donors before metric calculation.\n\n"
    report += "## Successful methods\n\n"
    report += summary.to_markdown(index=False) + "\n\n"
    report += "## Failures\n\n"
    fail = pd.DataFrame(failures)
    report += (fail.to_markdown(index=False) if not fail.empty else "No failures.") + "\n"
    (OUT / "EXT_P9_STRUCTURED_TRANSFER_REPORT.md").write_text(report)
    files = [p for p in sorted(OUT.rglob("*")) if p.is_file() and p.name != "MANIFEST_SHA256.tsv"]
    pd.DataFrame([{"path": str(p.relative_to(OUT)), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in files]).to_csv(
        OUT / "MANIFEST_SHA256.tsv", sep="\t", index=False
    )
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
