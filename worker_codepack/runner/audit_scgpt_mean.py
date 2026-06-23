#!/usr/bin/env python3
"""Strict audit of the migrated scgpt_mean formal permutation result."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/autodl-fs/data/rheumlens")
    ap.add_argument("--frozen-seed-table", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    root = Path(args.root)
    result_path = root / "results/P6_GSE174188_v1/permutation/GSE174188_CD4/scgpt_mean/permutation_auc.csv"
    summary_path = result_path.parent / "summary.json"
    method_summary_path = root / "results/P6_GSE174188_v1/final_tables/method_summary.csv"
    frame = pd.read_csv(result_path).rename(columns={"rep": "rep_id", "seed": "rep_seed", "auc": "auc"})
    seeds = pd.read_csv(args.frozen_seed_table)
    if list(frame.columns) != ["rep_id", "rep_seed", "auc"]:
        raise RuntimeError(f"unexpected columns: {list(frame.columns)}")
    if len(frame) != 1000 or frame.rep_id.duplicated().any():
        raise RuntimeError("scgpt permutation must have exactly 1000 unique rows")
    frame = frame.sort_values("rep_id").reset_index(drop=True)
    if not np.array_equal(frame.rep_id.to_numpy(), np.arange(1000)):
        raise RuntimeError("rep_id set is not 0..999")
    if not np.isfinite(frame.auc).all():
        raise RuntimeError("nonfinite null AUC")
    if not np.array_equal(frame.rep_seed.to_numpy(np.int64), seeds.rep_seed.iloc[:1000].to_numpy(np.int64)):
        raise RuntimeError("seed table mismatch")
    summary = json.loads(summary_path.read_text())
    if summary.get("method_id") != "scgpt_mean" or summary.get("n_requested") != 1000 or summary.get("n_finite") != 1000:
        raise RuntimeError("summary identity/count mismatch")
    methods = pd.read_csv(method_summary_path)
    row = methods[(methods.cohort == "GSE174188_CD4") & (methods.method_id == "scgpt_mean")]
    if len(row) != 1:
        raise RuntimeError("cannot recover unique observed scgpt_mean AUC")
    observed = float(row.iloc[0].roc_auc)
    p = float((1 + np.sum(frame.auc.to_numpy() >= observed)) / 1001)
    if not np.isclose(observed, float(summary["observed_auc"]), atol=1e-12, rtol=0):
        raise RuntimeError("observed AUC mismatch")
    if not np.isclose(p, float(summary["empirical_p"]), atol=1e-15, rtol=0):
        raise RuntimeError("empirical P mismatch")
    report = {
        "status": "ACCEPT_RECOMMENDED",
        "method_id": "scgpt_mean",
        "n_reps": 1000,
        "observed_auc": observed,
        "empirical_p": p,
        "result_path": str(result_path),
    }
    Path(args.output).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
