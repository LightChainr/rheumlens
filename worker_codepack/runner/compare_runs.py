#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load(run: Path) -> pd.DataFrame:
    path = run / "permutation_results.csv"
    frame = pd.read_csv(path)[["rep_id", "rep_seed", "auc"]].sort_values("rep_id").reset_index(drop=True)
    if frame.rep_id.duplicated().any() or not np.isfinite(frame.auc).all():
        raise RuntimeError(f"invalid result: {path}")
    return frame


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_a")
    ap.add_argument("run_b")
    args = ap.parse_args()
    a, b = load(Path(args.run_a)), load(Path(args.run_b))
    if not np.array_equal(a.rep_id.to_numpy(), b.rep_id.to_numpy()):
        raise RuntimeError("rep_id mismatch")
    if not np.array_equal(a.rep_seed.to_numpy(np.int64), b.rep_seed.to_numpy(np.int64)):
        raise RuntimeError("rep_seed mismatch")
    if not np.allclose(a.auc.to_numpy(), b.auc.to_numpy(), rtol=0, atol=1e-12):
        raise RuntimeError(f"AUC mismatch; max diff={np.max(np.abs(a.auc-b.auc))}")
    print(f"PASS_COMPARE rows={len(a)} max_abs_diff={np.max(np.abs(a.auc-b.auc))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
