"""Check that a built cohort satisfies the donor-level interface contract.

Run this after build_donor_level.py and before run_design_screen.py. It is the
gate that stops a malformed cohort from silently producing a plausible-looking
but wrong screen row.

    python scripts/cohorts/validate_cohort.py --all
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
INPUTS = REPO / "inputs"

REQUIRED_COV = ["y_true", "cells_per_donor", "mean_umi_per_cell",
                "mean_genes_per_cell", "log_cells_per_donor",
                "log_mean_umi_per_cell", "log_mean_genes_per_cell", "sex"]


def check(cohort: str) -> list[str]:
    problems: list[str] = []
    dl = INPUTS / "donor_level" / cohort
    cov_path = INPUTS / "design_metadata" / f"{cohort.lower()}_donor_covariates.tsv"

    for f in ["donor_labels.tsv", "donor_log1p_cpm.parquet"]:
        if not (dl / f).exists():
            problems.append(f"missing {dl / f}")
    if not cov_path.exists():
        problems.append(f"missing {cov_path}")
    if problems:
        return problems

    pb = pd.read_parquet(dl / "donor_log1p_cpm.parquet")
    pb.index = pb.index.astype(str)
    cov = pd.read_csv(cov_path, sep="\t", dtype={"donor_id": str}).set_index("donor_id")
    lab = pd.read_csv(dl / "donor_labels.tsv", sep="\t", dtype={"donor_id": str})

    for c in REQUIRED_COV:
        if c not in cov.columns:
            problems.append(f"covariates missing column '{c}'")

    if not set(cov.index).issubset(pb.index):
        problems.append(f"{len(set(cov.index) - set(pb.index))} donors lack pseudobulk rows")
    if len(lab) != len(cov):
        problems.append(f"donor_labels has {len(lab)} rows, covariates {len(cov)}")
    if cov.index.duplicated().any():
        problems.append("duplicate donor_id in covariates")

    if "y_true" in cov.columns:
        y = cov["y_true"]
        if set(y.unique()) - {0, 1}:
            problems.append(f"y_true is not binary: {sorted(y.unique())[:5]}")
        elif y.nunique() < 2:
            problems.append("y_true has a single class")
        else:
            minority = min(int(y.sum()), int((y == 0).sum()))
            if minority < 10:
                problems.append(f"minority class has only {minority} donors (need >=10)")

    if not np.isfinite(pb.to_numpy(dtype=np.float64)).all():
        problems.append("pseudobulk contains non-finite values")
    if (pb.to_numpy() < 0).any():
        problems.append("pseudobulk contains negative values (expected log1p CPM)")

    batch = [c for c in cov.columns if c.startswith("batch__")]
    if not batch and cov.get("assay", pd.Series(dtype=str)).nunique() <= 1:
        problems.append("WARN no batch-like design column and a single assay: "
                        "the design block will rest on QC + demographics only")
    return problems


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", nargs="*")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    cohorts = (sorted(p.name for p in (INPUTS / "donor_level").iterdir() if p.is_dir())
               if args.all else (args.cohorts or []))
    if not cohorts:
        raise SystemExit("pass --cohorts A B or --all")

    failed = False
    for c in cohorts:
        problems = check(c)
        hard = [p for p in problems if not p.startswith("WARN")]
        if not problems:
            print(f"[OK   ] {c}")
        else:
            print(f"[{'FAIL' if hard else 'WARN'}] {c}")
            for p in problems:
                print(f"         {p}")
            failed |= bool(hard)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
