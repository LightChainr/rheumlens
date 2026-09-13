"""How much does expression add over the metadata already recorded?

The screen reports the metadata-only AUC and the expression AUC side by side, which
answers "does the metadata predict the phenotype" and "does expression predict the
phenotype" but not the question a reader of a patient-level classifier actually has:
whether expression carries information the recorded metadata does not already have.

This script fits the third model. Inside each training fold it builds the metadata
matrix from the training donors (`pipeline_core.fold_design`), reduces expression to
principal components fitted on the same donors, and fits one classifier on the
concatenation. Every arm - metadata only, expression only, both - runs the identical
fixed-hyperparameter specification through the identical folds, so the differences are
differences of information and not of pipeline.

Two increments are reported per comparison:

    delta_over_metadata = AUC(D + X) - AUC(D)     what expression adds
    delta_over_expression = AUC(D + X) - AUC(X)   what the metadata adds

A permutation null for the first would answer "is the increment significant", but the
honest version of that test has to hold the metadata association fixed while permuting -
which is the collection-stratified permutation the screen already runs, on the statistic
it already runs it on. This script reports the increment as a descriptive quantity over
five split seeds and does not attach a p-value to it.

    python3 scripts/cohorts/run_incremental.py --all --seed 20260907 --out results/incremental/seed_20260907
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

try:
    from . import pipeline_core as core
except ImportError:
    import sys as _sys
    _here = str(Path(__file__).resolve().parent)
    if _here not in _sys.path:
        _sys.path.insert(0, _here)
    import pipeline_core as core

REPO = Path(__file__).resolve().parents[2]
INPUTS = REPO / "inputs"
SEED = 20260907


def joint_folds(X: np.ndarray, cov: pd.DataFrame, block: str, spec):
    """`get_fold` for [PCs of expression | metadata], both fitted on the fold.

    The expression side is reduced first, and to the same number of components the
    expression-only arm uses, so the two arms differ by the appended metadata columns
    and by nothing else. Reducing before concatenating also stops 20,000 gene columns
    from drowning a handful of one-hot columns in the standardiser.
    """
    def get(tr, te):
        Xtr, Xte = X[tr], X[te]
        if spec.n_top_var is not None and spec.n_top_var < Xtr.shape[1]:
            keep = np.argsort(Xtr.var(axis=0))[::-1][:spec.n_top_var]
            Xtr, Xte = Xtr[:, keep], Xte[:, keep]
        pre = StandardScaler().fit(Xtr)
        k = min(spec.n_pc, Xtr.shape[1], len(tr) - 2)
        pca = PCA(n_components=k, random_state=SEED).fit(pre.transform(Xtr))
        Ptr, Pte = pca.transform(pre.transform(Xtr)), pca.transform(pre.transform(Xte))
        Dtr, Dte = core.fold_design(cov, block, tr, te)
        return np.hstack([Ptr, Dtr]), np.hstack([Pte, Dte])
    return get


def one(cohort: str) -> dict | None:
    cov_path = INPUTS / "design_metadata" / f"{cohort.lower()}_donor_covariates.tsv"
    pb_path = INPUTS / "donor_level" / cohort / "donor_log1p_cpm.parquet"
    if not cov_path.exists() or not pb_path.exists():
        print(f"[{cohort}] SKIP")
        return None
    cov = pd.read_csv(cov_path, sep="\t", dtype={"donor_id": str}).set_index("donor_id")
    pb = pd.read_parquet(pb_path)
    pb.index = pb.index.astype(str)
    donors = [d for d in cov.index if d in pb.index]
    cov, pb = cov.loc[donors], pb.loc[donors]
    y = cov["y_true"].to_numpy(dtype=int)
    X = pb.to_numpy(dtype=np.float64)

    spec = core.SCREEN_FROZEN
    # The metadata-only arm applies no filter and no PCA: the matrix is already narrow,
    # and reducing it would change what "the recorded metadata" means.
    a_d = core.fold_contained_auc(core.design_folds(cov, "all"), y,
                                  core.SCREEN_METADATA_FROZEN, SEED)
    a_x = core.fold_contained_auc(core.array_folds(X), y, spec, SEED)
    a_dx = core.fold_contained_auc(joint_folds(X, cov, "all", spec), y,
                                   core.SCREEN_METADATA_FROZEN, SEED)
    rec = dict(cohort=cohort, seed=SEED, n_donor=len(y), n_case=int(y.sum()),
               auc_metadata_only=round(a_d, 4),
               auc_expression_only=round(a_x, 4),
               auc_joint=round(a_dx, 4),
               delta_over_metadata=round(a_dx - a_d, 4),
               delta_over_expression=round(a_dx - a_x, 4))
    print(f"[{cohort}] D={a_d:.4f}  X={a_x:.4f}  D+X={a_dx:.4f}  "
          f"+over_D={a_dx - a_d:+.4f}  +over_X={a_dx - a_x:+.4f}", flush=True)
    return rec


def main() -> None:
    global SEED
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    SEED = args.seed
    cohorts = (sorted(p.name for p in (INPUTS / "donor_level").iterdir() if p.is_dir())
               if args.all else args.cohorts)
    if not cohorts:
        raise SystemExit("pass --cohorts A B or --all")

    rows = [r for c in cohorts if (r := one(c)) is not None]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "incremental.tsv", sep="\t", index=False)
    (out / "config.json").write_text(json.dumps(
        {"seed": SEED, "cohorts": cohorts,
         "spec": core.SCREEN_FROZEN.describe(),
         "metadata_spec": core.SCREEN_METADATA_FROZEN.describe()}, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
