"""Two calibration studies asked for at pre-submission review.

A. How often the collection-stratified permutation rejects when the association
   that survives INSIDE the strata is not a disease effect.

   The stratified permutation conditions on collection strata only. The declared
   metadata matrix also holds sample-quality and demographic columns. If one of
   those is associated with the phenotype within a stratum, and drives expression,
   then a classifier can score above the stratified null with no disease effect
   present at all. This arm measures how often that happens, as a function of how
   strong the within-stratum association is.

   Note what this rejection rate is and is not. The generating model is
   QC -> Y and QC -> X, so conditioning on the collection strata leaves X and Y
   dependent: the conditional-exchangeability null the stratified test states is
   FALSE for gamma > 0, and rejecting it is not a type-I error. Only the gamma = 0
   row is a calibration measurement. The rows above it measure something more
   useful for a reader: passing this test does not certify disease biology,
   because a non-disease pathway outside the conditioning set can carry the
   classifier past the null.

B. Residualisation loss when the metadata matrix and the phenotype are
   independent by construction, as a function of matrix width.

   The CMV comparison loses AUC to residualisation while showing no detectable
   metadata-phenotype association. This arm tests the width attribution instead of
   asserting it: D is generated independently of Y, so any loss is a property of
   the adjustment, not of confounding.

   The two arms of B differ in the residualisation step and in nothing else.
   An earlier version fitted the unadjusted arm's PCA on all donors and the
   residualised arm's PCA inside the training fold, so the measured gap contained
   the change of representation as well as the cost of adjustment; the reported
   loss was negative at small widths, which is the signature of that asymmetry.

Both arms run `pipeline_core.SIMULATION`, which is the screen's frozen pipeline
with the same repeat budget, so "the same pipeline as the empirical screen" is
enforced by construction rather than asserted here.
"""
from __future__ import annotations
import hashlib, itertools, os, time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# Imported by name off sys.path, not exec'd from a path: joblib workers
# re-import the module the specs live in, and a spec-loaded module has no name
# they can import.
import sys as _sys
_CORE_DIR = str(Path(os.environ.get(
    "RHEUMLENS_REPO",
    Path(__file__).resolve().parents[2] / "20_repo_restructure_20260907"))
    / "scripts" / "cohorts")
if _CORE_DIR not in _sys.path:
    _sys.path.insert(0, _CORE_DIR)
import pipeline_core as core

OUT = Path(os.environ.get("RHEUMLENS_SIM_OUT",
                          Path(__file__).resolve().parent / "results"))
OUT.mkdir(exist_ok=True)

MASTER_SEED = 20260908
N_PERM      = 200   # resolution 1/201, ample for a 0.05 rejection rate
N_REP       = 300
C_FIXED     = 1.0
PERM_PCS    = 50
N_JOBS      = int(os.environ.get("RHEUMLENS_JOBS", "-1"))


# Level sizes of batch__batch_id x batch__pool_id in the CMV cohort, read off the
# released covariate table. Used to reproduce the shape of a real design matrix.
CMV_LEVEL_SIZES: list[int] = []


def load_cmv_level_sizes(path: Path | None = None) -> list[int]:
    """Locate the CMV covariate table without assuming any absolute path.

    Order: CMV_COVARIATES env var, then a copy sitting beside this script's
    parent, then the release repo's inputs directory.
    """
    global CMV_LEVEL_SIZES
    here = Path(__file__).resolve()
    for cand in ([Path(path)] if path else []) + [
            Path(os.environ["CMV_COVARIATES"]) if "CMV_COVARIATES" in os.environ else None,
            here.parent.parent / "cmv_hiha_donor_covariates.tsv",
            here.parents[2] / "20_repo_restructure_20260907" / "inputs" /
            "design_metadata" / "cmv_hiha_donor_covariates.tsv"]:
        if cand is not None and cand.exists():
            path = cand
            break
    else:
        return []
    if not path.exists():
        return []
    cov = pd.read_csv(path, sep="\t", dtype={"donor_id": str})
    key = cov[["batch__batch_id", "batch__pool_id"]].astype(str).agg("|".join, axis=1)
    CMV_LEVEL_SIZES = sorted(key.value_counts().tolist(), reverse=True)
    return CMV_LEVEL_SIZES


def cell_seed(*parts) -> int:
    key = ("|".join(str(x) for x in (MASTER_SEED,) + parts)).encode()
    return int.from_bytes(hashlib.blake2b(key, digest_size=4).digest(), "big") % (2**31)


def frozen_auc(X, yy, seed, residualise=None):
    """The screen's frozen pipeline, fitted entirely inside the training fold.

    `X` is the raw feature matrix, not a precomputed representation: the variance
    filter, the standardiser and the PCA are refitted on the training donors of
    every split, exactly as `run_design_screen.py` does. `residualise(tr, te)`
    returns the adjustment matrix for a split; the unadjusted and adjusted arms of
    arm B call this same function and differ only in that argument.
    """
    if len(np.unique(yy)) < 2 or min(np.bincount(yy)) < 3:
        return np.nan
    return core.fold_contained_auc(core.array_folds(X), yy, core.SIMULATION,
                                   seed, residualise=residualise)


# ------------------------------------------------------------------- arm A
def arm_a(gamma, rep, n=200, p=100, n_site=8):
    """No disease effect on expression. A QC variable predicts Y within site."""
    seed = cell_seed("A", round(gamma, 3), rep)
    rng = np.random.default_rng(seed)

    site = rng.integers(0, n_site, n)
    site_eff = rng.standard_normal(n_site)

    # qc is associated with the diagnosis WITHIN each site, at strength gamma
    qc = rng.standard_normal(n)
    lin = gamma * qc                      # no site term: the association is
    pr = 1.0 / (1.0 + np.exp(-lin))       # purely within-stratum
    y = (rng.random(n) < pr).astype(int)
    if min(np.bincount(y, minlength=2)) < 10:
        return None

    u_site = rng.standard_normal(p); u_site /= np.linalg.norm(u_site)
    u_qc   = rng.standard_normal(p); u_qc   /= np.linalg.norm(u_qc)

    # Expression depends on the site and on qc. It does NOT depend on y.
    X = (np.outer(site_eff[site], u_site) * 1.2
         + np.outer(qc, u_qc) * 1.2
         + rng.standard_normal((n, p)))

    # The raw matrix goes in, not a representation fitted on all n donors: the
    # standardiser and the PCA are refitted on the training donors of every fold,
    # for the observed statistic and for each permuted label alike.
    obs = frozen_auc(X, y, seed)
    if not np.isfinite(obs):
        return None

    rng2 = np.random.default_rng(seed + 1)
    hits_u = hits_c = 0
    usable = 0
    for _ in range(N_PERM):
        # unrestricted
        yu = rng2.permutation(y)
        a = frozen_auc(X, yu, seed)
        if np.isfinite(a) and a >= obs:
            hits_u += 1
        # within collection strata
        yc = y.copy()
        for s in np.unique(site):
            m = site == s
            yc[m] = rng2.permutation(y[m])
        if len(np.unique(yc)) < 2:
            continue
        usable += 1
        a = frozen_auc(X, yc, seed)
        if np.isfinite(a) and a >= obs:
            hits_c += 1

    # how much of the y-qc association survives inside strata
    within = float(np.mean([abs(np.corrcoef(qc[site == s], y[site == s])[0, 1])
                            for s in np.unique(site)
                            if len(np.unique(y[site == s])) > 1]))
    return dict(arm="A_within_stratum_non_disease_association",
                gamma=gamma, rep=rep, seed=seed,
                observed_auc=obs, n_case=int(y.sum()),
                within_stratum_abs_corr=within,
                p_unstratified=(hits_u + 1) / (N_PERM + 1),
                p_collection_preserving=((hits_c + 1) / (usable + 1)
                                         if usable else np.nan),
                n_usable=usable)


# ------------------------------------------------------------------- arm B
def arm_b(n_design_col, rep, n=108, p=4000, bio=5.5, layout="balanced"):
    """D generated independently of Y, so any loss is the adjustment, not confounding.

    Read the result as a property of ridge residualisation on a wide, ragged
    matrix, not as a measurement of how much confounding a matrix of that shape
    carries. Column count is not effective dimension: the real CMV collection
    matrix has 80 columns but a centred rank of 45, and the `cmv_ragged` layout is
    here because equal blocks of the same width span a different subspace.
    """
    seed = cell_seed("B", layout, n_design_col, rep)
    rng = np.random.default_rng(seed)

    y = np.zeros(n, dtype=int); y[: n // 2] = 1
    y = rng.permutation(y)

    # a real biological effect, so there is something for residualisation to lose
    # bio is calibrated so the unadjusted AUC lands near the CMV comparison's
    # expression AUC, otherwise a saturated 1.000 hides any loss.
    u_bio = rng.standard_normal(p); u_bio /= np.linalg.norm(u_bio)
    X = np.outer(y - y.mean(), u_bio) * bio + rng.standard_normal((n, p))

    # Design: one-hot blocks, assigned independently of y by construction.
    # Levels are dealt round-robin and then shuffled, so the requested width is
    # actually realised - drawing labels at random leaves empty levels and the
    # matrix comes out narrower than asked for.
    if layout == "cmv_ragged":
        # The real CMV design is not a set of equal blocks: batch_id has 36 levels
        # and pool_id 46, for 108 donors, so many levels hold one or two people.
        # A matrix of that shape spans far more of the sample than equal blocks of
        # the same column count. Here the observed level SIZES are reused and the
        # donors reshuffled into them, which keeps D independent of Y.
        sizes = CMV_LEVEL_SIZES
        lab = np.concatenate([[i] * k for i, k in enumerate(sizes)])[:n]
        if len(lab) < n:
            lab = np.concatenate([lab, np.full(n - len(lab), len(sizes))])
        lab = rng.permutation(lab)
    else:
        lab = rng.permutation(np.arange(n) % (n_design_col + 1))
    D = pd.get_dummies(pd.Series(lab).astype(str), drop_first=True).to_numpy(float)
    if D.shape[1] == 0:
        return None
    # let the design touch expression, as a real batch would
    for j in range(D.shape[1]):
        u = rng.standard_normal(p); u /= np.linalg.norm(u)
        X += np.outer(D[:, j], u) * 4.0

    # The two arms are the SAME call with one argument added. Folds, variance
    # filter, standardiser, PCA, classifier and scoring rule are shared, and each
    # is fitted on the training donors of the split; the only difference is
    # whether the training-fold ridge fit on D is subtracted first. That is what
    # makes the difference readable as the cost of adjustment.
    unadj = frozen_auc(X, y, seed)
    resid = frozen_auc(X, y, seed, residualise=core.array_folds(D))

    # metadata-only AUC, to confirm D really is independent of Y here
    d_auc = core.fold_contained_auc(core.array_folds(D), y,
                                    core.SCREEN_METADATA_FROZEN, seed) \
        if D.shape[1] >= 2 else np.nan

    return dict(arm="B_residualisation_width", layout=layout,
                n_design_col=int(D.shape[1]), rep=rep,
                seed=seed, n_donor=n, design_only_auc=d_auc,
                unadjusted_auc=unadj, residualised_auc=resid,
                residualisation_loss=unadj - resid)


def _aggregate_a(da: pd.DataFrame) -> pd.DataFrame:
    """Rejection rate by within-stratum association strength.

    Only the gamma = 0 row is a calibration measurement. Above it the
    conditional-exchangeability null is false by construction, so the rate is the
    chance of detecting a non-disease association that the conditioning set does
    not cover - not a type-I error rate. The column names say so.
    """
    return (da.assign(rej_u=da.p_unstratified <= 0.05,
                      rej_c=da.p_collection_preserving <= 0.05)
              .groupby("gamma")
              .agg(n=("rep", "size"),
                   mean_within_corr=("within_stratum_abs_corr", "mean"),
                   mean_observed_auc=("observed_auc", "mean"),
                   reject_unstratified=("rej_u", "mean"),
                   reject_collection_preserving=("rej_c", "mean"))
              .round(4))


def _aggregate_b(db: pd.DataFrame) -> pd.DataFrame:
    return (db.groupby(["layout", "n_design_col"])
              .agg(n=("rep", "size"),
                   design_only_auc=("design_only_auc", "mean"),
                   unadjusted=("unadjusted_auc", "mean"),
                   residualised=("residualised_auc", "mean"),
                   loss_mean=("residualisation_loss", "mean"),
                   loss_p025=("residualisation_loss", lambda s: s.quantile(.025)),
                   loss_p975=("residualisation_loss", lambda s: s.quantile(.975)))
              .round(4))


def _shard(cells, shard):
    """Take every n-th cell, so k machines cover the grid without coordination.

    Each cell is seeded from its own coordinates (`cell_seed`), so a sharded
    sweep and a single-process one produce the same rows in a different order,
    and the aggregate is identical after sorting.
    """
    if not shard:
        return cells
    i, n = (int(x) for x in shard.split("/"))
    return cells[i::n]


def main(shard: str = "", arms: str = "AB"):
    t0 = time.time()
    gammas = [0.0, 0.3, 0.6, 1.0, 1.5, 2.0]
    cells_a = _shard(list(itertools.product(gammas, range(N_REP))), shard)
    if "A" not in arms:
        cells_a = []
    print(f"arm A: {len(cells_a)} runs", flush=True)
    ra = Parallel(n_jobs=N_JOBS, verbose=5, batch_size=8)(
        delayed(arm_a)(g, r) for g, r in cells_a) if cells_a else []
    da = pd.DataFrame([x for x in ra if x])
    da.to_csv(OUT / "calibration_arm_a_raw.tsv", sep="\t", index=False)
    if shard or "A" not in arms:
        print("sharded run: raw rows only, aggregate with --merge-from", flush=True)
    if len(da):
        agg_a = _aggregate_a(da)
        agg_a.to_csv(OUT / "calibration_arm_a_summary.tsv", sep="\t")
        print(agg_a.to_string(), flush=True)

    sizes = load_cmv_level_sizes()
    print(f"CMV level sizes: {len(sizes)} levels, "
          f"largest {max(sizes) if sizes else 0}, singletons "
          f"{sum(1 for k in sizes if k == 1) if sizes else 0}", flush=True)
    widths = [1, 2, 4, 8, 16, 32, 48, 64, 80, 92]
    cells_b = [("balanced", w, r) for w in widths for r in range(60)]
    if sizes:
        cells_b += [("cmv_ragged", len(sizes) - 1, r) for r in range(60)]
    cells_b = _shard(cells_b, shard)
    if "B" not in arms:
        cells_b = []
    print(f"arm B: {len(cells_b)} runs", flush=True)
    rb = Parallel(n_jobs=N_JOBS, verbose=5, batch_size=4)(
        delayed(arm_b)(w, r, layout=lay) for lay, w, r in cells_b) if cells_b else []
    db = pd.DataFrame([x for x in rb if x])
    db.to_csv(OUT / "calibration_arm_b_raw.tsv", sep="\t", index=False)
    if len(db):
        aggb = _aggregate_b(db)
        aggb.to_csv(OUT / "calibration_arm_b_summary.tsv", sep="\t")
        print(aggb.to_string(), flush=True)
    print(f"done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


def merge(raw_dirs):
    """Concatenate sharded raw files and write the two summaries.

    The summary is a mean over replicates, so it must be formed once over the
    whole grid; forming it per shard and averaging would weight unequal shards
    equally.
    """
    import glob
    for arm, agg in (("a", _aggregate_a), ("b", _aggregate_b)):
        frames = []
        for d in raw_dirs:
            f = Path(d) / f"calibration_arm_{arm}_raw.tsv"
            if f.exists() and f.stat().st_size > 1:
                t = pd.read_csv(f, sep="\t")
                if len(t):
                    frames.append(t)
        if not frames:
            print(f"arm {arm}: nothing to merge")
            continue
        df = pd.concat(frames, ignore_index=True).drop_duplicates(
            subset=["seed"]).sort_values("seed").reset_index(drop=True)
        df.to_csv(OUT / f"calibration_arm_{arm}_raw.tsv", sep="\t", index=False)
        out = agg(df)
        out.to_csv(OUT / f"calibration_arm_{arm}_summary.tsv", sep="\t")
        print(f"=== arm {arm}: {len(df)} rows ===")
        print(out.to_string(), flush=True)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", default="", metavar="i/n",
                    help="run every n-th cell; results are identical to an "
                         "unsharded run because every cell seeds itself")
    ap.add_argument("--arms", default="AB")
    ap.add_argument("--merge-from", nargs="*", default=None,
                    help="directories of sharded raw files to aggregate")
    a = ap.parse_args()
    if a.merge_from:
        merge(a.merge_from)
    else:
        main(shard=a.shard, arms=a.arms)
