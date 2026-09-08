"""Two calibration studies asked for at pre-submission review.

A. Type-I error of the restricted (collection-preserving) permutation when a
   confounder survives INSIDE the strata and there is no biological signal.

   The restricted permutation conditions on collection strata only. The declared
   design matrix also holds sample-quality and demographic columns. If one of
   those is associated with the diagnosis within a stratum, and drives expression,
   then a classifier can score above the restricted null with no disease effect
   present at all. This arm measures how often that happens, as a function of how
   strong the within-stratum association is. It is the direct test of what the
   test does NOT protect against.

B. Residualisation loss when design and diagnosis are independent by
   construction, as a function of design-matrix width.

   The CMV cohort loses 0.284 AUC to residualisation while showing no
   design-diagnosis association. The manuscript attributed that to the width of
   the design matrix (many one-hot columns, few donors). This arm tests the
   attribution instead of asserting it: D is generated independently of Y, so any
   loss is width, not confounding.

Both arms use the same frozen pipeline as the empirical screen.
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


def frozen_auc(Xp, yy, seed):
    """The screen's frozen pipeline: fixed PCA, fixed C, one cross-fitting repeat."""
    if len(np.unique(yy)) < 2 or min(np.bincount(yy)) < 3:
        return np.nan
    oof = np.zeros(len(yy))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(Xp, yy):
        sc = StandardScaler().fit(Xp[tr])
        m = LogisticRegression(C=C_FIXED, max_iter=2000).fit(sc.transform(Xp[tr]), yy[tr])
        oof[te] = m.decision_function(sc.transform(Xp[te]))
    return roc_auc_score(yy, oof)


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

    n_pc = int(min(PERM_PCS, p, n - 2))
    Xp = PCA(n_components=n_pc, random_state=seed).fit_transform(
        StandardScaler().fit_transform(X))

    obs = frozen_auc(Xp, y, seed)
    if not np.isfinite(obs):
        return None

    rng2 = np.random.default_rng(seed + 1)
    hits_u = hits_c = 0
    usable = 0
    for _ in range(N_PERM):
        # unstratified
        yu = rng2.permutation(y)
        a = frozen_auc(Xp, yu, seed)
        if np.isfinite(a) and a >= obs:
            hits_u += 1
        # restricted to collection strata
        yc = y.copy()
        for s in np.unique(site):
            m = site == s
            yc[m] = rng2.permutation(y[m])
        if len(np.unique(yc)) < 2:
            continue
        usable += 1
        a = frozen_auc(Xp, yc, seed)
        if np.isfinite(a) and a >= obs:
            hits_c += 1

    # how much of the y-qc association survives inside strata
    within = float(np.mean([abs(np.corrcoef(qc[site == s], y[site == s])[0, 1])
                            for s in np.unique(site)
                            if len(np.unique(y[site == s])) > 1]))
    return dict(arm="A_within_stratum_confounder", gamma=gamma, rep=rep, seed=seed,
                observed_auc=obs, n_case=int(y.sum()),
                within_stratum_abs_corr=within,
                p_unstratified=(hits_u + 1) / (N_PERM + 1),
                p_collection_preserving=((hits_c + 1) / (usable + 1)
                                         if usable else np.nan),
                n_usable=usable)


# ------------------------------------------------------------------- arm B
def arm_b(n_design_col, rep, n=108, p=4000, bio=5.5, layout="balanced"):
    """D generated independently of Y. Any residualisation loss is width."""
    seed = cell_seed("B", layout, n_design_col, rep)
    rng = np.random.default_rng(seed)

    y = np.zeros(n, dtype=int); y[: n // 2] = 1
    y = rng.permutation(y)

    # a real biological effect, so there is something for residualisation to lose
    # bio is calibrated so the unadjusted AUC lands near the CMV cohort AUC of 0.845,
    # otherwise a saturated 1.000 hides any loss.
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

    var = X.var(axis=0)
    Xk = X[:, np.argsort(var)[::-1][:4000]]

    def auc_of(M):
        n_pc = int(min(PERM_PCS, M.shape[1], n - 2))
        Mp = PCA(n_components=n_pc, random_state=seed).fit_transform(
            StandardScaler().fit_transform(M))
        return frozen_auc(Mp, y, seed)

    unadj = auc_of(Xk)

    # fold-contained ridge residualisation, exactly as in the empirical analysis
    oof = np.zeros(n)
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(Xk, y):
        r = Ridge(alpha=1.0).fit(D[tr], Xk[tr])
        Rtr, Rte = Xk[tr] - r.predict(D[tr]), Xk[te] - r.predict(D[te])
        n_pc = int(min(PERM_PCS, Rtr.shape[1], len(tr) - 2))
        pca = PCA(n_components=n_pc, random_state=seed)
        sc = StandardScaler().fit(Rtr)
        Ptr = pca.fit_transform(sc.transform(Rtr))
        Pte = pca.transform(sc.transform(Rte))
        sc2 = StandardScaler().fit(Ptr)
        m = LogisticRegression(C=C_FIXED, max_iter=2000).fit(sc2.transform(Ptr), y[tr])
        oof[te] = m.decision_function(sc2.transform(Pte))
    resid = roc_auc_score(y, oof)

    # design-only AUC, to confirm D really is independent of Y here
    d_auc = frozen_auc(D, y, seed) if D.shape[1] >= 2 else np.nan

    return dict(arm="B_residualisation_width", layout=layout,
                n_design_col=int(D.shape[1]), rep=rep,
                seed=seed, n_donor=n, design_only_auc=d_auc,
                unadjusted_auc=unadj, residualised_auc=resid,
                residualisation_loss=unadj - resid)


def main():
    t0 = time.time()
    gammas = [0.0, 0.3, 0.6, 1.0, 1.5, 2.0]
    cells_a = list(itertools.product(gammas, range(N_REP)))
    print(f"arm A: {len(cells_a)} runs", flush=True)
    ra = Parallel(n_jobs=N_JOBS, verbose=5, batch_size=8)(
        delayed(arm_a)(g, r) for g, r in cells_a)
    da = pd.DataFrame([x for x in ra if x])
    da.to_csv(OUT / "calibration_arm_a_raw.tsv", sep="\t", index=False)
    agg_a = (da.assign(rej_u=da.p_unstratified <= 0.05,
                      rej_c=da.p_collection_preserving <= 0.05)
              .groupby("gamma")
              .agg(n=("rep", "size"),
                   mean_within_corr=("within_stratum_abs_corr", "mean"),
                   mean_observed_auc=("observed_auc", "mean"),
                   reject_unstratified=("rej_u", "mean"),
                   reject_collection_preserving=("rej_c", "mean"))
              .round(4))
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
    print(f"arm B: {len(cells_b)} runs", flush=True)
    rb = Parallel(n_jobs=N_JOBS, verbose=5, batch_size=4)(
        delayed(arm_b)(w, r, layout=lay) for lay, w, r in cells_b)
    db = pd.DataFrame([x for x in rb if x])
    db.to_csv(OUT / "calibration_arm_b_raw.tsv", sep="\t", index=False)
    aggb = (db.groupby(["layout", "n_design_col"])
              .agg(n=("rep", "size"),
                   design_only_auc=("design_only_auc", "mean"),
                   unadjusted=("unadjusted_auc", "mean"),
                   residualised=("residualised_auc", "mean"),
                   loss_mean=("residualisation_loss", "mean"),
                   loss_p025=("residualisation_loss", lambda s: s.quantile(.025)),
                   loss_p975=("residualisation_loss", lambda s: s.quantile(.975)))
              .round(4))
    aggb.to_csv(OUT / "calibration_arm_b_summary.tsv", sep="\t")
    print(aggb.to_string(), flush=True)
    print(f"done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
