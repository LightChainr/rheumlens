"""Extended donor-level simulation for the restructured manuscript (items B3/B4).

The previous simulation used one binary design variable, a continuous label, an
additive linear design effect and homoscedastic noise. Referees #1 and #4 asked
whether the conclusions survive when those four choices are relaxed. This script
crosses six data-generating regimes with seven levels of design-label association
and both label types, and reports the same quantities the empirical screen reports.

Nothing here is fitted to the empirical cohorts; the simulation is a check on the
estimator's behaviour, not a model of any dataset.
"""
from __future__ import annotations
import hashlib, itertools, json, os, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression, LinearRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.preprocessing import StandardScaler

OUT = Path(os.environ.get("RHEUMLENS_SIM_OUT",
                          Path(__file__).resolve().parent / "results"))
OUT.mkdir(exist_ok=True)

N_DONOR   = 200
N_FEATURE = 100
N_REP     = 200
RHOS      = [0.0, 0.25, 0.50, 0.75, 0.90, 0.98, 1.00]
ARMS      = ["base", "nonlinear_design", "heteroscedastic",
             "design_by_celltype", "imbalanced", "imbalanced_offset",
             "mechanism_shift"]
LABEL_TYPES = ["binary", "continuous"]
C_FIXED   = 1.0
MASTER_SEED = 20260908   # every cell seed is derived from this, reproducibly


# ----------------------------------------------------------------- generation
def make_cohort(rng, rho, arm, label_type, n=N_DONOR, p=N_FEATURE, shift=False,
                loadings=None):
    """Return X, y, y_cont, D, strata, loadings.

    `loadings` carries the biological and technical directions. An external target
    must reuse the SOURCE loadings: the biological effect is a property of the
    disease, not of the cohort. Regenerating them makes every transfer chance-level
    by construction and tests nothing.
    """
    # latent design/label pair with association rho
    z = rng.standard_normal(n)
    e = rng.standard_normal(n)
    y_lat = rho * z + np.sqrt(max(1.0 - rho**2, 0.0)) * e

    # Label threshold, and the site threshold that goes with it.
    #
    # The previous version moved the label to the 80th percentile for the
    # imbalanced arm but left the site indicator at z > 0. At rho = 1 the two
    # therefore stopped coinciding, so design-only AUC could not reach 1 even
    # under perfect confounding, and the resulting ceiling was reported as an
    # effect of prevalence. It is two effects at once. They are now separated:
    #   imbalanced         - both thresholds at the 80th percentile, so a rare
    #                        diagnosis is the ONLY thing that changed
    #   imbalanced_offset  - label at the 80th percentile, site at the median,
    #                        the realistic case where collection boundaries do
    #                        not track prevalence
    q = 0.80 if arm in ("imbalanced", "imbalanced_offset") else 0.50
    thresh = np.quantile(y_lat, q)
    site_q = 0.50 if arm == "imbalanced_offset" else q
    site = (z > np.quantile(z, site_q)).astype(float)
    qc   = 0.6 * z + 0.8 * rng.standard_normal(n)      # library-quality proxy
    D    = np.column_stack([site, qc])
    strata = site.astype(int)                           # collection stratum

    y = (y_lat > thresh).astype(int)
    if y.sum() < 5 or (n - y.sum()) < 5:                # guard degenerate draws
        y = (y_lat > np.median(y_lat)).astype(int)

    y_used = y_lat if label_type == "continuous" else y.astype(float)
    y_used = (y_used - y_used.mean()) / (y_used.std() + 1e-12)

    # loading vectors: biological and technical directions overlap partially
    if loadings is None:
        u_bio  = rng.standard_normal(p); u_bio  /= np.linalg.norm(u_bio)
        u_tech = rng.standard_normal(p); u_tech /= np.linalg.norm(u_tech)
        u_tech = 0.7 * u_tech + 0.3 * u_bio
        u_tech /= np.linalg.norm(u_tech)
    else:
        u_bio, u_tech = loadings
        if shift and arm == "mechanism_shift":
            # the target's technical direction is its own; the biology is shared
            u_tech = rng.standard_normal(p); u_tech /= np.linalg.norm(u_tech)

    beta, delta = 1.0, 1.0
    if shift and arm == "mechanism_shift":
        delta = 2.5                                      # target acquired differently
        zz = rng.standard_normal(n)
        site = (zz > np.quantile(zz, site_q)).astype(float)
        D = np.column_stack([site, qc])
        strata = site.astype(int)

    # design effect on X: linear in `site`, or a nonlinear function of qc
    if arm == "nonlinear_design":
        d_eff = np.tanh(2.0 * qc) + 0.5 * (qc ** 2 - 1.0)
    else:
        d_eff = site - site.mean()

    X = np.outer(beta * y_used, u_bio) + np.outer(delta * d_eff, u_tech)

    if arm == "design_by_celltype":
        # two cell-type blocks; the design shifts composition, the biological
        # effect lives only in block 2, so design and biology interact.
        half = p // 2
        w = 0.5 + 0.35 * (site - 0.5)                    # block-1 weight by site
        X[:, :half]  *= w[:, None]
        X[:, half:]  *= (1.0 - w)[:, None]
        X[:, half:]  += np.outer(0.8 * y_used * w, u_bio[half:])

    if arm == "heteroscedastic":
        sd = 1.0 + 0.8 * site                            # site-dependent noise
        X = X + rng.standard_normal((n, p)) * sd[:, None]
    else:
        X = X + rng.standard_normal((n, p))

    return X, y, y_lat, D, strata, (u_bio, u_tech)


# ------------------------------------------------------------------- metrics
def cindex(truth: np.ndarray, score: np.ndarray) -> float:
    """Concordance over comparable pairs.

    For a binary outcome this is exactly the ROC AUC, so binary and continuous
    arms are reported on one scale and can be compared. For a continuous outcome
    it is the fraction of pairs whose predicted order matches the true order,
    which is what "does the conclusion depend on outcome type" actually needs.
    """
    t = np.asarray(truth, dtype=float)
    sc = np.asarray(score, dtype=float)
    if len(np.unique(t)) < 2:
        return np.nan
    # rank-based, so it is O(n log n) rather than O(n^2)
    from scipy.stats import kendalltau
    tau = kendalltau(t, sc, variant="b").statistic
    return float((tau + 1.0) / 2.0) if np.isfinite(tau) else np.nan


def cv_score_continuous(X, y_cont, seed, alpha=1.0):
    """Cross-fitted ridge prediction of a CONTINUOUS outcome, scored by c-index."""
    oof = np.zeros(len(y_cont))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        sc = StandardScaler().fit(X[tr])
        m = Ridge(alpha=alpha).fit(sc.transform(X[tr]), y_cont[tr])
        oof[te] = m.predict(sc.transform(X[te]))
    return cindex(y_cont, oof)


def residualised_continuous(X, y_cont, D, seed, alpha=1.0):
    oof = np.zeros(len(y_cont))
    for tr, te in KFold(5, shuffle=True, random_state=seed).split(X):
        r = Ridge(alpha=1.0).fit(D[tr], X[tr])
        Xtr, Xte = X[tr] - r.predict(D[tr]), X[te] - r.predict(D[te])
        sc = StandardScaler().fit(Xtr)
        m = Ridge(alpha=alpha).fit(sc.transform(Xtr), y_cont[tr])
        oof[te] = m.predict(sc.transform(Xte))
    return cindex(y_cont, oof)


def restricted_continuous(X, y_cont, strata, seed):
    best, best_n = np.nan, 0
    for st in np.unique(strata):
        m = strata == st
        if m.sum() < 40:
            continue
        if m.sum() > best_n:
            best_n, best = m.sum(), cv_score_continuous(X[m], y_cont[m], seed)
    return best


def cv_auc(X, y, seed, C=C_FIXED):
    if len(np.unique(y)) < 2 or min(np.bincount(y)) < 3:
        return np.nan
    oof = np.zeros(len(y))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, y):
        sc = StandardScaler().fit(X[tr])
        m = LogisticRegression(C=C, max_iter=2000).fit(sc.transform(X[tr]), y[tr])
        oof[te] = m.decision_function(sc.transform(X[te]))
    return roc_auc_score(y, oof)


def id_contrast(y_c, D, seed, n_perm=100):
    """Cross-fitted V_D against a CROSS-FITTED permutation null.

    Observed statistic and null statistic must be the same quantity. The earlier
    version reported the cross-fitted value but tested the in-sample one, so the
    p-value did not belong to the number beside it.
    """
    n = len(y_c)
    folds = list(KFold(5, shuffle=True, random_state=seed).split(D))

    def vd(yy):
        o = np.zeros(n)
        for tr, te in folds:
            o[te] = LinearRegression().fit(D[tr], yy[tr]).predict(D[te])
        ss = ((yy - yy.mean()) ** 2).sum()
        r2 = 1.0 - ((yy - o) ** 2).sum() / ss if ss > 0 else 0.0
        return 1.0 - max(r2, 0.0)

    i_cv = vd(y_c)
    rng = np.random.default_rng(seed)
    null = np.array([vd(rng.permutation(y_c)) for _ in range(n_perm)])
    p = (np.sum(null <= i_cv) + 1) / (n_perm + 1)
    return i_cv, float(null.mean()), float(p)


def residualised_auc(X, y, D, seed):
    """Fold-contained ridge residualisation of X on D, then classify."""
    if len(np.unique(y)) < 2 or min(np.bincount(y)) < 3:
        return np.nan
    oof = np.zeros(len(y))
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, y):
        r = Ridge(alpha=1.0, fit_intercept=True).fit(D[tr], X[tr])
        Xtr, Xte = X[tr] - r.predict(D[tr]), X[te] - r.predict(D[te])
        sc = StandardScaler().fit(Xtr)
        m = LogisticRegression(C=C_FIXED, max_iter=2000).fit(sc.transform(Xtr), y[tr])
        oof[te] = m.decision_function(sc.transform(Xte))
    return roc_auc_score(y, oof)


def restricted_auc(X, y, strata, seed):
    """AUC inside the largest design stratum that contains both labels."""
    best, best_n = np.nan, 0
    for s in np.unique(strata):
        m = strata == s
        if len(np.unique(y[m])) < 2 or m.sum() < 40:
            continue
        if min(np.bincount(y[m])) < 5:
            continue
        if m.sum() > best_n:
            best_n, best = m.sum(), cv_auc(X[m], y[m], seed)
    return best


# ------------------------------------------------------------------ one cell
def cell_seed(arm, rho, label_type, rep) -> int:
    """Stable seed for one simulation cell.

    Python's built-in hash() is salted per interpreter process for str, so the
    previous `hash((arm, rho, label_type, rep))` produced different seeds on
    every fresh run unless PYTHONHASHSEED happened to be pinned. A blake2b digest
    of the same tuple is stable across processes, machines and versions.
    """
    key = f"{MASTER_SEED}|{arm}|{rho:.3f}|{label_type}|{rep}".encode()
    return int.from_bytes(hashlib.blake2b(key, digest_size=4).digest(), "big") % (2**31)


def one_rep(arm, rho, label_type, rep):
    seed = cell_seed(arm, rho, label_type, rep)
    rng = np.random.default_rng(seed)
    X, y, y_lat, D, strata, load = make_cohort(rng, rho, arm, label_type)

    # The outcome that is actually predicted. The previous version used the
    # continuous latent variable only to build the expression signal and then
    # classified the binary label in both arms, so "continuous" tested nothing
    # about outcome type. Each arm now predicts, and is scored on, its own
    # outcome, with the c-index as the common scale (it equals AUC when binary).
    if label_type == "continuous":
        target = y_lat
        des = cv_score_continuous(D, target, seed)
        dis = cv_score_continuous(X, target, seed)
        res = residualised_continuous(X, target, D, seed)
        rst = restricted_continuous(X, target, strata, seed)
        y_c = target - target.mean()
    else:
        target = y
        des = cv_auc(D, y, seed)
        dis = cv_auc(X, y, seed)
        res = residualised_auc(X, y, D, seed)
        rst = restricted_auc(X, y, strata, seed)
        y_c = (y.astype(float) - y.mean())
    i_cv, i_null, p_id = id_contrast(y_c, D, seed)

    # external target: independently generated, same or different mechanism
    Xt, yt, yt_lat, Dt, _, _ = make_cohort(np.random.default_rng(seed + 7), rho, arm,
                                      label_type, shift=(arm == "mechanism_shift"),
                                      loadings=load)
    ext = np.nan
    if label_type == "continuous":
        sc = StandardScaler().fit(X)
        m = Ridge(alpha=1.0).fit(sc.transform(X), y_lat)
        ext = cindex(yt_lat, m.predict(sc.transform(Xt)))
    elif len(np.unique(y)) == 2 and len(np.unique(yt)) == 2:
        sc = StandardScaler().fit(X)
        m = LogisticRegression(C=C_FIXED, max_iter=2000).fit(sc.transform(X), y)
        ext = roc_auc_score(yt, m.decision_function(sc.transform(Xt)))

    return dict(arm=arm, rho=rho, label_type=label_type, rep=rep, seed=seed,
                n_case=int(y.sum()), design_c=des, I_D_cv=i_cv,
                I_D_null=i_null, p_I_D=p_id, disease_c=dis,
                residualised_c=res, restricted_c=rst, external_c=ext)


def main():
    cells = list(itertools.product(ARMS, RHOS, LABEL_TYPES, range(N_REP)))
    print(f"{len(cells)} runs", flush=True)
    t0 = time.time()
    rows = Parallel(n_jobs=8, verbose=5, batch_size=16)(
        delayed(one_rep)(*c) for c in cells)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "extended_simulation_raw.tsv", sep="\t", index=False)
    agg = (df.groupby(["arm", "label_type", "rho"])
             .agg(["mean", "std"])[["design_c", "I_D_cv", "p_I_D", "disease_c",
                                    "residualised_c", "restricted_c",
                                    "external_c"]]
             .round(4))
    agg.to_csv(OUT / "extended_simulation_summary.tsv", sep="\t")
    print(f"done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
