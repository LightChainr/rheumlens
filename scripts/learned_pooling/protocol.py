"""Locked source-only transfer protocol: deterministic C selection, metrics, inference.

Reconstructed on 2026-07-26 on the local Mac after the cloud container that held the
original copy expired. Behaviour is specified by the handoff document
(`RheumLens_重构工作交接文档_20260726.md`, sections 2.2-2.4):

  * standardisation is fitted ONCE on the whole source cohort, then passed into the
    inner cross-validation (this is what the archived `LogisticRegressionCV` did --
    it is NOT fold-contained rescaling);
  * inner selection uses StratifiedKFold(5, shuffle=True, random_state=20260815) over
    Cs = logspace(-4, 4, 9), scoring mean fold ROC-AUC;
  * ties resolve to the FIRST maximum, i.e. the smaller C -- this makes the selector
    deterministic and independent of the scikit-learn version;
  * the target cohort contributes no labels and no fitted statistic.

Everything downstream (bootstrap CIs, paired DeLong, Benjamini-Hochberg) is implemented
here so the whole comparison chain lives in one auditable file.
"""
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

SEED = 20260815
CS = np.logspace(-4, 4, 9)


def _lr(C):
    return LogisticRegression(C=C, solver="liblinear", class_weight="balanced", max_iter=20000)


def roc_auc(y, s):
    """Rank-based ROC-AUC with mid-ranks for ties."""
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = np.arange(1, len(s) + 1)
    ss = s[order]
    i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    n1 = float((y == 1).sum())
    n0 = float((y == 0).sum())
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def pr_auc(y, s):
    """Average precision, identical convention to sklearn.average_precision_score."""
    y = np.asarray(y, dtype=float)
    s = np.asarray(s, dtype=float)
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y)
    prec = tp / np.arange(1, len(y) + 1)
    rec = tp / max(y.sum(), 1.0)
    drec = np.diff(np.concatenate(([0.0], rec)))
    return float((prec * drec).sum())


def brier(y, p):
    return float(np.mean((np.asarray(p, dtype=float) - np.asarray(y, dtype=float)) ** 2))


def ece_10bin(y, p, bins=10):
    """Expected calibration error over `bins` equal-width probability bins."""
    y = np.asarray(y, dtype=float)
    p = np.asarray(p, dtype=float)
    edges = np.linspace(0.0, 1.0, bins + 1)
    idx = np.clip(np.digitize(p, edges[1:-1], right=False), 0, bins - 1)
    e = 0.0
    for b in range(bins):
        m = idx == b
        if m.any():
            e += m.mean() * abs(p[m].mean() - y[m].mean())
    return float(e)


def select_C(X, y, Cs=CS, seed=SEED, folds=5):
    """Deterministic inner selection of C on source donors only.

    Returns (best_C, fold-mean AUC per C).
    """
    Z = StandardScaler().fit_transform(X)
    skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    splits = list(skf.split(Z, y))
    scores = []
    for C in Cs:
        aucs = []
        for tr, te in splits:
            m = _lr(C).fit(Z[tr], y[tr])
            aucs.append(roc_auc(y[te], m.predict_proba(Z[te])[:, 1]))
        scores.append(float(np.mean(aucs)))
    scores = np.asarray(scores)
    best = int(np.argmax(scores))  # first maximum -> smaller C on ties
    return float(Cs[best]), scores


def transfer(Xs, ys, Xt, C=None, Cs=CS, seed=SEED):
    """Fit on source, score target. Scaler and classifier are source-fitted only."""
    if C is None:
        C, _ = select_C(Xs, ys, Cs=Cs, seed=seed)
    sc = StandardScaler().fit(Xs)
    m = _lr(C).fit(sc.transform(Xs), ys)
    return m.predict_proba(sc.transform(Xt))[:, 1], float(C)


def evaluate(y, p):
    return dict(roc_auc=roc_auc(y, p), pr_auc=pr_auc(y, p), brier=brier(y, p), ece_10bin=ece_10bin(y, p))


def reproducible_at(Xs, ys, Xt, yt, target_auc, Cs=CS, tol=1e-9):
    """Which grid C, if any, reproduces an archived AUC exactly."""
    hits = []
    for C in Cs:
        p, _ = transfer(Xs, ys, Xt, C=C)
        a = roc_auc(yt, p)
        if abs(a - target_auc) <= tol:
            hits.append((float(C), a))
    return hits


# ---------------------------------------------------------------- inference


def _midrank(x):
    order = np.argsort(x, kind="mergesort")
    xs = x[order]
    n = len(x)
    r = np.empty(n, dtype=float)
    i = 0
    while i < n:
        j = i
        while j < n - 1 and xs[j + 1] == xs[i]:
            j += 1
        r[i:j + 1] = 0.5 * (i + j) + 1
        i = j + 1
    out = np.empty(n, dtype=float)
    out[order] = r
    return out


def delong_cov(y, scores):
    """Fast DeLong (Sun & Xu 2014). `scores` is (n_methods, n_donors)."""
    y = np.asarray(y)
    pos = scores[:, y == 1]
    neg = scores[:, y == 0]
    m, n = pos.shape[1], neg.shape[1]
    k = scores.shape[0]
    tx = np.array([_midrank(pos[r]) for r in range(k)])
    ty = np.array([_midrank(neg[r]) for r in range(k)])
    tz = np.array([_midrank(np.concatenate([pos[r], neg[r]])) for r in range(k)])
    aucs = (tz[:, :m].sum(axis=1) - m * (m + 1) / 2.0) / (m * n)
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    s01 = np.cov(v01)
    s10 = np.cov(v10)
    s01 = np.atleast_2d(s01)
    s10 = np.atleast_2d(s10)
    return aucs, s01 / m + s10 / n


def delong_paired_p(y, s_a, s_b):
    """Two-sided p value for AUC(a) - AUC(b) on the same donors."""
    from scipy import stats
    aucs, cov = delong_cov(y, np.vstack([s_a, s_b]))
    d = np.array([1.0, -1.0])
    var = float(d @ cov @ d)
    if var <= 0:
        return float(aucs[0] - aucs[1]), 1.0
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    return float(aucs[0] - aucs[1]), float(2.0 * stats.norm.sf(abs(z)))


def bh(pvals):
    """Benjamini-Hochberg adjusted p values, order preserved."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p, kind="mergesort")
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        prev = min(prev, p[i] * n / (rank + 1))
        adj[i] = prev
    return np.clip(adj, 0.0, 1.0)


def boot_auc_ci(y, p, n=5000, seed=SEED, alpha=0.05):
    """Case/control-stratified donor bootstrap percentile interval for one AUC."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    ipos = np.flatnonzero(y == 1)
    ineg = np.flatnonzero(y == 0)
    vals = np.empty(n)
    for b in range(n):
        idx = np.concatenate([rng.choice(ipos, len(ipos), replace=True),
                              rng.choice(ineg, len(ineg), replace=True)])
        vals[b] = roc_auc(y[idx], np.asarray(p)[idx])
    return float(np.quantile(vals, alpha / 2)), float(np.quantile(vals, 1 - alpha / 2))


def boot_auc_diff_ci(y, p_a, p_b, n=5000, seed=SEED, alpha=0.05):
    """Paired stratified bootstrap interval for AUC(a) - AUC(b) on identical donors."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    p_a = np.asarray(p_a)
    p_b = np.asarray(p_b)
    ipos = np.flatnonzero(y == 1)
    ineg = np.flatnonzero(y == 0)
    vals = np.empty(n)
    for b in range(n):
        idx = np.concatenate([rng.choice(ipos, len(ipos), replace=True),
                              rng.choice(ineg, len(ineg), replace=True)])
        vals[b] = roc_auc(y[idx], p_a[idx]) - roc_auc(y[idx], p_b[idx])
    return float(np.quantile(vals, alpha / 2)), float(np.quantile(vals, 1 - alpha / 2))
