"""Learned donor pooling under the locked reciprocal source-only transfer protocol.

Everything fitted -- cell standardisation, PCA basis, pooling-network weights,
hyperparameters -- is fitted on SOURCE donors only. The target cohort contributes
unlabelled cell embeddings and nothing else.
"""
import json, os, sys, time, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import poolnet

SEED = 20260815
PCA_K = 32
GRID = [dict(H=16, wd=1e-2), dict(H=16, wd=1e-1), dict(H=32, wd=1e-2), dict(H=32, wd=1e-1)]
EPOCHS = 150
LR = 0.02
ENSEMBLE_SEEDS = (1, 2, 3)
HEADS = 4
MODELS = ("deepsets", "gatedmil", "pma")

IN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

def load_cells(tag):
    if tag == "GSE285773":
        X = np.load(IN + "/GSE285773_cells_f16.npy")
        idx = np.load(IN + "/GSE285773_cell_donor_idx.npy")
    else:
        X = np.concatenate([np.load(IN + "/GSE174188_cells_grouped_f16_%02d.npy" % i) for i in range(10)], axis=0)
        idx = np.load(IN + "/GSE174188_cell_donor_idx_v2.npy")
    counts = np.bincount(idx.astype(np.int64))
    starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
    return X, starts.astype(np.int64), counts.astype(np.int64)

def load_donors(tag):
    f = IN + ("/GSE174188_donor_index_v2.csv" if tag == "GSE174188" else "/GSE285773_donor_index.csv")
    with open(f) as fh:
        return [r["donor_id"] for r in csv.DictReader(fh)]

LAB = {}
with open(IN + "/gse174188_final_donor_covariates.tsv") as fh:
    LAB["GSE174188"] = {r["donor_id"]: r["case_control"] for r in csv.DictReader(fh, delimiter="\t")}
with open(IN + "/labels_285773.tsv") as fh:
    LAB["GSE285773"] = {r["donor_id"]: r["case_control"] for r in csv.DictReader(fh, delimiter="\t")}

DATA = {}
for tag in ("GSE174188", "GSE285773"):
    X, starts, counts = load_cells(tag)
    donors = load_donors(tag)
    assert len(donors) == len(counts) and counts.sum() == len(X), (len(donors), len(counts), counts.sum(), len(X))
    y = np.array([1.0 if LAB[tag][d] == "case" else 0.0 for d in donors])
    DATA[tag] = dict(X=X, starts=starts, counts=counts, donors=donors, y=y)
    print(tag, X.shape, len(donors), "cases", int(y.sum()), flush=True)

def fit_projection(Xs, k=PCA_K):
    Xf = Xs.astype(np.float32)
    mu = Xf.mean(axis=0)
    sd = Xf.std(axis=0); sd[sd == 0] = 1.0
    Z = (Xf - mu) / sd
    cov = (Z.T @ Z) / (len(Z) - 1)
    w, V = np.linalg.eigh(cov.astype(np.float64))
    order = np.argsort(w)[::-1][:k]
    return dict(mu=mu, sd=sd, basis=V[:, order], scale=np.sqrt(np.maximum(w[order], 1e-12)))

def apply_projection(P, X):
    out = np.empty((len(X), P["basis"].shape[1]), dtype=np.float64)
    step = 20000
    for i in range(0, len(X), step):
        Z = (X[i:i+step].astype(np.float32) - P["mu"]) / P["sd"]
        out[i:i+step] = (Z @ P["basis"]) / P["scale"]
    return out

def donor_folds(y, k=5, seed=SEED):
    rng = np.random.default_rng(seed)
    folds = [[] for _ in range(k)]
    for cls in (0.0, 1.0):
        idx = np.flatnonzero(y == cls)
        idx = idx[rng.permutation(len(idx))]
        for j, i in enumerate(idx):
            folds[j % k].append(i)
    return [np.array(sorted(f)) for f in folds]

def subset(starts, counts, donor_idx):
    rows = np.concatenate([np.arange(starts[d], starts[d] + counts[d]) for d in donor_idx])
    c = counts[donor_idx]
    s = np.concatenate(([0], np.cumsum(c)[:-1]))
    return rows, s.astype(np.int64), c.astype(np.int64)

def auc(y, p):
    o = np.argsort(p, kind="mergesort")
    r = np.empty(len(p)); r[o] = np.arange(1, len(p) + 1)
    sp = p[o]; i = 0
    while i < len(sp):
        j = i
        while j + 1 < len(sp) and sp[j + 1] == sp[i]:
            j += 1
        if j > i:
            r[o[i:j+1]] = (i + j + 2) / 2.0
        i = j + 1
    n1 = float((y == 1).sum()); n0 = float((y == 0).sum())
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))

DIRECTIONS = [("GSE174188", "GSE285773"), ("GSE285773", "GSE174188")]
results, records = [], {}
for src, tgt in DIRECTIONS:
    t0 = time.time()
    S, T = DATA[src], DATA[tgt]
    print("\n=== %s -> %s" % (src, tgt), flush=True)
    P = fit_projection(S["X"])
    Zs = apply_projection(P, S["X"]); Zt = apply_projection(P, T["X"])
    print("  projection fitted on source cells: %s -> target %s (%.0fs)" % (Zs.shape, Zt.shape, time.time()-t0), flush=True)
    folds = donor_folds(S["y"])
    for kind in MODELS:
        cv_rows = []
        for cfg in GRID:
            oof = np.full(len(S["y"]), np.nan)
            for f in range(5):
                te = folds[f]
                tr = np.concatenate([folds[j] for j in range(5) if j != f])
                rtr, str_, ctr = subset(S["starts"], S["counts"], tr)
                rte, ste, cte = subset(S["starts"], S["counts"], te)
                p = poolnet.train(kind, Zs[rtr], str_, ctr, S["y"][tr], PCA_K, cfg["H"], EPOCHS, LR, cfg["wd"], seed=SEED + f, heads=HEADS)
                oof[te] = poolnet.predict(kind, p, Zs[rte], ste, cte, HEADS)
            cv_rows.append(dict(H=cfg["H"], wd=cfg["wd"], cv_auc=auc(S["y"], oof)))
            print("  %-9s H=%-3d wd=%-5g source-CV AUC=%.4f (%.0fs)" % (kind, cfg["H"], cfg["wd"], cv_rows[-1]["cv_auc"], time.time()-t0), flush=True)
        best = max(cv_rows, key=lambda r: r["cv_auc"])
        probs = []
        for s in ENSEMBLE_SEEDS:
            p = poolnet.train(kind, Zs, S["starts"], S["counts"], S["y"], PCA_K, best["H"], EPOCHS, LR, best["wd"], seed=SEED + 100 * s, heads=HEADS)
            probs.append(poolnet.predict(kind, p, Zt, T["starts"], T["counts"], HEADS))
        pt = np.mean(probs, axis=0)
        a = auc(T["y"], pt)
        print("  -> %s: best H=%d wd=%g | TARGET AUC=%.4f" % (kind, best["H"], best["wd"], a), flush=True)
        np.save(OUT + "/pred_%s_%s_to_%s.npy" % (kind, src, tgt), pt)
        results.append(dict(direction="%s_to_%s" % (src, tgt), method=kind, target_auc=a,
                            best_H=best["H"], best_wd=best["wd"], source_cv_auc=best["cv_auc"]))
        records["%s|%s_to_%s" % (kind, src, tgt)] = cv_rows
    ms = np.add.reduceat(Zs, S["starts"], axis=0) / S["counts"][:, None]
    mt = np.add.reduceat(Zt, T["starts"], axis=0) / T["counts"][:, None]
    np.save(OUT + "/feat_pca32_mean_source_%s_to_%s.npy" % (src, tgt), ms)
    np.save(OUT + "/feat_pca32_mean_target_%s_to_%s.npy" % (src, tgt), mt)
    np.save(OUT + "/y_source_%s_to_%s.npy" % (src, tgt), S["y"])
    np.save(OUT + "/y_target_%s_to_%s.npy" % (src, tgt), T["y"])

with open(OUT + "/learned_pooling_summary.json", "w") as fh:
    json.dump(dict(results=results, cv=records, config=dict(pca_k=PCA_K, grid=GRID, epochs=EPOCHS, lr=LR, seeds=list(ENSEMBLE_SEEDS), heads=HEADS, seed=SEED)), fh, indent=2)
print("\nDONE", flush=True)
for r in results:
    print(r, flush=True)
