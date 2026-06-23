#!/usr/bin/env python3
"""Auditable GSE135779 pediatric-to-adult donor transfer rerun."""

import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from fold_contained_benchmark import aggregate_matrix_by_donor, classifier, gene_names, select_hvg


ROOT = Path(__file__).resolve().parents[1]
SSD = Path("/Volumes/Mac Data/Research/论文资料库/04_来自Downloads_2026-06-17/P3_GPU服务器与MD/rheumlens_gpu_pull_20260617")
OUT = ROOT / "results" / "evidence_package" / "cross_age_transfer"


def digest(x):
    return hashlib.sha256(np.ascontiguousarray(x).view(np.uint8)).hexdigest()


def aggregate_adult():
    a = ad.read_h5ad(SSD / "gse135779_adult_embeddings.h5ad", backed="r")
    cell_donor = a.obs["donor_id"].astype(str).to_numpy()
    donors = np.asarray(sorted(pd.unique(cell_donor)))
    dmap = {d: i for i, d in enumerate(donors)}
    codes = np.asarray([dmap[d] for d in cell_donor], dtype=np.int32)
    ycell = a.obs["disease"].astype(str).str.upper().eq("SLE").astype(int).to_numpy()
    y = np.asarray([np.unique(ycell[codes == i]).item() for i in range(len(donors))])
    expression = aggregate_matrix_by_donor(a.X, codes, len(donors))
    embedding = aggregate_matrix_by_donor(a.obsm["X_scGPT"], codes, len(donors), chunk_size=16384)
    genes = gene_names(a)
    a.file.close()
    return donors, y, expression.astype(float), embedding.astype(float), genes


def bootstrap(y, p, n=10000):
    rng = np.random.default_rng(20260619)
    case, control = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    values = np.empty(n)
    for i in range(n):
        idx = np.concatenate([rng.choice(case, len(case), replace=True),
                              rng.choice(control, len(control), replace=True)])
        values[i] = roc_auc_score(y[idx], p[idx])
    return np.quantile(values, .025), np.quantile(values, .975), values


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    child = np.load(ROOT / "data" / "donor_features" / "GSE135779_donor_features.npz", allow_pickle=True)
    yc = child["y"].astype(int)
    dc, xc, ec, gc = child["donors"].astype(str), child["expression"].astype(float), child["scgpt"].astype(float), child["genes"].astype(str)
    da, ya, xa, ea, ga = aggregate_adult()
    pred = {"adult_donor_id": da, "label": ya}
    provenance = []

    for name, scaled in (("scgpt_source_scaled", True), ("scgpt_unscaled", False)):
        if scaled:
            scaler = StandardScaler().fit(ec)
            tr, te = scaler.transform(ec), scaler.transform(ea)
            scaler_hash = digest(np.concatenate([scaler.mean_, scaler.scale_]))
        else:
            tr, te, scaler_hash = ec, ea, "none"
        model = classifier().fit(tr, yc)
        p = model.predict_proba(te)[:, 1]
        pred[name] = p
        lo, hi, values = bootstrap(ya, p)
        np.save(OUT / f"{name}_bootstrap_auc.npy", values)
        provenance.append({"method": name, "auc": roc_auc_score(ya, p), "ci_low": lo, "ci_high": hi,
                           "coefficient_sha256": digest(model.coef_), "scaler_sha256": scaler_hash})

    child_lookup = {g.split("__dup", 1)[0]: i for i, g in enumerate(gc)}
    adult_lookup = {g.split("__dup", 1)[0]: i for i, g in enumerate(ga)}
    common = sorted(set(child_lookup) & set(adult_lookup))
    cidx = np.asarray([child_lookup[g] for g in common])
    aidx = np.asarray([adult_lookup[g] for g in common])
    hvg_local = select_hvg(xc[:, cidx], 2000)
    pca = PCA(n_components=25, svd_solver="full", random_state=42).fit(xc[:, cidx[hvg_local]])
    ctr, ate = pca.transform(xc[:, cidx[hvg_local]]), pca.transform(xa[:, aidx[hvg_local]])
    scaler = StandardScaler().fit(ctr)
    model = classifier().fit(scaler.transform(ctr), yc)
    p = model.predict_proba(scaler.transform(ate))[:, 1]
    pred["expression_pca_source_fitted"] = p
    lo, hi, values = bootstrap(ya, p)
    np.save(OUT / "expression_pca_source_fitted_bootstrap_auc.npy", values)
    provenance.append({"method": "expression_pca_source_fitted", "auc": roc_auc_score(ya, p),
                       "ci_low": lo, "ci_high": hi, "n_common_genes": len(common),
                       "n_hvg": len(hvg_local), "coefficient_sha256": digest(model.coef_),
                       "scaler_sha256": digest(np.concatenate([scaler.mean_, scaler.scale_]))})

    pd.DataFrame(pred).to_csv(OUT / "GSE135779_pediatric_to_adult_predictions.csv", index=False)
    pd.DataFrame(provenance).to_csv(OUT / "GSE135779_pediatric_to_adult_summary.csv", index=False)
    (OUT / "provenance.json").write_text(json.dumps({
        "source": "44 pediatric donors (33 SLE/11 HD)",
        "target": "12 archived adult donors (7 SLE/5 HD)",
        "preprocessing": "all fit operations estimated in pediatric source only",
        "models": provenance,
    }, indent=2), encoding="utf-8")
    print(pd.DataFrame(provenance).to_string(index=False))


if __name__ == "__main__":
    main()
