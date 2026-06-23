#!/usr/bin/env python3
"""Project-derived broad-lineage sensitivity for GSE135779.

This is explicitly not an author-provided per-cell annotation. Labels use a
fixed four-lineage marker-score rule and are suitable only as sensitivity data.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from fold_contained_benchmark import COHORTS, classifier, folds_from_json, gene_names


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "evidence_package" / "secondary_foldcontained"
DATA = ROOT / "data" / "secondary_features"
MARKERS = {
    "T": ["CD3D", "CD3E", "CD2"],
    "Monocyte": ["CD14", "LYZ", "S100A8"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
}


def resampling(y, p, folds, x, n_boot=10000, n_perm=1000):
    rng = np.random.default_rng(20260619)
    cases, controls = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.concatenate([rng.choice(cases, len(cases), replace=True),
                              rng.choice(controls, len(controls), replace=True)])
        boot[i] = roc_auc_score(y[idx], p[idx])
    observed = roc_auc_score(y, p)
    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = rng.permutation(y)
        pp = np.full(len(y), np.nan)
        if any(np.unique(yp[tr]).size < 2 for tr, _ in folds):
            null[i] = .5
            continue
        for tr, te in folds:
            model = make_pipeline(StandardScaler(), classifier()).fit(x[tr], yp[tr])
            pp[te] = model.predict_proba(x[te])[:, 1]
        null[i] = roc_auc_score(yp, pp)
    exceed = int((null >= observed).sum())
    return np.quantile(boot, .025), np.quantile(boot, .975), (exceed + 1) / (n_perm + 1), boot, null


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(COHORTS["GSE135779"]["h5ad"], backed="r")
    genes = gene_names(a)
    base = np.asarray([g.split("__dup", 1)[0] for g in genes])
    lookup = {g: i for i, g in enumerate(base)}
    marker_indices = {name: [lookup[g] for g in markers if g in lookup]
                      for name, markers in MARKERS.items()}
    donor_cell = a.obs["donor_id"].astype(str).to_numpy()
    donors = np.asarray(sorted(pd.unique(donor_cell)))
    dmap = {d: i for i, d in enumerate(donors)}
    codes = np.asarray([dmap[d] for d in donor_cell], dtype=np.int32)
    donor_disease = (a.obs.assign(_donor=donor_cell)
                     .drop_duplicates("_donor").set_index("_donor")["disease"]
                     .astype(str).str.upper())
    y = donor_disease.reindex(donors).eq("SLE").astype(int).to_numpy()
    folds = folds_from_json(Path(COHORTS["GSE135779"]["folds"]), donors, y)

    lineage_names = list(MARKERS)
    annotation = np.full(a.n_obs, -1, dtype=np.int8)
    chunk = 8192
    all_marker_indices = [i for name in lineage_names for i in marker_indices[name]]
    marker_slices = {}
    offset = 0
    for name in lineage_names:
        marker_slices[name] = slice(offset, offset + len(marker_indices[name]))
        offset += len(marker_indices[name])
    for start in range(0, a.n_obs, chunk):
        stop = min(start + chunk, a.n_obs)
        marker_x = a.X[start:stop, all_marker_indices]
        scores = np.zeros((stop - start, len(lineage_names)), dtype=np.float32)
        for j, name in enumerate(lineage_names):
            x = marker_x[:, marker_slices[name]]
            scores[:, j] = np.asarray(x.mean(axis=1)).ravel()
        best = scores.argmax(axis=1).astype(np.int8)
        best[scores.max(axis=1) <= 0] = -1
        annotation[start:stop] = best
        if stop == a.n_obs or stop % (chunk * 8) == 0:
            print(f"annotation rows {stop:,}/{a.n_obs:,}", flush=True)
    np.save(DATA / "GSE135779_project_marker_lineage.npy", annotation)

    counts = []
    embeddings = {"T": np.zeros((len(donors), 512)), "Monocyte": np.zeros((len(donors), 512))}
    n_lineage = {name: np.zeros(len(donors), dtype=np.int64) for name in lineage_names + ["Unassigned"]}
    for start in range(0, a.n_obs, chunk):
        stop = min(start + chunk, a.n_obs)
        c = codes[start:stop]
        lab = annotation[start:stop]
        emb = np.asarray(a.obsm["X_scGPT"][start:stop], dtype=np.float64)
        for j, name in enumerate(lineage_names):
            mask = lab == j
            np.add.at(n_lineage[name], c[mask], 1)
            if name in embeddings:
                np.add.at(embeddings[name], c[mask], emb[mask])
        mask = lab < 0
        np.add.at(n_lineage["Unassigned"], c[mask], 1)
    for name in embeddings:
        embeddings[name] /= n_lineage[name][:, None]

    count_frame = pd.DataFrame({"donor_id": donors, "label": y, **n_lineage})
    count_frame.to_csv(OUT / "GSE135779_project_marker_lineage_counts.csv", index=False)

    summaries = []
    predictions = {"donor_id": donors, "label": y}
    for name, x in embeddings.items():
        p = np.full(len(y), np.nan)
        for tr, te in folds:
            model = make_pipeline(StandardScaler(), classifier()).fit(x[tr], y[tr])
            p[te] = model.predict_proba(x[te])[:, 1]
        lo, hi, perm_p, boot, null = resampling(y, p, folds, x)
        predictions[f"{name.lower()}_scgpt_probability"] = p
        summaries.append({"cohort": "GSE135779", "lineage": name,
                          "annotation": "project_fixed_marker_score_not_author_cell_labels",
                          "auc": roc_auc_score(y, p), "ci_low": lo, "ci_high": hi,
                          "permutation_p": perm_p, "min_cells_per_donor": int(n_lineage[name].min()),
                          "median_cells_per_donor": float(np.median(n_lineage[name]))})
        np.savez_compressed(OUT / f"GSE135779_{name.lower()}_scgpt_resampling.npz",
                            bootstrap_auc=boot, permutation_null=null)
    pd.DataFrame(predictions).to_csv(OUT / "GSE135779_project_marker_lineage_oof_predictions.csv", index=False)
    pd.DataFrame(summaries).to_csv(OUT / "GSE135779_project_marker_lineage_summary.csv", index=False)

    # Raw-count donor sums for a transparent, optional T-lineage DE rerun.
    t_sums = np.zeros((len(donors), a.n_vars), dtype=np.float64)
    layer = a.layers["counts"] if "counts" in a.layers else a.raw.X
    for start in range(0, a.n_obs, chunk):
        stop = min(start + chunk, a.n_obs)
        mask = annotation[start:stop] == lineage_names.index("T")
        if not mask.any():
            continue
        x = layer[start:stop][mask]
        if not sparse.issparse(x):
            x = sparse.csr_matrix(x)
        cc = codes[start:stop][mask]
        selector = sparse.csr_matrix((np.ones(mask.sum()), (cc, np.arange(mask.sum()))),
                                     shape=(len(donors), mask.sum()))
        t_sums += (selector @ x).toarray()
    np.savez_compressed(DATA / "GSE135779_project_marker_T_raw_count_sums.npz",
                        donors=donors, y=y, counts=t_sums.astype(np.int64), genes=genes)
    pd.DataFrame(t_sums.astype(np.int64).T, index=genes, columns=donors).to_csv(
        DATA / "GSE135779_project_marker_T_raw_count_sums.csv.gz",
        index_label="gene", compression="gzip"
    )
    a.file.close()
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
