#!/usr/bin/env python3
"""Matched evaluation of valid Geneformer V2 and scGPT embeddings."""

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


BASE = Path("/autodl-fs/data/rheumlens_20260619")
RESULT = BASE / "results/geneformer_valid"
GF_CSV = RESULT / "GSE135779_v2_104M_all.csv"
SCGPT_NPZ = BASE / "data/geneformer_rerun/GSE135779_matched_scgpt_embeddings.npz"
FOLDS = BASE / "data/geneformer_rerun/gse135779_child_donor_folds.json"
SEED = 20260619


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def model():
    return make_pipeline(StandardScaler(), LogisticRegression(
        C=1.0, class_weight="balanced", solver="liblinear", max_iter=5000,
        random_state=42))


def fold_indices(donors):
    spec = json.loads(FOLDS.read_text())
    lookup = {str(d): i for i, d in enumerate(donors)}
    out = []
    for key in sorted(spec, key=int):
        f = spec[key]
        out.append((np.asarray([lookup[str(x)] for x in f["train_donors"]]),
                    np.asarray([lookup[str(x)] for x in f["test_donors"]])))
    return out


def oof(x, y, folds):
    p = np.full(len(y), np.nan)
    for tr, te in folds:
        m = model().fit(x[tr], y[tr])
        p[te] = m.predict_proba(x[te])[:, 1]
    return p


def resample(y, pred, x, folds, n_boot=10000, n_perm=1000):
    rng = np.random.default_rng(SEED)
    case, ctrl = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.r_[rng.choice(case, len(case), True), rng.choice(ctrl, len(ctrl), True)]
        boot[i] = roc_auc_score(y[idx], pred[idx])
    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = rng.permutation(y)
        null[i] = roc_auc_score(yp, oof(x, yp, folds))
    observed = roc_auc_score(y, pred)
    return boot, null, (int((null >= observed).sum()) + 1) / (n_perm + 1)


def main():
    gf = pd.read_csv(GF_CSV).drop(columns=["Unnamed: 0"], errors="ignore")
    emb_cols = [str(i) for i in range(768)]
    gf_x = gf[emb_cols].to_numpy(np.float32)
    source_idx = gf["source_row_index"].to_numpy(np.int64)

    z = np.load(SCGPT_NPZ, allow_pickle=True)
    pos = {int(v): i for i, v in enumerate(z["source_row_index"])}
    match = np.asarray([pos[int(v)] for v in source_idx])
    sc_x = z["embeddings"][match].astype(np.float32)
    assert np.array_equal(gf["donor_id"].astype(str), z["donor_id"][match].astype(str))
    assert np.array_equal(gf["marker_lineage"].astype(str), z["marker_lineage"][match].astype(str))

    labels = gf["marker_lineage"].astype(str).to_numpy()
    keep = labels != "Unassigned"
    ari_rows = []
    for representation, x in [("Geneformer_V2_104M_CLS", gf_x), ("scGPT", sc_x)]:
        for scaling in ["raw", "standardized"]:
            xx = x[keep] if scaling == "raw" else StandardScaler().fit_transform(x[keep])
            cluster = KMeans(n_clusters=4, n_init=20, random_state=SEED).fit_predict(xx)
            ari_rows.append({"representation": representation, "scaling": scaling,
                             "n_cells": int(keep.sum()),
                             "ari_vs_project_marker_lineage": adjusted_rand_score(labels[keep], cluster)})
    pd.DataFrame(ari_rows).to_csv(RESULT / "matched_lineage_ari.csv", index=False)

    donors = np.asarray(sorted(pd.unique(gf["donor_id"].astype(str))))
    disease = gf.drop_duplicates("donor_id").set_index("donor_id")["disease"].astype(str)
    y = disease.reindex(donors).str.upper().eq("SLE").astype(int).to_numpy()
    donor_gf = np.vstack([gf_x[gf["donor_id"].astype(str).to_numpy() == d].mean(0) for d in donors])
    donor_sc = np.vstack([sc_x[gf["donor_id"].astype(str).to_numpy() == d].mean(0) for d in donors])
    folds = fold_indices(donors)
    p_gf, p_sc = oof(donor_gf, y, folds), oof(donor_sc, y, folds)
    boot_gf, null_gf, pperm_gf = resample(y, p_gf, donor_gf, folds)
    boot_sc, null_sc, pperm_sc = resample(y, p_sc, donor_sc, folds)
    delta = boot_gf - boot_sc
    rows = []
    for name, pred, boot, pperm, x in [
        ("Geneformer_V2_104M_CLS", p_gf, boot_gf, pperm_gf, donor_gf),
        ("matched_500cell_scGPT", p_sc, boot_sc, pperm_sc, donor_sc)]:
        rows.append({"representation": name, "n_donors": len(y),
                     "embedding_dimensions": x.shape[1], "auc": roc_auc_score(y, pred),
                     "ci_low": np.quantile(boot, .025), "ci_high": np.quantile(boot, .975),
                     "permutation_p": pperm})
    pd.DataFrame(rows).to_csv(RESULT / "matched_donor_auc_summary.csv", index=False)
    pd.DataFrame({"donor_id": donors, "label": y,
                  "geneformer_probability": p_gf, "scgpt_probability": p_sc}).to_csv(
        RESULT / "matched_donor_oof_predictions.csv", index=False)
    pd.DataFrame([{"comparison": "Geneformer_minus_scGPT", "estimate": roc_auc_score(y,p_gf)-roc_auc_score(y,p_sc),
                   "ci_low": np.quantile(delta,.025), "ci_high": np.quantile(delta,.975)}]).to_csv(
        RESULT / "matched_auc_difference.csv", index=False)
    np.savez_compressed(RESULT / "matched_resampling.npz", gf_bootstrap=boot_gf,
                        scgpt_bootstrap=boot_sc, paired_delta=delta,
                        gf_permutation=null_gf, scgpt_permutation=null_sc)
    variance = pd.DataFrame({"representation": ["Geneformer_V2_104M_CLS", "scGPT"],
        "mean_dimension_variance": [gf_x.var(0).mean(), sc_x.var(0).mean()],
        "median_dimension_variance": [np.median(gf_x.var(0)), np.median(sc_x.var(0))],
        "constant_dimensions": [(gf_x.var(0) < 1e-12).sum(), (sc_x.var(0) < 1e-12).sum()]})
    variance.to_csv(RESULT / "matched_embedding_variance.csv", index=False)
    metadata = {"geneformer_csv_sha256": sha256(GF_CSV),
                "matched_scgpt_npz_sha256": sha256(SCGPT_NPZ),
                "n_cells": len(gf), "sampling": "500 cells per donor; seed 20260619",
                "lineage_labels": "project fixed marker score, not author annotations",
                "folds": str(FOLDS), "classifier": "StandardScaler + balanced L2 logistic C=1",
                "bootstrap": 10000, "permutations": 1000, "seed": SEED}
    (RESULT / "evaluation_metadata.json").write_text(json.dumps(metadata, indent=2))
    print(pd.DataFrame(ari_rows).to_string(index=False))
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
