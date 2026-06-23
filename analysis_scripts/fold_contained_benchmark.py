#!/usr/bin/env python3
"""Leakage-resistant donor-level benchmark for RheumLens.

Core guarantees:
- donor is the statistical unit;
- every fitted transform (HVG selection, PCA, scaling, IFN residualization)
  is estimated on training donors only;
- all methods share exactly the same held-out donors;
- AUC differences use paired, outcome-stratified donor bootstrap resampling.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "fold_contained"
FEATURE_DIR = ROOT / "data" / "donor_features"

SSD_PULL = Path(
    "/Volumes/Mac Data/Research/论文资料库/04_来自Downloads_2026-06-17/"
    "P3_GPU服务器与MD/rheumlens_gpu_pull_20260617"
)

COHORTS = {
    "GSE135779": {
        "h5ad": SSD_PULL / "gse135779_child_embeddings.h5ad",
        "embedding": "obsm:X_scGPT",
        "label_col": "disease",
        "positive": "SLE",
        "folds": ROOT / "data" / "folds" / "gse135779_child_donor_folds.json",
    },
    "GSE285773": {
        "h5ad": SSD_PULL / "gse285773_embeddings.h5ad",
        "embedding": "obsm:X_scGPT",
        "label_col": "disease",
        "positive": "SLE",
        "folds": ROOT / "data" / "folds" / "gse285773_donor_folds.json",
    },
    "GSE174188": {
        "h5ad": Path(
            "/Volumes/Mac Data/Research/论文项目总库/07_RheumLens_scGPT/"
            "05_原始数据/GSE174188_CELLxGENE_2025-11-08.h5ad"
        ),
        "embedding_h5ad": ROOT / "data" / "gse174188" / "gse174188_cd4_scgpt_embeddings.h5ad",
        "embedding": "X",
        "label_col": "disease_bin",
        "positive": 1,
        "folds": ROOT / "data" / "folds" / "gse174188_cd4_donor_folds.json",
    },
}

ISG_GENES = [
    "ISG15", "IFI6", "MX1", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT3",
    "IFI44", "IFI44L", "STAT1", "RSAD2", "IFITM1", "IFITM3", "HERC5",
]


def log(message: str) -> None:
    print(message, flush=True)


def normalize_label(values: pd.Series, positive) -> np.ndarray:
    if pd.api.types.is_categorical_dtype(values.dtype):
        values = values.astype(object)
    if isinstance(positive, str):
        return (values.astype(str).str.upper() == str(positive).upper()).astype(int).to_numpy()
    return (pd.to_numeric(values) == positive).astype(int).to_numpy()


def gene_names(adata: ad.AnnData) -> np.ndarray:
    if "gene_symbol" in adata.var.columns:
        names = adata.var["gene_symbol"].astype(str).to_numpy()
    elif "feature_name" in adata.var.columns:
        names = adata.var["feature_name"].astype(str).to_numpy()
    else:
        names = adata.var_names.astype(str).to_numpy()
    # Stable uniqueness is needed for auditable selected-HVG exports.
    seen: dict[str, int] = {}
    out = []
    for name in names:
        count = seen.get(name, 0)
        out.append(name if count == 0 else f"{name}__dup{count}")
        seen[name] = count + 1
    return np.asarray(out, dtype=object)


def aggregate_matrix_by_donor(matrix, donor_codes: np.ndarray, n_donors: int,
                              chunk_size: int = 8192) -> np.ndarray:
    n_obs, n_vars = matrix.shape
    sums = np.zeros((n_donors, n_vars), dtype=np.float64)
    counts = np.bincount(donor_codes, minlength=n_donors).astype(np.int64)
    for start in range(0, n_obs, chunk_size):
        stop = min(n_obs, start + chunk_size)
        chunk = matrix[start:stop]
        if not sparse.issparse(chunk):
            dense = np.asarray(chunk)
            if dense.dtype == np.float16:
                dense = dense.astype(np.float32)
            chunk = sparse.csr_matrix(dense)
        rows = donor_codes[start:stop]
        selector = sparse.csr_matrix(
            (np.ones(stop - start, dtype=np.float32), (rows, np.arange(stop - start))),
            shape=(n_donors, stop - start),
        )
        sums += (selector @ chunk).toarray()
        if start == 0 or stop == n_obs or (start // chunk_size) % 10 == 0:
            log(f"    aggregate rows {stop:,}/{n_obs:,}")
    if np.any(counts == 0):
        raise ValueError("At least one donor has zero cells")
    return (sums / counts[:, None]).astype(np.float32)


def aggregate_embedding(adata: ad.AnnData, location: str, donor_codes: np.ndarray,
                        n_donors: int) -> np.ndarray:
    if location == "X":
        matrix = adata.X
    elif location.startswith("obsm:"):
        matrix = adata.obsm[location.split(":", 1)[1]]
    else:
        raise ValueError(location)
    return aggregate_matrix_by_donor(matrix, donor_codes, n_donors, chunk_size=16384)


def prepare_gse174188_features(cfg: dict, out: Path) -> Path:
    """Recreate the archived Census preprocessing from direct CELLxGENE raw counts."""
    log(f"  opening official CELLxGENE object: {cfg['h5ad']}")
    a = ad.read_h5ad(cfg["h5ad"], backed="r")
    cd4_mask = a.obs["cell_type"].astype(str).str.contains("CD4", case=False, na=False).to_numpy()
    cd4_idx = np.flatnonzero(cd4_mask)
    raw_x = a.raw.X
    log(f"  CD4 before QC: {len(cd4_idx):,}")

    # Match the archived Scanpy filter_cells(min_genes=200).
    keep_cell_parts = []
    gene_detected = np.zeros(raw_x.shape[1], dtype=np.int64)
    chunk_size = 4096
    for start in range(0, len(cd4_idx), chunk_size):
        idx = cd4_idx[start:start + chunk_size]
        x = raw_x[idx].tocsr()
        keep = np.asarray(x.getnnz(axis=1)).ravel() >= 200
        kept_idx = idx[keep]
        keep_cell_parts.append(kept_idx)
        if keep.any():
            gene_detected += np.asarray((x[keep] > 0).sum(axis=0)).ravel().astype(np.int64)
        if start == 0 or start + chunk_size >= len(cd4_idx) or (start // chunk_size) % 20 == 0:
            log(f"    QC pass rows {min(start + chunk_size, len(cd4_idx)):,}/{len(cd4_idx):,}")
    kept_cells = np.concatenate(keep_cell_parts)
    keep_genes = np.flatnonzero(gene_detected >= 3)
    log(f"  after QC: cells={len(kept_cells):,}, genes={len(keep_genes):,}")

    obs = a.obs.iloc[kept_cells]
    donor_ids_cell = obs["donor_id"].astype(str).to_numpy()
    donors = np.asarray(sorted(pd.unique(donor_ids_cell)), dtype=object)
    donor_to_code = {d: i for i, d in enumerate(donors)}
    codes = np.asarray([donor_to_code[d] for d in donor_ids_cell], dtype=np.int32)
    labels_cell = (obs["disease"].astype(str) == "systemic lupus erythematosus").astype(int).to_numpy()
    labels = np.zeros(len(donors), dtype=np.int8)
    for i in range(len(donors)):
        unique = np.unique(labels_cell[codes == i])
        if len(unique) != 1:
            raise ValueError(f"donor {donors[i]} has inconsistent disease labels")
        labels[i] = unique[0]
    log(f"  donors={len(donors)} SLE={int(labels.sum())} HC={int((labels == 0).sum())}")

    # Per-cell normalize_total(1e4)+log1p, then donor mean-pool.
    sums = np.zeros((len(donors), len(keep_genes)), dtype=np.float64)
    counts = np.bincount(codes, minlength=len(donors)).astype(np.int64)
    for start in range(0, len(kept_cells), chunk_size):
        stop = min(start + chunk_size, len(kept_cells))
        idx = kept_cells[start:stop]
        x = raw_x[idx].tocsr().astype(np.float64)
        totals = np.asarray(x.sum(axis=1)).ravel()
        if np.any(totals <= 0):
            raise ValueError("zero-count cell survived min_genes QC")
        x = sparse.diags(1e4 / totals) @ x
        x.data = np.log1p(x.data)
        x = x[:, keep_genes]
        rows = codes[start:stop]
        selector = sparse.csr_matrix(
            (np.ones(stop - start), (rows, np.arange(stop - start))),
            shape=(len(donors), stop - start),
        )
        sums += (selector @ x).toarray()
        if start == 0 or stop == len(kept_cells) or (start // chunk_size) % 20 == 0:
            log(f"    normalize/aggregate rows {stop:,}/{len(kept_cells):,}")
    expression = (sums / counts[:, None]).astype(np.float32)

    raw_var = a.raw.var.iloc[keep_genes]
    if "feature_name" in raw_var.columns:
        genes = raw_var["feature_name"].astype(str).to_numpy()
    else:
        genes = raw_var.index.astype(str).to_numpy()
    # Stable duplicate suffixes.
    seen = {}
    genes_unique = []
    for name in genes:
        count = seen.get(name, 0)
        genes_unique.append(name if count == 0 else f"{name}__dup{count}")
        seen[name] = count + 1
    genes = np.asarray(genes_unique, dtype=object)
    a.file.close()

    e = ad.read_h5ad(cfg["embedding_h5ad"], backed="r")
    emb_donors_cell = e.obs["donor_id"].astype(str).to_numpy()
    if set(emb_donors_cell) != set(donors):
        raise ValueError("raw-expression and archived embedding donor sets differ")
    emb_codes = np.asarray([donor_to_code[d] for d in emb_donors_cell], dtype=np.int32)
    log(f"  aggregating archived scGPT embedding: {e.shape}")
    scgpt = aggregate_embedding(e, cfg["embedding"], emb_codes, len(donors))
    e.file.close()

    np.savez_compressed(out, donors=donors, y=labels, expression=expression,
                        scgpt=scgpt, genes=genes)
    log(f"  wrote {out}")
    return out


def prepare_features(cohort: str, force: bool = False) -> Path:
    cfg = COHORTS[cohort]
    FEATURE_DIR.mkdir(parents=True, exist_ok=True)
    out = FEATURE_DIR / f"{cohort}_donor_features.npz"
    if out.exists() and not force:
        log(f"  using cached donor features: {out}")
        return out

    if cohort == "GSE174188":
        return prepare_gse174188_features(cfg, out)

    if not Path(cfg["h5ad"]).exists():
        raise FileNotFoundError(cfg["h5ad"])
    log(f"  opening expression object: {cfg['h5ad']}")
    a = ad.read_h5ad(cfg["h5ad"], backed="r")
    donor_ids_cell = a.obs["donor_id"].astype(str).to_numpy()
    donors = np.asarray(sorted(pd.unique(donor_ids_cell)), dtype=object)
    donor_to_code = {d: i for i, d in enumerate(donors)}
    codes = np.asarray([donor_to_code[d] for d in donor_ids_cell], dtype=np.int32)
    labels_cell = normalize_label(a.obs[cfg["label_col"]], cfg["positive"])
    labels = np.zeros(len(donors), dtype=np.int8)
    for i in range(len(donors)):
        unique = np.unique(labels_cell[codes == i])
        if len(unique) != 1:
            raise ValueError(f"donor {donors[i]} has inconsistent labels: {unique}")
        labels[i] = unique[0]

    log(f"  donors={len(donors)} SLE={int(labels.sum())} HC={int((labels == 0).sum())}")
    log(f"  aggregating log-normalized expression {a.shape}")
    expression = aggregate_matrix_by_donor(a.X, codes, len(donors))
    genes = gene_names(a)

    if cohort == "GSE174188":
        emb_path = Path(cfg["embedding_h5ad"])
        e = ad.read_h5ad(emb_path, backed="r")
        emb_donor_cell = e.obs["donor_id"].astype(str).to_numpy()
        if set(emb_donor_cell) != set(donors):
            raise ValueError("GSE174188 raw-expression and embedding donor sets differ")
        emb_codes = np.asarray([donor_to_code[d] for d in emb_donor_cell], dtype=np.int32)
        log(f"  aggregating embedding {e.shape}")
        scgpt = aggregate_embedding(e, cfg["embedding"], emb_codes, len(donors))
        e.file.close()
    else:
        log("  aggregating scGPT embedding")
        scgpt = aggregate_embedding(a, cfg["embedding"], codes, len(donors))

    a.file.close()
    np.savez_compressed(
        out, donors=donors, y=labels, expression=expression, scgpt=scgpt,
        genes=genes,
    )
    log(f"  wrote {out}")
    return out


def folds_from_json(path: Path, donors: np.ndarray, y: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    donor_to_idx = {d: i for i, d in enumerate(donors)}
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        payload = {}
        for fold, (tr, te) in enumerate(skf.split(donors, y)):
            payload[str(fold)] = {
                "train_donors": donors[tr].tolist(),
                "test_donors": donors[te].tolist(),
                "source": "StratifiedKFold(n_splits=5, shuffle=True, random_state=42)",
            }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        log(f"  created deterministic folds: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    out = []
    all_test = []
    for key in sorted(payload, key=lambda x: int(x.replace("fold", "")) if x.replace("fold", "").isdigit() else x):
        fd = payload[key]
        test_names = fd.get("test_donors", fd) if isinstance(fd, dict) else fd
        te = np.asarray([donor_to_idx[str(d)] for d in test_names], dtype=int)
        if isinstance(fd, dict) and "train_donors" in fd:
            tr = np.asarray([donor_to_idx[str(d)] for d in fd["train_donors"]], dtype=int)
        else:
            tr = np.setdiff1d(np.arange(len(donors)), te)
        out.append((tr, te))
        all_test.extend(te.tolist())
    if sorted(all_test) != list(range(len(donors))):
        raise ValueError("Fold test sets do not partition donors exactly once")
    return out


def classifier() -> LogisticRegression:
    return LogisticRegression(
        max_iter=5000, class_weight="balanced", solver="liblinear",
        C=1.0, random_state=42,
    )


def select_hvg(x_train: np.ndarray, n_hvg: int) -> np.ndarray:
    variance = np.var(x_train, axis=0, ddof=1)
    # mergesort makes ties deterministic by retaining source-column order.
    return np.argsort(-variance, kind="mergesort")[: min(n_hvg, x_train.shape[1])]


def run_oof(expression: np.ndarray, scgpt: np.ndarray, y: np.ndarray,
            folds: list[tuple[np.ndarray, np.ndarray]], n_hvg: int = 2000,
            n_pcs: int = 25, isg_score: np.ndarray | None = None):
    predictions = {
        "scgpt": np.full(len(y), np.nan),
        "expression_pca": np.full(len(y), np.nan),
        "hvg_pseudobulk": np.full(len(y), np.nan),
    }
    if isg_score is not None:
        predictions["scgpt_ifn_residual"] = np.full(len(y), np.nan)
    selected = []
    for fold_id, (tr, te) in enumerate(folds):
        log(f"  fold {fold_id}: train={len(tr)} test={len(te)}")

        sc_pipe = make_pipeline(StandardScaler(), classifier())
        sc_pipe.fit(scgpt[tr], y[tr])
        predictions["scgpt"][te] = sc_pipe.predict_proba(scgpt[te])[:, 1]

        hvg = select_hvg(expression[tr], n_hvg)
        selected.append(hvg)

        pca = PCA(
            n_components=min(n_pcs, len(tr) - 1, len(hvg)),
            svd_solver="full", random_state=42,
        )
        train_pc = pca.fit_transform(expression[tr][:, hvg])
        test_pc = pca.transform(expression[te][:, hvg])
        pca_pipe = make_pipeline(StandardScaler(), classifier())
        pca_pipe.fit(train_pc, y[tr])
        predictions["expression_pca"][te] = pca_pipe.predict_proba(test_pc)[:, 1]

        hvg_pipe = make_pipeline(StandardScaler(), classifier())
        hvg_pipe.fit(expression[tr][:, hvg], y[tr])
        predictions["hvg_pseudobulk"][te] = hvg_pipe.predict_proba(expression[te][:, hvg])[:, 1]

        if isg_score is not None:
            resid = LinearRegression()
            resid.fit(isg_score[tr, None], scgpt[tr])
            train_resid = scgpt[tr] - resid.predict(isg_score[tr, None])
            test_resid = scgpt[te] - resid.predict(isg_score[te, None])
            resid_pipe = make_pipeline(StandardScaler(), classifier())
            resid_pipe.fit(train_resid, y[tr])
            predictions["scgpt_ifn_residual"][te] = resid_pipe.predict_proba(test_resid)[:, 1]

    for name, values in predictions.items():
        if np.isnan(values).any():
            raise ValueError(f"missing OOF predictions for {name}")
    return predictions, selected


def residual_permutation_test(scgpt: np.ndarray, isg_score: np.ndarray, y: np.ndarray,
                              folds: list[tuple[np.ndarray, np.ndarray]],
                              observed_auc: float, n_perm: int = 1000,
                              seed: int = 42):
    """Permutation test with residualization and scaling fitted in each train fold."""
    prepared = []
    for tr, te in folds:
        resid = LinearRegression().fit(isg_score[tr, None], scgpt[tr])
        x_tr = scgpt[tr] - resid.predict(isg_score[tr, None])
        x_te = scgpt[te] - resid.predict(isg_score[te, None])
        scaler = StandardScaler().fit(x_tr)
        prepared.append((tr, te, scaler.transform(x_tr), scaler.transform(x_te)))
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=np.float64)
    for i in range(n_perm):
        yp = rng.permutation(y)
        probs = np.full(len(y), np.nan)
        valid = True
        for tr, te, x_tr, x_te in prepared:
            if np.unique(yp[tr]).size < 2:
                valid = False
                break
            model = classifier().fit(x_tr, yp[tr])
            probs[te] = model.predict_proba(x_te)[:, 1]
        null[i] = roc_auc_score(yp, probs) if valid else 0.5
        if (i + 1) % 200 == 0:
            log(f"    residual permutation {i + 1}/{n_perm}")
    p = (1 + np.sum(null >= observed_auc)) / (n_perm + 1)
    return float(p), null


def stratified_bootstrap(y: np.ndarray, predictions: dict[str, np.ndarray],
                         n_boot: int = 10000, seed: int = 20260618):
    rng = np.random.default_rng(seed)
    cases = np.flatnonzero(y == 1)
    controls = np.flatnonzero(y == 0)
    names = list(predictions)
    samples = {name: np.empty(n_boot) for name in names}
    for i in range(n_boot):
        idx = np.concatenate([
            rng.choice(cases, size=len(cases), replace=True),
            rng.choice(controls, size=len(controls), replace=True),
        ])
        for name in names:
            samples[name][i] = roc_auc_score(y[idx], predictions[name][idx])
    rows = []
    for name in names:
        arr = samples[name]
        rows.append({
            "comparison": name,
            "estimate": roc_auc_score(y, predictions[name]),
            "ci_low": np.quantile(arr, 0.025),
            "ci_high": np.quantile(arr, 0.975),
            "bootstrap_replicates": n_boot,
        })
    for left, right in (
        ("expression_pca", "scgpt"),
        ("hvg_pseudobulk", "scgpt"),
        ("hvg_pseudobulk", "expression_pca"),
    ):
        delta = samples[left] - samples[right]
        rows.append({
            "comparison": f"{left}_minus_{right}",
            "estimate": roc_auc_score(y, predictions[left]) - roc_auc_score(y, predictions[right]),
            "ci_low": np.quantile(delta, 0.025),
            "ci_high": np.quantile(delta, 0.975),
            "bootstrap_replicates": n_boot,
        })
    return pd.DataFrame(rows), samples


def ifn_score(expression: np.ndarray, genes: np.ndarray):
    lookup = {g.split("__dup", 1)[0]: i for i, g in enumerate(genes)}
    available = [g for g in ISG_GENES if g in lookup]
    if len(available) < 8:
        raise ValueError(f"Only {len(available)} ISGs available: {available}")
    idx = [lookup[g] for g in available]
    return expression[:, idx].mean(axis=1), available


def evaluate(cohort: str, force_features: bool = False, n_boot: int = 10000):
    path = prepare_features(cohort, force=force_features)
    z = np.load(path, allow_pickle=True)
    donors = z["donors"].astype(str)
    y = z["y"].astype(int)
    # Float64 avoids BLAS overflow warnings observed with small-magnitude float32
    # donor matrices under recent NumPy/scikit-learn builds.
    expression = z["expression"].astype(np.float64)
    scgpt = z["scgpt"].astype(np.float64)
    genes = z["genes"].astype(str)
    folds = folds_from_json(Path(COHORTS[cohort]["folds"]), donors, y)
    isg = None
    available_isg = []
    if cohort == "GSE135779":
        isg, available_isg = ifn_score(expression, genes)
        log(f"  IFN residualization genes ({len(available_isg)}): {available_isg}")

    predictions, selected = run_oof(expression, scgpt, y, folds, isg_score=isg)
    cis, boot = stratified_bootstrap(y, predictions, n_boot=n_boot)
    cis.insert(0, "cohort", cohort)

    OUT.mkdir(parents=True, exist_ok=True)
    pred_df = pd.DataFrame({"cohort": cohort, "donor_id": donors, "label": y, **predictions})
    pred_df.to_csv(OUT / f"{cohort}_oof_predictions.csv", index=False)
    cis.to_csv(OUT / f"{cohort}_auc_bootstrap_ci.csv", index=False)
    np.savez_compressed(OUT / f"{cohort}_bootstrap_distributions.npz", **boot)

    residual_perm_p = None
    if isg is not None:
        observed = roc_auc_score(y, predictions["scgpt_ifn_residual"])
        residual_perm_p, residual_null = residual_permutation_test(
            scgpt, isg, y, folds, observed_auc=observed, n_perm=1000,
        )
        np.save(OUT / f"{cohort}_ifn_residual_permutation_null.npy", residual_null)
        pd.DataFrame([{
            "cohort": cohort, "observed_auc": observed,
            "permutation_p": residual_perm_p, "permutations": 1000,
            "null_mean": float(residual_null.mean()),
            "null_sd": float(residual_null.std(ddof=1)),
        }]).to_csv(OUT / f"{cohort}_ifn_residual_permutation.csv", index=False)

    hvg_rows = []
    for fold_id, hvg in enumerate(selected):
        hvg_rows.extend({"cohort": cohort, "fold": fold_id, "rank": rank + 1,
                         "gene": genes[idx], "column_index": int(idx)}
                        for rank, idx in enumerate(hvg))
    pd.DataFrame(hvg_rows).to_csv(OUT / f"{cohort}_fold_selected_hvgs.csv", index=False)

    meta = {
        "cohort": cohort,
        "n_donors": int(len(y)),
        "n_sle": int(y.sum()),
        "n_hc": int((y == 0).sum()),
        "n_genes": int(expression.shape[1]),
        "n_hvg_per_fold": 2000,
        "n_pcs": 25,
        "classifier": "StandardScaler -> LogisticRegression(L2,C=1,class_weight=balanced,solver=liblinear)",
        "pca": "training-fold HVG variance selection -> training-fold PCA -> training-fold PC scaling",
        "hvg_pseudobulk": "training-fold HVG variance selection -> training-fold gene scaling",
        "expression_representation": "donor mean of per-cell normalize_total(1e4)+log1p expression",
        "bootstrap": f"paired outcome-stratified donor bootstrap, {n_boot} replicates, seed=20260618",
        "ifn_genes_available": available_isg,
        "ifn_residual_permutation_p": residual_perm_p,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }
    (OUT / f"{cohort}_run_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log("\n" + cis.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    return cis


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", action="append", choices=COHORTS, required=True)
    parser.add_argument("--force-features", action="store_true")
    parser.add_argument("--bootstrap", type=int, default=10000)
    args = parser.parse_args()
    all_ci = []
    for cohort in args.cohort:
        log(f"\n=== {cohort} ===")
        all_ci.append(evaluate(cohort, force_features=args.force_features, n_boot=args.bootstrap))
    OUT.mkdir(parents=True, exist_ok=True)
    pd.concat(all_ci, ignore_index=True).to_csv(OUT / "all_cohorts_auc_bootstrap_ci.csv", index=False)


if __name__ == "__main__":
    main()
