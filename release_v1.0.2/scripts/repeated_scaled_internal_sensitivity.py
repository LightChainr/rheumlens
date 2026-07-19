#!/usr/bin/env python3
"""Repeat donor-level CV after fold-contained scaling of every feature block."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


PROJECT = Path(__file__).resolve().parents[2]
ARCHIVE = (
    PROJECT
    / "03_远程回传"
    / "final_gpu_figure_ready_20260709"
    / "extracted_plus"
    / "RheumLens_GPU_figure_ready_plus_20260709"
    / "extra"
)
OUT = PROJECT / "11_advanced_visualizations_20260717" / "source_data" / "scaled_internal_sensitivity"

DATASETS = ("SLE_GSE135779", "SLE_GSE174188_CD4", "SLE_GSE285773_CD4")
METHODS = ("Frozen Geneformer", "Mean-HVG pseudobulk", "PCA pseudobulk")
N_REPEATS = 20
N_SPLITS = 5
N_HVG = 4000
N_PC = 30
SEED_START = 20_260_801


def paths(dataset_id: str) -> tuple[Path, Path, Path]:
    method = "geneformer_v2_316m_cell_sample1000_clspool_logistic_maxlen4096_seed001"
    return (
        ARCHIVE / "04_models" / "Geneformer" / dataset_id / method / "donor_embedding.parquet",
        ARCHIVE / "pseudobulk" / dataset_id / "donor_log1p_cpm.parquet",
        ARCHIVE / "pseudobulk" / dataset_id / "donor_labels.tsv",
    )


def load_dataset(dataset_id: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    embedding_path, pseudobulk_path, label_path = paths(dataset_id)
    embedding = pd.read_parquet(embedding_path)
    pseudobulk = pd.read_parquet(pseudobulk_path)
    embedding.index = embedding.index.astype(str)
    pseudobulk.index = pseudobulk.index.astype(str)
    labels = pd.read_csv(label_path, sep="\t", dtype={"donor_id": str}).set_index("donor_id")["case_control"]
    donors = sorted(set(embedding.index) & set(pseudobulk.index) & set(labels.index))
    y = labels.loc[donors].eq("case").astype(int).to_numpy()
    return (
        embedding.loc[donors].to_numpy(dtype=np.float64),
        pseudobulk.loc[donors].to_numpy(dtype=np.float64),
        y,
        donors,
    )


def classifier(seed: int) -> LogisticRegression:
    return LogisticRegression(
        C=1.0,
        class_weight="balanced",
        max_iter=5000,
        random_state=seed,
        solver="liblinear",
    )


def tuned_classifier(seed: int) -> LogisticRegressionCV:
    inner = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    return LogisticRegressionCV(
        Cs=np.logspace(-4, 4, 9),
        class_weight="balanced",
        cv=inner,
        max_iter=20000,
        random_state=seed,
        scoring="roc_auc",
        solver="liblinear",
    )


def top_hvg(train_x: np.ndarray) -> np.ndarray:
    n_features = min(N_HVG, train_x.shape[1])
    variance = train_x.var(axis=0, dtype=np.float64)
    selected = np.argpartition(variance, -n_features)[-n_features:]
    return selected[np.argsort(variance[selected])[::-1]]


def scaled_fit_predict(
    train_x: np.ndarray,
    test_x: np.ndarray,
    y_train: np.ndarray,
    seed: int,
    tune_c: bool,
) -> np.ndarray:
    scaler = StandardScaler().fit(train_x)
    model = (tuned_classifier(seed) if tune_c else classifier(seed)).fit(
        scaler.transform(train_x), y_train
    )
    return model.predict_proba(scaler.transform(test_x))[:, 1]


def evaluate_repeat(
    dataset_id: str,
    embedding: np.ndarray,
    pseudobulk: np.ndarray,
    y: np.ndarray,
    donors: list[str],
    repeat: int,
    tune_c: bool,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    split_seed = SEED_START + repeat
    splitter = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=split_seed)
    probabilities = {method: np.zeros(len(y), dtype=float) for method in METHODS}
    fold_ids = np.zeros(len(y), dtype=int)

    for fold, (train, test) in enumerate(splitter.split(embedding, y), start=1):
        seed = split_seed * 10 + fold
        probabilities["Frozen Geneformer"][test] = scaled_fit_predict(
            embedding[train], embedding[test], y[train], seed, tune_c
        )

        hvg = top_hvg(pseudobulk[train])
        train_hvg = pseudobulk[train][:, hvg]
        test_hvg = pseudobulk[test][:, hvg]
        probabilities["Mean-HVG pseudobulk"][test] = scaled_fit_predict(
            train_hvg, test_hvg, y[train], seed, tune_c
        )

        gene_scaler = StandardScaler().fit(train_hvg)
        train_scaled = gene_scaler.transform(train_hvg)
        test_scaled = gene_scaler.transform(test_hvg)
        n_components = min(N_PC, len(train) - 1, train_scaled.shape[1])
        pca = PCA(n_components=n_components, random_state=seed, svd_solver="randomized")
        train_pc = pca.fit_transform(train_scaled)
        test_pc = pca.transform(test_scaled)
        probabilities["PCA pseudobulk"][test] = scaled_fit_predict(
            train_pc, test_pc, y[train], seed, tune_c
        )
        fold_ids[test] = fold

    metrics: list[dict[str, object]] = []
    predictions: list[dict[str, object]] = []
    for method, values in probabilities.items():
        metrics.append(
            {
                "dataset_id": dataset_id,
                "repeat": repeat + 1,
                "split_seed": split_seed,
                "method": method,
                "roc_auc": roc_auc_score(y, values),
                "pr_auc": average_precision_score(y, values),
                "brier": brier_score_loss(y, values),
            }
        )
        predictions.extend(
            {
                "dataset_id": dataset_id,
                "repeat": repeat + 1,
                "split_seed": split_seed,
                "method": method,
                "donor_id": donor,
                "fold": int(fold),
                "y_true": int(label),
                "prob_case": float(probability),
            }
            for donor, fold, label, probability in zip(donors, fold_ids, y, values)
        )
    return metrics, predictions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tune-c", action="store_true")
    args = parser.parse_args()
    output = OUT / ("tuned_c" if args.tune_c else "fixed_c1")
    output.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for dataset_id in DATASETS:
        embedding, pseudobulk, y, donors = load_dataset(dataset_id)
        results = Parallel(n_jobs=8, verbose=10)(
            delayed(evaluate_repeat)(dataset_id, embedding, pseudobulk, y, donors, repeat, args.tune_c)
            for repeat in range(N_REPEATS)
        )
        for metrics, predictions in results:
            metric_rows.extend(metrics)
            prediction_rows.extend(predictions)

    metrics = pd.DataFrame(metric_rows)
    predictions = pd.DataFrame(prediction_rows)
    summary = (
        metrics.groupby(["dataset_id", "method"], as_index=False)
        .agg(
            auc_mean=("roc_auc", "mean"),
            auc_sd=("roc_auc", "std"),
            auc_q025=("roc_auc", lambda x: x.quantile(0.025)),
            auc_q975=("roc_auc", lambda x: x.quantile(0.975)),
            pr_auc_mean=("pr_auc", "mean"),
            brier_mean=("brier", "mean"),
        )
    )
    metrics.to_csv(output / "repeat_level_metrics.csv", index=False)
    predictions.to_csv(output / "repeat_level_oof_predictions.csv", index=False)
    summary.to_csv(output / "summary.csv", index=False)
    (output / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "analysis": "20 repeated stratified five-fold donor CV with outer-training-only scaling of every final feature block",
                "n_repeats": N_REPEATS,
                "n_splits": N_SPLITS,
                "split_seeds": list(range(SEED_START, SEED_START + N_REPEATS)),
                "classifier": (
                    "balanced liblinear logistic regression; C selected from 1e-4 to 1e4 by five-fold inner CV"
                    if args.tune_c
                    else "balanced liblinear logistic regression, C=1.0"
                ),
                "pca": "HVG scaling, training-fold PCA, then training-fold scaling of PC scores",
                "purpose": "scale sensitivity for the archived unscaled fixed-C internal comparison; external source-only transfer is unchanged",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
