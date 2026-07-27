#!/usr/bin/env python3
"""Matched covariate sensitivity on the final GSE174188 donor universe."""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from common import WORKSPACE, embedding_path, label_path, load_embedding, load_labels, load_pseudobulk, pseudobulk_path, sha256_file


DATASET = "SLE_GSE174188_CD4"
OUT = WORKSPACE / "results" / "02_covariate_audit"
FIGURE = WORKSPACE / "figures" / "figure_covariate_audit"
N_REPEATS = int(os.environ.get("N_REPEATS", "20"))
N_JOBS = min(int(os.environ.get("N_JOBS", "6")), os.cpu_count() or 1)
N_SPLITS = 5
N_HVG = 4_000
N_PC = 30
SEED_START = 20_260_801
RIDGE_ALPHA = 1.0

METHODS = ("Frozen Geneformer", "Mean-HVG pseudobulk", "PCA pseudobulk")
BLOCKS = (
    "cell_yield",
    "qc_core",
    "acquisition_structure",
    "processing_cohort",
    "demographic",
    "combined",
)


def classifier(random_state: int) -> LogisticRegression:
    return LogisticRegression(
        C=1.0,
        class_weight="balanced",
        max_iter=5_000,
        random_state=random_state,
        solver="liblinear",
    )


def top_hvg(train_x: np.ndarray) -> np.ndarray:
    n_features = min(N_HVG, train_x.shape[1])
    variance = train_x.var(axis=0, dtype=np.float64)
    selected = np.argpartition(variance, -n_features)[-n_features:]
    return selected[np.argsort(variance[selected])[::-1]]


def covariate_columns(frame: pd.DataFrame, block: str) -> tuple[list[str], list[str]]:
    processing = sorted(column for column in frame.columns if column.startswith("processing_fraction_"))
    qc_core = [
        "log_cells_per_donor",
        "log_mean_umi_per_cell",
        "log_mean_genes_per_cell",
        "aggregate_pct_mito",
    ]
    acquisition = [
        "log_n_libraries",
        "log_n_samples",
        "log_n_suspensions",
        "n_processing_cohorts",
    ]
    if block == "cell_yield":
        return ["log_cells_per_donor"], []
    if block == "qc_core":
        return qc_core, []
    if block == "acquisition_structure":
        return acquisition, []
    if block == "processing_cohort":
        return processing, []
    if block == "demographic":
        return ["age_years"], ["sex"]
    if block == "combined":
        return [*qc_core, *acquisition, *processing, "age_years"], ["sex"]
    raise ValueError(block)


def make_design(
    covariates: pd.DataFrame, train_index: np.ndarray, test_index: np.ndarray, block: str
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    numeric, categorical = covariate_columns(covariates, block)
    transformers = []
    if numeric:
        transformers.append(("numeric", StandardScaler(), numeric))
    if categorical:
        transformers.append(("categorical", OneHotEncoder(handle_unknown="ignore", sparse_output=False), categorical))
    transformer = ColumnTransformer(transformers, remainder="drop", sparse_threshold=0.0)
    train = transformer.fit_transform(covariates.iloc[train_index])
    test = transformer.transform(covariates.iloc[test_index])
    names = list(transformer.get_feature_names_out())
    return np.asarray(train, dtype=float), np.asarray(test, dtype=float), names


def residualize(
    train_x: np.ndarray, test_x: np.ndarray, train_covariates: np.ndarray, test_covariates: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    model = Ridge(alpha=RIDGE_ALPHA, fit_intercept=True)
    model.fit(train_covariates, train_x)
    return train_x - model.predict(train_covariates), test_x - model.predict(test_covariates)


def fit_representation_fold(
    embedding: np.ndarray,
    pseudobulk: np.ndarray,
    y: np.ndarray,
    train_index: np.ndarray,
    test_index: np.ndarray,
    random_state: int,
) -> dict[str, np.ndarray]:
    output: dict[str, np.ndarray] = {}
    model = classifier(random_state)
    model.fit(embedding[train_index], y[train_index])
    output["Frozen Geneformer"] = model.predict_proba(embedding[test_index])[:, 1]

    hvg = top_hvg(pseudobulk[train_index])
    model = classifier(random_state)
    model.fit(pseudobulk[train_index][:, hvg], y[train_index])
    output["Mean-HVG pseudobulk"] = model.predict_proba(pseudobulk[test_index][:, hvg])[:, 1]

    n_components = min(N_PC, len(train_index) - 1, len(hvg))
    pca = PCA(n_components=n_components, random_state=random_state, svd_solver="randomized")
    train_pc = pca.fit_transform(pseudobulk[train_index][:, hvg])
    test_pc = pca.transform(pseudobulk[test_index][:, hvg])
    model = classifier(random_state)
    model.fit(train_pc, y[train_index])
    output["PCA pseudobulk"] = model.predict_proba(test_pc)[:, 1]
    return output


def fit_residual_fold(
    embedding: np.ndarray,
    pseudobulk: np.ndarray,
    y: np.ndarray,
    train_index: np.ndarray,
    test_index: np.ndarray,
    train_covariates: np.ndarray,
    test_covariates: np.ndarray,
    random_state: int,
) -> dict[str, np.ndarray]:
    output: dict[str, np.ndarray] = {}
    gf_train, gf_test = residualize(
        embedding[train_index], embedding[test_index], train_covariates, test_covariates
    )
    model = classifier(random_state)
    model.fit(gf_train, y[train_index])
    output["Frozen Geneformer"] = model.predict_proba(gf_test)[:, 1]

    pb_train, pb_test = residualize(
        pseudobulk[train_index], pseudobulk[test_index], train_covariates, test_covariates
    )
    hvg = top_hvg(pb_train)
    model = classifier(random_state)
    model.fit(pb_train[:, hvg], y[train_index])
    output["Mean-HVG pseudobulk"] = model.predict_proba(pb_test[:, hvg])[:, 1]

    n_components = min(N_PC, len(train_index) - 1, len(hvg))
    pca = PCA(n_components=n_components, random_state=random_state, svd_solver="randomized")
    train_pc = pca.fit_transform(pb_train[:, hvg])
    test_pc = pca.transform(pb_test[:, hvg])
    model = classifier(random_state)
    model.fit(train_pc, y[train_index])
    output["PCA pseudobulk"] = model.predict_proba(test_pc)[:, 1]
    return output


def evaluate_repeat(
    repeat: int,
    embedding: np.ndarray,
    pseudobulk: np.ndarray,
    y: np.ndarray,
    donors: list[str],
    covariates: pd.DataFrame,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, list[str]]]:
    split_seed = SEED_START + repeat - 1
    splitter = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=split_seed)
    probabilities: dict[tuple[str, str], np.ndarray] = {
        (method, "unadjusted"): np.zeros(len(y), dtype=float) for method in METHODS
    }
    for block in BLOCKS:
        probabilities[(f"Covariates only: {block}", block)] = np.zeros(len(y), dtype=float)
        for method in METHODS:
            probabilities[(method, block)] = np.zeros(len(y), dtype=float)
    folds = np.zeros(len(y), dtype=int)
    design_names: dict[str, list[str]] = {}

    for fold, (train_index, test_index) in enumerate(splitter.split(embedding, y), start=1):
        random_state = split_seed * 10 + fold
        unadjusted = fit_representation_fold(embedding, pseudobulk, y, train_index, test_index, random_state)
        for method, values in unadjusted.items():
            probabilities[(method, "unadjusted")][test_index] = values

        for block in BLOCKS:
            train_cov, test_cov, names = make_design(covariates, train_index, test_index, block)
            design_names.setdefault(block, names)
            cov_model = classifier(random_state)
            cov_model.fit(train_cov, y[train_index])
            probabilities[(f"Covariates only: {block}", block)][test_index] = cov_model.predict_proba(test_cov)[:, 1]
            adjusted = fit_residual_fold(
                embedding,
                pseudobulk,
                y,
                train_index,
                test_index,
                train_cov,
                test_cov,
                random_state,
            )
            for method, values in adjusted.items():
                probabilities[(method, block)][test_index] = values
        folds[test_index] = fold

    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    for (method, adjustment), values in probabilities.items():
        metric_rows.append(
            {
                "repeat": repeat,
                "split_seed": split_seed,
                "method": method,
                "adjustment": adjustment,
                "n_donors": len(y),
                "roc_auc": roc_auc_score(y, values),
                "pr_auc": average_precision_score(y, values),
                "brier": brier_score_loss(y, values),
            }
        )
        prediction_rows.extend(
            {
                "repeat": repeat,
                "split_seed": split_seed,
                "method": method,
                "adjustment": adjustment,
                "donor_id": donor,
                "fold": int(fold),
                "y_true": int(label),
                "prob_case": float(probability),
            }
            for donor, fold, label, probability in zip(donors, folds, y, values)
        )
    return metric_rows, prediction_rows, design_names


def overlap_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for repeat, group in predictions.groupby("repeat"):
        for support_block in BLOCKS:
            propensity = group.loc[
                group["method"].eq(f"Covariates only: {support_block}")
                & group["adjustment"].eq(support_block),
                ["donor_id", "y_true", "prob_case"],
            ].rename(columns={"prob_case": "propensity"})
            case = propensity.loc[propensity.y_true.eq(1), "propensity"]
            control = propensity.loc[propensity.y_true.eq(0), "propensity"]
            lower = max(case.quantile(0.05), control.quantile(0.05))
            upper = min(case.quantile(0.95), control.quantile(0.95))
            retained = propensity.loc[propensity.propensity.between(lower, upper), "donor_id"]
            selected = group.loc[
                group.donor_id.isin(retained)
                & group.adjustment.isin(("unadjusted", support_block))
                & group.method.isin(METHODS)
            ]
            for (method, adjustment), values in selected.groupby(["method", "adjustment"]):
                if values.y_true.nunique() < 2:
                    continue
                rows.append(
                    {
                        "repeat": repeat,
                        "support_block": support_block,
                        "method": method,
                        "adjustment": adjustment,
                        "n_overlap": len(values),
                        "overlap_lower": lower,
                        "overlap_upper": upper,
                        "roc_auc_overlap": roc_auc_score(values.y_true, values.prob_case),
                    }
                )
    return pd.DataFrame(rows)


def plot(metrics: pd.DataFrame) -> None:
    summary = metrics.groupby(["method", "adjustment"], as_index=False).agg(
        auc_mean=("roc_auc", "mean"), auc_low=("roc_auc", lambda x: x.quantile(0.025)), auc_high=("roc_auc", lambda x: x.quantile(0.975))
    )
    colors = {
        "Frozen Geneformer": "#2F78B7",
        "Mean-HVG pseudobulk": "#D76A39",
        "PCA pseudobulk": "#2FA38B",
    }
    order = ["unadjusted", *BLOCKS]
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), gridspec_kw={"width_ratios": [1.35, 1]})
    x = np.arange(len(order))
    for method in METHODS:
        sub = summary.loc[summary.method.eq(method)].set_index("adjustment").reindex(order)
        axes[0].plot(x, sub.auc_mean, marker="o", linewidth=2, color=colors[method], label=method)
        axes[0].fill_between(x, sub.auc_low, sub.auc_high, color=colors[method], alpha=0.12)
    axes[0].set_xticks(
        x,
        ["Unadjusted", "Cell yield", "QC core", "Acquisition", "Processing cohort", "Demographic", "Combined"],
        rotation=22,
        ha="right",
    )
    axes[0].set_ylabel("Repeated five-fold OOF ROC-AUC")
    axes[0].set_title("A  Final-representation covariate sensitivity", loc="left", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=8)

    cov = summary.loc[summary.method.str.startswith("Covariates only")].copy()
    cov["label"] = cov.adjustment.map(
        {
            "cell_yield": "Cell yield",
            "qc_core": "QC core",
            "acquisition_structure": "Acquisition",
            "processing_cohort": "Processing cohort",
            "demographic": "Demographic",
            "combined": "Combined",
        }
    )
    cov = cov.set_index("adjustment").reindex(BLOCKS)
    axes[1].barh(np.arange(len(cov)), cov.auc_mean, color=["#9DA8B3", "#6D7D8B", "#A68B6D", "#445D75"])
    axes[1].set_yticks(np.arange(len(cov)), cov.label)
    axes[1].invert_yaxis()
    axes[1].set_xlim(0.45, 1.0)
    axes[1].set_xlabel("Covariates-only ROC-AUC")
    axes[1].set_title("B  Measured covariates predict case status", loc="left", fontweight="bold")
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="y", color="#E1E8EE", linewidth=0.7)
        axis.set_axisbelow(True)
    fig.suptitle("Measured-covariate sensitivity in the final 261-donor CD4 analysis", fontweight="bold")
    fig.text(0.5, 0.005, "All covariate transforms and residualizers were fitted within each outer training fold. Attenuation is associational, not causal.", ha="center", fontsize=8, color="#52657A")
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    for suffix, kwargs in (("png", {"dpi": 320}), ("pdf", {})):
        fig.savefig(FIGURE.with_suffix(f".{suffix}"), bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    WORKSPACE.joinpath("figures").mkdir(parents=True, exist_ok=True)
    covariates = pd.read_csv(WORKSPACE / "results" / "01_inputs" / "gse174188_final_donor_covariates.tsv", sep="\t", dtype={"donor_id": str}).set_index("donor_id")
    labels = load_labels(DATASET)
    donors = labels.index.to_list()
    covariates = covariates.loc[donors].copy()
    covariates["log_cells_per_donor"] = np.log1p(covariates["cells_per_donor"])
    covariates["log_mean_umi_per_cell"] = np.log1p(covariates["mean_umi_per_cell"])
    covariates["log_mean_genes_per_cell"] = np.log1p(covariates["mean_genes_per_cell"])
    covariates["log_n_libraries"] = np.log1p(covariates["n_libraries"])
    covariates["log_n_samples"] = np.log1p(covariates["n_samples"])
    covariates["log_n_suspensions"] = np.log1p(covariates["n_suspensions"])
    if covariates[["age_years", "sex"]].isna().any().any():
        raise RuntimeError("Primary demographic covariates contain missing values")

    embedding = load_embedding(DATASET).loc[donors].to_numpy(dtype=np.float32)
    pseudobulk = load_pseudobulk(DATASET).loc[donors].to_numpy(dtype=np.float32)
    y = labels.eq("case").astype(int).to_numpy()
    results = Parallel(n_jobs=N_JOBS, backend="loky", verbose=10)(
        delayed(evaluate_repeat)(repeat, embedding, pseudobulk, y, donors, covariates)
        for repeat in range(1, N_REPEATS + 1)
    )
    metrics = pd.DataFrame([row for metric_rows, _, _ in results for row in metric_rows])
    predictions = pd.DataFrame([row for _, prediction_rows, _ in results for row in prediction_rows])
    metrics.to_csv(OUT / "repeat_level_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(OUT / "repeat_level_oof_predictions.tsv.gz", sep="\t", index=False, compression="gzip")
    summary = metrics.groupby(["method", "adjustment"], as_index=False).agg(
        auc_mean=("roc_auc", "mean"),
        auc_sd=("roc_auc", "std"),
        auc_q025=("roc_auc", lambda x: x.quantile(0.025)),
        auc_q500=("roc_auc", "median"),
        auc_q975=("roc_auc", lambda x: x.quantile(0.975)),
        pr_auc_mean=("pr_auc", "mean"),
        brier_mean=("brier", "mean"),
    )
    summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    paired = metrics.loc[metrics.method.isin(METHODS)].pivot(index=["repeat", "method"], columns="adjustment", values="roc_auc").reset_index()
    for block in BLOCKS:
        paired[f"delta_{block}_minus_unadjusted"] = paired[block] - paired["unadjusted"]
    paired.to_csv(OUT / "paired_auc_attenuation.tsv", sep="\t", index=False)
    overlap = overlap_metrics(predictions)
    overlap.to_csv(OUT / "common_support_metrics.tsv", sep="\t", index=False)
    plot(metrics)

    design_names = {block: sorted(set(name for _, _, names in results for name in names.get(block, []))) for block in BLOCKS}
    manifest = {
        "analysis": "20 repeated stratified five-fold final-representation covariate sensitivity",
        "n_repeats": N_REPEATS,
        "n_jobs": N_JOBS,
        "ridge_alpha": RIDGE_ALPHA,
        "blocks": design_names,
        "inputs": {
            "embedding": {"path": str(embedding_path(DATASET)), "sha256": sha256_file(embedding_path(DATASET))},
            "pseudobulk": {"path": str(pseudobulk_path(DATASET)), "sha256": sha256_file(pseudobulk_path(DATASET))},
            "labels": {"path": str(label_path(DATASET)), "sha256": sha256_file(label_path(DATASET))},
            "covariates": str(WORKSPACE / "results" / "01_inputs" / "gse174188_final_donor_covariates.tsv"),
        },
        "claim_boundary": "Sensitivity to measured covariates; no causal attribution or proof of biological independence.",
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
