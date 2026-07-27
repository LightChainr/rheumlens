#!/usr/bin/env python3
"""Recompute bidirectional CD4 transfer with a fully source-only feature path.

This release-facing analysis intentionally does not reuse the historical target
scores. For each direction it reconstructs the retained donor pseudobulk and
frozen Geneformer feature matrices, then fits every learned operation on source
donors only. It exports gene-universe, HVG, source-scaling, PCA, classifier,
and file-hash records so the transfer table can be independently audited.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
from path_config import DONOR_LEVEL_ROOT, RESULTS_ROOT

OUT = RESULTS_ROOT / "strict_source_only_transfer"

DATASETS = ("SLE_GSE174188_CD4", "SLE_GSE285773_CD4")
GENEFORMER_METHOD = "geneformer_v2_316m_cell_sample1000_clspool_logistic_maxlen4096_seed001"
METHODS = (
    ("frozen_geneformer", "Frozen Geneformer", "Geneformer V2-316M donor embedding"),
    ("source_hvg_pseudobulk", "Source-HVG pseudobulk", "4,000 source-selected log1p-CPM features"),
    ("source_pca_pseudobulk", "Source-PCA pseudobulk", "30 source-fitted PCs from source-selected log1p-CPM features"),
)
N_HVG = 4_000
N_PCS = 30
BOOTSTRAPS = 5_000
SEED = 20_260_815
CS = np.logspace(-4, 4, 9)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: list[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode("utf-8")).hexdigest()


def pseudobulk_path(dataset: str) -> Path:
    return DONOR_LEVEL_ROOT / dataset / "donor_log1p_cpm.parquet"


def label_path(dataset: str) -> Path:
    return DONOR_LEVEL_ROOT / dataset / "donor_labels.tsv"


def embedding_path(dataset: str) -> Path:
    return DONOR_LEVEL_ROOT / dataset / "donor_embedding.parquet"


def load_labels(dataset: str) -> pd.Series:
    labels = pd.read_csv(label_path(dataset), sep="\t", dtype={"donor_id": str}).set_index("donor_id")["case_control"]
    labels.index = labels.index.astype(str)
    if labels.index.has_duplicates or set(labels.astype(str)) != {"case", "control"}:
        raise RuntimeError(f"{dataset}: invalid donor label table")
    return labels.astype(str).sort_index()


def canonicalize_pseudobulk(frame: pd.DataFrame, dataset: str) -> pd.DataFrame:
    """Use retained identifiers exactly; reject ambiguity instead of silently merging genes."""

    result = frame.copy()
    result.index = result.index.astype(str)
    columns = pd.Index(result.columns.astype(str).str.strip())
    if result.index.has_duplicates:
        raise RuntimeError(f"{dataset}: duplicate donor identifiers")
    if (columns == "").any() or columns.has_duplicates:
        raise RuntimeError(f"{dataset}: empty or duplicate retained feature identifiers")
    if not np.isfinite(result.to_numpy(dtype=np.float32)).all():
        raise RuntimeError(f"{dataset}: non-finite retained log1p-CPM values")
    result.columns = columns
    return result.sort_index()


def load_pseudobulk(dataset: str, labels: pd.Series) -> pd.DataFrame:
    features = canonicalize_pseudobulk(pd.read_parquet(pseudobulk_path(dataset)), dataset)
    if set(features.index) != set(labels.index):
        raise RuntimeError(f"{dataset}: donor mismatch between pseudobulk and labels")
    return features.loc[labels.index]


def load_embedding(dataset: str, labels: pd.Series) -> pd.DataFrame:
    embedding = pd.read_parquet(embedding_path(dataset))
    embedding.index = embedding.index.astype(str)
    if embedding.index.has_duplicates or embedding.shape[1] != 1152:
        raise RuntimeError(f"{dataset}: unexpected Geneformer donor embedding shape {embedding.shape}")
    if set(embedding.index) != set(labels.index):
        raise RuntimeError(f"{dataset}: donor mismatch between Geneformer embedding and labels")
    array = embedding.to_numpy(dtype=np.float32)
    if not np.isfinite(array).all():
        raise RuntimeError(f"{dataset}: non-finite Geneformer embedding values")
    return embedding.loc[labels.index]


def bootstrap_auc_ci(y: np.ndarray, probability: np.ndarray, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    positive = np.flatnonzero(y == 1)
    negative = np.flatnonzero(y == 0)
    values = np.empty(BOOTSTRAPS, dtype=np.float64)
    for draw in range(BOOTSTRAPS):
        sampled = np.concatenate(
            (rng.choice(positive, size=len(positive), replace=True), rng.choice(negative, size=len(negative), replace=True))
        )
        values[draw] = roc_auc_score(y[sampled], probability[sampled])
    low, high = np.quantile(values, (0.025, 0.975))
    return float(low), float(high)


def expected_calibration_error(y: np.ndarray, probability: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    value = 0.0
    for lower, upper in zip(edges[:-1], edges[1:]):
        selected = (probability >= lower) & (probability < upper if upper < 1 else probability <= upper)
        if selected.any():
            value += float(selected.mean()) * abs(float(y[selected].mean()) - float(probability[selected].mean()))
    return float(value)


def placements(y: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    positive = score[y == 1]
    negative = score[y == 0]
    comparison = positive[:, None] - negative[None, :]
    values = (comparison > 0).astype(float) + 0.5 * (comparison == 0)
    return values.mean(axis=1), values.mean(axis=0), float(values.mean())


def paired_delong(y: np.ndarray, baseline: np.ndarray, geneformer: np.ndarray) -> dict[str, float]:
    base_v10, base_v01, base_auc = placements(y, baseline)
    gf_v10, gf_v01, gf_auc = placements(y, geneformer)
    covariance = np.cov(np.vstack((base_v10, gf_v10)), ddof=1) / len(base_v10)
    covariance += np.cov(np.vstack((base_v01, gf_v01)), ddof=1) / len(base_v01)
    variance = float(covariance[0, 0] + covariance[1, 1] - 2.0 * covariance[0, 1])
    standard_error = math.sqrt(max(variance, 0.0))
    delta = gf_auc - base_auc
    z_value = delta / standard_error if standard_error else 0.0
    p_value = 1.0 if standard_error == 0.0 and delta == 0.0 else math.erfc(abs(z_value) / math.sqrt(2.0))
    return {
        "baseline_auc": base_auc,
        "geneformer_auc": gf_auc,
        "delta_geneformer_minus_baseline": delta,
        "delong_standard_error": standard_error,
        "delong_z": z_value,
        "delong_p": p_value,
    }


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values), dtype=float)
    running = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = len(p_values) - reverse_rank + 1
        running = min(running, float(p_values[index]) * len(p_values) / rank)
        adjusted[index] = running
    return adjusted


def fit_classifier(source_x: np.ndarray, source_y: np.ndarray, target_x: np.ndarray, direction_seed: int) -> tuple[np.ndarray, float, int]:
    smallest_class = int(min(source_y.sum(), (1 - source_y).sum()))
    folds = max(2, min(5, smallest_class))
    inner_cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=direction_seed)
    classifier = LogisticRegressionCV(
        Cs=CS,
        cv=inner_cv,
        scoring="roc_auc",
        penalty="l2",
        solver="liblinear",
        class_weight="balanced",
        max_iter=20_000,
        n_jobs=1,
        refit=True,
    )
    classifier.fit(source_x, source_y)
    return classifier.predict_proba(target_x)[:, 1], float(classifier.C_[0]), folds


def make_pseudobulk_features(
    source: pd.DataFrame, target: pd.DataFrame, direction: str
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray]], dict[str, object], pd.DataFrame, pd.DataFrame]:
    common = sorted(set(source.columns).intersection(target.columns))
    if not common:
        raise RuntimeError(f"{direction}: no common retained features")
    source_common = source.loc[:, common]
    target_common = target.loc[:, common]
    variance = source_common.var(axis=0, ddof=0).sort_values(ascending=False, kind="mergesort")
    selected = variance.head(min(N_HVG, len(variance))).index.astype(str).tolist()
    source_hvg = source_common.loc[:, selected].to_numpy(dtype=np.float64)
    target_hvg = target_common.loc[:, selected].to_numpy(dtype=np.float64)
    scaler = StandardScaler().fit(source_hvg)
    source_scaled = scaler.transform(source_hvg)
    target_scaled = scaler.transform(target_hvg)
    if not np.isfinite(source_scaled).all() or not np.isfinite(target_scaled).all():
        raise RuntimeError(f"{direction}: non-finite values after source-only pseudobulk scaling")
    n_components = min(N_PCS, source_scaled.shape[0] - 2, source_scaled.shape[1])
    pca = PCA(n_components=n_components, svd_solver="full", random_state=SEED)
    source_pc = pca.fit_transform(source_scaled)
    target_pc = pca.transform(target_scaled)
    selected_frame = pd.DataFrame(
        {
            "feature_id": selected,
            "source_population_variance_log1p_cpm": variance.loc[selected].to_numpy(dtype=float),
            "source_mean_log1p_cpm": scaler.mean_,
            "source_scale_log1p_cpm": scaler.scale_,
        }
    )
    pca_frame = pd.DataFrame(
        {
            "component": np.arange(1, n_components + 1),
            "source_explained_variance": pca.explained_variance_,
            "source_explained_variance_ratio": pca.explained_variance_ratio_,
        }
    )
    metadata: dict[str, object] = {
        "identifier_harmonization": "stripped retained feature_name/10x var_names identifiers; exact case-sensitive intersection; no external identifier map",
        "duplicate_feature_policy": "fail on any duplicate retained identifier after whitespace stripping; no duplicate collapse was needed",
        "missing_feature_policy": "features absent from either cohort were excluded, never imputed as zero",
        "n_source_features_before_intersection": int(source.shape[1]),
        "n_target_features_before_intersection": int(target.shape[1]),
        "n_exact_common_features": len(common),
        "exact_common_feature_sha256": sha256_lines(common),
        "n_source_selected_hvgs": len(selected),
        "selected_hvg_sha256": sha256_lines(selected),
        "n_source_fitted_pcs": int(n_components),
        "source_hvg_rule": "top population variance across all source donors after donor-specific log1p CPM; stable lexical order breaks ties",
        "source_scaling": "StandardScaler mean and scale fit on source donor features only; target transformed with these source parameters",
        "pca": "PCA fit on source-scaled selected features only; target projected with source PCA loadings",
    }
    return {
        "source_hvg_pseudobulk": (source_scaled, target_scaled),
        "source_pca_pseudobulk": (source_pc, target_pc),
    }, metadata, selected_frame, pca_frame


def source_scale_embedding(source: pd.DataFrame, target: pd.DataFrame, direction: str) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    scaler = StandardScaler().fit(source.to_numpy(dtype=np.float64))
    source_scaled = scaler.transform(source.to_numpy(dtype=np.float64))
    target_scaled = scaler.transform(target.to_numpy(dtype=np.float64))
    if not np.isfinite(source_scaled).all() or not np.isfinite(target_scaled).all():
        raise RuntimeError(f"{direction}: non-finite values after source-only Geneformer scaling")
    return source_scaled, target_scaled, {
        "source_embedding_dimensions": int(source.shape[1]),
        "source_scaling": "StandardScaler fit on source frozen donor embeddings only; target transformed with source parameters",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    labels = {dataset: load_labels(dataset) for dataset in DATASETS}
    pseudobulk = {dataset: load_pseudobulk(dataset, labels[dataset]) for dataset in DATASETS}
    embedding = {dataset: load_embedding(dataset, labels[dataset]) for dataset in DATASETS}
    metric_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    input_paths = [pseudobulk_path(dataset) for dataset in DATASETS] + [label_path(dataset) for dataset in DATASETS] + [embedding_path(dataset) for dataset in DATASETS]

    for direction_index, (source_dataset, target_dataset) in enumerate(((DATASETS[0], DATASETS[1]), (DATASETS[1], DATASETS[0]))):
        direction = f"{source_dataset}_to_{target_dataset}"
        source_y = labels[source_dataset].eq("case").astype(int).to_numpy()
        target_y = labels[target_dataset].eq("case").astype(int).to_numpy()
        pseudobulk_features, pseudobulk_metadata, selected_frame, pca_frame = make_pseudobulk_features(
            pseudobulk[source_dataset], pseudobulk[target_dataset], direction
        )
        selected_frame.to_csv(OUT / f"{direction}_source_hvg_and_scaling.tsv", sep="\t", index=False)
        pca_frame.to_csv(OUT / f"{direction}_source_pca_variance.tsv", sep="\t", index=False)
        gf_source, gf_target, gf_metadata = source_scale_embedding(embedding[source_dataset], embedding[target_dataset], direction)
        feature_sets = {"frozen_geneformer": (gf_source, gf_target), **pseudobulk_features}
        predictions: dict[str, np.ndarray] = {}
        feature_rows.append({"direction": direction, **pseudobulk_metadata, **gf_metadata})

        for method_index, (method_id, label, description) in enumerate(METHODS):
            source_x, target_x = feature_sets[method_id]
            probability, selected_c, inner_folds = fit_classifier(
                source_x, source_y, target_x, direction_seed=SEED + direction_index * 100 + method_index
            )
            predictions[method_id] = probability
            ci_low, ci_high = bootstrap_auc_ci(target_y, probability, SEED + direction_index * 1000 + method_index)
            metric_rows.append(
                {
                    "direction": direction,
                    "source_dataset": source_dataset,
                    "target_dataset": target_dataset,
                    "method_id": method_id,
                    "method_label": label,
                    "representation": description,
                    "n_source": len(source_y),
                    "n_source_case": int(source_y.sum()),
                    "n_source_control": int((1 - source_y).sum()),
                    "n_target": len(target_y),
                    "n_target_case": int(target_y.sum()),
                    "n_target_control": int((1 - target_y).sum()),
                    "n_features_for_classifier": int(source_x.shape[1]),
                    "roc_auc": float(roc_auc_score(target_y, probability)),
                    "roc_auc_ci_low": ci_low,
                    "roc_auc_ci_high": ci_high,
                    "pr_auc": float(average_precision_score(target_y, probability)),
                    "brier": float(brier_score_loss(target_y, probability)),
                    "ece_10bin_supplementary": expected_calibration_error(target_y, probability),
                    "selected_c_source_internal_cv": selected_c,
                    "source_internal_cv_folds": inner_folds,
                    "source_internal_cv": "stratified, shuffled, ROC-AUC selection; fixed direction-method seed",
                    "target_adaptation": "none; no target labels, mean, standard deviation, PCA fit, feature selection, or classifier tuning used",
                }
            )
            prediction_rows.extend(
                {
                    "direction": direction,
                    "source_dataset": source_dataset,
                    "target_dataset": target_dataset,
                    "method_id": method_id,
                    "donor_id": donor_id,
                    "y_true": int(y_value),
                    "prob_case": float(probability_value),
                }
                for donor_id, y_value, probability_value in zip(labels[target_dataset].index, target_y, probability)
            )

        for baseline in ("source_hvg_pseudobulk", "source_pca_pseudobulk"):
            paired = paired_delong(target_y, predictions[baseline], predictions["frozen_geneformer"])
            comparison_rows.append(
                {
                    "direction": direction,
                    "target_dataset": target_dataset,
                    "baseline_method_id": baseline,
                    "comparison": f"Frozen Geneformer vs {baseline} on identical target donors",
                    "n_target": len(target_y),
                    **paired,
                }
            )

    comparisons = pd.DataFrame(comparison_rows)
    comparisons["delong_p_bh_across_four_tests"] = bh_adjust(comparisons["delong_p"].to_numpy(dtype=float))
    metrics = pd.DataFrame(metric_rows)
    predictions = pd.DataFrame(prediction_rows)
    features = pd.DataFrame(feature_rows)
    metrics.to_csv(OUT / "strict_source_only_transfer_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(OUT / "strict_source_only_transfer_predictions.tsv", sep="\t", index=False)
    comparisons.to_csv(OUT / "strict_source_only_transfer_paired_delong.tsv", sep="\t", index=False)
    features.to_csv(OUT / "strict_source_only_feature_harmonization.tsv", sep="\t", index=False)

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "analysis": "strict bidirectional source-only CD4 transfer recomputation",
        "datasets": list(DATASETS),
        "donor_pseudobulk": "donor-specific log1p(1e6 * summed raw counts / donor library size), supplied by retained matrix provenance",
        "identifier_harmonization": "strip whitespace, require unique retained identifiers, then exact case-sensitive intersection; no external annotation mapping",
        "missing_or_duplicate_features": "missing identifiers excluded from both cohorts; retained duplicates fail fast and are never summed or imputed",
        "hvg_selection": {"n_hvg": N_HVG, "source_only": True, "variance_ddof": 0},
        "scaling": "source-only StandardScaler for all classifier inputs",
        "pca": {"n_pcs_max": N_PCS, "source_only": True, "svd_solver": "full"},
        "classifier": {
            "class": "LogisticRegressionCV",
            "solver": "liblinear",
            "penalty": "l2",
            "class_weight": "balanced",
            "Cs": CS.tolist(),
            "selection_metric": "ROC-AUC",
            "internal_cv": "min(5, source minority class), stratified shuffled, source-only",
        },
        "inference": {
            "auc_ci": f"{BOOTSTRAPS} donor-stratified bootstrap resamples; percentile 2.5% and 97.5%",
            "paired_test": "paired DeLong on identical target donors; Benjamini-Hochberg adjustment across four prespecified comparisons",
        },
        "environment": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__, "scikit_learn": sklearn.__version__},
        "input_files": {
            str(path): {"bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in input_paths
        },
        "output_files": sorted(path.name for path in OUT.iterdir() if path.is_file()),
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    summary = ["# Strict Source-only Transfer Recalculation", "", "## Design", ""]
    summary.append("All learned preprocessing, feature selection, scaling, PCA, and classifier regularization selection was fitted on source donors only. Target cohorts contributed only already-normalized donor vectors and unlabeled feature rows.")
    summary.extend(["", "## Results", ""])
    for row in metrics.itertuples(index=False):
        summary.append(f"- {row.direction} / {row.method_label}: AUC {row.roc_auc:.3f} ({row.roc_auc_ci_low:.3f}-{row.roc_auc_ci_high:.3f}); Brier {row.brier:.3f}.")
    summary.extend(["", "## Paired comparisons", ""])
    for row in comparisons.itertuples(index=False):
        summary.append(f"- {row.direction} / {row.baseline_method_id}: Geneformer minus baseline AUC {row.delta_geneformer_minus_baseline:+.3f}; paired DeLong BH p={row.delong_p_bh_across_four_tests:.4g}.")
    summary.extend(["", "## Boundary", "", "This recomputation does not infer single-patient prospective deployment performance beyond the evaluated two-cohort transfer setting."])
    (OUT / "SUMMARY.md").write_text("\n".join(summary) + "\n")
    print(metrics.to_string(index=False))
    print(comparisons.to_string(index=False))


if __name__ == "__main__":
    main()
