#!/usr/bin/env python3
"""Re-evaluate cap500 pooling with shared source folds and fold-contained scaling."""

from __future__ import annotations

import argparse
import importlib.util
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SEED = 20_260_815
CS = np.logspace(-4, 4, 9)
warnings.filterwarnings("ignore", category=FutureWarning)


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("locked_pooling", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fit_strict(source: np.ndarray, y: np.ndarray, target: np.ndarray, seed: int):
    folds = max(2, min(5, int(min(y.sum(), (1 - y).sum()))))
    cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    pipeline = Pipeline(
        [
            ("scale", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    penalty="l2",
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=20_000,
                ),
            ),
        ]
    )
    search = GridSearchCV(
        pipeline,
        {"classifier__C": CS},
        scoring="roc_auc",
        cv=cv,
        refit=True,
        n_jobs=-1,
        return_train_score=False,
    )
    search.fit(source, y)
    return search.predict_proba(target)[:, 1], float(search.best_params_["classifier__C"]), folds


def fit_prefit_scaler(source: np.ndarray, y: np.ndarray, target: np.ndarray, seed: int):
    scaler = StandardScaler().fit(source)
    source_scaled = scaler.transform(source)
    target_scaled = scaler.transform(target)
    folds = max(2, min(5, int(min(y.sum(), (1 - y).sum()))))
    classifier = LogisticRegressionCV(
        Cs=CS,
        cv=StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed),
        scoring="roc_auc",
        penalty="l2",
        solver="liblinear",
        class_weight="balanced",
        max_iter=20_000,
        n_jobs=-1,
        refit=True,
    )
    classifier.fit(source_scaled, y)
    return classifier.predict_proba(target_scaled)[:, 1], float(classifier.C_[0]), folds


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-root", type=Path, required=True)
    parser.add_argument("--locked-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    locked = load_module(args.locked_script)

    embedding_root = args.archive_root / "outputs/cell_embeddings_cap500"
    reference_root = args.archive_root / "reference_inputs"
    datasets = locked.DATASETS
    labels = {dataset: locked.load_labels(reference_root, dataset) for dataset in datasets}
    metadata, cells, indices, aggregates = {}, {}, {}, {}
    for dataset in datasets:
        metadata[dataset], cells[dataset] = locked.load_cells(embedding_root / dataset)
        indices[dataset] = locked.donor_indices(metadata[dataset], labels[dataset].index)
        aggregates[dataset] = locked.coordinate_aggregates(
            cells[dataset], indices[dataset], labels[dataset].index
        )

    rows = []
    directions = ((datasets[0], datasets[1]), (datasets[1], datasets[0]))
    for direction_index, (source_id, target_id) in enumerate(directions):
        source_y = labels[source_id].eq("case").astype(int).to_numpy()
        target_y = labels[target_id].eq("case").astype(int).to_numpy()
        features = {
            method: (aggregates[source_id][method], aggregates[target_id][method])
            for method in ("mean", "median", "trimmed_mean_10", "mean_std")
        }
        features.update(
            locked.source_distribution_features(
                cells[source_id],
                cells[target_id],
                indices[source_id],
                indices[target_id],
                labels[source_id].index,
                labels[target_id].index,
                SEED + direction_index,
            )
        )
        for components in (8, 16):
            features[f"source_pca{components}_donor_mean"] = locked.donor_pca_features(
                aggregates[source_id]["mean"],
                aggregates[target_id]["mean"],
                components,
                SEED + direction_index,
            )

        shared_seed = SEED + direction_index
        for method_index, (method, (source_frame, target_frame)) in enumerate(features.items()):
            canonical_seed = SEED + direction_index * 100 + method_index
            protocols = (
                ("shared_folds_prefit_scaler", fit_prefit_scaler, shared_seed),
                ("canonical_folds_fold_contained_scaler", fit_strict, canonical_seed),
                ("shared_folds_fold_contained_scaler", fit_strict, shared_seed),
            )
            for protocol, fitter, protocol_seed in protocols:
                probability, selected_c, folds = fitter(
                    source_frame.to_numpy(dtype=float),
                    source_y,
                    target_frame.to_numpy(dtype=float),
                    protocol_seed,
                )
                rows.append(
                    {
                        "protocol": protocol,
                        "direction": f"{source_id}_to_{target_id}",
                        "source_dataset": source_id,
                        "target_dataset": target_id,
                        "method": method,
                        "n_source": len(source_y),
                        "n_target": len(target_y),
                        "n_features": source_frame.shape[1],
                        "roc_auc": roc_auc_score(target_y, probability),
                        "pr_auc": average_precision_score(target_y, probability),
                        "brier": brier_score_loss(target_y, probability),
                        "ece_10bin": locked.expected_calibration_error(target_y, probability),
                        "selected_c": selected_c,
                        "source_cv_seed": protocol_seed,
                        "source_internal_cv_folds": folds,
                        "fraction_extreme_probability": float(
                            ((probability < 0.01) | (probability > 0.99)).mean()
                        ),
                    }
                )

    strict = pd.DataFrame(rows)
    canonical = pd.read_csv(
        args.archive_root / "outputs/distributional_pooling_cap500/transfer_metrics.tsv", sep="\t"
    )
    comparison = strict.merge(
        canonical,
        on=["direction", "source_dataset", "target_dataset", "method", "n_source", "n_target", "n_features"],
        suffixes=("_strict", "_canonical"),
    )
    for metric in ("roc_auc", "pr_auc", "brier", "ece_10bin"):
        comparison[f"delta_{metric}_strict_minus_canonical"] = (
            comparison[f"{metric}_strict"] - comparison[f"{metric}_canonical"]
        )
    strict.to_csv(args.output / "strict_shared_fold_metrics.tsv", sep="\t", index=False)
    comparison.to_csv(args.output / "strict_vs_canonical.tsv", sep="\t", index=False)
    protocol_summary = (
        comparison.groupby("protocol")
        .agg(
            maximum_absolute_auc_change=(
                "delta_roc_auc_strict_minus_canonical",
                lambda values: float(values.abs().max()),
            ),
            maximum_absolute_brier_change=(
                "delta_brier_strict_minus_canonical",
                lambda values: float(values.abs().max()),
            ),
        )
        .reset_index()
    )
    protocol_summary.to_csv(args.output / "protocol_sensitivity_summary.tsv", sep="\t", index=False)
    summary = {
        "status": "complete",
        "scope": "Source-CV fold and scaler-containment sensitivity",
        "methods": int(strict["method"].nunique()),
        "directions": int(strict["direction"].nunique()),
        "protocol_summary": protocol_summary.to_dict(orient="records"),
        "mean_pooling_results": strict[strict["method"].eq("mean")].to_dict(orient="records"),
    }
    (args.output / "strict_shared_fold_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
