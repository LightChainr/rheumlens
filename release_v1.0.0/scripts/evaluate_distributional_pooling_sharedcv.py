#!/usr/bin/env python3
"""Paper-ready cap500 pooling sensitivity with shared source CV folds."""

from __future__ import annotations

import argparse
import importlib.util
import json
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


SEED = 20_260_815
CS = np.logspace(-4, 4, 9)
warnings.filterwarnings("ignore", category=FutureWarning)


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bh_adjust(values: pd.Series) -> np.ndarray:
    p_values = values.to_numpy(dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    output = np.empty_like(ranked)
    output[order] = np.minimum(ranked, 1.0)
    return output


def fit_classifier(source: np.ndarray, y: np.ndarray, target: np.ndarray, seed: int):
    folds = max(2, min(5, int(min(y.sum(), (1 - y).sum()))))
    search = GridSearchCV(
        Pipeline(
            [
                ("scale", StandardScaler()),
                (
                    "classifier",
                    LogisticRegression(
                        solver="liblinear",
                        class_weight="balanced",
                        max_iter=20_000,
                    ),
                ),
            ]
        ),
        {"classifier__C": CS},
        scoring="roc_auc",
        cv=StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed),
        refit=True,
        n_jobs=-1,
    )
    search.fit(source, y)
    return search.predict_proba(target)[:, 1], float(search.best_params_["classifier__C"]), folds


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-root", type=Path, required=True)
    parser.add_argument("--locked-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    locked = load_module(args.locked_script, "locked_pooling")
    started = time.time()

    datasets = locked.DATASETS
    embedding_root = args.archive_root / "outputs/cell_embeddings_cap500"
    reference_root = args.archive_root / "reference_inputs"
    labels = {dataset: locked.load_labels(reference_root, dataset) for dataset in datasets}
    metadata, cells, indices, aggregates = {}, {}, {}, {}
    for dataset in datasets:
        metadata[dataset], cells[dataset] = locked.load_cells(embedding_root / dataset)
        indices[dataset] = locked.donor_indices(metadata[dataset], labels[dataset].index)
        aggregates[dataset] = locked.coordinate_aggregates(
            cells[dataset], indices[dataset], labels[dataset].index
        )
    validation = locked.historical_validation(reference_root, labels, aggregates)
    validation.to_csv(args.output / "historical_mean_reproduction.tsv", sep="\t", index=False)

    metric_rows = []
    prediction_rows = []
    comparison_rows = []
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
        probabilities = {}
        for method, (source_frame, target_frame) in features.items():
            probability, selected_c, folds = fit_classifier(
                source_frame.to_numpy(dtype=float),
                source_y,
                target_frame.to_numpy(dtype=float),
                shared_seed,
            )
            probabilities[method] = probability
            metric_rows.append(
                {
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
                    "source_cv_seed": shared_seed,
                    "source_internal_cv_folds": folds,
                }
            )
            prediction_rows.extend(
                {
                    "direction": f"{source_id}_to_{target_id}",
                    "source_dataset": source_id,
                    "target_dataset": target_id,
                    "method": method,
                    "donor_id": donor,
                    "y_true": int(y),
                    "prob_case": float(p),
                }
                for donor, y, p in zip(labels[target_id].index, target_y, probability)
            )

        baseline = probabilities["mean"]
        for method_index, method in enumerate(features):
            if method == "mean":
                continue
            candidate = probabilities[method]
            low, median, high = locked.paired_bootstrap(
                target_y, baseline, candidate, SEED + direction_index * 10_000 + method_index
            )
            comparison_rows.append(
                {
                    "direction": f"{source_id}_to_{target_id}",
                    "target_dataset": target_id,
                    "method": method,
                    "baseline_auc": roc_auc_score(target_y, baseline),
                    "candidate_auc": roc_auc_score(target_y, candidate),
                    "delta_auc": roc_auc_score(target_y, candidate) - roc_auc_score(target_y, baseline),
                    "delta_bootstrap_q025": low,
                    "delta_bootstrap_q500": median,
                    "delta_bootstrap_q975": high,
                    **locked.paired_delong(target_y, baseline, candidate),
                }
            )

    metrics = pd.DataFrame(metric_rows)
    predictions = pd.DataFrame(prediction_rows)
    comparisons = pd.DataFrame(comparison_rows)
    comparisons["delong_bh_q_14_tests"] = bh_adjust(comparisons["delong_p"])
    metrics.to_csv(args.output / "transfer_metrics.tsv", sep="\t", index=False)
    predictions.to_parquet(args.output / "target_predictions.parquet", index=False)
    comparisons.to_csv(args.output / "paired_comparisons.tsv", sep="\t", index=False)

    gates = []
    for method in comparisons["method"].unique():
        method_comparison = comparisons[comparisons["method"].eq(method)].set_index("target_dataset")
        method_metrics = metrics[metrics["method"].eq(method)].set_index("target_dataset")
        baseline_metrics = metrics[metrics["method"].eq("mean")].set_index("target_dataset")
        large = method_comparison.loc["SLE_GSE174188_CD4"]
        small = method_comparison.loc["SLE_GSE285773_CD4"]
        calibration_ok = bool(
            ((method_metrics["brier"] - baseline_metrics["brier"]) <= 0.02).all()
            and ((method_metrics["ece_10bin"] - baseline_metrics["ece_10bin"]) <= 0.05).all()
        )
        first_gate = bool(large.delta_auc >= 0.010 and small.delta_auc >= 0 and calibration_ok)
        promotion = bool(first_gate and large.delta_bootstrap_q025 > 0 and small.delta_auc >= -0.020)
        gates.append(
            {
                "method": method,
                "delta_auc_target_261": large.delta_auc,
                "delta_ci_low_target_261": large.delta_bootstrap_q025,
                "delta_auc_target_26": small.delta_auc,
                "calibration_gate_passed": calibration_ok,
                "continue_to_1000_cells": first_gate,
                "promotion_gate_passed_at_500": promotion,
            }
        )
    gate_frame = pd.DataFrame(gates)
    gate_frame.to_csv(args.output / "gate_decisions.tsv", sep="\t", index=False)
    manifest = {
        "status": "complete",
        "analysis": "Cap500 Geneformer pooling with shared source CV folds and fold-contained scaling",
        "bootstrap_draws": locked.BOOTSTRAPS,
        "shared_source_cv_folds_across_methods": True,
        "classifier_scaling_fold_contained": True,
        "historical_mean_reproduction_passed": bool(validation["passed"].all()),
        "continue_to_1000_cells": gate_frame.loc[
            gate_frame["continue_to_1000_cells"], "method"
        ].tolist(),
        "promotion_gate_passed_at_500": gate_frame.loc[
            gate_frame["promotion_gate_passed_at_500"], "method"
        ].tolist(),
        "elapsed_seconds": time.time() - started,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
