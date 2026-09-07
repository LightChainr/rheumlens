#!/usr/bin/env python3
"""Locked same-cell-set CD4 composition audit for GSE174188."""

from __future__ import annotations

import os

import h5py
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from common import (
    BASELINE,
    H5AD_174188,
    SEEDS,
    WORKSPACE,
    effective_splits,
    predict_scaled,
    residualize,
    sha256,
    summarise_metrics,
    write_json,
)


OUT = WORKSPACE / "results" / "composition"
METADATA_PATH = (
    BASELINE
    / "inputs"
    / "design_metadata"
    / "gse174188_final_donor_covariates.tsv"
)
N_JOBS = min(int(os.environ.get("N_JOBS", "12")), os.cpu_count() or 1)
RIDGE_ALPHA = 1.0
ZERO_REPLACEMENT_COUNT = 0.5
PARTS = ("T4_naive", "T4_em", "T4_reg", "other_or_unassigned")


def read_categorical(group: h5py.Group, name: str) -> tuple[np.ndarray, np.ndarray]:
    item = group[name]
    categories = np.array(
        [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in item["categories"][:]
        ],
        dtype=object,
    )
    return categories, item["codes"][:]


def aggregate_features() -> pd.DataFrame:
    with h5py.File(H5AD_174188, "r") as handle:
        obs = handle["obs"]
        donor_categories, donor_codes = read_categorical(obs, "donor_id")
        cell_categories, cell_codes = read_categorical(obs, "cell_type")
        state_categories, state_codes = read_categorical(obs, "ct_cov")

    cd4_code = int(
        np.flatnonzero(
            cell_categories == "CD4-positive, alpha-beta T cell"
        )[0]
    )
    selected = (cell_codes == cd4_code) & (donor_codes >= 0)
    donors = donor_codes[selected]
    states = state_codes[selected]

    known_state_codes = {
        name: int(np.flatnonzero(state_categories == name)[0])
        for name in PARTS[:3]
    }
    part_index = np.full(len(states), 3, dtype=np.int8)
    for index, name in enumerate(PARTS[:3]):
        part_index[states == known_state_codes[name]] = index

    counts = np.zeros((len(donor_categories), len(PARTS)), dtype=np.int64)
    np.add.at(counts, (donors, part_index), 1)
    totals = counts.sum(axis=1)
    keep = totals > 0

    raw = counts[keep] / totals[keep, None]
    smoothed = counts[keep].astype(float) + ZERO_REPLACEMENT_COUNT
    closed = smoothed / smoothed.sum(axis=1, keepdims=True)
    clr = np.log(closed) - np.log(closed).mean(axis=1, keepdims=True)

    frame = pd.DataFrame(
        {
            "donor_id": donor_categories[keep].astype(str),
            "n_cd4_cells": totals[keep],
            **{
                f"count_{part}": counts[keep, index]
                for index, part in enumerate(PARTS)
            },
            **{
                f"prop_{part}": raw[:, index]
                for index, part in enumerate(PARTS)
            },
            **{
                f"clr_{part}": clr[:, index]
                for index, part in enumerate(PARTS)
            },
            "log1p_n_cd4_cells": np.log1p(totals[keep]),
        }
    )
    return frame.sort_values("donor_id").reset_index(drop=True)


def variants(data: pd.DataFrame) -> dict[str, np.ndarray]:
    raw_columns = [f"prop_{part}" for part in PARTS[:-1]]
    clr_columns = [f"clr_{part}" for part in PARTS[:-1]]
    return {
        "raw_proportions": data[raw_columns].to_numpy(dtype=float),
        "raw_plus_log_cells": data[
            [*raw_columns, "log1p_n_cd4_cells"]
        ].to_numpy(dtype=float),
        "clr_composition": data[clr_columns].to_numpy(dtype=float),
        "clr_plus_log_cells": data[
            [*clr_columns, "log1p_n_cd4_cells"]
        ].to_numpy(dtype=float),
    }


def metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(y, prediction)),
        "pr_auc": float(average_precision_score(y, prediction)),
        "brier": float(brier_score_loss(y, prediction)),
    }


def evaluate_repeat(
    name: str,
    x: np.ndarray,
    y: np.ndarray,
    wave: np.ndarray,
    repeat: int,
) -> tuple[list[dict], list[dict]]:
    seed = SEEDS[repeat]
    splitter = StratifiedKFold(
        effective_splits(y),
        shuffle=True,
        random_state=seed,
    )
    unadjusted = np.full(len(y), np.nan)
    residualised = np.full(len(y), np.nan)
    folds = np.zeros(len(y), dtype=int)
    for fold, (train, test) in enumerate(splitter.split(x, y), start=1):
        unadjusted[test], _ = predict_scaled(
            x[train], y[train], x[test], seed
        )
        train_residual, test_residual = residualize(
            x[train],
            x[test],
            wave[train],
            wave[test],
            alpha=RIDGE_ALPHA,
        )
        residualised[test], _ = predict_scaled(
            train_residual,
            y[train],
            test_residual,
            seed,
        )
        folds[test] = fold

    metric_rows = []
    prediction_rows = []
    for adjustment, values in (
        ("unadjusted", unadjusted),
        ("residual_processing_wave", residualised),
    ):
        metric_rows.append(
            {
                "variant": name,
                "adjustment": adjustment,
                "repeat": repeat + 1,
                "seed": seed,
                **metrics(y, values),
            }
        )
        prediction_rows.extend(
            {
                "variant": name,
                "adjustment": adjustment,
                "repeat": repeat + 1,
                "seed": seed,
                "row": int(index),
                "fold": int(fold),
                "y_true": int(label),
                "prob_case": float(probability),
            }
            for index, (fold, label, probability) in enumerate(
                zip(folds, y, values)
            )
        )
    return metric_rows, prediction_rows


def evaluate_subset(
    x: np.ndarray,
    y: np.ndarray,
    indices: np.ndarray,
    seed: int,
) -> float:
    labels = y[indices]
    splitter = StratifiedKFold(
        effective_splits(labels),
        shuffle=True,
        random_state=seed,
    )
    prediction = np.full(len(indices), np.nan)
    for train_local, test_local in splitter.split(x[indices], labels):
        prediction[test_local], _ = predict_scaled(
            x[indices][train_local],
            labels[train_local],
            x[indices][test_local],
            seed,
        )
    return float(roc_auc_score(labels, prediction))


def restriction_rows(
    name: str,
    x: np.ndarray,
    y: np.ndarray,
    observed: np.ndarray,
    stratum: str,
) -> list[dict]:
    rows = []
    for repeat, seed in enumerate(SEEDS, start=1):
        rows.append(
            {
                "variant": name,
                "stratum": stratum,
                "control_type": "observed",
                "draw": 0,
                "repeat": repeat,
                "seed": seed,
                "n_donors": len(observed),
                "n_cases": int(y[observed].sum()),
                "roc_auc": evaluate_subset(x, y, observed, seed),
            }
        )

    rng = np.random.default_rng(20260815)
    n_cases = int(y[observed].sum())
    n_controls = int(len(observed) - n_cases)
    cases = np.flatnonzero(y == 1)
    controls = np.flatnonzero(y == 0)
    for draw, seed in enumerate(SEEDS, start=1):
        matched = np.concatenate(
            [
                rng.choice(cases, n_cases, replace=False),
                rng.choice(controls, n_controls, replace=False),
            ]
        )
        rows.append(
            {
                "variant": name,
                "stratum": stratum,
                "control_type": "size_label_matched",
                "draw": draw,
                "repeat": draw,
                "seed": seed,
                "n_donors": len(matched),
                "n_cases": n_cases,
                "roc_auc": evaluate_subset(x, y, matched, seed),
            }
        )
    return rows


def wave_predictability(
    name: str,
    x: np.ndarray,
    dominant_wave: np.ndarray,
) -> list[dict]:
    rows = []
    for wave_value in sorted(np.unique(dominant_wave)):
        target = (dominant_wave == wave_value).astype(int)
        if np.bincount(target).min() < 2:
            continue
        for repeat, seed in enumerate(SEEDS[:5], start=1):
            splitter = StratifiedKFold(
                effective_splits(target),
                shuffle=True,
                random_state=seed,
            )
            prediction = np.full(len(target), np.nan)
            for train, test in splitter.split(x, target):
                prediction[test], _ = predict_scaled(
                    x[train],
                    target[train],
                    x[test],
                    seed,
                )
            rows.append(
                {
                    "variant": name,
                    "dominant_wave": str(wave_value),
                    "repeat": repeat,
                    "seed": seed,
                    "n_positive": int(target.sum()),
                    "roc_auc": float(roc_auc_score(target, prediction)),
                }
            )
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    features = aggregate_features()
    metadata = pd.read_csv(METADATA_PATH, sep="\t", dtype={"donor_id": str})
    data = metadata.merge(
        features,
        on="donor_id",
        how="inner",
        validate="one_to_one",
    ).sort_values("donor_id").reset_index(drop=True)
    if len(data) != len(metadata):
        raise RuntimeError(
            f"Composition donor mapping incomplete: {len(data)} of {len(metadata)}"
        )

    wave_columns = [
        column
        for column in data.columns
        if column.startswith("processing_fraction_")
    ]
    wave = data[wave_columns].to_numpy(dtype=float)
    dominant_wave = np.argmax(wave, axis=1) + 1
    pure = np.max(wave, axis=1) >= 0.99
    y = data["y_true"].to_numpy(dtype=int)
    feature_sets = variants(data)

    completed = Parallel(n_jobs=N_JOBS, prefer="processes", verbose=5)(
        delayed(evaluate_repeat)(name, x, y, wave, repeat)
        for name, x in feature_sets.items()
        for repeat in range(len(SEEDS))
    )
    metric_rows = []
    prediction_rows = []
    for metrics_part, predictions_part in completed:
        metric_rows.extend(metrics_part)
        prediction_rows.extend(predictions_part)
    metric_frame = pd.DataFrame(metric_rows)
    prediction_frame = pd.DataFrame(prediction_rows)
    metric_frame.to_csv(OUT / "composition_repeat_metrics.tsv", sep="\t", index=False)
    prediction_frame.to_csv(
        OUT / "composition_oof_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    summarise_metrics(
        metric_frame,
        ["variant", "adjustment"],
    ).to_csv(OUT / "composition_summary.tsv", sep="\t", index=False)

    restriction = []
    wave_rows = []
    for name, x in feature_sets.items():
        dominant_four = np.flatnonzero(dominant_wave == 4)
        pure_four = np.flatnonzero((dominant_wave == 4) & pure)
        restriction.extend(
            restriction_rows(
                name,
                x,
                y,
                dominant_four,
                "dominant_wave_4",
            )
        )
        restriction.extend(
            restriction_rows(
                name,
                x,
                y,
                pure_four,
                "pure_wave_4",
            )
        )
        wave_rows.extend(wave_predictability(name, x, dominant_wave))
    pd.DataFrame(restriction).to_csv(
        OUT / "composition_restriction_metrics.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(wave_rows).to_csv(
        OUT / "composition_wave_predictability.tsv",
        sep="\t",
        index=False,
    )
    data.to_csv(
        OUT / "gse174188_cd4_composition_features.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    write_json(
        OUT / "composition_manifest.json",
        {
            "h5ad_path": str(H5AD_174188),
            "h5ad_sha256": sha256(H5AD_174188),
            "metadata_path": str(METADATA_PATH),
            "metadata_sha256": sha256(METADATA_PATH),
            "cell_filter": "cell_type == CD4-positive, alpha-beta T cell",
            "parts": PARTS,
            "zero_replacement_count": ZERO_REPLACEMENT_COUNT,
            "clr_definition": (
                "add 0.5 to each donor-part count, close to proportions, "
                "take log, subtract donor mean log; omit final coordinate"
            ),
            "variants": {
                name: int(x.shape[1]) for name, x in feature_sets.items()
            },
            "seeds": SEEDS,
            "ridge_alpha": RIDGE_ALPHA,
            "classifier": (
                "balanced liblinear logistic; training-fold scaling; "
                "inner 5-fold C selection from 1e-4 through 1e4"
            ),
        },
    )
    print(
        summarise_metrics(metric_frame, ["variant", "adjustment"]).to_string(
            index=False
        )
    )
    print(f"Wrote locked CD4 composition audit to {OUT}")


if __name__ == "__main__":
    main()
