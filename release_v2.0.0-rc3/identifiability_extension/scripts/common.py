#!/usr/bin/env python3
"""Shared deterministic utilities for the identifiability upgrade."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


WORKSPACE = Path(__file__).resolve().parents[1]
PROJECT_ROOT = WORKSPACE.parent
BASELINE = (
    PROJECT_ROOT
    / "15_design_entanglement_20260726"
    / "release"
    / "v2.0.0-rc1"
)
H5AD_174188 = PROJECT_ROOT / "05_原始数据" / "GSE174188_CELLxGENE_2025-11-08.h5ad"

SEEDS = tuple(range(20260801, 20260821))
CS = np.logspace(-4, 4, 9)
N_SPLITS = 5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def effective_splits(y: np.ndarray, requested: int = N_SPLITS) -> int:
    counts = np.bincount(np.asarray(y, dtype=int))
    positive = counts[counts > 0]
    if len(positive) < 2:
        raise ValueError("Classification target has fewer than two classes")
    return int(min(requested, positive.min()))


def logistic(C: float = 1.0) -> LogisticRegression:
    return LogisticRegression(
        C=float(C),
        class_weight="balanced",
        solver="liblinear",
        max_iter=20_000,
    )


def select_c_prepared(x: np.ndarray, y: np.ndarray, seed: int) -> float:
    splitter = StratifiedKFold(
        effective_splits(y),
        shuffle=True,
        random_state=int(seed),
    )
    splits = list(splitter.split(x, y))
    scores: list[float] = []
    for C in CS:
        fold_scores = []
        for train, test in splits:
            model = logistic(C).fit(x[train], y[train])
            prediction = model.predict_proba(x[test])[:, 1]
            fold_scores.append(roc_auc_score(y[test], prediction))
        scores.append(float(np.mean(fold_scores)))
    return float(CS[int(np.argmax(scores))])


def predict_scaled(
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
    seed: int,
) -> tuple[np.ndarray, float]:
    scaler = StandardScaler().fit(train_x)
    train_z = scaler.transform(train_x)
    test_z = scaler.transform(test_x)
    C = select_c_prepared(train_z, train_y, seed)
    model = logistic(C).fit(train_z, train_y)
    return model.predict_proba(test_z)[:, 1], C


def repeated_oof(
    x: np.ndarray,
    y: np.ndarray,
    seeds: tuple[int, ...] = SEEDS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows = []
    prediction_rows = []
    for repeat, seed in enumerate(seeds, start=1):
        splitter = StratifiedKFold(
            effective_splits(y), shuffle=True, random_state=seed
        )
        prediction = np.full(len(y), np.nan)
        folds = np.zeros(len(y), dtype=int)
        selected_c = np.full(len(y), np.nan)
        for fold, (train, test) in enumerate(splitter.split(x, y), start=1):
            prediction[test], C = predict_scaled(
                x[train], y[train], x[test], seed
            )
            folds[test] = fold
            selected_c[test] = C
        metric_rows.append(
            {
                "repeat": repeat,
                "seed": seed,
                "roc_auc": roc_auc_score(y, prediction),
                "pr_auc": average_precision_score(y, prediction),
                "brier": brier_score_loss(y, prediction),
                "median_selected_c": float(np.median(selected_c)),
            }
        )
        prediction_rows.extend(
            {
                "repeat": repeat,
                "seed": seed,
                "row": int(index),
                "fold": int(fold),
                "y_true": int(label),
                "prob_case": float(probability),
            }
            for index, (fold, label, probability) in enumerate(
                zip(folds, y, prediction)
            )
        )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows)


def design_transformer(
    numeric: list[str],
    categorical: list[str],
) -> ColumnTransformer:
    transformers = []
    if numeric:
        transformers.append(
            (
                "numeric",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="median")),
                        ("scale", StandardScaler()),
                    ]
                ),
                numeric,
            )
        )
    if categorical:
        transformers.append(
            (
                "categorical",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="most_frequent")),
                        (
                            "onehot",
                            OneHotEncoder(
                                handle_unknown="ignore",
                                sparse_output=False,
                            ),
                        ),
                    ]
                ),
                categorical,
            )
        )
    return ColumnTransformer(
        transformers,
        remainder="drop",
        sparse_threshold=0.0,
    )


def encode_design(
    metadata: pd.DataFrame,
    train: np.ndarray,
    test: np.ndarray,
    numeric: list[str],
    categorical: list[str],
) -> tuple[np.ndarray, np.ndarray, ColumnTransformer]:
    transformer = design_transformer(numeric, categorical)
    train_x = transformer.fit_transform(metadata.iloc[train])
    test_x = transformer.transform(metadata.iloc[test])
    return (
        np.asarray(train_x, dtype=float),
        np.asarray(test_x, dtype=float),
        transformer,
    )


def residualize(
    train_x: np.ndarray,
    test_x: np.ndarray,
    train_design: np.ndarray,
    test_design: np.ndarray,
    alpha: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    model = Ridge(alpha=alpha, fit_intercept=True)
    model.fit(train_design, train_x)
    return (
        train_x - model.predict(train_design),
        test_x - model.predict(test_design),
    )


def summarise_metrics(
    metrics: pd.DataFrame,
    groups: list[str],
) -> pd.DataFrame:
    return (
        metrics.groupby(groups, as_index=False)
        .agg(
            n_repeats=("repeat", "nunique"),
            auc_mean=("roc_auc", "mean"),
            auc_sd=("roc_auc", "std"),
            auc_p025=("roc_auc", lambda x: np.quantile(x, 0.025)),
            auc_p975=("roc_auc", lambda x: np.quantile(x, 0.975)),
            pr_auc_mean=("pr_auc", "mean"),
            brier_mean=("brier", "mean"),
        )
        .sort_values(groups)
    )


def information_fraction(y: np.ndarray, design: np.ndarray) -> dict[str, float]:
    y_float = np.asarray(y, dtype=float)
    centered = y_float - y_float.mean()
    denominator = float(centered @ centered)
    augmented = np.column_stack([np.ones(len(y_float)), design])
    fitted = augmented @ np.linalg.lstsq(augmented, y_float, rcond=None)[0]
    residual = y_float - fitted
    numerator = float(residual @ residual)
    fraction = numerator / denominator if denominator > 0 else np.nan
    return {
        "r2_y_on_design": float(1.0 - fraction),
        "information_fraction": float(fraction),
        "variance_inflation": float(1.0 / fraction) if fraction > 0 else np.inf,
        "design_rank": int(np.linalg.matrix_rank(augmented)),
        "n_donors": int(len(y_float)),
    }
