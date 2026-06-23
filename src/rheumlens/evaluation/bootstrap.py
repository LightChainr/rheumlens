from __future__ import annotations

import numpy as np
from sklearn.metrics import roc_auc_score


def stratified_bootstrap_indices(y: np.ndarray, n_reps: int, seed: int = 0) -> list[np.ndarray]:
    y = np.asarray(y).astype(int)
    rng = np.random.default_rng(seed)
    classes = [np.flatnonzero(y == value) for value in np.unique(y)]
    return [np.concatenate([rng.choice(idx, len(idx), replace=True) for idx in classes]) for _ in range(n_reps)]


def paired_auc_bootstrap(
    y: np.ndarray,
    score_a: np.ndarray,
    score_b: np.ndarray,
    n_reps: int = 10_000,
    seed: int = 0,
) -> dict[str, float | list[float]]:
    y = np.asarray(y).astype(int)
    score_a = np.asarray(score_a)
    score_b = np.asarray(score_b)
    observed = roc_auc_score(y, score_a) - roc_auc_score(y, score_b)
    values = []
    for idx in stratified_bootstrap_indices(y, n_reps, seed):
        values.append(roc_auc_score(y[idx], score_a[idx]) - roc_auc_score(y[idx], score_b[idx]))
    array = np.asarray(values)
    return {
        "observed_delta_auc": float(observed),
        "bootstrap_mean": float(array.mean()),
        "ci_low": float(np.quantile(array, 0.025)),
        "ci_high": float(np.quantile(array, 0.975)),
        "superiority_frequency": float(np.mean(array > 0)),
        "values": array.tolist(),
    }
