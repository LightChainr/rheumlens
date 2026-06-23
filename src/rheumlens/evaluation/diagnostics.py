from __future__ import annotations

import numpy as np
from sklearn.metrics import roc_auc_score


def feature_cellcount_correlations(X: np.ndarray, n_cells: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    n_cells = np.asarray(n_cells, dtype=float)
    values = []
    for col in X.T:
        if np.std(col) <= 1e-12 or np.std(n_cells) <= 1e-12:
            values.append(np.nan)
        else:
            values.append(np.corrcoef(col, n_cells)[0, 1])
    return np.asarray(values)


def uncertainty_only_auc(y: np.ndarray, uncertainty: np.ndarray) -> float:
    score = np.mean(uncertainty, axis=1)
    if len(np.unique(y)) < 2:
        return np.nan
    auc = roc_auc_score(y, score)
    return float(max(auc, 1.0 - auc))
