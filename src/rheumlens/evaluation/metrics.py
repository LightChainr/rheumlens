from __future__ import annotations

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score


def binary_metrics(y: np.ndarray, score: np.ndarray) -> dict[str, float]:
    y = np.asarray(y).astype(int)
    score = np.asarray(score).astype(float)
    if len(np.unique(y)) < 2:
        return {"roc_auc": np.nan, "pr_auc": np.nan, "brier": np.nan}
    return {
        "roc_auc": float(roc_auc_score(y, score)),
        "pr_auc": float(average_precision_score(y, score)),
        "brier": float(brier_score_loss(y, score)),
    }
