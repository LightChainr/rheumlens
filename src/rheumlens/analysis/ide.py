from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


@dataclass
class InvariantEvidenceSelector:
    """Select features with stable sign across training environments only."""

    min_frequency: float = 0.7
    C: float = 0.1
    max_cohort_auc: float = 0.8
    random_state: int = 0

    def __post_init__(self) -> None:
        self.selected_: np.ndarray | None = None
        self.signs_: np.ndarray | None = None
        self.scaler_: StandardScaler | None = None

    def fit(self, X: np.ndarray, y: np.ndarray, environments: np.ndarray) -> "InvariantEvidenceSelector":
        environments = np.asarray(environments).astype(str)
        coefficients = []
        for env in np.unique(environments):
            mask = environments != env
            if len(np.unique(y[mask])) < 2:
                continue
            scaler = StandardScaler().fit(X[mask])
            model = LogisticRegression(
                penalty="l1", solver="liblinear", C=self.C, class_weight="balanced", random_state=self.random_state
            ).fit(scaler.transform(X[mask]), y[mask])
            coefficients.append(model.coef_[0])
        if not coefficients:
            raise ValueError("insufficient training environments for IDE")
        coef = np.vstack(coefficients)
        nonzero_frequency = np.mean(np.abs(coef) > 1e-12, axis=0)
        positive = np.mean(coef > 0, axis=0)
        negative = np.mean(coef < 0, axis=0)
        sign_consistency = np.maximum(positive, negative)
        self.selected_ = np.flatnonzero((nonzero_frequency >= self.min_frequency) & (sign_consistency >= self.min_frequency))
        if len(self.selected_) == 0:
            # Deterministic conservative fallback: highest median absolute coefficients.
            rank = np.argsort(np.median(np.abs(coef), axis=0))[::-1]
            self.selected_ = rank[: min(10, X.shape[1])]
        self.signs_ = np.sign(np.median(coef[:, self.selected_], axis=0))
        self.scaler_ = StandardScaler().fit(X[:, self.selected_])
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self.selected_ is None or self.scaler_ is None:
            raise RuntimeError("IDE selector is not fitted")
        return self.scaler_.transform(X[:, self.selected_])
