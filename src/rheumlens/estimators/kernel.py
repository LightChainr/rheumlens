from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from rheumlens.estimators.base import DonorEstimator


@dataclass
class RBFKernelSVMEstimator(DonorEstimator):
    C: float = 1.0
    gamma: str | float = "scale"
    class_weight: str | dict | None = "balanced"
    random_state: int = 0

    def __post_init__(self) -> None:
        self.model_: Pipeline | None = None

    def fit(self, X: np.ndarray, y: np.ndarray, sample_weight: np.ndarray | None = None) -> "RBFKernelSVMEstimator":
        self.model_ = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                (
                    "svc",
                    SVC(
                        C=self.C,
                        kernel="rbf",
                        gamma=self.gamma,
                        probability=True,
                        class_weight=self.class_weight,
                        random_state=self.random_state,
                    ),
                ),
            ]
        )
        fit_params = {}
        if sample_weight is not None:
            fit_params["svc__sample_weight"] = sample_weight
        self.model_.fit(X, y, **fit_params)
        return self

    def predict_score(self, X: np.ndarray) -> np.ndarray:
        if self.model_ is None:
            raise RuntimeError("estimator is not fitted")
        return self.model_.predict_proba(X)[:, 1]
