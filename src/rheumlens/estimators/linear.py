from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from rheumlens.estimators.base import DonorEstimator


@dataclass
class LogisticL2Estimator(DonorEstimator):
    C: float = 1.0
    class_weight: str | dict | None = "balanced"
    max_iter: int = 5000
    random_state: int = 0

    def __post_init__(self) -> None:
        self.model_: Pipeline | None = None

    def fit(self, X: np.ndarray, y: np.ndarray, sample_weight: np.ndarray | None = None) -> "LogisticL2Estimator":
        self.model_ = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median", add_indicator=False)),
                ("scaler", StandardScaler()),
                (
                    "logistic",
                    LogisticRegression(
                        C=self.C,
                        penalty="l2",
                        solver="liblinear",
                        class_weight=self.class_weight,
                        max_iter=self.max_iter,
                        random_state=self.random_state,
                    ),
                ),
            ]
        )
        fit_params = {}
        if sample_weight is not None:
            fit_params["logistic__sample_weight"] = sample_weight
        self.model_.fit(X, y, **fit_params)
        return self

    def predict_score(self, X: np.ndarray) -> np.ndarray:
        if self.model_ is None:
            raise RuntimeError("estimator is not fitted")
        return self.model_.predict_proba(X)[:, 1]


@dataclass
class ElasticNetEstimator(DonorEstimator):
    C: float = 1.0
    l1_ratio: float = 0.5
    class_weight: str | dict | None = "balanced"
    max_iter: int = 5000
    random_state: int = 0

    def __post_init__(self) -> None:
        self.model_: Pipeline | None = None

    def fit(self, X: np.ndarray, y: np.ndarray, sample_weight: np.ndarray | None = None) -> "ElasticNetEstimator":
        self.model_ = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                (
                    "logistic",
                    LogisticRegression(
                        C=self.C,
                        penalty="elasticnet",
                        solver="saga",
                        l1_ratio=self.l1_ratio,
                        class_weight=self.class_weight,
                        max_iter=self.max_iter,
                        random_state=self.random_state,
                    ),
                ),
            ]
        )
        fit_params = {}
        if sample_weight is not None:
            fit_params["logistic__sample_weight"] = sample_weight
        self.model_.fit(X, y, **fit_params)
        return self

    def predict_score(self, X: np.ndarray) -> np.ndarray:
        if self.model_ is None:
            raise RuntimeError("estimator is not fitted")
        return self.model_.predict_proba(X)[:, 1]
