from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from rheumlens.postprocessors.base import DonorPostprocessor
from rheumlens.types import DonorFeatures, FitContext


@dataclass
class CovariateResidualizer(DonorPostprocessor):
    """Fold-contained residualization of donor features against available covariates.

    `FitContext.donor_covariates` must be a DataFrame with a `donor_id` column. Numeric and
    categorical columns are inferred unless explicitly supplied. The transformation is fitted
    only on outer-training donors and is therefore suitable for sensitivity analyses, not a
    causal deconfounding claim.
    """

    alpha: float = 1.0
    numeric_columns: tuple[str, ...] | None = None
    categorical_columns: tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        self.encoder_: ColumnTransformer | None = None
        self.regression_: Ridge | None = None
        self.columns_: list[str] = []

    @staticmethod
    def _aligned_covariates(features: DonorFeatures, context: FitContext) -> pd.DataFrame:
        if context.donor_covariates is None:
            raise ValueError("CovariateResidualizer requires FitContext.donor_covariates")
        frame = context.donor_covariates.copy()
        if "donor_id" not in frame.columns:
            raise ValueError("donor covariates require a donor_id column")
        frame["donor_id"] = frame["donor_id"].astype(str)
        if frame["donor_id"].duplicated().any():
            raise ValueError("donor covariates contain duplicate donor_id rows")
        aligned = frame.set_index("donor_id").reindex(features.donor_ids.astype(str))
        if aligned.isna().all(axis=1).any():
            missing = aligned.index[aligned.isna().all(axis=1)].tolist()
            raise ValueError(f"missing covariate rows for donors: {missing[:10]}")
        return aligned.reset_index(drop=True)

    def fit(self, features: DonorFeatures, context: FitContext) -> "CovariateResidualizer":
        cov = self._aligned_covariates(features, context)
        numeric = list(self.numeric_columns or tuple(cov.select_dtypes(include=[np.number, "bool"]).columns))
        categorical = list(self.categorical_columns or tuple(c for c in cov.columns if c not in numeric))
        # Avoid including accidental labels in the sensitivity model.
        forbidden = {"disease", "label", "y", "case_control"}
        numeric = [c for c in numeric if c.lower() not in forbidden]
        categorical = [c for c in categorical if c.lower() not in forbidden]
        if not numeric and not categorical:
            raise ValueError("no usable covariate columns")
        transformers = []
        if numeric:
            transformers.append(
                (
                    "num",
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
                    "cat",
                    Pipeline(
                        [
                            ("impute", SimpleImputer(strategy="most_frequent")),
                            ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
                        ]
                    ),
                    categorical,
                )
            )
        self.encoder_ = ColumnTransformer(transformers, remainder="drop", sparse_threshold=0.0)
        design = np.asarray(self.encoder_.fit_transform(cov), dtype=np.float64)
        self.regression_ = Ridge(alpha=self.alpha, fit_intercept=True).fit(design, features.X)
        self.columns_ = numeric + categorical
        return self

    def transform(self, features: DonorFeatures, context: FitContext) -> DonorFeatures:
        if self.encoder_ is None or self.regression_ is None:
            raise RuntimeError("CovariateResidualizer is not fitted")
        cov = self._aligned_covariates(features, context)
        design = np.asarray(self.encoder_.transform(cov), dtype=np.float64)
        fitted = self.regression_.predict(design)
        residual = np.asarray(features.X, dtype=np.float64) - fitted
        diagnostics = dict(features.diagnostics)
        diagnostics["covariate_residualized"] = True
        diagnostics["covariate_columns"] = self.columns_
        return DonorFeatures(
            residual.astype(np.float32),
            features.donor_ids.copy(),
            features.feature_names.copy(),
            uncertainty=features.uncertainty,
            diagnostics=diagnostics,
        )
