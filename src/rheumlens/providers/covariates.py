from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from rheumlens.providers.base import RepresentationProvider
from rheumlens.types import CellDataset, DonorFeatures, FitContext


@dataclass
class DonorCovariateProvider(RepresentationProvider):
    """Construct a fold-contained covariates-only donor representation."""

    numeric_columns: tuple[str, ...] | None = None
    categorical_columns: tuple[str, ...] | None = None
    output_level: str = "donor"

    def __post_init__(self) -> None:
        self.encoder_: ColumnTransformer | None = None
        self.feature_names_: list[str] = []

    @staticmethod
    def _table(context: FitContext, donors: list[str]) -> pd.DataFrame:
        if context.donor_covariates is None:
            raise ValueError("DonorCovariateProvider requires donor_covariates")
        frame = context.donor_covariates.copy()
        if "donor_id" not in frame:
            raise ValueError("donor covariates require donor_id")
        frame["donor_id"] = frame["donor_id"].astype(str)
        if frame.donor_id.duplicated().any():
            raise ValueError("duplicate donor covariate rows")
        aligned = frame.set_index("donor_id").reindex(list(map(str, donors)))
        if aligned.isna().all(axis=1).any():
            missing = aligned.index[aligned.isna().all(axis=1)].tolist()
            raise ValueError(f"missing donor covariates: {missing[:10]}")
        return aligned.reset_index(drop=True)

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "DonorCovariateProvider":
        table = self._table(context, train_donors)
        numeric = list(self.numeric_columns or tuple(table.select_dtypes(include=[np.number, "bool"]).columns))
        categorical = list(self.categorical_columns or tuple(c for c in table.columns if c not in numeric))
        forbidden = {"disease", "label", "y", "case_control"}
        numeric = [c for c in numeric if c.lower() not in forbidden]
        categorical = [c for c in categorical if c.lower() not in forbidden]
        transformers = []
        if numeric:
            transformers.append(
                (
                    "num",
                    Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]),
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
        if not transformers:
            raise ValueError("no usable covariate columns")
        self.encoder_ = ColumnTransformer(transformers, sparse_threshold=0.0)
        self.encoder_.fit(table)
        self.feature_names_ = self.encoder_.get_feature_names_out().astype(str).tolist()
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.encoder_ is None:
            raise RuntimeError("provider is not fitted")
        table = self._table(context, donors)
        X = np.asarray(self.encoder_.transform(table), dtype=np.float32)
        return DonorFeatures(X, np.asarray(donors, dtype=str), self.feature_names_)
