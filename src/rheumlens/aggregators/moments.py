from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import skew
from sklearn.covariance import LedoitWolf

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import as_dense, grouped_indices


def _blocks(data: CellDataset, donors: list[str]) -> tuple[CellDataset, np.ndarray, dict[str, np.ndarray]]:
    subset = data.subset_donors(donors)
    groups = grouped_indices(subset.donor_ids)
    order = np.asarray(list(groups), dtype=str)
    return subset, order, groups


@dataclass
class MeanVarianceAggregator(DonorAggregator):
    ddof: int = 1

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "MeanVarianceAggregator":
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        subset, order, groups = _blocks(data, donors)
        rows = []
        for donor in order:
            X = as_dense(subset.X[groups[donor]])
            rows.append(np.concatenate([X.mean(axis=0), X.var(axis=0, ddof=min(self.ddof, max(0, X.shape[0] - 1))) ]))
        names = [f"mean_{x}" for x in subset.feature_names] + [f"var_{x}" for x in subset.feature_names]
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, names)


@dataclass
class MedianMADAggregator(DonorAggregator):
    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "MedianMADAggregator":
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        subset, order, groups = _blocks(data, donors)
        rows = []
        for donor in order:
            X = as_dense(subset.X[groups[donor]])
            med = np.median(X, axis=0)
            mad = np.median(np.abs(X - med), axis=0)
            rows.append(np.concatenate([med, mad]))
        names = [f"median_{x}" for x in subset.feature_names] + [f"mad_{x}" for x in subset.feature_names]
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, names)


@dataclass
class MeanSkewAggregator(DonorAggregator):
    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "MeanSkewAggregator":
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        subset, order, groups = _blocks(data, donors)
        rows = []
        for donor in order:
            X = as_dense(subset.X[groups[donor]])
            values = skew(X, axis=0, bias=False, nan_policy="omit")
            rows.append(np.concatenate([X.mean(axis=0), np.nan_to_num(values)]))
        names = [f"mean_{x}" for x in subset.feature_names] + [f"skew_{x}" for x in subset.feature_names]
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, names)


@dataclass
class QuantileAggregator(DonorAggregator):
    quantiles: tuple[float, ...] = (0.1, 0.5, 0.9, 0.95, 0.99)
    include_mean: bool = True

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "QuantileAggregator":
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        subset, order, groups = _blocks(data, donors)
        rows = []
        for donor in order:
            X = as_dense(subset.X[groups[donor]])
            parts = [X.mean(axis=0)] if self.include_mean else []
            parts.extend(np.quantile(X, q, axis=0) for q in self.quantiles)
            rows.append(np.concatenate(parts))
        names = []
        if self.include_mean:
            names.extend(f"mean_{x}" for x in subset.feature_names)
        for q in self.quantiles:
            names.extend(f"q{q:g}_{x}" for x in subset.feature_names)
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, list(names))


@dataclass
class TailFractionAggregator(DonorAggregator):
    quantiles: tuple[float, ...] = (0.95, 0.99)
    controls_only: bool = True

    def __post_init__(self) -> None:
        self.thresholds_: dict[float, np.ndarray] = {}

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "TailFractionAggregator":
        subset = data.subset_donors(train_donors)
        if self.controls_only:
            labels = np.asarray([context.y_by_donor[str(d)] for d in subset.donor_ids])
            mask = labels == 0
            if not mask.any():
                raise ValueError("no training controls for tail thresholds")
            X = as_dense(subset.X[mask])
        else:
            X = as_dense(subset.X)
        self.thresholds_ = {q: np.quantile(X, q, axis=0) for q in self.quantiles}
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if not self.thresholds_:
            raise RuntimeError("aggregator is not fitted")
        subset, order, groups = _blocks(data, donors)
        rows = []
        for donor in order:
            X = as_dense(subset.X[groups[donor]])
            rows.append(np.concatenate([(X > self.thresholds_[q]).mean(axis=0) for q in self.quantiles]))
        names = [f"tail_q{q:g}_{x}" for q in self.quantiles for x in subset.feature_names]
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, names)


@dataclass
class ShrinkageCovarianceAggregator(DonorAggregator):
    include_mean: bool = True
    diagonal_only: bool = False

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "ShrinkageCovarianceAggregator":
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        subset, order, groups = _blocks(data, donors)
        rows = []
        p = subset.X.shape[1]
        tri = np.triu_indices(p)
        for donor in order:
            X = as_dense(subset.X[groups[donor]])
            if X.shape[0] < 2:
                covariance = np.zeros((p, p), dtype=np.float64)
            else:
                covariance = LedoitWolf().fit(X).covariance_
            cov_values = np.diag(covariance) if self.diagonal_only else covariance[tri]
            parts = [X.mean(axis=0), cov_values] if self.include_mean else [cov_values]
            rows.append(np.concatenate(parts))
        names: list[str] = []
        if self.include_mean:
            names.extend(f"mean_{x}" for x in subset.feature_names)
        if self.diagonal_only:
            names.extend(f"covdiag_{x}" for x in subset.feature_names)
        else:
            names.extend(f"cov_{i}_{j}" for i, j in zip(*tri, strict=True))
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, names)
