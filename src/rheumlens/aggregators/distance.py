from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import sqrtm
from sklearn.cluster import MiniBatchKMeans
from sklearn.covariance import LedoitWolf
from sklearn.metrics import pairwise_distances

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import as_dense, grouped_indices, seed_sequence


def _psd_sqrt(matrix: np.ndarray) -> np.ndarray:
    value = sqrtm((matrix + matrix.T) / 2)
    return np.real_if_close(value).astype(np.float64)


def bures_distance_squared(mean_a: np.ndarray, cov_a: np.ndarray, mean_b: np.ndarray, cov_b: np.ndarray) -> float:
    root_a = _psd_sqrt(cov_a)
    middle = root_a @ cov_b @ root_a
    transport = _psd_sqrt(middle)
    covariance_term = np.trace(cov_a + cov_b - 2 * transport)
    return float(np.square(mean_a - mean_b).sum() + max(float(covariance_term), 0.0))


@dataclass
class GaussianBuresReferenceAggregator(DonorAggregator):
    """Distances to training donor Gaussian references using the true Bures formula.

    Direct donor-donor kernels are handled separately; this fixed-vector variant uses distances
    to a controlled number of training reference donors, selected without test information.
    """

    n_references: int = 16
    seed: int = 0

    def __post_init__(self) -> None:
        self.references_: list[tuple[str, np.ndarray, np.ndarray]] = []

    @staticmethod
    def _gaussian(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if X.shape[0] < 2:
            return X.mean(axis=0), np.eye(X.shape[1]) * 1e-6
        return X.mean(axis=0), LedoitWolf().fit(X).covariance_

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "GaussianBuresReferenceAggregator":
        subset = data.subset_donors(train_donors)
        groups = grouped_indices(subset.donor_ids)
        rng = np.random.default_rng(self.seed)
        candidates = np.asarray(list(groups), dtype=str)
        if len(candidates) > self.n_references:
            candidates = rng.choice(candidates, self.n_references, replace=False)
        self.references_ = []
        for donor in candidates:
            self.references_.append((donor, *self._gaussian(as_dense(subset.X[groups[donor]]))))
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if not self.references_:
            raise RuntimeError("aggregator is not fitted")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rows = []
        for donor in order:
            mean, cov = self._gaussian(as_dense(subset.X[groups[donor]]))
            rows.append([bures_distance_squared(mean, cov, ref_mean, ref_cov) for _, ref_mean, ref_cov in self.references_])
        return DonorFeatures(np.asarray(rows, dtype=np.float32), order, [f"bures_to_{d}" for d, _, _ in self.references_])


@dataclass
class SlicedWassersteinReferenceAggregator(DonorAggregator):
    n_directions: int = 64
    n_quantiles: int = 64
    n_references: int = 16
    seed: int = 0

    def __post_init__(self) -> None:
        self.directions_: np.ndarray | None = None
        self.references_: dict[str, np.ndarray] = {}

    def _signature(self, X: np.ndarray) -> np.ndarray:
        assert self.directions_ is not None
        projections = X @ self.directions_.T
        q = np.linspace(0.0, 1.0, self.n_quantiles)
        return np.quantile(projections, q, axis=0).T

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "SlicedWassersteinReferenceAggregator":
        subset = data.subset_donors(train_donors)
        rng = np.random.default_rng(self.seed)
        directions = rng.normal(size=(self.n_directions, subset.X.shape[1]))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True) + 1e-12
        self.directions_ = directions
        groups = grouped_indices(subset.donor_ids)
        candidates = np.asarray(list(groups), dtype=str)
        if len(candidates) > self.n_references:
            candidates = rng.choice(candidates, self.n_references, replace=False)
        self.references_ = {d: self._signature(as_dense(subset.X[groups[d]])) for d in candidates}
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if not self.references_:
            raise RuntimeError("aggregator is not fitted")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        keys = list(self.references_)
        rows = []
        for donor in order:
            signature = self._signature(as_dense(subset.X[groups[donor]]))
            rows.append([float(np.square(signature - self.references_[key]).mean()) for key in keys])
        return DonorFeatures(np.asarray(rows, dtype=np.float32), order, [f"sw2_to_{x}" for x in keys])


@dataclass
class PrototypeAggregator(DonorAggregator):
    n_prototypes: int = 32
    max_fit_cells: int = 200_000
    temperature: float = 1.0
    seed: int = 0

    def __post_init__(self) -> None:
        self.kmeans_: MiniBatchKMeans | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "PrototypeAggregator":
        subset = data.subset_donors(train_donors)
        X = as_dense(subset.X)
        rng = np.random.default_rng(self.seed)
        if X.shape[0] > self.max_fit_cells:
            X = X[rng.choice(X.shape[0], self.max_fit_cells, replace=False)]
        self.kmeans_ = MiniBatchKMeans(
            n_clusters=min(self.n_prototypes, X.shape[0]),
            random_state=self.seed,
            batch_size=min(4096, max(256, X.shape[0])),
            n_init=10,
        ).fit(X)
        return self

    def _features(self, X: np.ndarray) -> np.ndarray:
        if self.kmeans_ is None:
            raise RuntimeError("aggregator is not fitted")
        distances = pairwise_distances(X, self.kmeans_.cluster_centers_)
        logits = -distances / max(self.temperature, 1e-6)
        logits -= logits.max(axis=1, keepdims=True)
        weights = np.exp(logits)
        weights /= weights.sum(axis=1, keepdims=True)
        occupancy = weights.mean(axis=0)
        mean_distance = (weights * distances).sum(axis=0) / np.maximum(weights.sum(axis=0), 1e-8)
        q90 = np.quantile(distances, 0.9, axis=0)
        entropy = -np.sum(weights * np.log(weights + 1e-12), axis=1).mean(keepdims=True)
        return np.concatenate([occupancy, mean_distance, q90, entropy])

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rows = [self._features(as_dense(subset.X[groups[d]])) for d in order]
        k = self.kmeans_.n_clusters if self.kmeans_ is not None else self.n_prototypes
        names = (
            [f"proto_occupancy_{i}" for i in range(k)]
            + [f"proto_mean_distance_{i}" for i in range(k)]
            + [f"proto_q90_distance_{i}" for i in range(k)]
            + ["soft_assignment_entropy"]
        )
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, names)
