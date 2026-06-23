from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import pairwise_distances

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import as_dense, grouped_indices


@dataclass
class REDAggregator(DonorAggregator):
    n_states: int = 24
    max_fit_cells: int = 200_000
    temperature: float = 1.0
    tail_quantile: float = 0.95
    seed: int = 0

    def __post_init__(self) -> None:
        self.kmeans_: MiniBatchKMeans | None = None
        self.control_occupancy_: np.ndarray | None = None
        self.control_dispersion_: np.ndarray | None = None
        self.tail_thresholds_: np.ndarray | None = None

    def _soft_weights(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.kmeans_ is None:
            raise RuntimeError("aggregator is not fitted")
        distances = pairwise_distances(X, self.kmeans_.cluster_centers_)
        logits = -distances / max(self.temperature, 1e-6)
        logits -= logits.max(axis=1, keepdims=True)
        weights = np.exp(logits)
        weights /= weights.sum(axis=1, keepdims=True)
        return weights, distances

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "REDAggregator":
        subset = data.subset_donors(train_donors)
        control_mask = np.asarray([context.y_by_donor[str(d)] == 0 for d in subset.donor_ids])
        controls = as_dense(subset.X[control_mask])
        if controls.shape[0] < self.n_states:
            raise ValueError("not enough control cells for RED states")
        rng = np.random.default_rng(self.seed)
        fit_X = controls
        if fit_X.shape[0] > self.max_fit_cells:
            fit_X = fit_X[rng.choice(fit_X.shape[0], self.max_fit_cells, replace=False)]
        self.kmeans_ = MiniBatchKMeans(
            n_clusters=self.n_states,
            random_state=self.seed,
            batch_size=min(4096, fit_X.shape[0]),
            n_init=10,
        ).fit(fit_X)
        weights, distances = self._soft_weights(controls)
        self.control_occupancy_ = weights.mean(axis=0)
        self.control_dispersion_ = (weights * distances).sum(axis=0) / np.maximum(weights.sum(axis=0), 1e-12)
        self.tail_thresholds_ = np.quantile(distances, self.tail_quantile, axis=0)
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.control_occupancy_ is None or self.control_dispersion_ is None or self.tail_thresholds_ is None:
            raise RuntimeError("aggregator is not fitted")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rows = []
        for donor in order:
            X = as_dense(subset.X[groups[donor]])
            weights, distances = self._soft_weights(X)
            occupancy = weights.mean(axis=0)
            dispersion = (weights * distances).sum(axis=0) / np.maximum(weights.sum(axis=0), 1e-12)
            tail = (distances > self.tail_thresholds_[None, :]).mean(axis=0)
            missing = (occupancy < 1.0 / max(len(X), 1)).astype(float)
            rows.append(
                np.concatenate(
                    [
                        occupancy - self.control_occupancy_,
                        dispersion - self.control_dispersion_,
                        tail,
                        missing,
                    ]
                )
            )
        names = (
            [f"abundance_deviation_{i}" for i in range(self.n_states)]
            + [f"dispersion_deviation_{i}" for i in range(self.n_states)]
            + [f"tail_burden_{i}" for i in range(self.n_states)]
            + [f"missing_state_{i}" for i in range(self.n_states)]
        )
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, names)
