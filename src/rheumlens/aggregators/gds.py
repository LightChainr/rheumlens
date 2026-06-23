from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import pairwise_distances

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import as_dense, grouped_indices


@dataclass
class GraphDonorSignatureAggregator(DonorAggregator):
    n_nodes: int = 32
    graph_k: int = 5
    n_spectral: int = 8
    temperature: float = 1.0
    max_fit_cells: int = 200_000
    seed: int = 0

    def __post_init__(self) -> None:
        self.kmeans_: MiniBatchKMeans | None = None
        self.eigenvectors_: np.ndarray | None = None
        self.node_coords_: np.ndarray | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "GraphDonorSignatureAggregator":
        subset = data.subset_donors(train_donors)
        X = as_dense(subset.X)
        rng = np.random.default_rng(self.seed)
        if X.shape[0] > self.max_fit_cells:
            X = X[rng.choice(X.shape[0], self.max_fit_cells, replace=False)]
        self.kmeans_ = MiniBatchKMeans(
            n_clusters=min(self.n_nodes, X.shape[0]),
            random_state=self.seed,
            batch_size=min(4096, X.shape[0]),
            n_init=10,
        ).fit(X)
        centers = self.kmeans_.cluster_centers_
        distances = pairwise_distances(centers)
        adjacency = np.zeros_like(distances)
        for i in range(len(centers)):
            neighbors = np.argsort(distances[i])[1 : self.graph_k + 1]
            scale = np.median(distances[i, neighbors]) + 1e-12
            adjacency[i, neighbors] = np.exp(-np.square(distances[i, neighbors]) / (2 * scale**2))
        adjacency = np.maximum(adjacency, adjacency.T)
        degree = adjacency.sum(axis=1)
        inv_sqrt = np.diag(1.0 / np.sqrt(np.maximum(degree, 1e-12)))
        laplacian = np.eye(len(centers)) - inv_sqrt @ adjacency @ inv_sqrt
        _, eigenvectors = eigh(laplacian)
        self.eigenvectors_ = eigenvectors[:, 1 : 1 + min(self.n_spectral, len(centers) - 1)]
        self.node_coords_ = centers
        return self

    def _soft_mass(self, X: np.ndarray) -> np.ndarray:
        assert self.node_coords_ is not None
        distances = pairwise_distances(X, self.node_coords_)
        logits = -distances / max(self.temperature, 1e-6)
        logits -= logits.max(axis=1, keepdims=True)
        weights = np.exp(logits)
        weights /= weights.sum(axis=1, keepdims=True)
        return weights.mean(axis=0)

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.eigenvectors_ is None:
            raise RuntimeError("aggregator is not fitted")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rows = []
        for donor in order:
            mass = self._soft_mass(as_dense(subset.X[groups[donor]]))
            entropy = -np.sum(mass * np.log(mass + 1e-12))
            spectral = mass @ self.eigenvectors_
            center = mass @ self.node_coords_
            rows.append(np.concatenate([mass, spectral, center, [entropy]]))
        names = (
            [f"graph_mass_{i}" for i in range(self.eigenvectors_.shape[0])]
            + [f"graph_spectral_{i}" for i in range(self.eigenvectors_.shape[1])]
            + [f"graph_center_{i}" for i in range(self.node_coords_.shape[1])]
            + ["graph_entropy"]
        )
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, names)
