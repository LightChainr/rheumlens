from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler

from rheumlens.providers.base import RepresentationProvider
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import as_dense, grouped_indices


def _donor_means(data: CellDataset, donors: list[str]) -> tuple[np.ndarray, np.ndarray]:
    subset = data.subset_donors(donors)
    groups = grouped_indices(subset.donor_ids)
    order = np.asarray(list(groups), dtype=str)
    rows = []
    for donor in order:
        block = subset.X[groups[donor]]
        rows.append(np.asarray(block.mean(axis=0)).ravel())
    return np.vstack(rows).astype(np.float32), order


@dataclass
class DonorMeanHVGProvider(RepresentationProvider):
    n_hvg: int = 2000
    output_level: str = "donor"

    def __post_init__(self) -> None:
        self.selected_: np.ndarray | None = None
        self.feature_names_: list[str] | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "DonorMeanHVGProvider":
        means, _ = _donor_means(data, train_donors)
        variances = np.var(means, axis=0, ddof=1)
        n = min(self.n_hvg, means.shape[1])
        self.selected_ = np.argsort(variances)[-n:]
        self.feature_names_ = data.feature_names[self.selected_].astype(str).tolist()
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.selected_ is None or self.feature_names_ is None:
            raise RuntimeError("provider is not fitted")
        means, order = _donor_means(data, donors)
        return DonorFeatures(means[:, self.selected_], order, self.feature_names_)


@dataclass
class DonorExpressionPCAProvider(RepresentationProvider):
    n_hvg: int = 2000
    n_components: int = 25
    output_level: str = "donor"

    def __post_init__(self) -> None:
        self.selected_: np.ndarray | None = None
        self.scaler_: StandardScaler | None = None
        self.pca_: PCA | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "DonorExpressionPCAProvider":
        means, _ = _donor_means(data, train_donors)
        variances = np.var(means, axis=0, ddof=1)
        n = min(self.n_hvg, means.shape[1])
        self.selected_ = np.argsort(variances)[-n:]
        selected = means[:, self.selected_]
        self.scaler_ = StandardScaler().fit(selected)
        scaled = self.scaler_.transform(selected)
        n_comp = min(self.n_components, scaled.shape[0] - 1, scaled.shape[1])
        if n_comp < 1:
            raise ValueError("insufficient training donors for PCA")
        self.pca_ = PCA(n_components=n_comp, random_state=context.random_state).fit(scaled)
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.selected_ is None or self.scaler_ is None or self.pca_ is None:
            raise RuntimeError("provider is not fitted")
        means, order = _donor_means(data, donors)
        X = self.pca_.transform(self.scaler_.transform(means[:, self.selected_])).astype(np.float32)
        names = [f"donor_pc_{i+1}" for i in range(X.shape[1])]
        return DonorFeatures(X, order, names, diagnostics={"explained_variance": self.pca_.explained_variance_ratio_.tolist()})


@dataclass
class CellPCAProvider(RepresentationProvider):
    n_components: int = 50
    max_fit_cells: int = 200_000
    standardize: bool = True
    output_level: str = "cell"

    def __post_init__(self) -> None:
        self.scaler_: StandardScaler | None = None
        self.pca_: PCA | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "CellPCAProvider":
        subset = data.subset_donors(train_donors)
        rng = np.random.default_rng(context.random_state)
        indices = np.arange(subset.X.shape[0])
        if len(indices) > self.max_fit_cells:
            indices = rng.choice(indices, self.max_fit_cells, replace=False)
        X = as_dense(subset.X[indices])
        if self.standardize:
            self.scaler_ = StandardScaler().fit(X)
            X = self.scaler_.transform(X)
        n_comp = min(self.n_components, X.shape[0] - 1, X.shape[1])
        self.pca_ = PCA(n_components=n_comp, random_state=context.random_state).fit(X)
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> CellDataset:
        if self.pca_ is None:
            raise RuntimeError("provider is not fitted")
        subset = data.subset_donors(donors)
        X = as_dense(subset.X)
        if self.scaler_ is not None:
            X = self.scaler_.transform(X)
        Z = self.pca_.transform(X).astype(np.float32)
        return CellDataset(
            X=Z,
            cell_ids=subset.cell_ids,
            donor_ids=subset.donor_ids,
            y=subset.y,
            feature_names=np.asarray([f"cell_pc_{i+1}" for i in range(Z.shape[1])]),
            cell_types=subset.cell_types,
            cohorts=subset.cohorts,
            metadata=subset.metadata,
            name=f"{subset.name}_cell_pca",
        )


@dataclass
class RawPseudobulkProvider(RepresentationProvider):
    min_total_count: int = 10
    n_features: int = 2000
    output_level: str = "donor"

    def __post_init__(self) -> None:
        self.selected_: np.ndarray | None = None
        self.feature_names_: list[str] | None = None

    @staticmethod
    def _aggregate(data: CellDataset, donors: list[str]) -> tuple[np.ndarray, np.ndarray]:
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rows = []
        for donor in order:
            rows.append(np.asarray(subset.X[groups[donor]].sum(axis=0)).ravel())
        return np.vstack(rows).astype(np.float64), order

    @staticmethod
    def _logcpm(counts: np.ndarray) -> np.ndarray:
        library = counts.sum(axis=1, keepdims=True)
        return np.log1p(counts / np.maximum(library, 1.0) * 1e6)

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "RawPseudobulkProvider":
        counts, _ = self._aggregate(data, train_donors)
        keep = counts.sum(axis=0) >= self.min_total_count
        candidate = np.flatnonzero(keep)
        transformed = self._logcpm(counts[:, candidate])
        variances = transformed.var(axis=0, ddof=1)
        n = min(self.n_features, len(candidate))
        self.selected_ = candidate[np.argsort(variances)[-n:]]
        self.feature_names_ = data.feature_names[self.selected_].astype(str).tolist()
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.selected_ is None or self.feature_names_ is None:
            raise RuntimeError("provider is not fitted")
        counts, order = self._aggregate(data, donors)
        X = self._logcpm(counts[:, self.selected_]).astype(np.float32)
        return DonorFeatures(X, order, self.feature_names_)


@dataclass
class ISGProvider(RepresentationProvider):
    genes: tuple[str, ...]
    vector: bool = False
    output_level: str = "donor"

    def __post_init__(self) -> None:
        self.indices_: np.ndarray | None = None
        self.names_: list[str] | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "ISGProvider":
        lookup = {str(g).upper(): i for i, g in enumerate(data.feature_names)}
        present = [(gene, lookup[gene.upper()]) for gene in self.genes if gene.upper() in lookup]
        if not present:
            raise ValueError("none of the requested ISG genes are present")
        self.names_ = [x[0] for x in present]
        self.indices_ = np.asarray([x[1] for x in present], dtype=int)
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.indices_ is None or self.names_ is None:
            raise RuntimeError("provider is not fitted")
        means, order = _donor_means(data, donors)
        values = means[:, self.indices_]
        if self.vector:
            return DonorFeatures(values.astype(np.float32), order, self.names_)
        scalar = values.mean(axis=1, keepdims=True).astype(np.float32)
        return DonorFeatures(scalar, order, ["isg_mean_score"])
