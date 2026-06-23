from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import as_dense, grouped_indices, seed_sequence


def _median_bandwidth(X: np.ndarray, max_points: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    if X.shape[0] > max_points:
        X = X[rng.choice(X.shape[0], max_points, replace=False)]
    distances = pairwise_distances(X, metric="euclidean")
    upper = distances[np.triu_indices_from(distances, k=1)]
    positive = upper[upper > 0]
    if positive.size == 0:
        return 1.0
    return float(np.median(positive))


class RFFMap:
    def __init__(self, dim: int, bandwidth: float, seed: int) -> None:
        self.dim = dim
        self.bandwidth = max(float(bandwidth), 1e-8)
        self.seed = seed
        self.W: np.ndarray | None = None
        self.b: np.ndarray | None = None

    def fit(self, input_dim: int) -> "RFFMap":
        rng = np.random.default_rng(self.seed)
        self.W = rng.normal(0.0, 1.0 / self.bandwidth, size=(input_dim, self.dim))
        self.b = rng.uniform(0.0, 2 * np.pi, size=self.dim)
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self.W is None or self.b is None:
            raise RuntimeError("RFFMap is not fitted")
        return np.sqrt(2.0 / self.dim) * np.cos(X @ self.W + self.b)


@dataclass
class MultiScaleRFFKMEAggregator(DonorAggregator):
    rff_dim: int = 256
    scales: tuple[float, ...] = (0.5, 1.0, 2.0)
    seed: int = 0
    max_bandwidth_points: int = 4000
    max_diagnostic_points: int = 512
    standardize_cells: bool = True
    include_linear_mean: bool = False

    def __post_init__(self) -> None:
        self.scaler_: StandardScaler | None = None
        self.bandwidth_: float | None = None
        self.maps_: list[RFFMap] = []
        self.train_reference_: np.ndarray | None = None
        self.diagnostics_: dict[str, object] = {}

    def _prepare_fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> np.ndarray:
        X = as_dense(data.subset_donors(train_donors).X)
        if self.standardize_cells:
            self.scaler_ = StandardScaler().fit(X)
            X = self.scaler_.transform(X)
        self.bandwidth_ = _median_bandwidth(X, self.max_bandwidth_points, context.random_state)
        seeds = seed_sequence(self.seed, len(self.scales))
        self.maps_ = [
            RFFMap(self.rff_dim, self.bandwidth_ * scale, seed).fit(X.shape[1])
            for scale, seed in zip(self.scales, seeds, strict=True)
        ]
        rng = np.random.default_rng(context.random_state + 9173)
        if len(X) > self.max_diagnostic_points:
            sample = X[rng.choice(len(X), self.max_diagnostic_points, replace=False)]
        else:
            sample = X
        gram = np.exp(-pairwise_distances(sample, squared=True) / (2 * self.bandwidth_**2))
        eig = np.linalg.eigvalsh((gram + gram.T) / 2)
        effective_rank = float((eig.sum() ** 2) / np.maximum(np.square(eig).sum(), 1e-12))
        offdiag = gram[np.triu_indices_from(gram, k=1)]
        self.diagnostics_ = {
            "bandwidth": self.bandwidth_,
            "effective_rank": effective_rank,
            "offdiag_q01_q50_q99": np.quantile(offdiag, [0.01, 0.5, 0.99]).tolist() if offdiag.size else [],
        }
        return X

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "MultiScaleRFFKMEAggregator":
        self._prepare_fit(data, train_donors, context)
        reference = self.transform(data, train_donors, context)
        self.train_reference_ = reference.X.mean(axis=0)
        return self

    def _transform_cells(self, X: np.ndarray) -> np.ndarray:
        if self.scaler_ is not None:
            X = self.scaler_.transform(X)
        return X

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if not self.maps_:
            raise RuntimeError("aggregator is not fitted")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rows = []
        for donor in order:
            X = self._transform_cells(as_dense(subset.X[groups[donor]]))
            parts: list[np.ndarray] = []
            if self.include_linear_mean:
                parts.append(X.mean(axis=0))
            parts.extend(mapping.transform(X).mean(axis=0) for mapping in self.maps_)
            rows.append(np.concatenate(parts))
        names = []
        if self.include_linear_mean:
            names.extend(f"linear_mean_{x}" for x in subset.feature_names)
        for scale in self.scales:
            names.extend(f"rff_s{scale:g}_{i}" for i in range(self.rff_dim))
        X_donor = np.vstack(rows).astype(np.float32)
        diagnostics = dict(self.diagnostics_)
        norms = np.linalg.norm(X_donor, axis=1)
        counts = np.asarray([len(groups[d]) for d in order], dtype=float)
        diagnostics["feature_norm_cellcount_corr"] = (
            float(np.corrcoef(norms, counts)[0, 1])
            if len(order) > 2 and np.std(norms) > 1e-12 and np.std(counts) > 1e-12
            else np.nan
        )
        return DonorFeatures(X_donor, order, names, diagnostics=diagnostics)


@dataclass
class RobustMedianOfMeansKMEAggregator(MultiScaleRFFKMEAggregator):
    n_blocks: int = 5

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if not self.maps_:
            raise RuntimeError("aggregator is not fitted")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rng = np.random.default_rng(context.random_state)
        rows = []
        for donor in order:
            X = self._transform_cells(as_dense(subset.X[groups[donor]]))
            permutation = rng.permutation(X.shape[0])
            blocks = np.array_split(permutation, min(self.n_blocks, X.shape[0]))
            parts = []
            if self.include_linear_mean:
                parts.append(np.median(np.vstack([X[idx].mean(axis=0) for idx in blocks]), axis=0))
            for mapping in self.maps_:
                block_means = np.vstack([mapping.transform(X[idx]).mean(axis=0) for idx in blocks])
                parts.append(np.median(block_means, axis=0))
            rows.append(np.concatenate(parts))
        base = super().transform(data, donors, context)
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, base.feature_names, diagnostics=base.diagnostics)


@dataclass
class MMDReferenceAggregator(MultiScaleRFFKMEAggregator):
    cell_type_conditional: bool = False

    def __post_init__(self) -> None:
        super().__post_init__()
        self.control_references_: dict[str, np.ndarray] = {}

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "MMDReferenceAggregator":
        self._prepare_fit(data, train_donors, context)
        subset = data.subset_donors(train_donors)
        control_mask = np.asarray([context.y_by_donor[str(d)] == 0 for d in subset.donor_ids])
        if not control_mask.any():
            raise ValueError("MMD reference requires training controls")
        types = ["global"]
        if self.cell_type_conditional:
            if subset.cell_types is None:
                raise ValueError("cell-type conditional MMD requires cell_types")
            types += sorted(map(str, np.unique(subset.cell_types)))
        for cell_type in types:
            mask = control_mask.copy()
            if cell_type != "global":
                mask &= subset.cell_types.astype(str) == cell_type
            X = self._transform_cells(as_dense(subset.X[mask]))
            if X.shape[0] == 0:
                continue
            self.control_references_[cell_type] = np.concatenate([m.transform(X).mean(axis=0) for m in self.maps_])
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if not self.control_references_:
            raise RuntimeError("aggregator is not fitted")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        keys = list(self.control_references_)
        rows = []
        for donor in order:
            donor_idx = groups[donor]
            values = []
            for key in keys:
                idx = donor_idx
                if key != "global":
                    if subset.cell_types is None:
                        values.append(np.nan)
                        continue
                    idx = donor_idx[subset.cell_types[donor_idx].astype(str) == key]
                if len(idx) == 0:
                    values.append(np.nan)
                    continue
                X = self._transform_cells(as_dense(subset.X[idx]))
                rep = np.concatenate([m.transform(X).mean(axis=0) for m in self.maps_])
                values.append(float(np.square(rep - self.control_references_[key]).sum()))
            rows.append(values)
        matrix = np.asarray(rows, dtype=np.float32)
        col_medians = np.nanmedian(matrix, axis=0)
        inds = np.where(np.isnan(matrix))
        matrix[inds] = col_medians[inds[1]]
        return DonorFeatures(matrix, order, [f"mmd2_to_control_{k}" for k in keys], diagnostics=self.diagnostics_)


@dataclass
class CCKMEUAggregator(MultiScaleRFFKMEAggregator):
    subsample_size: int = 500
    n_subsamples: int = 20
    tau: float = 1.0
    include_uncertainty_features: bool = True

    def __post_init__(self) -> None:
        super().__post_init__()
        self.kme_reference_: np.ndarray | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "CCKMEUAggregator":
        self._prepare_fit(data, train_donors, context)
        raw = self._subsample_features(data, train_donors, context)
        kme_start = data.X.shape[1] if self.include_linear_mean else 0
        kme_end = kme_start + len(self.scales) * self.rff_dim
        self.kme_reference_ = raw.X[:, kme_start:kme_end].mean(axis=0)
        return self

    def _subsample_features(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rng = np.random.default_rng(context.random_state)
        rows, uncertainties = [], []
        for donor in order:
            X = self._transform_cells(as_dense(subset.X[groups[donor]]))
            reps = []
            n_draw = min(self.subsample_size, X.shape[0])
            for _ in range(self.n_subsamples):
                idx = rng.choice(X.shape[0], n_draw, replace=X.shape[0] < self.subsample_size)
                block = X[idx]
                parts = [block.mean(axis=0)] if self.include_linear_mean else []
                parts.extend(mapping.transform(block).mean(axis=0) for mapping in self.maps_)
                reps.append(np.concatenate(parts))
            reps_arr = np.vstack(reps)
            rows.append(reps_arr.mean(axis=0))
            uncertainties.append(reps_arr.var(axis=0, ddof=1) if len(reps) > 1 else np.zeros(reps_arr.shape[1]))
        matrix = np.vstack(rows)
        uncertainty = np.vstack(uncertainties)
        names = []
        if self.include_linear_mean:
            names.extend(f"linear_mean_{x}" for x in subset.feature_names)
        for scale in self.scales:
            names.extend(f"rff_s{scale:g}_{i}" for i in range(self.rff_dim))
        return DonorFeatures(matrix.astype(np.float32), order, names, uncertainty=uncertainty.astype(np.float32))

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.kme_reference_ is None:
            raise RuntimeError("aggregator is not fitted")
        raw = self._subsample_features(data, donors, context)
        matrix = raw.X.copy()
        uncertainty = raw.uncertainty
        assert uncertainty is not None
        linear_dim = data.X.shape[1] if self.include_linear_mean else 0
        kme_slice = slice(linear_dim, linear_dim + len(self.scales) * self.rff_dim)
        mean_u = uncertainty[:, kme_slice].mean(axis=1)
        weights = 1.0 / (1.0 + self.tau * mean_u)
        matrix[:, kme_slice] = (
            weights[:, None] * matrix[:, kme_slice]
            + (1.0 - weights[:, None]) * self.kme_reference_[None, :]
        )
        names = raw.feature_names.copy()
        if self.include_uncertainty_features:
            matrix = np.hstack([matrix, np.log1p(uncertainty)])
            names += [f"uncertainty_{x}" for x in raw.feature_names]
        diagnostics = dict(self.diagnostics_)
        diagnostics["reliability_weights"] = weights.tolist()
        return DonorFeatures(matrix.astype(np.float32), raw.donor_ids, names, uncertainty=uncertainty, diagnostics=diagnostics)
