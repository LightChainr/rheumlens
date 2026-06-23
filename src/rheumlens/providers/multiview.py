from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rheumlens.providers.base import RepresentationProvider
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import grouped_indices


@dataclass
class CellTypePseudobulkProvider(RepresentationProvider):
    n_features_per_type: int = 250
    min_cells: int = 5
    output_level: str = "donor"

    def __post_init__(self) -> None:
        self.types_: list[str] = []
        self.selected_: dict[str, np.ndarray] = {}
        self.names_: list[str] = []

    @staticmethod
    def _logcpm(counts: np.ndarray) -> np.ndarray:
        return np.log1p(counts / np.maximum(counts.sum(axis=1, keepdims=True), 1.0) * 1e6)

    def _matrix(self, data: CellDataset, donors: list[str], cell_type: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if data.cell_types is None:
            raise ValueError("cell-type pseudobulk requires cell_types")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rows, present = [], []
        for donor in order:
            idx = groups[donor]
            idx = idx[subset.cell_types[idx].astype(str) == cell_type]
            present.append(len(idx) >= self.min_cells)
            if len(idx) == 0:
                rows.append(np.zeros(subset.X.shape[1]))
            else:
                rows.append(np.asarray(subset.X[idx].sum(axis=0)).ravel())
        return np.vstack(rows), order, np.asarray(present, dtype=float)

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "CellTypePseudobulkProvider":
        if data.cell_types is None:
            raise ValueError("cell-type pseudobulk requires cell_types")
        self.types_ = sorted(map(str, np.unique(data.subset_donors(train_donors).cell_types)))
        self.selected_.clear()
        self.names_.clear()
        for cell_type in self.types_:
            counts, _, present = self._matrix(data, train_donors, cell_type)
            transformed = self._logcpm(counts)
            variance = transformed.var(axis=0, ddof=1)
            n = min(self.n_features_per_type, transformed.shape[1])
            selected = np.argsort(variance)[-n:]
            self.selected_[cell_type] = selected
            self.names_.extend(f"{cell_type}__{data.feature_names[i]}" for i in selected)
            self.names_.append(f"{cell_type}__present")
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if not self.types_:
            raise RuntimeError("provider is not fitted")
        blocks = []
        order = None
        for cell_type in self.types_:
            counts, current_order, present = self._matrix(data, donors, cell_type)
            order = current_order if order is None else order
            transformed = self._logcpm(counts)[:, self.selected_[cell_type]]
            blocks.append(np.hstack([transformed, present[:, None]]))
        return DonorFeatures(np.hstack(blocks).astype(np.float32), order, self.names_)


@dataclass
class ScFeaturesStyleProvider(RepresentationProvider):
    n_genes_per_type: int = 50
    output_level: str = "donor"

    def __post_init__(self) -> None:
        self.pb_ = CellTypePseudobulkProvider(n_features_per_type=self.n_genes_per_type)
        self.types_: list[str] = []

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "ScFeaturesStyleProvider":
        self.pb_.fit(data, train_donors, context)
        self.types_ = self.pb_.types_
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        pb = self.pb_.transform(data, donors, context)
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        composition = np.zeros((len(pb.donor_ids), len(self.types_)), dtype=float)
        lookup = {x: i for i, x in enumerate(self.types_)}
        for row, donor in enumerate(pb.donor_ids):
            idx = groups[str(donor)]
            values, counts = np.unique(subset.cell_types[idx].astype(str), return_counts=True)
            for value, count in zip(values, counts, strict=True):
                if value in lookup:
                    composition[row, lookup[value]] = count / len(idx)
        entropy = -np.sum(composition * np.log(composition + 1e-12), axis=1, keepdims=True)
        X = np.hstack([composition, entropy, pb.X])
        names = [f"composition_{x}" for x in self.types_] + ["composition_entropy"] + pb.feature_names
        return DonorFeatures(X.astype(np.float32), pb.donor_ids, names)
