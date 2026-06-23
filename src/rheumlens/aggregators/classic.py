from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.preprocessing import normalize

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import as_dense, grouped_indices


@dataclass
class MeanAggregator(DonorAggregator):
    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "MeanAggregator":
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        rows = [np.asarray(subset.X[groups[d]].mean(axis=0)).ravel() for d in order]
        return DonorFeatures(np.vstack(rows).astype(np.float32), order, subset.feature_names.astype(str).tolist())


@dataclass
class CompositionCLRAggregator(DonorAggregator):
    pseudocount: float = 0.5

    def __post_init__(self) -> None:
        self.types_: list[str] | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "CompositionCLRAggregator":
        if data.cell_types is None:
            raise ValueError("composition requires cell_types")
        subset = data.subset_donors(train_donors)
        self.types_ = sorted(map(str, np.unique(subset.cell_types)))
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.types_ is None:
            raise RuntimeError("aggregator is not fitted")
        if data.cell_types is None:
            raise ValueError("composition requires cell_types")
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        order = np.asarray(list(groups), dtype=str)
        matrix = np.zeros((len(order), len(self.types_)), dtype=np.float64)
        type_to_idx = {t: i for i, t in enumerate(self.types_)}
        for row, donor in enumerate(order):
            values, counts = np.unique(subset.cell_types[groups[donor]].astype(str), return_counts=True)
            for value, count in zip(values, counts, strict=True):
                if value in type_to_idx:
                    matrix[row, type_to_idx[value]] = count
        matrix += self.pseudocount
        proportions = matrix / matrix.sum(axis=1, keepdims=True)
        logp = np.log(proportions)
        clr = logp - logp.mean(axis=1, keepdims=True)
        return DonorFeatures(clr.astype(np.float32), order, [f"clr_{t}" for t in self.types_])
