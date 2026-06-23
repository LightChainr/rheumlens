from __future__ import annotations

import copy
from dataclasses import dataclass

import numpy as np

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.types import CellDataset, DonorFeatures, FitContext
from rheumlens.utils.arrays import grouped_indices


@dataclass
class UDERWrapper(DonorAggregator):
    """Repeated cell-subsampling wrapper for any fixed donor aggregator."""

    base: DonorAggregator
    subsample_size: int = 500
    n_subsamples: int = 20
    append_log_variance: bool = True
    seed: int = 0

    def __post_init__(self) -> None:
        self.base_: DonorAggregator | None = None

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "UDERWrapper":
        self.base_ = copy.deepcopy(self.base).fit(data, train_donors, context)
        return self

    def _sample_dataset(self, data: CellDataset, donors: list[str], rng: np.random.Generator) -> CellDataset:
        subset = data.subset_donors(donors)
        groups = grouped_indices(subset.donor_ids)
        keep = []
        for donor, idx in groups.items():
            n = min(self.subsample_size, len(idx))
            # CellDataset uses a boolean subset, so duplicated indices would collapse.
            # For donors below the target count we therefore retain all available cells;
            # uncertainty is estimated from donors that can actually be resampled.
            keep.extend(rng.choice(idx, n, replace=False).tolist())
        mask = np.zeros(subset.X.shape[0], dtype=bool)
        mask[np.asarray(keep, dtype=int)] = True
        return subset.subset_cells(mask)

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        if self.base_ is None:
            raise RuntimeError("wrapper is not fitted")
        rng = np.random.default_rng(self.seed + context.random_state)
        reps = []
        donor_order = None
        names = None
        for _ in range(self.n_subsamples):
            sampled = self._sample_dataset(data, donors, rng)
            result = self.base_.transform(sampled, donors, context)
            donor_order = result.donor_ids
            names = result.feature_names
            reps.append(result.X)
        stack = np.stack(reps, axis=0)
        mean = stack.mean(axis=0)
        variance = stack.var(axis=0, ddof=1) if len(reps) > 1 else np.zeros_like(mean)
        output_names = list(names or [])
        if self.append_log_variance:
            mean = np.hstack([mean, np.log1p(variance)])
            output_names += [f"uder_variance_{x}" for x in output_names]
        return DonorFeatures(mean.astype(np.float32), donor_order, output_names, uncertainty=variance.astype(np.float32))
