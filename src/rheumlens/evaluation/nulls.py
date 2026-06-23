from __future__ import annotations

import numpy as np
from sklearn.covariance import LedoitWolf

from rheumlens.types import CellDataset
from rheumlens.utils.arrays import as_dense, grouped_indices


def mean_sufficient_surrogate(data: CellDataset, seed: int = 0) -> CellDataset:
    """Preserve donor/cell-type counts and means, remove label-linked higher moments.

    A label-agnostic pooled residual covariance is estimated for each cell type, then Gaussian
    residuals are added to each donor/type mean. This is a calibration surrogate, not a
    biological generative model.
    """

    rng = np.random.default_rng(seed)
    X = as_dense(data.X, dtype=np.float64)
    cell_types = data.cell_types if data.cell_types is not None else np.repeat("global", len(data.cell_ids))
    generated = np.zeros_like(X)
    for cell_type in np.unique(cell_types.astype(str)):
        type_idx = np.flatnonzero(cell_types.astype(str) == cell_type)
        donor_groups = grouped_indices(data.donor_ids[type_idx])
        residuals = []
        donor_means: dict[str, np.ndarray] = {}
        for donor, local_idx in donor_groups.items():
            idx = type_idx[local_idx]
            mean = X[idx].mean(axis=0)
            donor_means[donor] = mean
            residuals.append(X[idx] - mean)
        residual_matrix = np.vstack(residuals)
        covariance = LedoitWolf().fit(residual_matrix).covariance_ if len(residual_matrix) > 1 else np.eye(X.shape[1])
        for donor, local_idx in donor_groups.items():
            idx = type_idx[local_idx]
            generated[idx] = rng.multivariate_normal(donor_means[donor], covariance, size=len(idx), method="svd")
    return CellDataset(
        X=generated.astype(np.float32),
        cell_ids=data.cell_ids,
        donor_ids=data.donor_ids,
        y=data.y,
        feature_names=data.feature_names,
        cell_types=data.cell_types,
        cohorts=data.cohorts,
        metadata=data.metadata,
        name=f"{data.name}_mean_sufficient_null",
    )


def moment2_sufficient_surrogate(data: CellDataset, seed: int = 0) -> CellDataset:
    """Preserve each donor's mean and shrinkage covariance; destroy non-Gaussian shape."""

    rng = np.random.default_rng(seed)
    X = as_dense(data.X, dtype=np.float64)
    generated = np.zeros_like(X)
    for donor, idx in grouped_indices(data.donor_ids).items():
        block = X[idx]
        mean = block.mean(axis=0)
        covariance = LedoitWolf().fit(block).covariance_ if len(block) > 1 else np.eye(X.shape[1]) * 1e-6
        generated[idx] = rng.multivariate_normal(mean, covariance, size=len(idx), method="svd")
    return CellDataset(
        X=generated.astype(np.float32),
        cell_ids=data.cell_ids,
        donor_ids=data.donor_ids,
        y=data.y,
        feature_names=data.feature_names,
        cell_types=data.cell_types,
        cohorts=data.cohorts,
        metadata=data.metadata,
        name=f"{data.name}_moment2_null",
    )
