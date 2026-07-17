from __future__ import annotations

import numpy as np

from rheumlens.types import CellDataset


def _select_features(dataset: CellDataset, indices: np.ndarray, names: np.ndarray) -> CellDataset:
    return CellDataset(
        X=dataset.X[:, indices],
        cell_ids=dataset.cell_ids.copy(),
        donor_ids=dataset.donor_ids.copy(),
        y=dataset.y.copy(),
        feature_names=names.copy(),
        cell_types=None if dataset.cell_types is None else dataset.cell_types.copy(),
        cohorts=None if dataset.cohorts is None else dataset.cohorts.copy(),
        metadata=None if dataset.metadata is None else dataset.metadata.copy(),
        name=dataset.name,
    )


def align_features(
    source: CellDataset, target: CellDataset
) -> tuple[CellDataset, CellDataset]:
    """Align two datasets to shared features while preserving source order.

    Feature selection is defined entirely by source names. Target labels and
    values do not influence the shared feature space.
    """

    source_names = source.feature_names.astype(str)
    target_names = target.feature_names.astype(str)
    if len(set(source_names)) != len(source_names):
        raise ValueError("source feature_names must be unique")
    if len(set(target_names)) != len(target_names):
        raise ValueError("target feature_names must be unique")

    target_lookup = {name: index for index, name in enumerate(target_names)}
    source_indices = np.asarray(
        [index for index, name in enumerate(source_names) if name in target_lookup], dtype=int
    )
    if source_indices.size == 0:
        raise ValueError("source and target have no shared features")
    shared_names = source_names[source_indices]
    target_indices = np.asarray([target_lookup[name] for name in shared_names], dtype=int)
    return (
        _select_features(source, source_indices, shared_names),
        _select_features(target, target_indices, shared_names),
    )
