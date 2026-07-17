from __future__ import annotations

from typing import Any

import numpy as np
from scipy import sparse

from rheumlens.types import CellDataset


def _matrix_values(dataset: CellDataset) -> np.ndarray:
    if sparse.issparse(dataset.X):
        return np.asarray(dataset.X.data)
    return np.asarray(dataset.X)


def _base_report(dataset: CellDataset) -> dict[str, Any]:
    values = _matrix_values(dataset)
    donor_labels = dataset.donor_label_map()
    labels = np.asarray(list(donor_labels.values()))
    unique_labels, label_counts = np.unique(labels, return_counts=True)
    return {
        "name": dataset.name,
        "n_cells": int(dataset.X.shape[0]),
        "n_features": int(dataset.X.shape[1]),
        "n_donors": int(len(donor_labels)),
        "label_counts": {
            str(label): int(count) for label, count in zip(unique_labels, label_counts)
        },
        "sparse": bool(sparse.issparse(dataset.X)),
        "finite": bool(np.isfinite(values).all()),
        "n_missing_values": int((~np.isfinite(values)).sum()),
    }


def _validate_common(dataset: CellDataset) -> dict[str, Any]:
    report = _base_report(dataset)
    if dataset.X.ndim != 2:
        raise ValueError("X must be two-dimensional")
    if dataset.X.shape[0] == 0 or dataset.X.shape[1] == 0:
        raise ValueError("X must contain at least one cell and one feature")
    if not report["finite"]:
        raise ValueError("X contains NaN or infinite values")
    if len(set(dataset.feature_names.astype(str))) != len(dataset.feature_names):
        raise ValueError("feature_names must be unique")
    if len(report["label_counts"]) < 2:
        report["warning"] = "dataset contains fewer than two donor-level classes"
    return report


def validate_embedding(dataset: CellDataset) -> dict[str, Any]:
    """Validate a finite cell-level expression or embedding matrix."""

    report = _validate_common(dataset)
    values = _matrix_values(dataset)
    report["dtype"] = str(dataset.X.dtype)
    report["value_min"] = float(values.min()) if values.size else 0.0
    report["value_max"] = float(values.max()) if values.size else 0.0
    return report


def validate_raw_counts(dataset: CellDataset) -> dict[str, Any]:
    """Validate a non-negative, integer-valued raw-count matrix."""

    report = _validate_common(dataset)
    values = _matrix_values(dataset)
    nonnegative = bool((values >= 0).all())
    integer_like = bool(np.allclose(values, np.rint(values), rtol=0.0, atol=1e-6))
    report.update(
        {
            "dtype": str(dataset.X.dtype),
            "nonnegative": nonnegative,
            "integer_like": integer_like,
            "total_counts": float(values.sum()),
        }
    )
    if not nonnegative:
        raise ValueError("raw counts contain negative values")
    if not integer_like:
        raise ValueError("raw counts are not integer-valued")
    return report
