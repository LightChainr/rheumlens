from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

from rheumlens.types import CellDataset


_FORMAT_VERSION = 1


def _portable_array(values: np.ndarray) -> np.ndarray:
    """Return an array that NumPy can load without enabling pickle."""

    array = np.asarray(values)
    if array.dtype.kind == "O":
        return array.astype(str)
    return array


def save_npz_dataset(path: str | Path, dataset: CellDataset) -> None:
    """Save a :class:`CellDataset` as a compressed, pickle-free NPZ archive."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "format_version": np.asarray(_FORMAT_VERSION, dtype=np.int16),
        "cell_ids": _portable_array(dataset.cell_ids),
        "donor_ids": _portable_array(dataset.donor_ids),
        "y": _portable_array(dataset.y),
        "feature_names": _portable_array(dataset.feature_names),
        "name": np.asarray(dataset.name),
        "has_cell_types": np.asarray(dataset.cell_types is not None),
        "has_cohorts": np.asarray(dataset.cohorts is not None),
        "has_metadata": np.asarray(dataset.metadata is not None),
    }
    if dataset.cell_types is not None:
        payload["cell_types"] = _portable_array(dataset.cell_types)
    if dataset.cohorts is not None:
        payload["cohorts"] = _portable_array(dataset.cohorts)
    if dataset.metadata is not None:
        if len(dataset.metadata) != dataset.X.shape[0]:
            raise ValueError("metadata row count must match the number of cells")
        payload["metadata_json"] = np.asarray(
            dataset.metadata.to_json(orient="table", index=False)
        )

    if sparse.issparse(dataset.X):
        matrix = sparse.csr_matrix(dataset.X)
        payload.update(
            {
                "matrix_format": np.asarray("csr"),
                "X_data": matrix.data,
                "X_indices": matrix.indices,
                "X_indptr": matrix.indptr,
                "X_shape": np.asarray(matrix.shape, dtype=np.int64),
            }
        )
    else:
        payload["matrix_format"] = np.asarray("dense")
        payload["X"] = np.asarray(dataset.X)
    np.savez_compressed(destination, **payload)


def _scalar(archive: np.lib.npyio.NpzFile, key: str, default: Any = None) -> Any:
    if key not in archive.files:
        return default
    return archive[key].item()


def load_npz_dataset(path: str | Path) -> CellDataset:
    """Load a dataset written by :func:`save_npz_dataset`.

    Archives from the early project layout that contain a plain ``X`` key are
    accepted as dense datasets even when ``matrix_format`` is absent.
    """

    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(source)
    with np.load(source, allow_pickle=False) as archive:
        required = {"cell_ids", "donor_ids", "y", "feature_names"}
        missing = sorted(required - set(archive.files))
        if missing:
            raise ValueError(f"{source} is missing required arrays: {', '.join(missing)}")

        matrix_format = str(_scalar(archive, "matrix_format", "dense"))
        if matrix_format == "csr":
            sparse_keys = {"X_data", "X_indices", "X_indptr", "X_shape"}
            missing_sparse = sorted(sparse_keys - set(archive.files))
            if missing_sparse:
                raise ValueError(
                    f"{source} is missing CSR arrays: {', '.join(missing_sparse)}"
                )
            shape = tuple(int(value) for value in archive["X_shape"])
            X = sparse.csr_matrix(
                (archive["X_data"], archive["X_indices"], archive["X_indptr"]),
                shape=shape,
            )
        elif matrix_format == "dense" and "X" in archive.files:
            X = archive["X"]
        else:
            raise ValueError(f"unsupported or incomplete matrix format: {matrix_format}")

        metadata = None
        if bool(_scalar(archive, "has_metadata", "metadata_json" in archive.files)):
            metadata_json = str(_scalar(archive, "metadata_json"))
            metadata = pd.read_json(metadata_json, orient="table")

        return CellDataset(
            X=X,
            cell_ids=archive["cell_ids"],
            donor_ids=archive["donor_ids"],
            y=archive["y"],
            feature_names=archive["feature_names"],
            cell_types=(
                archive["cell_types"]
                if bool(_scalar(archive, "has_cell_types", "cell_types" in archive.files))
                else None
            ),
            cohorts=(
                archive["cohorts"]
                if bool(_scalar(archive, "has_cohorts", "cohorts" in archive.files))
                else None
            ),
            metadata=metadata,
            name=str(_scalar(archive, "name", source.stem)),
        )
