from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

ArrayLike2D = np.ndarray | sparse.spmatrix


@dataclass
class CellDataset:
    """A cell-aligned matrix plus immutable donor metadata.

    `y` is stored at cell level for convenient slicing, but validation enforces a single label
    per donor. Cells are never treated as independent supervised disease samples.
    """

    X: ArrayLike2D
    cell_ids: np.ndarray
    donor_ids: np.ndarray
    y: np.ndarray
    feature_names: np.ndarray
    cell_types: np.ndarray | None = None
    cohorts: np.ndarray | None = None
    metadata: pd.DataFrame | None = None
    name: str = "dataset"

    def __post_init__(self) -> None:
        n = self.X.shape[0]
        for label, value in {
            "cell_ids": self.cell_ids,
            "donor_ids": self.donor_ids,
            "y": self.y,
        }.items():
            if len(value) != n:
                raise ValueError(f"{label} length {len(value)} != number of cells {n}")
        if self.X.shape[1] != len(self.feature_names):
            raise ValueError("feature_names do not match X columns")
        if self.cell_types is not None and len(self.cell_types) != n:
            raise ValueError("cell_types length mismatch")
        if self.cohorts is not None and len(self.cohorts) != n:
            raise ValueError("cohorts length mismatch")
        if len(np.unique(self.cell_ids)) != n:
            raise ValueError("cell_ids must be unique")
        self.validate_donor_labels()

    def validate_donor_labels(self) -> None:
        frame = pd.DataFrame({"donor": self.donor_ids, "y": self.y})
        counts = frame.groupby("donor", observed=True)["y"].nunique()
        bad = counts[counts != 1]
        if not bad.empty:
            raise ValueError(f"donors with non-constant labels: {bad.index.tolist()[:10]}")

    @property
    def donors(self) -> np.ndarray:
        return pd.unique(self.donor_ids)

    def donor_label_map(self) -> dict[str, int]:
        frame = pd.DataFrame({"donor": self.donor_ids.astype(str), "y": self.y.astype(int)})
        return frame.drop_duplicates("donor").set_index("donor")["y"].to_dict()

    def subset_donors(self, donors: Sequence[str]) -> "CellDataset":
        donor_set = set(map(str, donors))
        mask = np.fromiter((str(x) in donor_set for x in self.donor_ids), dtype=bool)
        return self.subset_cells(mask)

    def subset_cells(self, mask: np.ndarray) -> "CellDataset":
        if mask.dtype != bool or len(mask) != self.X.shape[0]:
            raise ValueError("mask must be a cell-length boolean vector")
        meta = None
        if self.metadata is not None:
            meta = self.metadata.loc[mask].reset_index(drop=True)
        return CellDataset(
            X=self.X[mask],
            cell_ids=self.cell_ids[mask],
            donor_ids=self.donor_ids[mask],
            y=self.y[mask],
            feature_names=self.feature_names.copy(),
            cell_types=None if self.cell_types is None else self.cell_types[mask],
            cohorts=None if self.cohorts is None else self.cohorts[mask],
            metadata=meta,
            name=self.name,
        )

    def to_dense(self, dtype: np.dtype = np.float32) -> np.ndarray:
        if sparse.issparse(self.X):
            return self.X.toarray().astype(dtype, copy=False)
        return np.asarray(self.X, dtype=dtype)


@dataclass
class DonorFeatures:
    X: np.ndarray
    donor_ids: np.ndarray
    feature_names: list[str]
    uncertainty: np.ndarray | None = None
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.X.ndim != 2:
            raise ValueError("DonorFeatures.X must be 2D")
        if self.X.shape[0] != len(self.donor_ids):
            raise ValueError("donor_ids length mismatch")
        if self.X.shape[1] != len(self.feature_names):
            raise ValueError("feature_names length mismatch")
        if self.uncertainty is not None and self.uncertainty.shape[0] != self.X.shape[0]:
            raise ValueError("uncertainty row mismatch")


@dataclass(frozen=True)
class OuterFold:
    split_id: str
    fold: int
    train_donors: tuple[str, ...]
    test_donors: tuple[str, ...]


@dataclass
class FitContext:
    y_by_donor: Mapping[str, int]
    random_state: int = 0
    expression: CellDataset | None = None
    donor_covariates: pd.DataFrame | None = None
    extras: dict[str, Any] = field(default_factory=dict)
