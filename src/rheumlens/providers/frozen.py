from __future__ import annotations

from dataclasses import dataclass

from rheumlens.providers.base import RepresentationProvider
from rheumlens.types import CellDataset, FitContext


@dataclass
class FrozenCellProvider(RepresentationProvider):
    """Pass-through provider for precomputed cell embeddings or fixed expression matrices."""

    output_level: str = "cell"

    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "FrozenCellProvider":
        return self

    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> CellDataset:
        return data.subset_donors(donors)
