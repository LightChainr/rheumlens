from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Literal

from rheumlens.types import CellDataset, DonorFeatures, FitContext


class RepresentationProvider(ABC):
    output_level: Literal["cell", "donor"] = "cell"

    @abstractmethod
    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "RepresentationProvider":
        raise NotImplementedError

    @abstractmethod
    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> CellDataset | DonorFeatures:
        raise NotImplementedError
