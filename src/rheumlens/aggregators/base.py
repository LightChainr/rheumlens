from __future__ import annotations

from abc import ABC, abstractmethod

from rheumlens.types import CellDataset, DonorFeatures, FitContext


class DonorAggregator(ABC):
    @abstractmethod
    def fit(self, data: CellDataset, train_donors: list[str], context: FitContext) -> "DonorAggregator":
        raise NotImplementedError

    @abstractmethod
    def transform(self, data: CellDataset, donors: list[str], context: FitContext) -> DonorFeatures:
        raise NotImplementedError
