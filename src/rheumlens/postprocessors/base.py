from __future__ import annotations

from abc import ABC, abstractmethod

from rheumlens.types import DonorFeatures, FitContext


class DonorPostprocessor(ABC):
    @abstractmethod
    def fit(self, features: DonorFeatures, context: FitContext) -> "DonorPostprocessor":
        raise NotImplementedError

    @abstractmethod
    def transform(self, features: DonorFeatures, context: FitContext) -> DonorFeatures:
        raise NotImplementedError
