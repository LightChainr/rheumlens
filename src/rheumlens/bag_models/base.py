from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class DonorBagModel(ABC):
    @abstractmethod
    def fit(self, bags: list[np.ndarray], y: np.ndarray) -> "DonorBagModel":
        raise NotImplementedError

    @abstractmethod
    def predict_score(self, bags: list[np.ndarray]) -> np.ndarray:
        raise NotImplementedError
