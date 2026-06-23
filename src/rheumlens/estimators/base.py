from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np


class DonorEstimator(ABC):
    @abstractmethod
    def fit(self, X: np.ndarray, y: np.ndarray, sample_weight: np.ndarray | None = None) -> "DonorEstimator":
        raise NotImplementedError

    @abstractmethod
    def predict_score(self, X: np.ndarray) -> np.ndarray:
        raise NotImplementedError
