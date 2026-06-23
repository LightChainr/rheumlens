from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rheumlens.estimators.linear import LogisticL2Estimator


@dataclass
class ReliabilityWeightedLogistic(LogisticL2Estimator):
    floor: float = 1e-3
    use_uncertainty_weights: bool = True

    @staticmethod
    def weights_from_uncertainty(uncertainty: np.ndarray) -> np.ndarray:
        score = np.mean(uncertainty, axis=1)
        weights = 1.0 / np.maximum(score, 1e-8)
        weights /= np.mean(weights)
        return np.clip(weights, 0.1, 10.0)
