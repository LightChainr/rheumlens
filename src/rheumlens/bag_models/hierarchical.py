from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rheumlens.bag_models.base import DonorBagModel
from rheumlens.estimators.linear import LogisticL2Estimator


@dataclass
class HierarchicalCellTypeMIL(DonorBagModel):
    """Low-capacity hierarchical MIL summary using cell-type conditioned moments.

    The bag entries are expected to include the cell-type label in a parallel list supplied
    through `fit_with_types`; the plain `fit` method is intentionally unsupported.
    """

    seed: int = 0

    def __post_init__(self) -> None:
        self.types_: list[str] = []
        self.estimator_: LogisticL2Estimator | None = None

    @staticmethod
    def _summary(bag: np.ndarray, types: np.ndarray, vocab: list[str]) -> np.ndarray:
        parts = []
        for value in vocab:
            block = bag[types.astype(str) == value]
            if len(block) == 0:
                parts.extend([np.zeros(bag.shape[1]), np.zeros(bag.shape[1]), np.asarray([0.0])])
            else:
                parts.extend([block.mean(axis=0), block.var(axis=0), np.asarray([len(block) / len(bag)])])
        return np.concatenate(parts)

    def fit_with_types(self, bags: list[np.ndarray], bag_types: list[np.ndarray], y: np.ndarray) -> "HierarchicalCellTypeMIL":
        self.types_ = sorted({str(x) for values in bag_types for x in values})
        X = np.vstack([self._summary(bag, types, self.types_) for bag, types in zip(bags, bag_types, strict=True)])
        self.estimator_ = LogisticL2Estimator(random_state=self.seed).fit(X, y)
        return self

    def fit(self, bags: list[np.ndarray], y: np.ndarray) -> "HierarchicalCellTypeMIL":
        raise ValueError("use fit_with_types for HierarchicalCellTypeMIL")

    def predict_score_with_types(self, bags: list[np.ndarray], bag_types: list[np.ndarray]) -> np.ndarray:
        if self.estimator_ is None:
            raise RuntimeError("model is not fitted")
        X = np.vstack([self._summary(bag, types, self.types_) for bag, types in zip(bags, bag_types, strict=True)])
        return self.estimator_.predict_score(X)

    def predict_score(self, bags: list[np.ndarray]) -> np.ndarray:
        raise ValueError("use predict_score_with_types for HierarchicalCellTypeMIL")
