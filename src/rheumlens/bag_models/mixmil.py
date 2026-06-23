from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.mixture import GaussianMixture

from rheumlens.bag_models.base import DonorBagModel
from rheumlens.estimators.linear import LogisticL2Estimator


@dataclass
class GaussianMixMIL(DonorBagModel):
    """A lightweight mixture-based patient model.

    It fits a shared cell-state GMM on training cells, then classifies donors using component
    occupancy and responsibility-weighted state means. This is an internal reproducible
    approximation, not a claim of source-code equivalence to an external MixMIL package.
    """

    n_components: int = 16
    max_fit_cells: int = 200_000
    seed: int = 0

    def __post_init__(self) -> None:
        self.gmm_: GaussianMixture | None = None
        self.estimator_: LogisticL2Estimator | None = None

    def _features(self, bag: np.ndarray) -> np.ndarray:
        assert self.gmm_ is not None
        resp = self.gmm_.predict_proba(bag)
        occupancy = resp.mean(axis=0)
        weighted = (resp.T @ bag) / np.maximum(resp.sum(axis=0)[:, None], 1e-12)
        return np.concatenate([occupancy, weighted.ravel()])

    def fit(self, bags: list[np.ndarray], y: np.ndarray) -> "GaussianMixMIL":
        rng = np.random.default_rng(self.seed)
        cells = np.vstack(bags)
        if cells.shape[0] > self.max_fit_cells:
            cells = cells[rng.choice(cells.shape[0], self.max_fit_cells, replace=False)]
        self.gmm_ = GaussianMixture(
            n_components=min(self.n_components, len(cells)),
            covariance_type="diag",
            random_state=self.seed,
            reg_covar=1e-5,
            max_iter=300,
        ).fit(cells)
        X = np.vstack([self._features(bag) for bag in bags])
        self.estimator_ = LogisticL2Estimator(random_state=self.seed).fit(X, y)
        return self

    def predict_score(self, bags: list[np.ndarray]) -> np.ndarray:
        if self.estimator_ is None:
            raise RuntimeError("model is not fitted")
        X = np.vstack([self._features(bag) for bag in bags])
        return self.estimator_.predict_score(X)
