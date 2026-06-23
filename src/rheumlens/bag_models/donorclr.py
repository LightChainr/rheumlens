from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rheumlens.bag_models.base import DonorBagModel
from rheumlens.bag_models.torch_models import _require_torch
from rheumlens.estimators.linear import LogisticL2Estimator


@dataclass
class DonorCLR(DonorBagModel):
    """Low-capacity contrastive donor encoder trained from independent cell subsamples."""

    hidden: int = 64
    cells_per_view: int = 256
    epochs: int = 100
    temperature: float = 0.2
    learning_rate: float = 1e-3
    seed: int = 0
    device: str = "auto"

    def __post_init__(self) -> None:
        self.encoder_ = None
        self.classifier_: LogisticL2Estimator | None = None
        self.device_ = None

    def _summary(self, bag: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        n = min(self.cells_per_view, len(bag))
        block = bag[rng.choice(len(bag), n, replace=len(bag) < n)]
        return np.concatenate([block.mean(axis=0), block.var(axis=0)])

    def fit(self, bags: list[np.ndarray], y: np.ndarray) -> "DonorCLR":
        torch, nn, F = _require_torch()
        self.device_ = torch.device("cuda" if self.device == "auto" and torch.cuda.is_available() else "cpu")
        torch.manual_seed(self.seed)
        rng = np.random.default_rng(self.seed)
        input_dim = bags[0].shape[1] * 2
        self.encoder_ = nn.Sequential(nn.Linear(input_dim, self.hidden), nn.ReLU(), nn.Linear(self.hidden, self.hidden)).to(self.device_)
        optimizer = torch.optim.AdamW(self.encoder_.parameters(), lr=self.learning_rate)
        for _ in range(self.epochs):
            view1 = torch.from_numpy(np.vstack([self._summary(b, rng) for b in bags]).astype(np.float32)).to(self.device_)
            view2 = torch.from_numpy(np.vstack([self._summary(b, rng) for b in bags]).astype(np.float32)).to(self.device_)
            z1 = F.normalize(self.encoder_(view1), dim=1)
            z2 = F.normalize(self.encoder_(view2), dim=1)
            logits = z1 @ z2.T / self.temperature
            labels = torch.arange(len(bags), device=self.device_)
            loss = 0.5 * (F.cross_entropy(logits, labels) + F.cross_entropy(logits.T, labels))
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        features = self._encode(bags, rng)
        self.classifier_ = LogisticL2Estimator(random_state=self.seed).fit(features, y)
        return self

    def _encode(self, bags: list[np.ndarray], rng: np.random.Generator) -> np.ndarray:
        torch, _, _ = _require_torch()
        summaries = np.vstack([self._summary(b, rng) for b in bags]).astype(np.float32)
        with torch.no_grad():
            return self.encoder_(torch.from_numpy(summaries).to(self.device_)).cpu().numpy()

    def predict_score(self, bags: list[np.ndarray]) -> np.ndarray:
        if self.classifier_ is None:
            raise RuntimeError("model is not fitted")
        features = self._encode(bags, np.random.default_rng(self.seed + 10000))
        return self.classifier_.predict_score(features)
