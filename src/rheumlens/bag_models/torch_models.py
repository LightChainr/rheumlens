from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from rheumlens.bag_models.base import DonorBagModel


def _require_torch():
    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Install rheumlens[gpu] for bag models") from exc
    return torch, nn, F


class _DeepSetNet:
    @staticmethod
    def build(input_dim: int, hidden: int):
        torch, nn, _ = _require_torch()

        class Net(nn.Module):
            def __init__(self) -> None:
                super().__init__()
                self.phi = nn.Sequential(nn.Linear(input_dim, hidden), nn.ReLU(), nn.Linear(hidden, hidden), nn.ReLU())
                self.head = nn.Linear(hidden, 1)

            def forward(self, x, mask):
                z = self.phi(x)
                z = (z * mask.unsqueeze(-1)).sum(1) / mask.sum(1, keepdim=True).clamp_min(1.0)
                return self.head(z).squeeze(-1)

        return Net()


class _AttentionNet:
    @staticmethod
    def build(input_dim: int, hidden: int, gated: bool = True, topk_fraction: float | None = None):
        torch, nn, _ = _require_torch()

        class Net(nn.Module):
            def __init__(self) -> None:
                super().__init__()
                self.embed = nn.Sequential(nn.Linear(input_dim, hidden), nn.Tanh())
                self.v = nn.Linear(hidden, hidden)
                self.u = nn.Linear(hidden, hidden) if gated else None
                self.a = nn.Linear(hidden, 1)
                self.head = nn.Linear(hidden, 1)

            def forward(self, x, mask):
                h = self.embed(x)
                score = torch.tanh(self.v(h))
                if self.u is not None:
                    score = score * torch.sigmoid(self.u(h))
                score = self.a(score).squeeze(-1)
                score = score.masked_fill(~mask.bool(), float("-inf"))
                if topk_fraction is not None:
                    keep_mask = torch.zeros_like(mask, dtype=torch.bool)
                    for row in range(score.shape[0]):
                        n_valid = int(mask[row].sum().item())
                        k = max(1, int(np.ceil(topk_fraction * n_valid)))
                        idx = torch.topk(score[row], k=k).indices
                        keep_mask[row, idx] = True
                    score = score.masked_fill(~keep_mask, float("-inf"))
                weights = torch.softmax(score, dim=1)
                pooled = torch.sum(h * weights.unsqueeze(-1), dim=1)
                return self.head(pooled).squeeze(-1)

        return Net()


class _SetTransformerNet:
    @staticmethod
    def build(input_dim: int, hidden: int, n_heads: int):
        torch, nn, _ = _require_torch()

        class Net(nn.Module):
            def __init__(self) -> None:
                super().__init__()
                self.input = nn.Linear(input_dim, hidden)
                self.attn1 = nn.MultiheadAttention(hidden, n_heads, batch_first=True)
                self.attn2 = nn.MultiheadAttention(hidden, n_heads, batch_first=True)
                self.seed = nn.Parameter(torch.randn(1, 1, hidden) * 0.02)
                self.pool = nn.MultiheadAttention(hidden, n_heads, batch_first=True)
                self.head = nn.Linear(hidden, 1)

            def forward(self, x, mask):
                h = self.input(x)
                key_padding = ~mask.bool()
                z, _ = self.attn1(h, h, h, key_padding_mask=key_padding)
                h = h + z
                z, _ = self.attn2(h, h, h, key_padding_mask=key_padding)
                h = h + z
                query = self.seed.expand(h.shape[0], -1, -1)
                pooled, _ = self.pool(query, h, h, key_padding_mask=key_padding)
                return self.head(pooled[:, 0]).squeeze(-1)

        return Net()


@dataclass
class TorchBagClassifier(DonorBagModel):
    architecture: Literal["deepsets", "attention", "topk_attention", "set_transformer"] = "attention"
    hidden: int = 64
    n_heads: int = 4
    topk_fraction: float = 0.05
    cells_per_bag: int = 512
    epochs: int = 100
    patience: int = 15
    learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    batch_size: int = 8
    seed: int = 0
    device: str = "auto"

    def __post_init__(self) -> None:
        self.model_ = None
        self.device_ = None
        self.input_dim_: int | None = None

    def _resolve_device(self):
        torch, _, _ = _require_torch()
        if self.device == "auto":
            return torch.device("cuda" if torch.cuda.is_available() else "cpu")
        return torch.device(self.device)

    def _build(self, input_dim: int):
        if self.architecture == "deepsets":
            return _DeepSetNet.build(input_dim, self.hidden)
        if self.architecture == "attention":
            return _AttentionNet.build(input_dim, self.hidden, gated=True)
        if self.architecture == "topk_attention":
            return _AttentionNet.build(input_dim, self.hidden, gated=True, topk_fraction=self.topk_fraction)
        if self.architecture == "set_transformer":
            return _SetTransformerNet.build(input_dim, self.hidden, self.n_heads)
        raise ValueError(f"unknown architecture {self.architecture}")

    def _sample_batch(self, bags: list[np.ndarray], indices: np.ndarray, rng: np.random.Generator):
        torch, _, _ = _require_torch()
        p = bags[0].shape[1]
        max_len = min(self.cells_per_bag, max(len(bags[i]) for i in indices))
        X = np.zeros((len(indices), max_len, p), dtype=np.float32)
        mask = np.zeros((len(indices), max_len), dtype=bool)
        for row, idx in enumerate(indices):
            bag = bags[int(idx)]
            n = min(len(bag), self.cells_per_bag)
            choice = rng.choice(len(bag), n, replace=len(bag) < n)
            X[row, :n] = bag[choice]
            mask[row, :n] = True
        return torch.from_numpy(X), torch.from_numpy(mask)

    def fit(self, bags: list[np.ndarray], y: np.ndarray) -> "TorchBagClassifier":
        torch, nn, _ = _require_torch()
        if not bags:
            raise ValueError("empty bag list")
        self.input_dim_ = bags[0].shape[1]
        self.device_ = self._resolve_device()
        torch.manual_seed(self.seed)
        np.random.seed(self.seed)
        self.model_ = self._build(self.input_dim_).to(self.device_)
        optimizer = torch.optim.AdamW(self.model_.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay)
        y = np.asarray(y, dtype=np.float32)
        rng = np.random.default_rng(self.seed)
        order = np.arange(len(bags))
        # Donor-level internal validation; cells are never split across train/validation.
        val_size = max(2, int(round(0.2 * len(order)))) if len(order) >= 10 else max(1, len(order) // 5)
        for _ in range(100):
            rng.shuffle(order)
            val_idx = order[:val_size]
            train_idx = order[val_size:]
            if len(np.unique(y[val_idx])) == 2 and len(np.unique(y[train_idx])) == 2:
                break
        else:
            val_idx = order[:val_size]
            train_idx = order[val_size:]
        loss_fn = nn.BCEWithLogitsLoss()
        best_state, best_loss, stale = None, float("inf"), 0
        for _epoch in range(self.epochs):
            self.model_.train()
            rng.shuffle(train_idx)
            for start in range(0, len(train_idx), self.batch_size):
                idx = train_idx[start : start + self.batch_size]
                xb, mb = self._sample_batch(bags, idx, rng)
                target = torch.from_numpy(y[idx]).to(self.device_)
                logits = self.model_(xb.to(self.device_), mb.to(self.device_))
                loss = loss_fn(logits, target)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
            self.model_.eval()
            with torch.no_grad():
                xb, mb = self._sample_batch(bags, val_idx, rng)
                target = torch.from_numpy(y[val_idx]).to(self.device_)
                val_loss = float(loss_fn(self.model_(xb.to(self.device_), mb.to(self.device_)), target).item())
            if val_loss < best_loss - 1e-5:
                best_loss = val_loss
                best_state = {k: v.detach().cpu().clone() for k, v in self.model_.state_dict().items()}
                stale = 0
            else:
                stale += 1
                if stale >= self.patience:
                    break
        if best_state is not None:
            self.model_.load_state_dict(best_state)
        return self

    def predict_score(self, bags: list[np.ndarray]) -> np.ndarray:
        torch, _, _ = _require_torch()
        if self.model_ is None or self.device_ is None:
            raise RuntimeError("model is not fitted")
        self.model_.eval()
        rng = np.random.default_rng(self.seed + 10_000)
        scores = []
        with torch.no_grad():
            for start in range(0, len(bags), self.batch_size):
                idx = np.arange(start, min(start + self.batch_size, len(bags)))
                xb, mb = self._sample_batch(bags, idx, rng)
                logits = self.model_(xb.to(self.device_), mb.to(self.device_))
                scores.extend(torch.sigmoid(logits).cpu().numpy().tolist())
        return np.asarray(scores, dtype=float)
