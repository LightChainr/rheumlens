from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rheumlens.bag_models.base import DonorBagModel
from rheumlens.bag_models.torch_models import _require_torch


@dataclass
class FocusAdapterBagModel(DonorBagModel):
    """Low-rank supervised adaptation of a frozen query bank.

    The model receives fixed initial query vectors. Only a low-rank query update and a small
    donor classifier are trained. scGPT/Geneformer embeddings remain frozen.
    """

    initial_queries: np.ndarray
    rank: int = 4
    cells_per_bag: int = 512
    epochs: int = 100
    patience: int = 15
    learning_rate: float = 1e-3
    topk_fraction: float = 0.05
    seed: int = 0
    device: str = "auto"

    def __post_init__(self) -> None:
        self.model_ = None
        self.device_ = None

    def _sample(self, bags: list[np.ndarray], rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        max_len = min(self.cells_per_bag, max(len(x) for x in bags))
        p = bags[0].shape[1]
        X = np.zeros((len(bags), max_len, p), dtype=np.float32)
        mask = np.zeros((len(bags), max_len), dtype=bool)
        for i, bag in enumerate(bags):
            n = min(len(bag), self.cells_per_bag)
            idx = rng.choice(len(bag), n, replace=len(bag) < n)
            X[i, :n] = bag[idx]
            mask[i, :n] = True
        return X, mask

    def fit(self, bags: list[np.ndarray], y: np.ndarray) -> "FocusAdapterBagModel":
        torch, nn, F = _require_torch()
        self.device_ = torch.device("cuda" if self.device == "auto" and torch.cuda.is_available() else "cpu")
        q0 = np.asarray(self.initial_queries, dtype=np.float32)
        if q0.ndim != 2 or q0.shape[1] != bags[0].shape[1]:
            raise ValueError("initial_queries must be n_queries × embedding_dim")
        torch.manual_seed(self.seed)
        rng = np.random.default_rng(self.seed)
        n_queries, p = q0.shape

        class Net(nn.Module):
            def __init__(self) -> None:
                super().__init__()
                self.register_buffer("q0", torch.from_numpy(q0))
                self.A = nn.Parameter(torch.zeros(n_queries, self.rank))
                self.B = nn.Parameter(torch.randn(self.rank, p) * 0.01)
                self.head = nn.Linear(n_queries * 4, 1)

            def forward(self, x, mask):
                q = F.normalize(self.q0 + self.A @ self.B, dim=1)
                xn = F.normalize(x, dim=-1)
                score = torch.einsum("bnp,qp->bnq", xn, q)
                score = score.masked_fill(~mask.unsqueeze(-1), float("nan"))
                mean = torch.nanmean(score, dim=1)
                var = torch.nanmean((score - mean.unsqueeze(1)) ** 2, dim=1)
                q95 = torch.nanquantile(score, 0.95, dim=1)
                topk_values = []
                for b in range(score.shape[0]):
                    n = int(mask[b].sum().item())
                    k = max(1, int(np.ceil(self.topk_fraction * n)))
                    topk_values.append(torch.topk(score[b, :n], k=k, dim=0).values.mean(dim=0))
                topk = torch.stack(topk_values)
                features = torch.cat([mean, var, q95, topk], dim=1)
                return self.head(features).squeeze(-1)

        model = Net().to(self.device_)
        optimizer = torch.optim.AdamW(model.parameters(), lr=self.learning_rate, weight_decay=1e-4)
        loss_fn = nn.BCEWithLogitsLoss()
        y = np.asarray(y, dtype=np.float32)
        indices = np.arange(len(bags))
        val_size = max(2, int(round(0.2 * len(indices)))) if len(indices) >= 10 else max(1, len(indices) // 5)
        for _ in range(100):
            rng.shuffle(indices)
            val_idx, train_idx = indices[:val_size], indices[val_size:]
            if len(np.unique(y[val_idx])) == 2 and len(np.unique(y[train_idx])) == 2:
                break
        best_state, best_loss, stale = None, float("inf"), 0
        for _epoch in range(self.epochs):
            model.train()
            X, M = self._sample([bags[i] for i in train_idx], rng)
            logits = model(torch.from_numpy(X).to(self.device_), torch.from_numpy(M).to(self.device_))
            target = torch.from_numpy(y[train_idx]).to(self.device_)
            loss = loss_fn(logits, target)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            model.eval()
            with torch.no_grad():
                Xv, Mv = self._sample([bags[i] for i in val_idx], rng)
                val_logits = model(torch.from_numpy(Xv).to(self.device_), torch.from_numpy(Mv).to(self.device_))
                val_loss = float(loss_fn(val_logits, torch.from_numpy(y[val_idx]).to(self.device_)).item())
            if val_loss < best_loss - 1e-5:
                best_loss = val_loss
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
                stale = 0
            else:
                stale += 1
                if stale >= self.patience:
                    break
        if best_state is not None:
            model.load_state_dict(best_state)
        self.model_ = model
        return self

    def predict_score(self, bags: list[np.ndarray]) -> np.ndarray:
        torch, _, _ = _require_torch()
        if self.model_ is None:
            raise RuntimeError("model is not fitted")
        rng = np.random.default_rng(self.seed + 10000)
        X, M = self._sample(bags, rng)
        self.model_.eval()
        with torch.no_grad():
            logits = self.model_(torch.from_numpy(X).to(self.device_), torch.from_numpy(M).to(self.device_))
            return torch.sigmoid(logits).cpu().numpy()
