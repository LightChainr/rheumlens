from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy import sparse


def stable_hash(payload: object) -> str:
    raw = json.dumps(payload, sort_keys=True, default=str, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def row_sums(X: np.ndarray | sparse.spmatrix) -> np.ndarray:
    return np.asarray(X.sum(axis=1)).ravel()


def column_means(X: np.ndarray | sparse.spmatrix) -> np.ndarray:
    return np.asarray(X.mean(axis=0)).ravel()


def as_dense(X: np.ndarray | sparse.spmatrix, dtype=np.float32) -> np.ndarray:
    if sparse.issparse(X):
        return X.toarray().astype(dtype, copy=False)
    return np.asarray(X, dtype=dtype)


def grouped_indices(groups: np.ndarray) -> dict[str, np.ndarray]:
    result: dict[str, list[int]] = {}
    for idx, value in enumerate(groups):
        result.setdefault(str(value), []).append(idx)
    return {key: np.asarray(value, dtype=int) for key, value in result.items()}


def safe_quantile(values: np.ndarray, q: float, axis: int = 0) -> np.ndarray:
    if values.shape[0] == 0:
        raise ValueError("cannot compute quantile of empty array")
    return np.quantile(values, q, axis=axis)


def ensure_dir(path: str | Path) -> Path:
    result = Path(path)
    result.mkdir(parents=True, exist_ok=True)
    return result


def finite_or_raise(X: np.ndarray, name: str) -> None:
    if not np.isfinite(X).all():
        bad = np.argwhere(~np.isfinite(X))[:10]
        raise ValueError(f"{name} contains non-finite values at {bad.tolist()}")


def seed_sequence(seed: int, n: int) -> list[int]:
    return [int(x.generate_state(1)[0]) for x in np.random.SeedSequence(seed).spawn(n)]
