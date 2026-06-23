from __future__ import annotations

import os
import sys
from pathlib import Path

# Prevent BLAS thread oversubscription on 144-core A800 hosts and CI containers.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import numpy as np
import pytest

from rheumlens.types import CellDataset


@pytest.fixture
def toy_data() -> tuple[CellDataset, CellDataset, CellDataset]:
    rng = np.random.default_rng(1)
    n_donors, cells, p = 20, 40, 12
    W = rng.normal(size=(p, 8))
    X, E, C, donors, labels, cell_ids, cell_types = [], [], [], [], [], [], []
    genes = np.asarray(["ISG15", "IFI6", "MX1", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT3", "IFI44", "IFI44L", "STAT1", "HERC5"])
    for d in range(n_donors):
        y = d % 2
        for i in range(cells):
            expr = rng.normal(size=p)
            if y and i < 4:
                expr[:6] += 2.0
            emb = np.tanh(expr @ W)
            raw = np.maximum(0, np.rint(np.exp(np.clip(expr, -2, 3)))).astype(int)
            X.append(emb)
            E.append(expr)
            C.append(raw)
            donors.append(f"D{d:02d}")
            labels.append(y)
            cell_ids.append(f"D{d:02d}_C{i:03d}")
            cell_types.append("CD4_T" if i % 2 else "monocyte")
    common = dict(
        cell_ids=np.asarray(cell_ids), donor_ids=np.asarray(donors), y=np.asarray(labels),
        cell_types=np.asarray(cell_types)
    )
    emb = CellDataset(np.asarray(X, dtype=np.float32), feature_names=np.asarray([f"e{i}" for i in range(8)]), name="emb", **common)
    expr = CellDataset(np.asarray(E, dtype=np.float32), feature_names=genes, name="expr", **common)
    raw = CellDataset(np.asarray(C, dtype=np.int32), feature_names=genes, name="raw", **common)
    return emb, expr, raw
