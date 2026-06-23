from __future__ import annotations

import copy
from typing import Callable

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from rheumlens.types import CellDataset


def permute_donor_labels(data: CellDataset, seed: int) -> CellDataset:
    rng = np.random.default_rng(seed)
    mapping = data.donor_label_map()
    donors = np.asarray(list(mapping), dtype=str)
    values = np.asarray([mapping[d] for d in donors], dtype=int)
    permuted = values[rng.permutation(len(values))]
    new_mapping = dict(zip(donors, permuted, strict=True))
    new_y = np.asarray([new_mapping[str(d)] for d in data.donor_ids], dtype=int)
    return CellDataset(
        X=data.X,
        cell_ids=data.cell_ids,
        donor_ids=data.donor_ids,
        y=new_y,
        feature_names=data.feature_names,
        cell_types=data.cell_types,
        cohorts=data.cohorts,
        metadata=data.metadata,
        name=f"{data.name}_permuted",
    )


def full_pipeline_permutation(
    data: CellDataset,
    runner: Callable[[CellDataset, int], pd.DataFrame],
    observed_auc: float,
    n_reps: int = 1000,
    seed: int = 0,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    null = []
    for rep in range(n_reps):
        permuted = permute_donor_labels(data, int(rng.integers(0, 2**31 - 1)))
        frame = runner(permuted, rep)
        good = frame.status.eq("SUCCESS") & frame.score.notna()
        if good.sum() == 0 or frame.loc[good, "y_true"].nunique() < 2:
            null.append(np.nan)
        else:
            null.append(roc_auc_score(frame.loc[good, "y_true"], frame.loc[good, "score"]))
    values = np.asarray(null, dtype=float)
    finite = values[np.isfinite(values)]
    p = (1 + np.sum(finite >= observed_auc)) / (1 + len(finite))
    return {"observed_auc": observed_auc, "empirical_p": float(p), "null_auc": values.tolist()}
