from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

from rheumlens.types import OuterFold


def make_stratified_folds(
    donor_ids: np.ndarray,
    y: np.ndarray,
    n_splits: int = 5,
    seed: int = 42,
    split_id: str | None = None,
) -> list[OuterFold]:
    donor_ids = np.asarray(donor_ids).astype(str)
    y = np.asarray(y).astype(int)
    if len(donor_ids) != len(y):
        raise ValueError("donor_ids/y length mismatch")
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    name = split_id or f"stratified_{n_splits}fold_seed{seed}"
    folds = []
    for fold, (train_idx, test_idx) in enumerate(splitter.split(donor_ids, y)):
        folds.append(
            OuterFold(
                split_id=name,
                fold=fold,
                train_donors=tuple(donor_ids[train_idx].tolist()),
                test_donors=tuple(donor_ids[test_idx].tolist()),
            )
        )
    return folds


def save_folds(path: str | Path, folds: list[OuterFold]) -> None:
    rows = []
    for item in folds:
        rows.extend(
            {"split_id": item.split_id, "fold": item.fold, "role": "train", "donor_id": donor}
            for donor in item.train_donors
        )
        rows.extend(
            {"split_id": item.split_id, "fold": item.fold, "role": "test", "donor_id": donor}
            for donor in item.test_donors
        )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(target, index=False)


def load_folds(path: str | Path) -> list[OuterFold]:
    frame = pd.read_csv(path, dtype={"donor_id": str})
    folds = []
    for (split_id, fold), group in frame.groupby(["split_id", "fold"], sort=True):
        train = tuple(group.loc[group.role == "train", "donor_id"].astype(str))
        test = tuple(group.loc[group.role == "test", "donor_id"].astype(str))
        folds.append(OuterFold(str(split_id), int(fold), train, test))
    return folds
