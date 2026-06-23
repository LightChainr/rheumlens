from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from rheumlens.evaluation.bootstrap import paired_auc_bootstrap
from rheumlens.evaluation.metrics import binary_metrics


def summarize_oof(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    good = frame.loc[frame.status.eq("SUCCESS") & frame.score.notna()].copy()
    for keys, group in good.groupby(["cohort", "estimand", "split_id", "method_id"], dropna=False):
        metrics = binary_metrics(group.y_true.to_numpy(), group.score.to_numpy())
        rows.append(
            {
                "cohort": keys[0],
                "estimand": keys[1],
                "split_id": keys[2],
                "method_id": keys[3],
                "n_donors": group.donor_id.nunique(),
                **metrics,
            }
        )
    return pd.DataFrame(rows)


def paired_method_table(frame: pd.DataFrame, reference: str, n_reps: int = 10_000, seed: int = 0) -> pd.DataFrame:
    rows = []
    for (cohort, estimand, split_id), subset in frame.groupby(["cohort", "estimand", "split_id"]):
        ref = subset.loc[subset.method_id.eq(reference) & subset.status.eq("SUCCESS"), ["donor_id", "y_true", "score"]]
        if ref.empty:
            continue
        ref = ref.rename(columns={"score": "score_ref"})
        for method_id, current in subset.loc[subset.status.eq("SUCCESS")].groupby("method_id"):
            if method_id == reference:
                continue
            merged = ref.merge(current[["donor_id", "score"]], on="donor_id", how="inner")
            if merged.y_true.nunique() < 2:
                continue
            result = paired_auc_bootstrap(
                merged.y_true.to_numpy(), merged.score.to_numpy(), merged.score_ref.to_numpy(), n_reps, seed
            )
            result.pop("values", None)
            rows.append(
                {
                    "cohort": cohort,
                    "estimand": estimand,
                    "split_id": split_id,
                    "method_id": method_id,
                    "reference": reference,
                    **result,
                }
            )
    return pd.DataFrame(rows)


def load_result_tree(root: str | Path) -> pd.DataFrame:
    paths = list(Path(root).rglob("oof.csv"))
    if not paths:
        return pd.DataFrame()
    return pd.concat([pd.read_csv(path) for path in paths], ignore_index=True)
