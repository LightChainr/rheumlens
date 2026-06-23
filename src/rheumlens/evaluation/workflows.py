from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from tqdm import trange

from rheumlens.evaluation.engine import BagMethod, FixedMethod, run_bag_oof, run_fixed_oof
from rheumlens.evaluation.nulls import mean_sufficient_surrogate, moment2_sufficient_surrogate
from rheumlens.evaluation.permutation import permute_donor_labels
from rheumlens.evaluation.splits import make_stratified_folds
from rheumlens.types import CellDataset, OuterFold
from rheumlens.utils.arrays import ensure_dir
from rheumlens.utils.manifests import write_json


def successful_auc(frame: pd.DataFrame) -> float:
    good = frame.status.eq("SUCCESS") & frame.score.notna()
    subset = frame.loc[good]
    if subset.empty or subset.y_true.nunique() < 2:
        return float("nan")
    return float(roc_auc_score(subset.y_true, subset.score))


def run_registered_method(
    data: CellDataset,
    folds: list[OuterFold],
    method: FixedMethod | BagMethod,
    cohort: str,
    estimand: str,
    expression_data: CellDataset | None = None,
    donor_covariates: pd.DataFrame | None = None,
    random_state: int = 0,
) -> pd.DataFrame:
    if isinstance(method, FixedMethod):
        frame, _ = run_fixed_oof(
            data,
            folds,
            method,
            cohort=cohort,
            estimand=estimand,
            expression_data=expression_data,
            donor_covariates=donor_covariates,
            random_state=random_state,
        )
        return frame
    return run_bag_oof(
        data,
        folds,
        method,
        cohort=cohort,
        estimand=estimand,
        random_state=random_state,
    )


def permutation_workflow(
    data: CellDataset,
    folds: list[OuterFold],
    method: FixedMethod | BagMethod,
    cohort: str,
    estimand: str,
    observed_auc: float,
    n_reps: int,
    seed: int,
    expression_data: CellDataset | None = None,
    donor_covariates: pd.DataFrame | None = None,
    rebuild_stratified_folds: bool = False,
    output_dir: str | Path | None = None,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    values = np.full(n_reps, np.nan, dtype=float)
    rep_seeds = rng.integers(0, 2**31 - 1, size=n_reps, dtype=np.int64)
    for rep in trange(n_reps, desc=f"Permutation {cohort}/{method.id}"):
        permuted = permute_donor_labels(data, int(rep_seeds[rep]))
        rep_folds = folds
        if rebuild_stratified_folds:
            mapping = permuted.donor_label_map()
            donors = np.asarray(list(mapping), dtype=str)
            y = np.asarray([mapping[d] for d in donors], dtype=int)
            rep_folds = make_stratified_folds(
                donors,
                y,
                n_splits=len(folds),
                seed=int(rep_seeds[rep]),
                split_id=f"perm_{rep:04d}",
            )
        frame = run_registered_method(
            permuted,
            rep_folds,
            method,
            cohort,
            f"{estimand}_permutation",
            expression_data=expression_data,
            donor_covariates=donor_covariates,
            random_state=int(rep_seeds[rep]),
        )
        values[rep] = successful_auc(frame)
    finite = values[np.isfinite(values)]
    empirical_p = float((1 + np.sum(finite >= observed_auc)) / (1 + len(finite)))
    result: dict[str, object] = {
        "cohort": cohort,
        "method_id": method.id,
        "observed_auc": float(observed_auc),
        "n_requested": int(n_reps),
        "n_finite": int(len(finite)),
        "empirical_p": empirical_p,
        "null_mean": float(np.mean(finite)) if len(finite) else np.nan,
        "null_q025": float(np.quantile(finite, 0.025)) if len(finite) else np.nan,
        "null_q975": float(np.quantile(finite, 0.975)) if len(finite) else np.nan,
    }
    if output_dir is not None:
        out = ensure_dir(output_dir)
        np.save(out / "null_auc.npy", values)
        pd.DataFrame({"rep": np.arange(n_reps), "seed": rep_seeds, "auc": values}).to_csv(
            out / "permutation_auc.csv", index=False
        )
        write_json(out / "summary.json", result)
    return result


def repeated_cv_workflow(
    data: CellDataset,
    method: FixedMethod | BagMethod,
    cohort: str,
    estimand: str,
    seeds: list[int],
    n_splits: int = 5,
    expression_data: CellDataset | None = None,
    donor_covariates: pd.DataFrame | None = None,
    output_dir: str | Path | None = None,
) -> pd.DataFrame:
    mapping = data.donor_label_map()
    donors = np.asarray(list(mapping), dtype=str)
    y = np.asarray([mapping[d] for d in donors], dtype=int)
    rows: list[dict[str, object]] = []
    out = ensure_dir(output_dir) if output_dir is not None else None
    for repeat_index, seed in enumerate(seeds):
        folds = make_stratified_folds(
            donors,
            y,
            n_splits=n_splits,
            seed=int(seed),
            split_id=f"repeated_seed_{seed}",
        )
        frame = run_registered_method(
            data,
            folds,
            method,
            cohort,
            f"{estimand}_repeated_cv",
            expression_data=expression_data,
            donor_covariates=donor_covariates,
            random_state=int(seed),
        )
        auc = successful_auc(frame)
        rows.append(
            {
                "cohort": cohort,
                "method_id": method.id,
                "repeat_index": repeat_index,
                "seed": int(seed),
                "auc": auc,
                "n_success": int((frame.status == "SUCCESS").sum()),
            }
        )
        if out is not None:
            frame.to_csv(out / f"oof_seed_{seed}.csv", index=False)
    result = pd.DataFrame(rows)
    if out is not None:
        result.to_csv(out / "repeated_cv_summary.csv", index=False)
    return result


def surrogate_delta_workflow(
    data: CellDataset,
    folds: list[OuterFold],
    target_method: FixedMethod | BagMethod,
    reference_method: FixedMethod | BagMethod,
    cohort: str,
    estimand: str,
    null_kind: str,
    n_reps: int,
    seed: int,
    donor_covariates: pd.DataFrame | None = None,
    output_dir: str | Path | None = None,
) -> dict[str, object]:
    if null_kind not in {"mean_sufficient", "moment2_sufficient"}:
        raise ValueError("null_kind must be mean_sufficient or moment2_sufficient")
    generator: Callable[[CellDataset, int], CellDataset] = (
        mean_sufficient_surrogate if null_kind == "mean_sufficient" else moment2_sufficient_surrogate
    )
    rng = np.random.default_rng(seed)
    values = np.full(n_reps, np.nan, dtype=float)
    rep_seeds = rng.integers(0, 2**31 - 1, size=n_reps, dtype=np.int64)
    for rep in trange(n_reps, desc=f"Surrogate {null_kind} {cohort}"):
        surrogate = generator(data, int(rep_seeds[rep]))
        target = run_registered_method(
            surrogate,
            folds,
            target_method,
            cohort,
            f"{estimand}_{null_kind}",
            donor_covariates=donor_covariates,
            random_state=int(rep_seeds[rep]),
        )
        reference = run_registered_method(
            surrogate,
            folds,
            reference_method,
            cohort,
            f"{estimand}_{null_kind}",
            donor_covariates=donor_covariates,
            random_state=int(rep_seeds[rep]),
        )
        target_auc = successful_auc(target)
        reference_auc = successful_auc(reference)
        if np.isfinite(target_auc) and np.isfinite(reference_auc):
            values[rep] = target_auc - reference_auc
    finite = values[np.isfinite(values)]
    result: dict[str, object] = {
        "cohort": cohort,
        "target_method": target_method.id,
        "reference_method": reference_method.id,
        "null_kind": null_kind,
        "n_requested": int(n_reps),
        "n_finite": int(len(finite)),
        "delta_mean": float(np.mean(finite)) if len(finite) else np.nan,
        "delta_q025": float(np.quantile(finite, 0.025)) if len(finite) else np.nan,
        "delta_q95": float(np.quantile(finite, 0.95)) if len(finite) else np.nan,
        "delta_q975": float(np.quantile(finite, 0.975)) if len(finite) else np.nan,
    }
    if output_dir is not None:
        out = ensure_dir(output_dir)
        np.save(out / "delta_auc_null.npy", values)
        pd.DataFrame({"rep": np.arange(n_reps), "seed": rep_seeds, "delta_auc": values}).to_csv(
            out / "surrogate_delta_auc.csv", index=False
        )
        write_json(out / "summary.json", result)
    return result
