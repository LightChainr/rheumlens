from __future__ import annotations

import copy
import resource
import time
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from rheumlens.aggregators.base import DonorAggregator
from rheumlens.bag_models.base import DonorBagModel
from rheumlens.estimators.base import DonorEstimator
from rheumlens.postprocessors.base import DonorPostprocessor
from rheumlens.types import CellDataset, DonorFeatures, FitContext, OuterFold
from rheumlens.utils.arrays import grouped_indices, stable_hash


@dataclass
class FixedMethod:
    id: str
    provider: Any
    aggregator: DonorAggregator | None
    estimator: DonorEstimator
    postprocessor: DonorPostprocessor | None = None
    params: dict[str, Any] = field(default_factory=dict)


@dataclass
class BagMethod:
    id: str
    provider: Any
    model: DonorBagModel
    params: dict[str, Any] = field(default_factory=dict)


def _labels(order: np.ndarray, mapping: dict[str, int]) -> np.ndarray:
    return np.asarray([mapping[str(d)] for d in order], dtype=int)


def _n_cells(data: CellDataset, donors: list[str]) -> dict[str, int]:
    subset = data.subset_donors(donors)
    groups = grouped_indices(subset.donor_ids)
    return {d: len(idx) for d, idx in groups.items()}


def run_fixed_oof(
    data: CellDataset,
    folds: list[OuterFold],
    method: FixedMethod,
    cohort: str,
    estimand: str,
    expression_data: CellDataset | None = None,
    donor_covariates: pd.DataFrame | None = None,
    random_state: int = 0,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    y_map = data.donor_label_map()
    all_rows: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for fold in folds:
        started = time.perf_counter()
        provider = copy.deepcopy(method.provider)
        aggregator = copy.deepcopy(method.aggregator)
        estimator = copy.deepcopy(method.estimator)
        postprocessor = copy.deepcopy(method.postprocessor)
        context = FitContext(
            y_by_donor=y_map,
            random_state=random_state + fold.fold,
            expression=expression_data,
            donor_covariates=donor_covariates,
            extras={"cohort": cohort, "estimand": estimand, "split_id": fold.split_id},
        )
        status = "SUCCESS"
        try:
            provider.fit(data, list(fold.train_donors), context)
            train_rep = provider.transform(data, list(fold.train_donors), context)
            test_rep = provider.transform(data, list(fold.test_donors), context)
            if isinstance(train_rep, DonorFeatures):
                if aggregator is not None:
                    raise ValueError("donor-level provider must not be followed by cell aggregator")
                train_features, test_features = train_rep, test_rep
            else:
                if aggregator is None:
                    raise ValueError("cell-level provider requires an aggregator")
                aggregator.fit(train_rep, list(fold.train_donors), context)
                train_features = aggregator.transform(train_rep, list(fold.train_donors), context)
                test_features = aggregator.transform(test_rep, list(fold.test_donors), context)
            if postprocessor is not None:
                postprocessor.fit(train_features, context)
                train_features = postprocessor.transform(train_features, context)
                test_features = postprocessor.transform(test_features, context)
            y_train = _labels(train_features.donor_ids, y_map)
            y_test = _labels(test_features.donor_ids, y_map)
            sample_weight = None
            if train_features.uncertainty is not None and getattr(estimator, "use_uncertainty_weights", False):
                sample_weight = estimator.weights_from_uncertainty(train_features.uncertainty)
            estimator.fit(train_features.X, y_train, sample_weight=sample_weight)
            score = estimator.predict_score(test_features.X)
            run_hash = stable_hash(
                {
                    "method": method.id,
                    "params": method.params,
                    "split": fold.split_id,
                    "fold": fold.fold,
                    "train": sorted(fold.train_donors),
                    "test": sorted(fold.test_donors),
                }
            )
            cell_counts = _n_cells(data, list(test_features.donor_ids))
            runtime = time.perf_counter() - started
            peak_ram = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
            for donor, y_value, pred in zip(test_features.donor_ids, y_test, score, strict=True):
                all_rows.append(
                    {
                        "cohort": cohort,
                        "estimand": estimand,
                        "split_id": fold.split_id,
                        "fold": fold.fold,
                        "method_id": method.id,
                        "donor_id": str(donor),
                        "y_true": int(y_value),
                        "score": float(pred),
                        "n_cells": int(cell_counts[str(donor)]),
                        "status": status,
                        "runtime_sec": float(runtime),
                        "peak_ram_mb": float(peak_ram),
                        "run_hash": run_hash,
                    }
                )
            diagnostics.append(
                {
                    "method_id": method.id,
                    "split_id": fold.split_id,
                    "fold": fold.fold,
                    "provider": getattr(train_rep, "diagnostics", {}),
                    "aggregator": getattr(test_features, "diagnostics", {}),
                    "postprocessor": None if postprocessor is None else postprocessor.__class__.__name__,
                }
            )
        except Exception as exc:
            status = "NUMERIC_FAIL" if isinstance(exc, (FloatingPointError, np.linalg.LinAlgError)) else "INVALID_INPUT"
            runtime = time.perf_counter() - started
            for donor in fold.test_donors:
                all_rows.append(
                    {
                        "cohort": cohort,
                        "estimand": estimand,
                        "split_id": fold.split_id,
                        "fold": fold.fold,
                        "method_id": method.id,
                        "donor_id": str(donor),
                        "y_true": int(y_map[str(donor)]),
                        "score": np.nan,
                        "n_cells": _n_cells(data, [donor]).get(str(donor), 0),
                        "status": status,
                        "runtime_sec": float(runtime),
                        "run_hash": stable_hash({"method": method.id, "fold": fold.fold, "error": repr(exc)}),
                        "notes": repr(exc),
                    }
                )
    return pd.DataFrame(all_rows), diagnostics


def _bags(data: CellDataset, donors: list[str]) -> tuple[list[np.ndarray], np.ndarray]:
    subset = data.subset_donors(donors)
    groups = grouped_indices(subset.donor_ids)
    order = np.asarray(list(groups), dtype=str)
    dense = subset.to_dense()
    return [dense[groups[d]] for d in order], order


def run_bag_oof(
    data: CellDataset,
    folds: list[OuterFold],
    method: BagMethod,
    cohort: str,
    estimand: str,
    random_state: int = 0,
) -> pd.DataFrame:
    y_map = data.donor_label_map()
    rows: list[dict[str, Any]] = []
    for fold in folds:
        started = time.perf_counter()
        provider = copy.deepcopy(method.provider)
        model = copy.deepcopy(method.model)
        context = FitContext(y_by_donor=y_map, random_state=random_state + fold.fold)
        try:
            provider.fit(data, list(fold.train_donors), context)
            train_rep = provider.transform(data, list(fold.train_donors), context)
            test_rep = provider.transform(data, list(fold.test_donors), context)
            if not isinstance(train_rep, CellDataset) or not isinstance(test_rep, CellDataset):
                raise ValueError("bag model requires cell-level provider")
            train_bags, train_order = _bags(train_rep, list(fold.train_donors))
            test_bags, test_order = _bags(test_rep, list(fold.test_donors))
            y_train = _labels(train_order, y_map)
            y_test = _labels(test_order, y_map)
            model.fit(train_bags, y_train)
            score = model.predict_score(test_bags)
            runtime = time.perf_counter() - started
            counts = _n_cells(data, list(test_order))
            run_hash = stable_hash({"method": method.id, "fold": fold.fold, "params": method.params})
            for donor, y_value, pred in zip(test_order, y_test, score, strict=True):
                rows.append(
                    {
                        "cohort": cohort,
                        "estimand": estimand,
                        "split_id": fold.split_id,
                        "fold": fold.fold,
                        "method_id": method.id,
                        "donor_id": str(donor),
                        "y_true": int(y_value),
                        "score": float(pred),
                        "n_cells": counts[str(donor)],
                        "status": "SUCCESS",
                        "runtime_sec": runtime,
                        "run_hash": run_hash,
                    }
                )
        except Exception as exc:
            for donor in fold.test_donors:
                rows.append(
                    {
                        "cohort": cohort,
                        "estimand": estimand,
                        "split_id": fold.split_id,
                        "fold": fold.fold,
                        "method_id": method.id,
                        "donor_id": str(donor),
                        "y_true": int(y_map[str(donor)]),
                        "score": np.nan,
                        "n_cells": _n_cells(data, [donor]).get(str(donor), 0),
                        "status": "INVALID_INPUT",
                        "runtime_sec": time.perf_counter() - started,
                        "run_hash": stable_hash({"method": method.id, "fold": fold.fold, "error": repr(exc)}),
                        "notes": repr(exc),
                    }
                )
    return pd.DataFrame(rows)
