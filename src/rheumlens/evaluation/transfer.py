from __future__ import annotations

import copy

import numpy as np
import pandas as pd

from rheumlens.data.operations import align_features
from rheumlens.evaluation.engine import FixedMethod
from rheumlens.types import CellDataset, DonorFeatures, FitContext


def _maybe_align(source: CellDataset, target: CellDataset) -> tuple[CellDataset, CellDataset]:
    if np.array_equal(source.feature_names.astype(str), target.feature_names.astype(str)):
        return source, target
    # Frozen embeddings should have the same named dimensions. Expression/count datasets can
    # be aligned to shared genes using source ordering.
    source_names = source.feature_names.astype(str)
    target_names = target.feature_names.astype(str)
    looks_like_embedding = all(name.startswith(("dim_", "cell_pc_", "scgpt_", "geneformer_")) for name in source_names[: min(10, len(source_names))])
    if looks_like_embedding:
        raise ValueError("source/target embedding dimensions do not match")
    return align_features(source, target)


def run_source_target(
    source: CellDataset,
    target: CellDataset,
    method: FixedMethod,
    expression_source: CellDataset | None = None,
    expression_target: CellDataset | None = None,
    donor_covariates_source: pd.DataFrame | None = None,
    donor_covariates_target: pd.DataFrame | None = None,
    random_state: int = 0,
) -> pd.DataFrame:
    source, target = _maybe_align(source, target)
    if expression_source is not None and expression_target is not None:
        expression_source, expression_target = _maybe_align(expression_source, expression_target)
    source_y = source.donor_label_map()
    target_y = target.donor_label_map()
    provider = copy.deepcopy(method.provider)
    aggregator = copy.deepcopy(method.aggregator)
    estimator = copy.deepcopy(method.estimator)
    postprocessor = copy.deepcopy(method.postprocessor)
    source_context = FitContext(
        source_y,
        random_state=random_state,
        expression=expression_source,
        donor_covariates=donor_covariates_source,
        extras={"transfer_role": "source"},
    )
    target_context = FitContext(
        source_y,
        random_state=random_state,
        expression=expression_target,
        donor_covariates=donor_covariates_target,
        extras={"transfer_role": "target"},
    )
    provider.fit(source, list(source.donors.astype(str)), source_context)
    source_rep = provider.transform(source, list(source.donors.astype(str)), source_context)
    target_rep = provider.transform(target, list(target.donors.astype(str)), target_context)
    if isinstance(source_rep, DonorFeatures):
        if not isinstance(target_rep, DonorFeatures):
            raise TypeError("provider returned inconsistent levels")
        source_features, target_features = source_rep, target_rep
    else:
        if aggregator is None:
            raise ValueError("cell provider requires aggregator")
        aggregator.fit(source_rep, list(source.donors.astype(str)), source_context)
        source_features = aggregator.transform(source_rep, list(source.donors.astype(str)), source_context)
        target_features = aggregator.transform(target_rep, list(target.donors.astype(str)), target_context)
    if postprocessor is not None:
        postprocessor.fit(source_features, source_context)
        source_features = postprocessor.transform(source_features, source_context)
        target_features = postprocessor.transform(target_features, target_context)
    y_train = np.asarray([source_y[str(x)] for x in source_features.donor_ids])
    estimator.fit(source_features.X, y_train)
    score = estimator.predict_score(target_features.X)
    return pd.DataFrame(
        {
            "donor_id": target_features.donor_ids,
            "y_true": [target_y[str(x)] for x in target_features.donor_ids],
            "score": score,
            "source": source.name,
            "target": target.name,
            "method_id": method.id,
        }
    )
