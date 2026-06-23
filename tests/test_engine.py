from __future__ import annotations

import numpy as np

from rheumlens.aggregators.classic import MeanAggregator
from rheumlens.estimators.linear import LogisticL2Estimator
from rheumlens.evaluation.engine import FixedMethod, run_fixed_oof
from rheumlens.evaluation.metrics import binary_metrics
from rheumlens.evaluation.splits import make_stratified_folds
from rheumlens.providers.frozen import FrozenCellProvider


def test_oof_engine(toy_data):
    emb, _, _ = toy_data
    mapping = emb.donor_label_map()
    donors = np.asarray(list(mapping))
    y = np.asarray([mapping[d] for d in donors])
    folds = make_stratified_folds(donors, y, n_splits=4, seed=42)
    method = FixedMethod("mean", FrozenCellProvider(), MeanAggregator(), LogisticL2Estimator())
    frame, diagnostics = run_fixed_oof(emb, folds, method, "toy", "test")
    assert frame.donor_id.nunique() == len(donors)
    assert frame.status.eq("SUCCESS").all()
    metrics = binary_metrics(frame.y_true.to_numpy(), frame.score.to_numpy())
    assert 0 <= metrics["roc_auc"] <= 1
    assert len(diagnostics) == 4
