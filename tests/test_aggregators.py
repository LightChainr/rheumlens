from __future__ import annotations

import numpy as np

from rheumlens.aggregators.focus import FocusLiteAggregator, FocusQuery
from rheumlens.aggregators.kme import CCKMEUAggregator, MultiScaleRFFKMEAggregator
from rheumlens.aggregators.moments import QuantileAggregator, TailFractionAggregator
from rheumlens.aggregators.red import REDAggregator
from rheumlens.types import FitContext


def test_tail_threshold_is_train_fitted(toy_data):
    emb, _, _ = toy_data
    donors = list(emb.donors.astype(str))
    context = FitContext(emb.donor_label_map())
    agg = TailFractionAggregator().fit(emb, donors[:14], context)
    result = agg.transform(emb, donors[14:], context)
    assert result.X.shape[0] == 6
    assert not np.allclose(result.X.mean(axis=1), 0.05)


def test_kme_and_cckme(toy_data):
    emb, _, _ = toy_data
    donors = list(emb.donors.astype(str))
    context = FitContext(emb.donor_label_map(), random_state=2)
    kme = MultiScaleRFFKMEAggregator(rff_dim=16, scales=(1.0,), include_linear_mean=True)
    kme.fit(emb, donors[:14], context)
    result = kme.transform(emb, donors[14:], context)
    assert result.X.shape == (6, 8 + 16)
    cc = CCKMEUAggregator(rff_dim=8, scales=(1.0,), include_linear_mean=True, subsample_size=20, n_subsamples=4)
    cc.fit(emb, donors[:14], context)
    cc_result = cc.transform(emb, donors[14:], context)
    assert cc_result.uncertainty is not None
    assert np.isfinite(cc_result.X).all()


def test_focus_and_red(toy_data):
    emb, expr, _ = toy_data
    donors = list(emb.donors.astype(str))
    context = FitContext(emb.donor_label_map(), expression=expr)
    focus = FocusLiteAggregator((FocusQuery("ifn", ("ISG15", "IFI6", "MX1"), "global"),), random_queries=1)
    focus.fit(emb, donors[:14], context)
    result = focus.transform(emb, donors[14:], context)
    assert result.X.shape[1] == 12
    assert np.std(result.X[:, 4]) > 0
    red = REDAggregator(n_states=4, max_fit_cells=1000)
    red.fit(emb, donors[:14], context)
    red_result = red.transform(emb, donors[14:], context)
    assert red_result.X.shape == (6, 16)
