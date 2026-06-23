from __future__ import annotations

import numpy as np

from rheumlens.evaluation.nulls import mean_sufficient_surrogate, moment2_sufficient_surrogate
from rheumlens.utils.arrays import grouped_indices


def test_surrogates_preserve_donor_means(toy_data):
    emb, _, _ = toy_data
    null = mean_sufficient_surrogate(emb, seed=3)
    for donor, idx in grouped_indices(emb.donor_ids).items():
        assert np.allclose(emb.X[idx].mean(axis=0), null.X[idx].mean(axis=0), atol=0.6)
    null2 = moment2_sufficient_surrogate(emb, seed=4)
    assert null2.X.shape == emb.X.shape
