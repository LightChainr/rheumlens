from __future__ import annotations

import numpy as np

from rheumlens.data.io import load_npz_dataset, save_npz_dataset
from rheumlens.data.validation import validate_raw_counts


def test_npz_roundtrip(tmp_path, toy_data):
    emb, _, raw = toy_data
    path = tmp_path / "emb.npz"
    save_npz_dataset(path, emb)
    loaded = load_npz_dataset(path)
    assert loaded.X.shape == emb.X.shape
    assert np.array_equal(loaded.cell_ids, emb.cell_ids)
    report = validate_raw_counts(raw)
    assert report["n_donors"] == 20
