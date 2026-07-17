from __future__ import annotations

import numpy as np
from scipy import sparse

from rheumlens.data.io import load_npz_dataset, save_npz_dataset
from rheumlens.data.operations import align_features
from rheumlens.data.validation import validate_raw_counts
from rheumlens.types import CellDataset


def test_npz_roundtrip(tmp_path, toy_data):
    emb, _, raw = toy_data
    path = tmp_path / "emb.npz"
    save_npz_dataset(path, emb)
    loaded = load_npz_dataset(path)
    assert loaded.X.shape == emb.X.shape
    assert np.array_equal(loaded.cell_ids, emb.cell_ids)
    report = validate_raw_counts(raw)
    assert report["n_donors"] == 20


def test_sparse_npz_roundtrip(tmp_path, toy_data):
    _, _, raw = toy_data
    sparse_raw = CellDataset(
        X=sparse.csr_matrix(raw.X),
        cell_ids=raw.cell_ids,
        donor_ids=raw.donor_ids,
        y=raw.y,
        feature_names=raw.feature_names,
        cell_types=raw.cell_types,
        name="sparse_raw",
    )
    path = tmp_path / "sparse_raw.npz"
    save_npz_dataset(path, sparse_raw)
    loaded = load_npz_dataset(path)
    assert sparse.isspmatrix_csr(loaded.X)
    assert np.array_equal(loaded.X.toarray(), raw.X)
    assert validate_raw_counts(loaded)["integer_like"]


def test_align_features_preserves_source_order(toy_data):
    _, expression, _ = toy_data
    source = expression
    target_order = np.asarray(["OAS1", "ISG15", "IFI6", "UNSHARED"])
    target = CellDataset(
        X=np.column_stack(
            [expression.X[:, 3], expression.X[:, 0], expression.X[:, 1], expression.X[:, 2]]
        ),
        cell_ids=expression.cell_ids,
        donor_ids=expression.donor_ids,
        y=expression.y,
        feature_names=target_order,
        cell_types=expression.cell_types,
        name="target",
    )
    aligned_source, aligned_target = align_features(source, target)
    assert aligned_source.feature_names.tolist() == ["ISG15", "IFI6", "OAS1"]
    assert np.array_equal(aligned_source.feature_names, aligned_target.feature_names)
    assert np.allclose(aligned_source.X, aligned_target.X)
