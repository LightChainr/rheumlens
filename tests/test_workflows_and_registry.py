from __future__ import annotations

import numpy as np
import pandas as pd

from rheumlens.evaluation.engine import FixedMethod
from rheumlens.evaluation.splits import make_stratified_folds
from rheumlens.evaluation.workflows import repeated_cv_workflow, successful_auc
from rheumlens.registry import build_method


def test_method_suffix_and_repeated_cv(toy_data, tmp_path):
    emb, expr, _ = toy_data
    method = build_method(
        "moments_mean_var@scgpt",
        defaults={"cell_pca_dim": 4, "estimator_C": 1.0},
        seed=1,
    )
    assert method.data_key == "scgpt"
    assert isinstance(method.method, FixedMethod)
    frame = repeated_cv_workflow(
        emb,
        method.method,
        cohort="toy",
        estimand="test",
        seeds=[1, 2],
        n_splits=2,
        output_dir=tmp_path,
    )
    assert len(frame) == 2
    assert frame.auc.notna().all()


def test_covariate_residualized_method(toy_data):
    emb, _, _ = toy_data
    mapping = emb.donor_label_map()
    donors = np.asarray(list(mapping), dtype=str)
    y = np.asarray([mapping[d] for d in donors], dtype=int)
    cov = pd.DataFrame({"donor_id": donors, "cells": 40, "batch": ["a", "b"] * (len(donors) // 2)})
    folds = make_stratified_folds(donors, y, n_splits=2, seed=3)
    method = build_method("scgpt_mean_covres", seed=3)
    from rheumlens.evaluation.engine import run_fixed_oof

    frame, _ = run_fixed_oof(emb, folds, method.method, "toy", "covres", donor_covariates=cov)
    assert successful_auc(frame) >= 0.0


def test_uder_methods_are_registered() -> None:
    from rheumlens.registry import build_method

    for method_id in ("uder_meanvar@scgpt", "uder_kme_weighted@geneformer"):
        registered = build_method(
            method_id,
            {
                "cell_pca_dim": 3,
                "matched_cells_per_donor": 4,
                "uder_subsamples": 2,
                "kme_rff_dim": 4,
                "kme_scales": [1.0],
                "kme_max_bandwidth_points": 32,
                "kme_max_diagnostic_points": 16,
            },
        )
        assert registered.data_key in {"scgpt", "geneformer"}
