from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from common import information_fraction


def test_information_fraction_is_one_minus_r_squared() -> None:
    y = np.array([0, 0, 0, 1, 1, 1], dtype=float)
    design = np.array([[0], [0], [1], [0], [1], [1]], dtype=float)
    result = information_fraction(y, design)

    augmented = np.column_stack([np.ones(len(y)), design])
    fitted = augmented @ np.linalg.lstsq(augmented, y, rcond=None)[0]
    expected_r2 = 1 - np.sum((y - fitted) ** 2) / np.sum((y - y.mean()) ** 2)

    assert np.isclose(result["r2_y_on_design"], expected_r2)
    assert np.isclose(result["information_fraction"], 1 - expected_r2)


def test_information_vanishes_under_perfect_collinearity() -> None:
    y = np.array([0, 1] * 20, dtype=float)
    result = information_fraction(y, y[:, None])

    assert result["information_fraction"] < 1e-12
    assert result["r2_y_on_design"] > 1 - 1e-12


def test_observational_equivalence_under_perfect_collinearity() -> None:
    rng = np.random.default_rng(17)
    y = rng.integers(0, 2, size=100)
    d = y.copy()
    noise = rng.normal(size=100)

    biological_only = 1.7 * y + noise
    technical_only = 1.7 * d + noise

    assert np.array_equal(biological_only, technical_only)


def test_information_is_near_one_under_independence() -> None:
    # Exact balanced 2x2 table avoids a finite-sample random association.
    y = np.array([0, 0, 1, 1] * 25, dtype=float)
    design = np.array([0, 1, 0, 1] * 25, dtype=float)[:, None]
    result = information_fraction(y, design)

    assert np.isclose(result["information_fraction"], 1.0)
    assert np.isclose(result["variance_inflation"], 1.0)
