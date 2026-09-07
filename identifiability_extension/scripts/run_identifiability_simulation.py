#!/usr/bin/env python3
"""Ground-truth simulation for design-label non-identifiability."""

from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from common import WORKSPACE, information_fraction, write_json


OUT = WORKSPACE / "results" / "simulation"
BETAS = (0.0, 0.5, 1.0, 2.0)
DELTAS = (0.0, 0.5, 1.0, 2.0)
RHOS = (0.0, 0.25, 0.5, 0.75, 0.9, 0.98, 1.0)
N_REPLICATES = int(os.environ.get("SIM_REPLICATES", "100"))
N_DONORS = int(os.environ.get("SIM_N_DONORS", "240"))
N_TARGET = int(os.environ.get("SIM_N_TARGET", "2000"))
N_FEATURES = 20
BASE_SEED = 2026072701
N_JOBS = min(int(os.environ.get("N_JOBS", "12")), os.cpu_count() or 1)


@dataclass(frozen=True)
class Simulated:
    x: np.ndarray
    y: np.ndarray
    d: np.ndarray
    biological_direction: np.ndarray


def labels_and_design(
    rng: np.random.Generator,
    n: int,
    rho: float,
) -> tuple[np.ndarray, np.ndarray]:
    y = rng.integers(0, 2, size=n, dtype=np.int8)
    if rho >= 1.0:
        return y, y.copy()
    align_probability = (1.0 + rho) / 2.0
    align = rng.random(n) < align_probability
    d = np.where(align, y, 1 - y).astype(np.int8)
    return y, d


def generate(
    rng: np.random.Generator,
    n: int,
    beta: float,
    delta: float,
    rho: float,
) -> Simulated:
    y, d = labels_and_design(rng, n, rho)
    biological = np.zeros(N_FEATURES)
    technical = np.zeros(N_FEATURES)
    biological[0] = 1.0
    technical[1] = 1.0
    y_centered = 2.0 * y - 1.0
    d_centered = 2.0 * d - 1.0
    noise = rng.normal(size=(n, N_FEATURES))
    x = (
        beta * y_centered[:, None] * biological[None, :]
        + delta * d_centered[:, None] * technical[None, :]
        + noise
    )
    return Simulated(x=x, y=y, d=d, biological_direction=biological)


def classifier() -> LogisticRegression:
    return LogisticRegression(
        C=1.0,
        class_weight="balanced",
        solver="liblinear",
        max_iter=10_000,
    )


def internal_predictions(data: Simulated, residualised: bool) -> np.ndarray:
    splitter = StratifiedKFold(5, shuffle=True, random_state=71)
    prediction = np.full(len(data.y), np.nan)
    for train, test in splitter.split(data.x, data.y):
        train_x = data.x[train]
        test_x = data.x[test]
        if residualised:
            design_train = data.d[train, None].astype(float)
            design_test = data.d[test, None].astype(float)
            residualiser = Ridge(alpha=1.0).fit(design_train, train_x)
            train_x = train_x - residualiser.predict(design_train)
            test_x = test_x - residualiser.predict(design_test)
        scaler = StandardScaler().fit(train_x)
        model = classifier().fit(scaler.transform(train_x), data.y[train])
        prediction[test] = model.predict_proba(scaler.transform(test_x))[:, 1]
    return prediction


def restricted_auc(data: Simulated) -> tuple[float, int]:
    candidates = []
    for design_value in (0, 1):
        indices = np.flatnonzero(data.d == design_value)
        labels = data.y[indices]
        if len(np.unique(labels)) < 2 or np.bincount(labels).min() < 2:
            continue
        candidates.append(
            (
                int(np.bincount(labels).min()),
                len(indices),
                -design_value,
                indices,
            )
        )
    if not candidates:
        return np.nan, 0

    # Evaluate one observed overlap stratum. Probabilities from independently
    # fitted stratum-specific models are not commensurate and must not be pooled.
    indices = max(candidates)[-1]
    labels = data.y[indices]
    prediction = np.full(len(indices), np.nan)
    folds = min(5, int(np.bincount(labels).min()))
    splitter = StratifiedKFold(folds, shuffle=True, random_state=73)
    for train_local, test_local in splitter.split(data.x[indices], labels):
        train = indices[train_local]
        test = indices[test_local]
        scaler = StandardScaler().fit(data.x[train])
        model = classifier().fit(
            scaler.transform(data.x[train]), data.y[train]
        )
        prediction[test_local] = model.predict_proba(
            scaler.transform(data.x[test])
        )[:, 1]
    return float(roc_auc_score(labels, prediction)), len(indices)


def external_auc(source: Simulated, target: Simulated) -> float:
    scaler = StandardScaler().fit(source.x)
    model = classifier().fit(scaler.transform(source.x), source.y)
    prediction = model.predict_proba(scaler.transform(target.x))[:, 1]
    return float(roc_auc_score(target.y, prediction))


def coefficient_error(data: Simulated, beta: float) -> float:
    y_centered = 2.0 * data.y - 1.0
    d_centered = 2.0 * data.d - 1.0
    design = np.column_stack(
        [np.ones(len(data.y)), y_centered, d_centered]
    )
    if np.linalg.matrix_rank(design) < design.shape[1]:
        return np.nan
    coefficients = np.linalg.pinv(design) @ data.x
    estimated = coefficients[1]
    truth = beta * data.biological_direction
    return float(np.linalg.norm(estimated - truth))


def run_one(beta: float, delta: float, rho: float, replicate: int) -> dict:
    seed = (
        BASE_SEED
        + int(beta * 1000) * 1_000_000
        + int(delta * 1000) * 10_000
        + int(rho * 1000) * 100
        + replicate
    )
    rng = np.random.default_rng(seed)
    source = generate(rng, N_DONORS, beta, delta, rho)
    target = generate(rng, N_TARGET, beta, delta, 0.0)
    unadjusted = internal_predictions(source, residualised=False)
    residualised = internal_predictions(source, residualised=True)
    restricted, retained = restricted_auc(source)
    info = information_fraction(source.y, source.d[:, None].astype(float))
    design_auc = roc_auc_score(source.y, source.d)
    design_auc = max(float(design_auc), float(1.0 - design_auc))
    return {
        "beta": beta,
        "delta": delta,
        "rho": rho,
        "replicate": replicate,
        "seed": seed,
        "n_source": N_DONORS,
        "n_restricted": retained,
        "design_only_auc": design_auc,
        **info,
        "auc_internal_unadjusted": roc_auc_score(source.y, unadjusted),
        "auc_internal_residualised": roc_auc_score(source.y, residualised),
        "auc_internal_restricted": restricted,
        "auc_external_source_only": external_auc(source, target),
        "coefficient_l2_error": coefficient_error(source, beta),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [
        (beta, delta, rho, replicate)
        for beta in BETAS
        for delta in DELTAS
        for rho in RHOS
        for replicate in range(1, N_REPLICATES + 1)
    ]
    rows = Parallel(n_jobs=N_JOBS, prefer="processes", verbose=5)(
        delayed(run_one)(*job) for job in jobs
    )
    raw = pd.DataFrame(rows)
    raw.to_csv(
        OUT / "simulation_replicates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    measures = [
        "design_only_auc",
        "information_fraction",
        "variance_inflation",
        "n_restricted",
        "auc_internal_unadjusted",
        "auc_internal_residualised",
        "auc_internal_restricted",
        "auc_external_source_only",
        "coefficient_l2_error",
    ]
    summary = (
        raw.groupby(["beta", "delta", "rho"], as_index=False)[measures]
        .agg(["mean", "std", "median"])
    )
    summary.columns = [
        "_".join(part for part in column if part)
        for column in summary.columns.to_flat_index()
    ]
    summary.to_csv(OUT / "simulation_summary.tsv", sep="\t", index=False)

    # A source-observational equivalence demonstration at D = Y.
    rng = np.random.default_rng(BASE_SEED + 999)
    y, d = labels_and_design(rng, 10_000, 1.0)
    shared = rng.normal(size=(10_000, 1))
    model_biological = 1.5 * (2 * y - 1)[:, None] + shared
    model_technical = 1.5 * (2 * d - 1)[:, None] + shared
    equivalence = {
        "max_absolute_difference": float(
            np.max(np.abs(model_biological - model_technical))
        ),
        "statement": (
            "At D=Y, beta=1.5/delta=0 and beta=0/delta=1.5 "
            "produce identical observed source data when the effect direction "
            "and noise are shared."
        ),
    }
    write_json(OUT / "observational_equivalence.json", equivalence)
    write_json(
        OUT / "simulation_manifest.json",
        {
            "betas": BETAS,
            "deltas": DELTAS,
            "rhos": RHOS,
            "n_replicates": N_REPLICATES,
            "n_source": N_DONORS,
            "n_target": N_TARGET,
            "n_features": N_FEATURES,
            "base_seed": BASE_SEED,
            "classifier": "balanced liblinear logistic regression, C=1",
            "residualiser": "ridge alpha=1 fitted in each outer training fold",
            "target_design_label_rho": 0.0,
        },
    )
    print(summary.head(20).to_string(index=False))
    print(f"Wrote {len(raw):,} simulation replicates to {OUT}")


if __name__ == "__main__":
    main()
