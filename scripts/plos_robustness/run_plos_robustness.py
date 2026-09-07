#!/usr/bin/env python3
"""PLOS-focused robustness analyses for the design-identifiability study."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    RandomForestClassifier,
    RandomForestRegressor,
)
from sklearn.linear_model import Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "design_validity"))
sys.path.insert(0, str(ROOT / "identifiability_extension" / "scripts"))

import locked_validity_audit as audit  # noqa: E402
import run_identifiability_simulation as simulation  # noqa: E402


OUT = ROOT / "results" / "plos_robustness"
N_JOBS = min(int(os.environ.get("N_JOBS", "4")), os.cpu_count() or 1)
FOREST_TREES = int(os.environ.get("FOREST_TREES", "96"))
BOOTSTRAPS = int(os.environ.get("BOOTSTRAPS", "5000"))
LEAKAGE_REPLICATES = int(os.environ.get("LEAKAGE_REPLICATES", "100"))
BASE_SEED = 2026090100


def class_weights(y: np.ndarray) -> np.ndarray:
    counts = np.bincount(y.astype(int), minlength=2)
    return np.asarray([len(y) / (2.0 * counts[value]) for value in y])


def select_features(
    rep: str,
    train_x: np.ndarray,
    test_x: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if rep == "geneformer":
        return train_x, test_x
    n_keep = min(audit.N_HVG, train_x.shape[1])
    keep = np.argsort(train_x.var(axis=0), kind="stable")[::-1][:n_keep]
    return train_x[:, keep], test_x[:, keep]


def ridge_residuals(
    train_x: np.ndarray,
    test_x: np.ndarray,
    train_design: np.ndarray,
    test_design: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    model = Ridge(alpha=1.0, fit_intercept=True).fit(train_design, train_x)
    return train_x - model.predict(train_design), test_x - model.predict(test_design)


def forest_residuals(
    train_x: np.ndarray,
    test_x: np.ndarray,
    train_design: np.ndarray,
    test_design: np.ndarray,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    model = RandomForestRegressor(
        n_estimators=FOREST_TREES,
        max_depth=4,
        min_samples_leaf=3,
        max_features=0.8,
        bootstrap=True,
        random_state=seed,
        n_jobs=1,
    ).fit(train_design, train_x)
    return train_x - model.predict(train_design), test_x - model.predict(test_design)


def batch_centered(
    train_x: np.ndarray,
    test_x: np.ndarray,
    train_batch: np.ndarray,
    test_batch: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Training-only location adjustment analogous to the location arm of ComBat."""

    grand = train_x.mean(axis=0)
    effects = {
        value: train_x[train_batch == value].mean(axis=0) - grand
        for value in np.unique(train_batch)
    }
    train_effect = np.vstack([effects[value] for value in train_batch])
    zero = np.zeros(train_x.shape[1], dtype=float)
    test_effect = np.vstack([effects.get(value, zero) for value in test_batch])
    return train_x - train_effect, test_x - test_effect


def overlap_weighted_prediction(
    rep: str,
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
    train_design: np.ndarray,
    seed: int,
) -> np.ndarray:
    train_z, test_z = audit.prepare_representation(rep, train_x, test_x)
    propensity_c = audit.select_c_prepared(train_design, train_y, seed)
    propensity = (
        audit.lr(propensity_c)
        .fit(train_design, train_y)
        .predict_proba(train_design)[:, 1]
    )
    propensity = np.clip(propensity, 0.025, 0.975)
    weights = np.where(train_y == 1, 1.0 - propensity, propensity)
    weights /= weights.mean()
    classifier_c = audit.select_c_prepared(train_z, train_y, seed)
    model = audit.lr(classifier_c).fit(
        train_z,
        train_y,
        sample_weight=weights,
    )
    return model.predict_proba(test_z)[:, 1]


def evaluate_residualisation_repeat(repeat: int) -> tuple[list[dict], list[dict]]:
    cohort = audit.load_cohort_135779()
    seed = audit.SEEDS[repeat]
    splitter = StratifiedKFold(
        audit.N_SPLITS,
        shuffle=True,
        random_state=seed,
    )
    methods = (
        "unadjusted",
        "ridge_full_design",
        "random_forest_full_design",
        "batch_location_adjustment",
        "overlap_weighted_full_design",
    )
    predictions = {
        (rep, method): np.full(len(cohort.y), np.nan)
        for rep in audit.REPS
        for method in methods
    }
    folds = np.zeros(len(cohort.y), dtype=int)
    raw = {
        "geneformer": cohort.geneformer,
        "hvg_pseudobulk": cohort.pseudobulk,
        "pca_pseudobulk": cohort.pseudobulk,
    }

    for fold, (train, test) in enumerate(
        splitter.split(cohort.geneformer, cohort.y),
        start=1,
    ):
        full_train, full_test = audit.make_design(
            cohort.metadata,
            train,
            test,
            *audit.blocks_for(cohort)["all"],
        )
        train_batch = cohort.metadata.iloc[train]["batch"].astype(str).to_numpy()
        test_batch = cohort.metadata.iloc[test]["batch"].astype(str).to_numpy()

        selected: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for rep in audit.REPS:
            selected[rep] = select_features(
                rep,
                raw[rep][train],
                raw[rep][test],
            )

        residual_cache: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]] = {}
        for source in ("geneformer", "pseudobulk"):
            rep = "geneformer" if source == "geneformer" else "hvg_pseudobulk"
            train_x, test_x = selected[rep]
            residual_cache[(source, "ridge_full_design")] = ridge_residuals(
                train_x,
                test_x,
                full_train,
                full_test,
            )
            residual_cache[(source, "random_forest_full_design")] = forest_residuals(
                train_x,
                test_x,
                full_train,
                full_test,
                seed + fold,
            )
            residual_cache[(source, "batch_location_adjustment")] = batch_centered(
                train_x,
                test_x,
                train_batch,
                test_batch,
            )

        for rep in audit.REPS:
            train_x, test_x = selected[rep]
            predictions[(rep, "unadjusted")][test] = audit.predict_representation(
                rep,
                train_x,
                cohort.y[train],
                test_x,
                seed,
            )
            source = "geneformer" if rep == "geneformer" else "pseudobulk"
            for method in (
                "ridge_full_design",
                "random_forest_full_design",
                "batch_location_adjustment",
            ):
                adjusted_train, adjusted_test = residual_cache[(source, method)]
                predictions[(rep, method)][test] = audit.predict_representation(
                    rep,
                    adjusted_train,
                    cohort.y[train],
                    adjusted_test,
                    seed,
                )
            predictions[(rep, "overlap_weighted_full_design")][test] = (
                overlap_weighted_prediction(
                    rep,
                    train_x,
                    cohort.y[train],
                    test_x,
                    full_train,
                    seed,
                )
            )
        folds[test] = fold

    metric_rows = []
    prediction_rows = []
    for (rep, method), probability in predictions.items():
        metric_rows.append(
            {
                "repeat": repeat + 1,
                "seed": seed,
                "representation": rep,
                "adjustment": method,
                "n_donors": len(cohort.y),
                "n_cases": int(cohort.y.sum()),
                "roc_auc": float(roc_auc_score(cohort.y, probability)),
            }
        )
        prediction_rows.extend(
            {
                "repeat": repeat + 1,
                "seed": seed,
                "representation": rep,
                "adjustment": method,
                "donor_id": donor,
                "fold": int(fold),
                "y_true": int(label),
                "prob_case": float(score),
            }
            for donor, fold, label, score in zip(
                cohort.donors,
                folds,
                cohort.y,
                probability,
            )
        )
    return metric_rows, prediction_rows


def nonlinear_design_only() -> tuple[pd.DataFrame, pd.DataFrame]:
    cohort = audit.load_cohort_135779()
    numeric, categorical = audit.blocks_for(cohort)["all"]
    metric_rows = []
    prediction_rows = []
    for repeat, seed in enumerate(audit.SEEDS, start=1):
        splitter = StratifiedKFold(
            audit.N_SPLITS,
            shuffle=True,
            random_state=seed,
        )
        predictions = {
            name: np.full(len(cohort.y), np.nan)
            for name in ("logistic", "random_forest", "gradient_boosting")
        }
        folds = np.zeros(len(cohort.y), dtype=int)
        for fold, (train, test) in enumerate(
            splitter.split(cohort.geneformer, cohort.y),
            start=1,
        ):
            train_x, test_x = audit.make_design(
                cohort.metadata,
                train,
                test,
                numeric,
                categorical,
            )
            C = audit.select_c_prepared(train_x, cohort.y[train], seed)
            predictions["logistic"][test] = (
                audit.lr(C)
                .fit(train_x, cohort.y[train])
                .predict_proba(test_x)[:, 1]
            )
            forest = RandomForestClassifier(
                n_estimators=256,
                max_depth=4,
                min_samples_leaf=3,
                max_features="sqrt",
                class_weight="balanced",
                random_state=seed + fold,
                n_jobs=1,
            ).fit(train_x, cohort.y[train])
            predictions["random_forest"][test] = forest.predict_proba(test_x)[:, 1]
            boosting = HistGradientBoostingClassifier(
                learning_rate=0.05,
                max_iter=150,
                max_leaf_nodes=7,
                min_samples_leaf=5,
                l2_regularization=1.0,
                random_state=seed + fold,
            ).fit(
                train_x,
                cohort.y[train],
                sample_weight=class_weights(cohort.y[train]),
            )
            predictions["gradient_boosting"][test] = boosting.predict_proba(
                test_x
            )[:, 1]
            folds[test] = fold

        for model, probability in predictions.items():
            metric_rows.append(
                {
                    "repeat": repeat,
                    "seed": seed,
                    "model": model,
                    "roc_auc": float(roc_auc_score(cohort.y, probability)),
                }
            )
            prediction_rows.extend(
                {
                    "repeat": repeat,
                    "seed": seed,
                    "model": model,
                    "donor_id": donor,
                    "fold": int(fold),
                    "y_true": int(label),
                    "prob_case": float(score),
                }
                for donor, fold, label, score in zip(
                    cohort.donors,
                    folds,
                    cohort.y,
                    probability,
                )
            )
    return pd.DataFrame(metric_rows), pd.DataFrame(prediction_rows)


def batch_exposure_negative_control() -> pd.DataFrame:
    design = pd.read_csv(
        ROOT / "results" / "locked_validity_audit" / "summary.tsv",
        sep="\t",
    )
    batch = pd.read_csv(
        ROOT
        / "results"
        / "locked_validity_audit"
        / "representation_batch_predictability_summary.tsv",
        sep="\t",
    )
    batch = batch[batch["cohort"].eq("GSE135779")]
    rows = []
    for rep in audit.REPS:
        unadjusted = design[
            design["cohort"].eq("GSE135779")
            & design["method"].eq(rep)
            & design["adjustment"].eq("unadjusted")
        ]["auc_mean"].iloc[0]
        residual = design[
            design["cohort"].eq("GSE135779")
            & design["method"].eq(rep)
            & design["adjustment"].eq("residual_batch")
        ]["auc_mean"].iloc[0]
        rows.append(
            {
                "representation": rep,
                "best_batch_auc": float(
                    batch[batch["method"].eq(rep)]["auc_mean"].max()
                ),
                "batch_only_disease_auc": float(
                    design[
                        design["cohort"].eq("GSE135779")
                        & design["method"].eq("covariates_batch")
                    ]["auc_mean"].iloc[0]
                ),
                "unadjusted_disease_auc": float(unadjusted),
                "batch_residualised_disease_auc": float(residual),
                "delta_batch_residualised_minus_unadjusted": float(
                    residual - unadjusted
                ),
            }
        )
    return pd.DataFrame(rows)


def paired_split_bootstrap() -> pd.DataFrame:
    metrics = pd.read_csv(
        ROOT / "results" / "locked_validity_audit" / "repeat_metrics.tsv",
        sep="\t",
    )
    headline = pd.read_csv(ROOT / "results" / "HEADLINE_stratification_table.tsv", sep="\t")
    rng = np.random.default_rng(BASE_SEED + 500)
    rows = []
    for rep in audit.REPS:
        selected = metrics[
            metrics["cohort"].eq("GSE174188_CD4")
            & metrics["method"].eq(rep)
            & metrics["adjustment"].isin(["unadjusted", "residual_batch"])
        ].pivot(index="repeat", columns="adjustment", values="roc_auc")
        residual_loss = (
            selected["unadjusted"] - selected["residual_batch"]
        ).to_numpy()
        for stratum, column in (
            ("dominant_wave_4", "delta_n89"),
            ("pure_wave_4", "delta_n66"),
        ):
            restriction_loss = -float(
                headline.loc[
                    headline["representation"].eq(rep),
                    column,
                ].iloc[0]
            )
            draws = np.empty(BOOTSTRAPS)
            for index in range(BOOTSTRAPS):
                sample = rng.choice(residual_loss, len(residual_loss), replace=True)
                draws[index] = sample.mean() - restriction_loss
            estimate = float(residual_loss.mean() - restriction_loss)
            rows.append(
                {
                    "representation": rep,
                    "restriction_stratum": stratum,
                    "residualisation_loss_mean": float(residual_loss.mean()),
                    "matched_restriction_loss": restriction_loss,
                    "loss_difference": estimate,
                    "loss_difference_bootstrap_p025": float(
                        np.quantile(draws, 0.025)
                    ),
                    "loss_difference_bootstrap_p975": float(
                        np.quantile(draws, 0.975)
                    ),
                    "loss_ratio": float(
                        residual_loss.mean() / max(restriction_loss, 1e-12)
                    ),
                    "bootstrap_unit": "paired repeated-split residualisation losses",
                }
            )
    return pd.DataFrame(rows)


def positive_control_protocol() -> pd.DataFrame:
    raw = pd.read_csv(
        ROOT
        / "identifiability_extension"
        / "results"
        / "simulation"
        / "simulation_replicates.tsv.gz",
        sep="\t",
    )
    selected = raw[
        raw["beta"].eq(1.0)
        & raw["delta"].eq(0.5)
        & raw["rho"].eq(0.25)
    ]
    means = selected.mean(numeric_only=True)
    rng = np.random.default_rng(BASE_SEED + 750)
    positive = simulation.generate(
        rng,
        n=240,
        beta=1.0,
        delta=0.5,
        rho=0.25,
    )
    observed_prediction = simulation.internal_predictions(
        positive,
        residualised=False,
    )
    observed_auc = roc_auc_score(positive.y, observed_prediction)
    null_aucs = []
    for permutation in range(250):
        permuted = simulation.Simulated(
            x=positive.x,
            y=rng.permutation(positive.y),
            d=positive.d,
            biological_direction=positive.biological_direction,
        )
        null_aucs.append(
            roc_auc_score(
                permuted.y,
                simulation.internal_predictions(
                    permuted,
                    residualised=False,
                ),
            )
        )
    permutation_p = (1 + np.sum(np.asarray(null_aucs) >= observed_auc)) / 251
    checks = [
        (
            "design_label_estimability",
            float(means["information_fraction"]),
            "information_fraction >= 0.80",
            means["information_fraction"] >= 0.80,
        ),
        (
            "complete_pipeline_permutation_null",
            float(permutation_p),
            "empirical upper-tail p <= 0.01",
            permutation_p <= 0.01,
        ),
        (
            "representation_to_design_exposure",
            float(means["design_only_auc"]),
            "design_only_auc <= 0.70",
            means["design_only_auc"] <= 0.70,
        ),
        (
            "adjustment_and_overlap_stability",
            float(
                max(
                    means["auc_internal_unadjusted"]
                    - means["auc_internal_residualised"],
                    means["auc_internal_unadjusted"]
                    - means["auc_internal_restricted"],
                )
            ),
            "maximum AUC loss <= 0.05",
            max(
                means["auc_internal_unadjusted"]
                - means["auc_internal_residualised"],
                means["auc_internal_unadjusted"]
                - means["auc_internal_restricted"],
            )
            <= 0.05,
        ),
        (
            "source_only_external_stability",
            float(
                means["auc_internal_unadjusted"]
                - means["auc_external_source_only"]
            ),
            "internal minus external AUC <= 0.05",
            (
                means["auc_internal_unadjusted"]
                - means["auc_external_source_only"]
            )
            <= 0.05,
        ),
    ]
    return pd.DataFrame(
        {
            "simulation_setting": "beta=1; delta=0.5; rho=0.25; 100 replicates",
            "check": [item[0] for item in checks],
            "value": [item[1] for item in checks],
            "illustrative_pass_rule": [item[2] for item in checks],
            "passed": [bool(item[3]) for item in checks],
        }
    )


def leakage_one(replicate: int) -> list[dict]:
    seed = BASE_SEED + 10_000 + replicate
    rng = np.random.default_rng(seed)
    data = simulation.generate(
        rng,
        n=240,
        beta=0.0,
        delta=1.0,
        rho=0.75,
    )
    splitter = StratifiedKFold(5, shuffle=True, random_state=seed)
    predictions = {
        "raw_random_forest": np.full(len(data.y), np.nan),
        "global_residualisation_random_forest": np.full(len(data.y), np.nan),
        "fold_contained_residualisation_random_forest": np.full(len(data.y), np.nan),
    }
    global_model = Ridge(alpha=1.0).fit(data.d[:, None], data.x)
    global_residual = data.x - global_model.predict(data.d[:, None])
    for fold, (train, test) in enumerate(splitter.split(data.x, data.y), start=1):
        fold_model = Ridge(alpha=1.0).fit(data.d[train, None], data.x[train])
        fold_train = data.x[train] - fold_model.predict(data.d[train, None])
        fold_test = data.x[test] - fold_model.predict(data.d[test, None])
        variants = {
            "raw_random_forest": (data.x[train], data.x[test]),
            "global_residualisation_random_forest": (
                global_residual[train],
                global_residual[test],
            ),
            "fold_contained_residualisation_random_forest": (
                fold_train,
                fold_test,
            ),
        }
        for name, (train_x, test_x) in variants.items():
            model = RandomForestClassifier(
                n_estimators=160,
                max_depth=5,
                min_samples_leaf=3,
                max_features="sqrt",
                class_weight="balanced",
                random_state=seed + fold,
                n_jobs=1,
            ).fit(train_x, data.y[train])
            predictions[name][test] = model.predict_proba(test_x)[:, 1]
    return [
        {
            "replicate": replicate,
            "seed": seed,
            "scenario": name,
            "roc_auc": float(roc_auc_score(data.y, probability)),
        }
        for name, probability in predictions.items()
    ]


def summarise(
    frame: pd.DataFrame,
    groups: list[str],
    measure: str = "roc_auc",
) -> pd.DataFrame:
    return (
        frame.groupby(groups, as_index=False)[measure]
        .agg(
            mean="mean",
            sd="std",
            p025=lambda values: np.quantile(values, 0.025),
            p975=lambda values: np.quantile(values, 0.975),
            n="count",
        )
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    completed = Parallel(n_jobs=N_JOBS, prefer="processes", verbose=5)(
        delayed(evaluate_residualisation_repeat)(repeat)
        for repeat in range(len(audit.SEEDS))
    )
    residual_metrics = pd.DataFrame(
        [row for metrics, _ in completed for row in metrics]
    )
    residual_predictions = pd.DataFrame(
        [row for _, predictions in completed for row in predictions]
    )
    residual_metrics.to_csv(
        OUT / "residualisation_sensitivity_metrics.tsv",
        sep="\t",
        index=False,
    )
    residual_predictions.to_csv(
        OUT / "residualisation_sensitivity_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    residual_summary = summarise(
        residual_metrics,
        ["representation", "adjustment"],
    )
    residual_summary.to_csv(
        OUT / "residualisation_sensitivity_summary.tsv",
        sep="\t",
        index=False,
    )

    design_metrics, design_predictions = nonlinear_design_only()
    design_metrics.to_csv(
        OUT / "nonlinear_design_only_metrics.tsv",
        sep="\t",
        index=False,
    )
    design_predictions.to_csv(
        OUT / "nonlinear_design_only_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    summarise(design_metrics, ["model"]).to_csv(
        OUT / "nonlinear_design_only_summary.tsv",
        sep="\t",
        index=False,
    )

    batch_exposure_negative_control().to_csv(
        OUT / "batch_exposure_negative_control.tsv",
        sep="\t",
        index=False,
    )
    paired_split_bootstrap().to_csv(
        OUT / "attenuation_difference_bootstrap.tsv",
        sep="\t",
        index=False,
    )
    positive_control_protocol().to_csv(
        OUT / "synthetic_positive_control.tsv",
        sep="\t",
        index=False,
    )

    leakage_rows = Parallel(n_jobs=N_JOBS, prefer="processes", verbose=5)(
        delayed(leakage_one)(replicate)
        for replicate in range(1, LEAKAGE_REPLICATES + 1)
    )
    leakage = pd.DataFrame([row for group in leakage_rows for row in group])
    leakage.to_csv(
        OUT / "confound_leakage_replicates.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    summarise(leakage, ["scenario"]).to_csv(
        OUT / "confound_leakage_summary.tsv",
        sep="\t",
        index=False,
    )

    manifest = {
        "analysis": "PLOS-focused robustness extension",
        "seeds": list(audit.SEEDS),
        "outer_folds": audit.N_SPLITS,
        "forest_trees": FOREST_TREES,
        "bootstrap_replicates": BOOTSTRAPS,
        "confound_leakage_replicates": LEAKAGE_REPLICATES,
        "nonlinear_residualiser": (
            "RandomForestRegressor fitted only on outer-training donors; "
            "max_depth=4, min_samples_leaf=3."
        ),
        "batch_location_adjustment": (
            "Training-only batch location centring; no empirical-Bayes scale "
            "shrinkage and therefore labelled ComBat-style rather than ComBat."
        ),
        "propensity_adjustment": (
            "Training-only logistic propensity model with clipped overlap weights."
        ),
        "bootstrap_scope": (
            "Paired resampling of repeat-level residualisation losses; matched "
            "restriction discrepancies are locked point estimates. Intervals are "
            "split-sensitivity intervals, not population confidence intervals."
        ),
    }
    (OUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    print("\nResidualisation sensitivity")
    print(residual_summary.to_string(index=False))
    print("\nNonlinear design-only prediction")
    print(
        pd.read_csv(OUT / "nonlinear_design_only_summary.tsv", sep="\t").to_string(
            index=False
        )
    )
    print("\nAttenuation differences")
    print(
        pd.read_csv(
            OUT / "attenuation_difference_bootstrap.tsv",
            sep="\t",
        ).to_string(index=False)
    )
    print(f"\nWrote PLOS robustness outputs to {OUT}")


if __name__ == "__main__":
    main()
