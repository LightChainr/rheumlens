#!/usr/bin/env python3
"""Locked two-cohort validity audit for donor-level SLE classification.

This script puts unadjusted representation performance, design-covariate
predictability, fold-contained residualisation, and matched design restriction
under one implementation. It covers GSE174188 and the restored GSE135779
childhood cohort. All outer test donors remain untouched by feature selection,
scaling, PCA, covariate encoding, residualisation, and classifier fitting.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


HERE = Path(__file__).resolve()
WORKSPACE = HERE.parents[1]
PROJECT = HERE.parents[2]
EXTRA = (
    PROJECT
    / "03_远程回传"
    / "final_gpu_figure_ready_20260709"
    / "extracted_plus"
    / "RheumLens_GPU_figure_ready_plus_20260709"
    / "extra"
)
OUT = WORKSPACE / "results" / "locked_validity_audit"
META135 = (
    WORKSPACE
    / "results"
    / "gse135779_metadata"
    / "gse135779_donor_metadata_restored.tsv"
)
META174 = (
    PROJECT
    / "05_patient_classifier_generalization_20260716"
    / "results"
    / "01_inputs"
    / "gse174188_final_donor_covariates.tsv"
)

SEEDS = list(range(20260801, 20260821))
BATCH_SEEDS = SEEDS[:5]
CS = np.logspace(-4, 4, 9)
N_SPLITS = 5
N_HVG = 4000
N_PC = 30
RIDGE_ALPHA = 1.0
N_JOBS = min(int(os.environ.get("N_JOBS", "6")), os.cpu_count() or 1)
REPS = ("geneformer", "hvg_pseudobulk", "pca_pseudobulk")

sys.path.insert(
    0,
    str(PROJECT / "14_learned_pooling_20260726" / "scripts"),
)
import protocol  # noqa: E402


@dataclass
class Cohort:
    name: str
    donors: list[str]
    y: np.ndarray
    geneformer: np.ndarray
    pseudobulk: np.ndarray
    metadata: pd.DataFrame


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def lr(C: float) -> LogisticRegression:
    return LogisticRegression(
        C=C,
        solver="liblinear",
        class_weight="balanced",
        max_iter=20_000,
    )


def effective_splits(y: np.ndarray, requested: int = N_SPLITS) -> int:
    counts = np.bincount(np.asarray(y, dtype=int))
    positive = counts[counts > 0]
    if len(positive) < 2:
        raise ValueError("Classification target has fewer than two classes")
    return int(min(requested, positive.min()))


def select_c_prepared(x: np.ndarray, y: np.ndarray, seed: int) -> float:
    """Select C inside the outer training set after outer-fold preprocessing.

    Scaling and any HVG/PCA construction are fit on the full outer training set,
    never on outer test donors. The inner resampling selects only C conditional
    on that outer-training representation.
    """

    splitter = StratifiedKFold(
        effective_splits(y), shuffle=True, random_state=seed
    )
    splits = list(splitter.split(x, y))
    means = []
    for C in CS:
        fold_auc = []
        for train, test in splits:
            model = lr(float(C)).fit(x[train], y[train])
            fold_auc.append(
                protocol.roc_auc(y[test], model.predict_proba(x[test])[:, 1])
            )
        means.append(float(np.mean(fold_auc)))
    return float(CS[int(np.argmax(means))])


def prepare_representation(
    rep: str,
    train_x: np.ndarray,
    test_x: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if rep != "geneformer":
        variance = train_x.var(axis=0, dtype=np.float64)
        n_keep = min(N_HVG, train_x.shape[1])
        keep = np.argsort(variance, kind="stable")[::-1][:n_keep]
        train_x = train_x[:, keep]
        test_x = test_x[:, keep]

    scaler = StandardScaler().fit(train_x)
    train_z = scaler.transform(train_x)
    test_z = scaler.transform(test_x)

    if rep == "pca_pseudobulk":
        n_components = min(N_PC, train_z.shape[0] - 1, train_z.shape[1])
        pca = PCA(n_components=n_components, random_state=0).fit(train_z)
        train_z = pca.transform(train_z)
        test_z = pca.transform(test_z)
    return np.asarray(train_z), np.asarray(test_z)


def predict_representation(
    rep: str,
    train_x: np.ndarray,
    train_y: np.ndarray,
    test_x: np.ndarray,
    seed: int,
) -> np.ndarray:
    train_z, test_z = prepare_representation(rep, train_x, test_x)
    C = select_c_prepared(train_z, train_y, seed)
    return lr(C).fit(train_z, train_y).predict_proba(test_z)[:, 1]


def design_transformer(
    numeric: list[str], categorical: list[str]
) -> ColumnTransformer:
    transformers = []
    if numeric:
        transformers.append(
            (
                "numeric",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="median")),
                        ("scale", StandardScaler()),
                    ]
                ),
                numeric,
            )
        )
    if categorical:
        transformers.append(
            (
                "categorical",
                Pipeline(
                    [
                        ("impute", SimpleImputer(strategy="most_frequent")),
                        (
                            "onehot",
                            OneHotEncoder(
                                handle_unknown="ignore", sparse_output=False
                            ),
                        ),
                    ]
                ),
                categorical,
            )
        )
    return ColumnTransformer(transformers, remainder="drop", sparse_threshold=0.0)


def make_design(
    metadata: pd.DataFrame,
    train: np.ndarray,
    test: np.ndarray,
    numeric: list[str],
    categorical: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    transformer = design_transformer(numeric, categorical)
    train_z = transformer.fit_transform(metadata.iloc[train])
    test_z = transformer.transform(metadata.iloc[test])
    return np.asarray(train_z, dtype=float), np.asarray(test_z, dtype=float)


def residualize(
    train_x: np.ndarray,
    test_x: np.ndarray,
    train_covariates: np.ndarray,
    test_covariates: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    model = Ridge(alpha=RIDGE_ALPHA, fit_intercept=True)
    model.fit(train_covariates, train_x)
    return (
        train_x - model.predict(train_covariates),
        test_x - model.predict(test_covariates),
    )


def load_cohort_135779() -> Cohort:
    metadata = pd.read_csv(META135, sep="\t", dtype={"donor_id": str})
    gf_path = (
        EXTRA
        / "04_models"
        / "Geneformer"
        / "SLE_GSE135779"
        / "geneformer_v2_316m_cell_sample500_clspool_logistic_maxlen4096_seed001"
        / "donor_embedding.parquet"
    )
    pb_path = EXTRA / "pseudobulk" / "SLE_GSE135779" / "donor_log1p_cpm.parquet"
    gf = pd.read_parquet(gf_path)
    pb = pd.read_parquet(pb_path)
    gf.index = gf.index.astype(str)
    pb.index = pb.index.astype(str)
    donors = metadata["donor_id"].astype(str).tolist()
    if not set(donors).issubset(gf.index) or not set(donors).issubset(pb.index):
        raise RuntimeError("GSE135779 representation donor mapping is incomplete")
    metadata = metadata.set_index("donor_id").loc[donors].reset_index()
    return Cohort(
        name="GSE135779",
        donors=donors,
        y=metadata["y_true"].to_numpy(dtype=int),
        geneformer=gf.loc[donors].to_numpy(dtype=np.float64),
        pseudobulk=pb.loc[donors].to_numpy(dtype=np.float32),
        metadata=metadata,
    )


def load_cohort_174188() -> Cohort:
    metadata = pd.read_csv(META174, sep="\t", dtype={"donor_id": str})
    gf_path = (
        EXTRA
        / "04_models"
        / "Geneformer"
        / "SLE_GSE174188_CD4"
        / "geneformer_v2_316m_cell_sample1000_clspool_logistic_maxlen4096_seed001"
        / "donor_embedding.parquet"
    )
    pb_path = (
        EXTRA
        / "pseudobulk"
        / "SLE_GSE174188_CD4"
        / "donor_log1p_cpm.parquet"
    )
    gf = pd.read_parquet(gf_path)
    pb = pd.read_parquet(pb_path)
    gf.index = gf.index.astype(str)
    pb.index = pb.index.astype(str)
    donors = metadata["donor_id"].astype(str).tolist()
    if not set(donors).issubset(gf.index) or not set(donors).issubset(pb.index):
        raise RuntimeError("GSE174188 representation donor mapping is incomplete")
    metadata = metadata.set_index("donor_id").loc[donors].reset_index()
    processing = [c for c in metadata if c.startswith("processing_fraction_")]
    metadata["batch"] = (
        metadata[processing].to_numpy().argmax(axis=1) + 1
    ).astype(str)
    metadata["log_cells_per_donor"] = np.log1p(metadata["cells_per_donor"])
    metadata["log_mean_umi_per_cell"] = np.log1p(
        metadata["mean_umi_per_cell"]
    )
    metadata["log_mean_genes_per_cell"] = np.log1p(
        metadata["mean_genes_per_cell"]
    )
    return Cohort(
        name="GSE174188_CD4",
        donors=donors,
        y=metadata["y_true"].to_numpy(dtype=int),
        geneformer=gf.loc[donors].to_numpy(dtype=np.float64),
        pseudobulk=pb.loc[donors].to_numpy(dtype=np.float32),
        metadata=metadata,
    )


def blocks_for(cohort: Cohort) -> dict[str, tuple[list[str], list[str]]]:
    if cohort.name == "GSE135779":
        qc = [
            "log_estimated_cells",
            "log_mean_reads_per_cell",
            "log_median_genes_per_cell",
            "sequencing_saturation",
            "fraction_reads_in_cells",
            "log_median_umi_per_cell",
        ]
        demographic_numeric = ["age_years"]
        demographic_categorical = ["sex", "race", "ethnicity"]
        return {
            "batch": ([], ["batch"]),
            "collection_year": ([], ["collection_year"]),
            "qc": (qc, []),
            "demographic": (demographic_numeric, demographic_categorical),
            "all": (
                qc + demographic_numeric,
                ["batch", "collection_year", *demographic_categorical],
            ),
        }

    processing = [
        c for c in cohort.metadata if c.startswith("processing_fraction_")
    ]
    qc = [
        "log_cells_per_donor",
        "log_mean_umi_per_cell",
        "log_mean_genes_per_cell",
        "aggregate_pct_mito",
    ]
    return {
        "batch": (processing, []),
        "qc": (qc, []),
        "demographic": (["age_years"], ["sex"]),
        "all": (processing + qc + ["age_years"], ["sex"]),
    }


def evaluate_repeat(
    cohort: Cohort,
    repeat: int,
    residual_blocks: tuple[str, ...],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    seed = SEEDS[repeat]
    splitter = StratifiedKFold(N_SPLITS, shuffle=True, random_state=seed)
    probabilities: dict[tuple[str, str], np.ndarray] = {}
    for rep in REPS:
        probabilities[(rep, "unadjusted")] = np.full(len(cohort.y), np.nan)
        for block in residual_blocks:
            probabilities[(rep, f"residual_{block}")] = np.full(
                len(cohort.y), np.nan
            )
    for block in blocks_for(cohort):
        probabilities[(f"covariates_{block}", "covariates_only")] = np.full(
            len(cohort.y), np.nan
        )
    folds = np.zeros(len(cohort.y), dtype=int)

    raw = {
        "geneformer": cohort.geneformer,
        "hvg_pseudobulk": cohort.pseudobulk,
        "pca_pseudobulk": cohort.pseudobulk,
    }
    for fold, (train, test) in enumerate(
        splitter.split(cohort.geneformer, cohort.y), start=1
    ):
        for rep in REPS:
            probabilities[(rep, "unadjusted")][test] = predict_representation(
                rep,
                raw[rep][train],
                cohort.y[train],
                raw[rep][test],
                seed,
            )

        designs: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for block, (numeric, categorical) in blocks_for(cohort).items():
            train_cov, test_cov = make_design(
                cohort.metadata, train, test, numeric, categorical
            )
            designs[block] = (train_cov, test_cov)
            C = select_c_prepared(train_cov, cohort.y[train], seed)
            probabilities[(f"covariates_{block}", "covariates_only")][test] = (
                lr(C)
                .fit(train_cov, cohort.y[train])
                .predict_proba(test_cov)[:, 1]
            )

        for block in residual_blocks:
            train_cov, test_cov = designs[block]
            gf_train, gf_test = residualize(
                cohort.geneformer[train],
                cohort.geneformer[test],
                train_cov,
                test_cov,
            )
            pb_train, pb_test = residualize(
                cohort.pseudobulk[train],
                cohort.pseudobulk[test],
                train_cov,
                test_cov,
            )
            adjusted = {
                "geneformer": (gf_train, gf_test),
                "hvg_pseudobulk": (pb_train, pb_test),
                "pca_pseudobulk": (pb_train, pb_test),
            }
            for rep in REPS:
                train_x, test_x = adjusted[rep]
                probabilities[(rep, f"residual_{block}")][test] = (
                    predict_representation(
                        rep, train_x, cohort.y[train], test_x, seed
                    )
                )
        folds[test] = fold

    metrics = []
    predictions = []
    for (method, adjustment), values in probabilities.items():
        metrics.append(
            {
                "cohort": cohort.name,
                "repeat": repeat + 1,
                "seed": seed,
                "method": method,
                "adjustment": adjustment,
                "n_donors": len(cohort.y),
                "n_cases": int(cohort.y.sum()),
                "roc_auc": protocol.roc_auc(cohort.y, values),
                "pr_auc": protocol.pr_auc(cohort.y, values),
                "brier": protocol.brier(cohort.y, values),
                "ece_10bin": protocol.ece_10bin(cohort.y, values),
            }
        )
        predictions.extend(
            {
                "cohort": cohort.name,
                "repeat": repeat + 1,
                "seed": seed,
                "method": method,
                "adjustment": adjustment,
                "donor_id": donor,
                "fold": int(fold),
                "y_true": int(label),
                "prob_case": float(probability),
            }
            for donor, fold, label, probability in zip(
                cohort.donors, folds, cohort.y, values
            )
        )
    return metrics, predictions


def repeated_subset(
    cohort: Cohort,
    indices: np.ndarray,
    stratum: str,
) -> list[dict[str, object]]:
    rows = []
    for rep in REPS:
        values = []
        raw = (
            cohort.geneformer
            if rep == "geneformer"
            else cohort.pseudobulk
        )
        for seed in SEEDS:
            splitter = StratifiedKFold(
                effective_splits(cohort.y[indices]),
                shuffle=True,
                random_state=seed,
            )
            oof = np.full(len(indices), np.nan)
            for train_local, test_local in splitter.split(
                raw[indices], cohort.y[indices]
            ):
                oof[test_local] = predict_representation(
                    rep,
                    raw[indices][train_local],
                    cohort.y[indices][train_local],
                    raw[indices][test_local],
                    seed,
                )
            values.append(protocol.roc_auc(cohort.y[indices], oof))
        for repeat, (seed, auc) in enumerate(zip(SEEDS, values), start=1):
            rows.append(
                {
                    "cohort": cohort.name,
                    "stratum": stratum,
                    "draw": 0,
                    "repeat": repeat,
                    "seed": seed,
                    "method": rep,
                    "n_donors": len(indices),
                    "n_cases": int(cohort.y[indices].sum()),
                    "roc_auc": auc,
                }
            )
    return rows


def matched_subset_controls(
    cohort: Cohort,
    observed: np.ndarray,
    stratum: str,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(20260815)
    n_case = int(cohort.y[observed].sum())
    n_control = int(len(observed) - n_case)
    cases = np.flatnonzero(cohort.y == 1)
    controls = np.flatnonzero(cohort.y == 0)
    rows = []
    for draw, seed in enumerate(SEEDS, start=1):
        indices = np.concatenate(
            [
                rng.choice(cases, n_case, replace=False),
                rng.choice(controls, n_control, replace=False),
            ]
        )
        for rep in REPS:
            raw = (
                cohort.geneformer
                if rep == "geneformer"
                else cohort.pseudobulk
            )
            splitter = StratifiedKFold(
                effective_splits(cohort.y[indices]),
                shuffle=True,
                random_state=seed,
            )
            oof = np.full(len(indices), np.nan)
            for train_local, test_local in splitter.split(
                raw[indices], cohort.y[indices]
            ):
                oof[test_local] = predict_representation(
                    rep,
                    raw[indices][train_local],
                    cohort.y[indices][train_local],
                    raw[indices][test_local],
                    seed,
                )
            rows.append(
                {
                    "cohort": cohort.name,
                    "stratum": f"matched_random_for_{stratum}",
                    "draw": draw,
                    "repeat": draw,
                    "seed": seed,
                    "method": rep,
                    "n_donors": len(indices),
                    "n_cases": n_case,
                    "roc_auc": protocol.roc_auc(cohort.y[indices], oof),
                }
            )
    return rows


def batch_predictability(cohort: Cohort) -> list[dict[str, object]]:
    rows = []
    for batch in sorted(cohort.metadata["batch"].astype(str).unique()):
        target = cohort.metadata["batch"].astype(str).eq(batch).to_numpy(dtype=int)
        if min(target.sum(), len(target) - target.sum()) < 2:
            continue
        for rep in REPS:
            raw = (
                cohort.geneformer
                if rep == "geneformer"
                else cohort.pseudobulk
            )
            for repeat, seed in enumerate(BATCH_SEEDS, start=1):
                splitter = StratifiedKFold(
                    effective_splits(target), shuffle=True, random_state=seed
                )
                oof = np.full(len(target), np.nan)
                for train, test in splitter.split(raw, target):
                    oof[test] = predict_representation(
                        rep, raw[train], target[train], raw[test], seed
                    )
                rows.append(
                    {
                        "cohort": cohort.name,
                        "batch": batch,
                        "repeat": repeat,
                        "seed": seed,
                        "method": rep,
                        "n_positive": int(target.sum()),
                        "roc_auc": protocol.roc_auc(target, oof),
                    }
                )
    return rows


def summarise(metrics: pd.DataFrame) -> pd.DataFrame:
    return (
        metrics.groupby(["cohort", "method", "adjustment"], as_index=False)
        .agg(
            n_donors=("n_donors", "first"),
            n_cases=("n_cases", "first"),
            n_repeats=("repeat", "nunique"),
            auc_mean=("roc_auc", "mean"),
            auc_sd=("roc_auc", "std"),
            auc_p025=("roc_auc", lambda x: np.quantile(x, 0.025)),
            auc_p975=("roc_auc", lambda x: np.quantile(x, 0.975)),
            pr_auc_mean=("pr_auc", "mean"),
            brier_mean=("brier", "mean"),
            ece_10bin_mean=("ece_10bin", "mean"),
        )
        .sort_values(["cohort", "adjustment", "method"])
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cohorts = [load_cohort_174188(), load_cohort_135779()]
    all_metrics = []
    all_predictions = []

    for cohort in cohorts:
        residual_blocks = ("batch",)
        print(
            f"{cohort.name}: n={len(cohort.y)}, cases={cohort.y.sum()}, "
            f"GF={cohort.geneformer.shape}, PB={cohort.pseudobulk.shape}",
            flush=True,
        )
        completed = Parallel(
            n_jobs=N_JOBS,
            prefer="processes",
            max_nbytes="10M",
            mmap_mode="r",
            verbose=5,
        )(
            delayed(evaluate_repeat)(cohort, repeat, residual_blocks)
            for repeat in range(len(SEEDS))
        )
        for metrics, predictions in completed:
            all_metrics.extend(metrics)
            all_predictions.extend(predictions)

    metrics = pd.DataFrame(all_metrics)
    predictions = pd.DataFrame(all_predictions)
    summary = summarise(metrics)
    metrics.to_csv(OUT / "repeat_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(OUT / "oof_predictions.tsv.gz", sep="\t", index=False)
    summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)

    subset_rows = []
    cohort135 = cohorts[1]
    mixed = np.flatnonzero(~cohort135.metadata["batch"].eq("B1").to_numpy())
    subset_rows.extend(repeated_subset(cohort135, mixed, "mixed_batches_B2_B6"))
    subset_rows.extend(
        matched_subset_controls(cohort135, mixed, "mixed_batches_B2_B6")
    )
    subset = pd.DataFrame(subset_rows)
    subset.to_csv(OUT / "gse135779_mixed_batch_matched.tsv", sep="\t", index=False)

    batch = pd.DataFrame(batch_predictability(cohort135))
    batch.to_csv(
        OUT / "gse135779_representation_batch_predictability.tsv",
        sep="\t",
        index=False,
    )
    batch_summary = (
        batch.groupby(["cohort", "batch", "method"], as_index=False)
        .agg(
            n_positive=("n_positive", "first"),
            n_repeats=("repeat", "nunique"),
            auc_mean=("roc_auc", "mean"),
            auc_sd=("roc_auc", "std"),
        )
    )
    existing_wave = pd.read_csv(
        WORKSPACE / "results" / "wave_predictability.tsv", sep="\t"
    )
    existing_wave = existing_wave[
        existing_wave["target"].str.match(r"wave\d+_vs_rest")
    ].copy()
    existing_wave["cohort"] = "GSE174188_CD4"
    existing_wave["batch"] = existing_wave["target"].str.extract(
        r"wave(\d+)"
    )
    existing_wave = existing_wave.rename(
        columns={
            "representation": "method",
            "n_repeats": "n_repeats",
        }
    )[
        [
            "cohort",
            "batch",
            "method",
            "n_positive",
            "n_repeats",
            "auc_mean",
            "auc_sd",
        ]
    ]
    pd.concat([existing_wave, batch_summary], ignore_index=True).to_csv(
        OUT / "representation_batch_predictability_summary.tsv",
        sep="\t",
        index=False,
    )

    input_paths = [
        META135,
        META174,
        EXTRA
        / "04_models"
        / "Geneformer"
        / "SLE_GSE135779"
        / "geneformer_v2_316m_cell_sample500_clspool_logistic_maxlen4096_seed001"
        / "donor_embedding.parquet",
        EXTRA / "pseudobulk" / "SLE_GSE135779" / "donor_log1p_cpm.parquet",
        EXTRA
        / "04_models"
        / "Geneformer"
        / "SLE_GSE174188_CD4"
        / "geneformer_v2_316m_cell_sample1000_clspool_logistic_maxlen4096_seed001"
        / "donor_embedding.parquet",
        EXTRA
        / "pseudobulk"
        / "SLE_GSE174188_CD4"
        / "donor_log1p_cpm.parquet",
    ]
    manifest = {
        "analysis": "locked two-cohort donor-level validity audit",
        "seeds": SEEDS,
        "batch_predictability_seeds": BATCH_SEEDS,
        "outer_folds": N_SPLITS,
        "inner_folds": N_SPLITS,
        "C_grid": CS.tolist(),
        "n_hvg": N_HVG,
        "n_pc": N_PC,
        "ridge_alpha": RIDGE_ALPHA,
        "standardisation": (
            "Fit once on each outer-training fold; inner resampling selects C "
            "conditional on the outer-training representation."
        ),
        "residualisation": (
            "Ridge map from outer-training design matrix to each raw representation "
            "coordinate; train and test residuals passed through the identical "
            "outer-fold representation and classifier pipeline."
        ),
        "input_sha256": {str(path): sha256(path) for path in input_paths},
        "claim_boundary": (
            "Residualisation and design restriction estimate different quantities "
            "under label-design collinearity; neither alone identifies a causal "
            "technical effect."
        ),
    }
    (OUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print("\nSummary")
    print(summary.to_string(index=False))
    print(f"\nWrote locked audit to {OUT}")


if __name__ == "__main__":
    main()
