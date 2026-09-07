#!/usr/bin/env python3
"""Full-pipeline null and degrees-of-freedom audit for GSE135779 design models."""

from __future__ import annotations

import itertools
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from common import (
    BASELINE,
    SEEDS,
    WORKSPACE,
    effective_splits,
    encode_design,
    information_fraction,
    logistic,
    select_c_prepared,
    sha256,
    write_json,
)


OUT = WORKSPACE / "results" / "gse135779_null"
METADATA_PATH = (
    BASELINE
    / "results"
    / "gse135779_metadata"
    / "gse135779_donor_metadata_restored.tsv"
)
RESTRICTION_PATH = (
    BASELINE
    / "results"
    / "locked_validity_audit"
    / "gse135779_mixed_batch_matched.tsv"
)
N_LABEL_PERMUTATIONS = int(os.environ.get("N_LABEL_PERMUTATIONS", "1000"))
N_RANDOM_DESIGNS = int(os.environ.get("N_RANDOM_DESIGNS", "250"))
N_JOBS = min(int(os.environ.get("N_JOBS", "12")), os.cpu_count() or 1)
BASE_SEED = 2026072800


QC = [
    "log_estimated_cells",
    "log_mean_reads_per_cell",
    "log_median_genes_per_cell",
    "sequencing_saturation",
    "fraction_reads_in_cells",
    "log_median_umi_per_cell",
]
DEMOGRAPHIC_NUMERIC = ["age_years"]
DEMOGRAPHIC_CATEGORICAL = ["sex", "race", "ethnicity"]
BLOCKS: dict[str, tuple[list[str], list[str]]] = {
    "batch": ([], ["batch"]),
    "collection_year": ([], ["collection_year"]),
    "qc": (QC, []),
    "demographic": (DEMOGRAPHIC_NUMERIC, DEMOGRAPHIC_CATEGORICAL),
    "all": (
        QC + DEMOGRAPHIC_NUMERIC,
        ["batch", "collection_year", *DEMOGRAPHIC_CATEGORICAL],
    ),
}


def effective_df(x: np.ndarray, probability: np.ndarray, C: float) -> float:
    """Penalised-logistic effective-df approximation.

    Liblinear penalises the intercept, so an intercept column is included in the
    same ridge-style trace approximation. This is a diagnostic, not an exact
    finite-sample degree of freedom for the fitted classifier.
    """

    augmented = np.column_stack([np.ones(len(x)), x])
    weights = np.clip(probability * (1.0 - probability), 1e-8, None)
    gram = augmented.T @ (weights[:, None] * augmented)
    penalty = np.eye(gram.shape[0]) / float(C)
    return float(np.trace(gram @ np.linalg.pinv(gram + penalty)))


def evaluate_design(
    metadata: pd.DataFrame,
    y: np.ndarray,
    block: str,
    seed: int,
    include_diagnostics: bool = False,
) -> dict[str, float]:
    numeric, categorical = BLOCKS[block]
    splitter = StratifiedKFold(
        effective_splits(y),
        shuffle=True,
        random_state=seed,
    )
    prediction = np.full(len(y), np.nan)
    feature_counts = []
    selected_cs = []
    effective_dfs = []
    for train, test in splitter.split(metadata, y):
        train_x, test_x, _ = encode_design(
            metadata,
            train,
            test,
            numeric,
            categorical,
        )
        C = select_c_prepared(train_x, y[train], seed)
        model = logistic(C).fit(train_x, y[train])
        prediction[test] = model.predict_proba(test_x)[:, 1]
        if include_diagnostics:
            feature_counts.append(train_x.shape[1])
            selected_cs.append(C)
            effective_dfs.append(
                effective_df(
                    train_x,
                    model.predict_proba(train_x)[:, 1],
                    C,
                )
            )
    result = {"roc_auc": float(roc_auc_score(y, prediction))}
    if include_diagnostics:
        result.update(
            {
                "encoded_features_mean": float(np.mean(feature_counts)),
                "encoded_features_min": int(np.min(feature_counts)),
                "encoded_features_max": int(np.max(feature_counts)),
                "selected_c_median": float(np.median(selected_cs)),
                "effective_df_mean": float(np.mean(effective_dfs)),
            }
        )
    return result


def label_permutation(
    metadata: pd.DataFrame,
    y: np.ndarray,
    block: str,
    permutation: int,
) -> dict:
    seed = BASE_SEED + permutation
    rng = np.random.default_rng(seed)
    permuted = rng.permutation(y)
    result = evaluate_design(metadata, permuted, block, seed)
    return {
        "null_type": "label_permutation",
        "block": block,
        "iteration": permutation,
        "seed": seed,
        **result,
    }


def shuffled_design(
    metadata: pd.DataFrame,
    y: np.ndarray,
    block: str,
    iteration: int,
) -> dict:
    seed = BASE_SEED + 1_000_000 + iteration
    rng = np.random.default_rng(seed)
    numeric, categorical = BLOCKS[block]
    shuffled = metadata.copy()
    for column in [*numeric, *categorical]:
        shuffled[column] = rng.permutation(shuffled[column].to_numpy())
    result = evaluate_design(shuffled, y, block, seed)
    return {
        "null_type": "independent_column_shuffle",
        "block": block,
        "iteration": iteration,
        "seed": seed,
        **result,
    }


def information_table(metadata: pd.DataFrame, y: np.ndarray) -> pd.DataFrame:
    rows = []
    all_indices = np.arange(len(y))
    for block, (numeric, categorical) in BLOCKS.items():
        design, _, _ = encode_design(
            metadata,
            all_indices,
            all_indices,
            numeric,
            categorical,
        )
        rows.append({"block": block, **information_fraction(y, design)})
    return pd.DataFrame(rows)


def overlap_table(metadata: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for column in ("batch", "collection_year"):
        table = pd.crosstab(metadata[column], metadata["y_true"])
        mixed_levels = table.index[(table > 0).all(axis=1)]
        in_mixed = metadata[column].isin(mixed_levels)
        rows.append(
            {
                "design_variable": column,
                "n_levels": int(len(table)),
                "n_mixed_levels": int(len(mixed_levels)),
                "n_donors_in_mixed_levels": int(in_mixed.sum()),
                "fraction_donors_in_mixed_levels": float(in_mixed.mean()),
                "mixed_levels": ",".join(map(str, mixed_levels)),
            }
        )
    return pd.DataFrame(rows)


def exact_sign_flip_p(differences: np.ndarray) -> float:
    observed = abs(float(np.mean(differences)))
    n = len(differences)
    exceed = 0
    total = 0
    for signs in itertools.product((-1.0, 1.0), repeat=n):
        statistic = abs(float(np.mean(differences * np.asarray(signs))))
        exceed += int(statistic >= observed - 1e-15)
        total += 1
    return exceed / total


def restriction_tests() -> pd.DataFrame:
    data = pd.read_csv(RESTRICTION_PATH, sep="\t")
    observed = data[data["stratum"].eq("mixed_batches_B2_B6")]
    matched = data[
        data["stratum"].eq("matched_random_for_mixed_batches_B2_B6")
    ]
    rows = []
    for method in sorted(observed["method"].unique()):
        left = (
            observed[observed["method"].eq(method)]
            .set_index("seed")["roc_auc"]
            .sort_index()
        )
        right = (
            matched[matched["method"].eq(method)]
            .set_index("seed")["roc_auc"]
            .sort_index()
        )
        paired = left.to_frame("observed").join(
            right.to_frame("matched"),
            how="inner",
        )
        differences = (paired["observed"] - paired["matched"]).to_numpy()
        rows.append(
            {
                "method": method,
                "n_pairs": len(differences),
                "observed_auc_mean": float(paired["observed"].mean()),
                "matched_auc_mean": float(paired["matched"].mean()),
                "delta_observed_minus_matched": float(np.mean(differences)),
                "exact_two_sided_sign_flip_p": exact_sign_flip_p(differences),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(METADATA_PATH, sep="\t")
    y = metadata["y_true"].to_numpy(dtype=int)

    observed_rows = []
    for block in BLOCKS:
        for repeat, seed in enumerate(SEEDS, start=1):
            observed_rows.append(
                {
                    "block": block,
                    "repeat": repeat,
                    "seed": seed,
                    **evaluate_design(
                        metadata,
                        y,
                        block,
                        seed,
                        include_diagnostics=True,
                    ),
                }
            )
    observed = pd.DataFrame(observed_rows)
    observed.to_csv(OUT / "observed_design_metrics.tsv", sep="\t", index=False)

    jobs = [
        ("label", block, iteration)
        for block in BLOCKS
        for iteration in range(1, N_LABEL_PERMUTATIONS + 1)
    ] + [
        ("design", block, iteration)
        for block in BLOCKS
        for iteration in range(1, N_RANDOM_DESIGNS + 1)
    ]

    def dispatch(kind: str, block: str, iteration: int) -> dict:
        if kind == "label":
            return label_permutation(metadata, y, block, iteration)
        return shuffled_design(metadata, y, block, iteration)

    null_rows = Parallel(n_jobs=N_JOBS, prefer="processes", verbose=5)(
        delayed(dispatch)(*job) for job in jobs
    )
    nulls = pd.DataFrame(null_rows)
    nulls.to_csv(
        OUT / "design_null_distributions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    summary_rows = []
    for block in BLOCKS:
        observed_mean = float(
            observed.loc[observed["block"].eq(block), "roc_auc"].mean()
        )
        for null_type in nulls["null_type"].unique():
            values = nulls.loc[
                nulls["block"].eq(block)
                & nulls["null_type"].eq(null_type),
                "roc_auc",
            ].to_numpy()
            summary_rows.append(
                {
                    "block": block,
                    "null_type": null_type,
                    "observed_auc_mean": observed_mean,
                    "null_n": len(values),
                    "null_auc_mean": float(np.mean(values)),
                    "null_auc_p950": float(np.quantile(values, 0.95)),
                    "null_auc_p990": float(np.quantile(values, 0.99)),
                    "empirical_upper_p": float(
                        (1 + np.sum(values >= observed_mean))
                        / (1 + len(values))
                    ),
                }
            )
    pd.DataFrame(summary_rows).to_csv(
        OUT / "design_null_summary.tsv",
        sep="\t",
        index=False,
    )
    information_table(metadata, y).to_csv(
        OUT / "design_information_fraction.tsv",
        sep="\t",
        index=False,
    )
    overlap_table(metadata).to_csv(
        OUT / "discrete_overlap_summary.tsv",
        sep="\t",
        index=False,
    )
    restriction_tests().to_csv(
        OUT / "restriction_sign_flip_tests.tsv",
        sep="\t",
        index=False,
    )
    write_json(
        OUT / "null_audit_manifest.json",
        {
            "metadata_path": str(METADATA_PATH),
            "metadata_sha256": sha256(METADATA_PATH),
            "restriction_path": str(RESTRICTION_PATH),
            "restriction_sha256": sha256(RESTRICTION_PATH),
            "blocks": BLOCKS,
            "observed_seeds": SEEDS,
            "n_label_permutations_per_block": N_LABEL_PERMUTATIONS,
            "n_random_designs_per_block": N_RANDOM_DESIGNS,
            "base_seed": BASE_SEED,
            "classifier": (
                "fold-contained design encoding; balanced liblinear logistic; "
                "inner 5-fold selection of C from 1e-4 through 1e4"
            ),
        },
    )
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print(f"Wrote GSE135779 null audit to {OUT}")


if __name__ == "__main__":
    main()
