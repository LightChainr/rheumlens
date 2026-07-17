#!/usr/bin/env python3
"""Audit and summarize the locked RTX 5090 cap500 analysis archive."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA


DATASETS = ("SLE_GSE174188_CD4", "SLE_GSE285773_CD4")
METHOD_LABELS = {
    "mean": "Coordinate mean",
    "median": "Coordinate median",
    "trimmed_mean_10": "10% trimmed mean",
    "mean_std": "Mean + SD",
    "source_pca32_quantiles": "Source-PCA32 quantiles",
    "source_pca16_moments": "Source-PCA16 moments",
    "source_pca8_donor_mean": "Source-PCA8 donor mean",
    "source_pca16_donor_mean": "Source-PCA16 donor mean",
}
BLUE = "#246B8E"
GOLD = "#C58A1C"
INK = "#20252B"
MID = "#69727D"
LIGHT = "#DCE2E7"
PINK = "#C45A7A"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bh_adjust(p_values: pd.Series) -> np.ndarray:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(ranked)
    result[order] = np.minimum(ranked, 1.0)
    return result


def effective_rank(matrix: np.ndarray) -> tuple[float, float, float]:
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    singular = np.linalg.svd(centered, compute_uv=False)
    eigen = singular**2
    probability = eigen[eigen > 0] / eigen.sum()
    entropy_rank = float(np.exp(-(probability * np.log(probability)).sum()))
    participation_rank = float(1.0 / np.square(probability).sum())
    pc1_fraction = float(probability[0])
    return entropy_rank, participation_rank, pc1_fraction


def rms_dispersion(matrix: np.ndarray, labels: np.ndarray) -> float:
    residuals = []
    for label in np.unique(labels):
        subset = matrix[labels == label]
        residuals.append(subset - subset.mean(axis=0, keepdims=True))
    merged = np.concatenate(residuals, axis=0)
    return float(np.sqrt(np.mean(np.square(merged).sum(axis=1))))


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "axes.edgecolor": "#AAB2BA",
            "axes.linewidth": 0.8,
            "axes.titleweight": "bold",
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "axes.labelcolor": INK,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def audit_archive(root: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    rows = []
    manifest = json.loads((root / "archive_manifest.json").read_text())
    for item in manifest["files"]:
        path = root / item["path"]
        state = "missing" if not path.exists() else "verified" if sha256(path) == item["sha256"] else "mismatch"
        rows.append(
            {
                "scope": "archive_manifest",
                "path": item["path"],
                "status": state,
                "expected_sha256": item["sha256"],
                "actual_sha256": sha256(path) if path.exists() else "",
            }
        )

    shard_count = 0
    shard_failures = 0
    for dataset in DATASETS:
        shard_root = root / "outputs/cell_embeddings_cap500" / dataset / "shards"
        for record_path in sorted(shard_root.glob("*.json")):
            shard_count += 1
            record = json.loads(record_path.read_text())
            shard_id = record["shard_id"]
            meta = shard_root / f"{shard_id}.metadata.parquet"
            emb = shard_root / f"{shard_id}.embeddings.npy"
            valid = (
                meta.exists()
                and emb.exists()
                and sha256(meta) == record["metadata_sha256"]
                and sha256(emb) == record["embeddings_sha256"]
            )
            shard_failures += int(not valid)
        success = root / "outputs/cell_embeddings_cap500" / dataset / "_SUCCESS"
        rows.append(
            {
                "scope": "completion_marker",
                "path": str(success.relative_to(root)),
                "status": "present" if success.exists() else "missing",
                "expected_sha256": "",
                "actual_sha256": sha256(success) if success.exists() else "",
            }
        )

    invalidated = root / "outputs/distributional_pooling_cap500_sampling_order_mismatch"
    rows.append(
        {
            "scope": "invalidated_output",
            "path": str(invalidated.relative_to(root)),
            "status": "present_exclude_from_analysis" if invalidated.exists() else "absent",
            "expected_sha256": "",
            "actual_sha256": "",
        }
    )
    summary = {
        "archive_status": manifest.get("status"),
        "listed_files": len(manifest["files"]),
        "verified_listed_files": sum(row["status"] == "verified" for row in rows),
        "missing_listed_files": [row["path"] for row in rows if row["status"] == "missing"],
        "shard_records_checked": shard_count,
        "shard_checksum_failures": shard_failures,
        "canonical_analysis": "outputs/distributional_pooling_cap500",
        "excluded_analysis": "outputs/distributional_pooling_cap500_sampling_order_mismatch",
    }
    return pd.DataFrame(rows), summary


def load_donor_means(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames = []
    geometry_rows = []
    for dataset in DATASETS:
        reference = root / "reference_inputs" / dataset
        means = pd.read_parquet(reference / "historical_cap500_donor_embedding.parquet")
        means.index = means.index.astype(str)
        labels = pd.read_csv(reference / "donor_labels.tsv", sep="\t", dtype={"donor_id": str})
        labels = labels.set_index("donor_id").loc[means.index]
        matrix = means.to_numpy(dtype=float)
        entropy_rank, participation_rank, pc1_fraction = effective_rank(matrix)
        y = labels["case_control"].astype(str).to_numpy()
        disease_distance = np.linalg.norm(matrix[y == "case"].mean(axis=0) - matrix[y == "control"].mean(axis=0))
        within = rms_dispersion(matrix, y)
        geometry_rows.append(
            {
                "dataset": dataset,
                "n_donors": len(matrix),
                "embedding_dim": matrix.shape[1],
                "maximum_centered_rank": min(matrix.shape[0] - 1, matrix.shape[1]),
                "entropy_effective_rank": entropy_rank,
                "participation_effective_rank": participation_rank,
                "pc1_variance_fraction": pc1_fraction,
                "case_control_centroid_distance": float(disease_distance),
                "within_class_rms_dispersion": within,
                "normalized_disease_separation": float(disease_distance / within),
            }
        )
        frame = means.copy()
        frame.insert(0, "case_control", labels["case_control"].to_numpy())
        frame.insert(0, "donor_id", means.index)
        frame.insert(0, "dataset", dataset)
        frames.append(frame.reset_index(drop=True))

    combined = pd.concat(frames, ignore_index=True)
    feature_columns = [column for column in combined.columns if column not in {"dataset", "donor_id", "case_control"}]
    matrix = combined[feature_columns].to_numpy(dtype=float)
    cohort = combined["dataset"].to_numpy()
    cohort_distance = np.linalg.norm(
        matrix[cohort == DATASETS[0]].mean(axis=0) - matrix[cohort == DATASETS[1]].mean(axis=0)
    )
    cohort_within = rms_dispersion(matrix, cohort)
    geometry_rows.append(
        {
            "dataset": "pooled_cohort_shift",
            "n_donors": len(matrix),
            "embedding_dim": matrix.shape[1],
            "maximum_centered_rank": min(matrix.shape[0] - 1, matrix.shape[1]),
            "entropy_effective_rank": np.nan,
            "participation_effective_rank": np.nan,
            "pc1_variance_fraction": np.nan,
            "case_control_centroid_distance": np.nan,
            "within_class_rms_dispersion": cohort_within,
            "normalized_disease_separation": float(cohort_distance / cohort_within),
        }
    )
    pca = PCA(n_components=min(20, len(matrix) - 1), svd_solver="full")
    coordinates = pca.fit_transform(matrix)
    combined["PC1"] = coordinates[:, 0]
    combined["PC2"] = coordinates[:, 1]
    combined.attrs["explained_variance_ratio"] = pca.explained_variance_ratio_
    return combined, pd.DataFrame(geometry_rows)


def sample_cell_geometry(root: Path, sample_size: int = 10_000) -> pd.DataFrame:
    rng = np.random.default_rng(20_260_717)
    rows = []
    for dataset in DATASETS:
        shard_root = root / "outputs/cell_embeddings_cap500" / dataset / "shards"
        arrays = [np.load(path, mmap_mode="r") for path in sorted(shard_root.glob("*.embeddings.npy"))]
        n_cells = sum(len(array) for array in arrays)
        take = min(sample_size, n_cells)
        selected = np.sort(rng.choice(n_cells, take, replace=False))
        sample = np.empty((take, 1152), dtype=np.float32)
        offset = 0
        output_offset = 0
        for array in arrays:
            local = selected[(selected >= offset) & (selected < offset + len(array))] - offset
            sample[output_offset : output_offset + len(local)] = array[local]
            output_offset += len(local)
            offset += len(array)
        components = min(256, take - 1, sample.shape[1])
        pca = PCA(n_components=components, svd_solver="randomized", random_state=20_260_717).fit(sample)
        explained = pca.explained_variance_ratio_
        probability = explained / explained.sum()
        rows.append(
            {
                "dataset": dataset,
                "n_cells_total": n_cells,
                "n_cells_sampled": take,
                "sample_seed": 20_260_717,
                "spectrum_components": components,
                "spectrum_variance_captured": float(explained.sum()),
                "pc1_variance_fraction": float(explained[0]),
                "entropy_effective_rank_top256": float(
                    np.exp(-(probability * np.log(probability)).sum())
                ),
                "participation_effective_rank_top256": float(1.0 / np.square(probability).sum()),
            }
        )
    return pd.DataFrame(rows)


def summarize_predictions(canonical: Path, metrics: pd.DataFrame) -> pd.DataFrame:
    predictions = pd.read_parquet(canonical / "target_predictions.parquet")
    rows = []
    for (direction, method), frame in predictions.groupby(["direction", "method"], sort=False):
        probability = frame["prob_case"].to_numpy(float)
        y = frame["y_true"].to_numpy(int)
        metric = metrics.loc[metrics["direction"].eq(direction) & metrics["method"].eq(method)].iloc[0]
        rows.append(
            {
                "direction": direction,
                "target_dataset": metric["target_dataset"],
                "method": method,
                "method_label": METHOD_LABELS[method],
                "n_target": len(frame),
                "case_fraction": float(y.mean()),
                "roc_auc": metric["roc_auc"],
                "pr_auc": metric["pr_auc"],
                "brier": metric["brier"],
                "ece_10bin": metric["ece_10bin"],
                "probability_min": float(probability.min()),
                "probability_max": float(probability.max()),
                "fraction_probability_lt_0.01": float((probability < 0.01).mean()),
                "fraction_probability_gt_0.99": float((probability > 0.99).mean()),
                "fraction_extreme_probability": float(((probability < 0.01) | (probability > 0.99)).mean()),
                "mean_probability_controls": float(probability[y == 0].mean()),
                "mean_probability_cases": float(probability[y == 1].mean()),
            }
        )
    return pd.DataFrame(rows)


def plot_pooling_summary(metrics: pd.DataFrame, comparisons: pd.DataFrame, output: Path) -> None:
    setup_style()
    order = metrics.groupby("method")["roc_auc"].mean().sort_values(ascending=True).index.tolist()
    directions = metrics["direction"].drop_duplicates().tolist()
    target_labels = {
        directions[0]: "Target: GSE285773 (n=26)",
        directions[1]: "Target: GSE174188 (n=261)",
    }
    fig = plt.figure(figsize=(12.4, 9.2))
    grid = fig.add_gridspec(
        2,
        2,
        height_ratios=(1.15, 1.0),
        left=0.16,
        right=0.985,
        bottom=0.105,
        top=0.885,
        hspace=0.30,
        wspace=0.22,
    )
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    for axis, direction in zip(axes, directions):
        subset = metrics[metrics["direction"].eq(direction)].set_index("method").loc[order]
        y = np.arange(len(order))
        colors = [GOLD if method == "mean" else BLUE for method in order]
        fills = [GOLD if method == "mean" else "white" for method in order]
        axis.hlines(y, 0.75, subset["roc_auc"], color=LIGHT, linewidth=2.0, zorder=1)
        for index, (value, edge, fill) in enumerate(zip(subset["roc_auc"], colors, fills)):
            axis.scatter(value, index, s=54, edgecolor=edge, facecolor=fill, linewidth=1.6, zorder=3)
            axis.text(value + 0.006, index, f"{value:.3f}", va="center", fontsize=8)
        baseline = float(subset.loc["mean", "roc_auc"])
        axis.axvline(baseline, color=GOLD, linestyle=(0, (3, 3)), linewidth=1.0)
        axis.set_yticks(y, [METHOD_LABELS[method] for method in order])
        axis.set_xlim(0.75, 0.985)
        axis.set_xlabel("External ROC-AUC")
        axis.set_title(target_labels[direction], loc="left")
        axis.grid(axis="x", color="#E8ECEF", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
    axes[1].set_yticklabels([])
    axes[0].text(-0.17, 1.05, "A", transform=axes[0].transAxes, fontsize=14, fontweight="bold")

    delta_axis = fig.add_subplot(grid[1, 0])
    comparison_order = [method for method in order if method != "mean"]
    pivot = comparisons.pivot(index="method", columns="target_dataset", values="delta_auc").loc[comparison_order]
    y = np.arange(len(comparison_order))
    small = pivot["SLE_GSE285773_CD4"].to_numpy()
    large = pivot["SLE_GSE174188_CD4"].to_numpy()
    delta_axis.hlines(y, np.minimum(small, large), np.maximum(small, large), color=MID, linewidth=1.2)
    delta_axis.scatter(small, y, color=GOLD, edgecolor="white", linewidth=0.8, s=48, zorder=3)
    delta_axis.scatter(large, y, color=BLUE, edgecolor="white", linewidth=0.8, s=48, zorder=3)
    delta_axis.axvline(0, color=INK, linewidth=0.9)
    delta_axis.set_yticks(y, [METHOD_LABELS[method] for method in comparison_order])
    delta_axis.set_xlabel("AUC difference from coordinate mean")
    delta_axis.set_title("Direction-specific changes relative to mean pooling", loc="left")
    delta_axis.grid(axis="x", color="#E8ECEF", linewidth=0.7)
    delta_axis.spines[["top", "right"]].set_visible(False)
    delta_axis.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=GOLD, label="Target n=26"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, label="Target n=261"),
        ],
        frameon=False,
        loc="lower right",
    )
    delta_axis.text(-0.17, 1.06, "B", transform=delta_axis.transAxes, fontsize=14, fontweight="bold")

    calibration_axis = fig.add_subplot(grid[1, 1])
    for direction, color, marker in zip(directions, (GOLD, BLUE), ("o", "s")):
        subset = metrics[metrics["direction"].eq(direction)]
        calibration_axis.scatter(
            subset["brier"],
            subset["roc_auc"],
            s=42 + subset["ece_10bin"].to_numpy() * 150,
            color=color,
            marker=marker,
            alpha=0.78,
            edgecolor="white",
            linewidth=0.7,
            label=target_labels[direction],
        )
        for row in subset.itertuples():
            if row.method in {"mean", "source_pca8_donor_mean", "source_pca16_moments"}:
                calibration_axis.annotate(
                    METHOD_LABELS[row.method],
                    (row.brier, row.roc_auc),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7,
                )
    calibration_axis.set_xlabel("Brier score (lower is better)")
    calibration_axis.set_ylabel("External ROC-AUC (higher is better)")
    calibration_axis.set_title("Discrimination and probability error", loc="left")
    calibration_axis.grid(color="#E8ECEF", linewidth=0.7)
    calibration_axis.spines[["top", "right"]].set_visible(False)
    calibration_axis.legend(frameon=False, fontsize=8, loc="lower left")
    calibration_axis.text(-0.17, 1.06, "C", transform=calibration_axis.transAxes, fontsize=14, fontweight="bold")

    fig.suptitle(
        "Fixed Geneformer pooling under reciprocal source-only transfer",
        fontsize=15,
        fontweight="bold",
        x=0.02,
        y=0.975,
        ha="left",
    )
    fig.text(
        0.01,
        0.025,
        "Cap500 cell-level embeddings; logistic-regression tuning used source donors only. Point size in panel C reflects 10-bin ECE.",
        fontsize=8,
        color=MID,
    )
    for suffix in ("pdf", "png"):
        fig.savefig(output / f"FIG_5090_cap500_pooling_stress_test.{suffix}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def plot_geometry(
    combined: pd.DataFrame, geometry: pd.DataFrame, cell_geometry: pd.DataFrame, output: Path
) -> None:
    setup_style()
    feature_columns = [column for column in combined.columns if column not in {"dataset", "donor_id", "case_control", "PC1", "PC2"}]
    fig, (scatter, scree, rank_axis) = plt.subplots(1, 3, figsize=(16.0, 5.1))
    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.17, top=0.82, wspace=0.23)
    cohort_style = {
        DATASETS[0]: (BLUE, "GSE174188"),
        DATASETS[1]: (GOLD, "GSE285773"),
    }
    marker_style = {"case": "o", "control": "^"}
    for dataset, (color, label) in cohort_style.items():
        for disease, marker in marker_style.items():
            subset = combined[combined["dataset"].eq(dataset) & combined["case_control"].eq(disease)]
            scatter.scatter(
                subset["PC1"],
                subset["PC2"],
                s=28 if dataset == DATASETS[0] else 46,
                marker=marker,
                facecolor=color if disease == "case" else "white",
                edgecolor=color,
                linewidth=0.9,
                alpha=0.78,
                label=f"{label}, {disease}",
            )
    variance = combined.attrs["explained_variance_ratio"]
    scatter.set_xlabel(f"Pooled donor PC1 ({variance[0] * 100:.1f}% variance)")
    scatter.set_ylabel(f"Pooled donor PC2 ({variance[1] * 100:.1f}% variance)")
    scatter.set_title("Donor-mean Geneformer geometry", loc="left")
    scatter.grid(color="#E8ECEF", linewidth=0.7)
    scatter.spines[["top", "right"]].set_visible(False)
    scatter.legend(frameon=False, fontsize=8, loc="best")
    scatter.text(-0.14, 1.04, "A", transform=scatter.transAxes, fontsize=14, fontweight="bold")

    for dataset, (color, label) in cohort_style.items():
        subset = combined[combined["dataset"].eq(dataset)]
        matrix = subset[feature_columns].to_numpy(dtype=float)
        pca = PCA(n_components=min(20, len(matrix) - 1), svd_solver="full").fit(matrix)
        cumulative = np.cumsum(pca.explained_variance_ratio_)
        scree.plot(np.arange(1, len(cumulative) + 1), cumulative, color=color, marker="o", markersize=3, label=label)
    scree.set_xlim(1, 20)
    scree.set_ylim(0, 1.02)
    scree.set_xlabel("Number of donor-level PCs")
    scree.set_ylabel("Cumulative variance explained")
    scree.set_title("Low-dimensional concentration", loc="left")
    scree.grid(color="#E8ECEF", linewidth=0.7)
    scree.spines[["top", "right"]].set_visible(False)
    scree.legend(frameon=False, loc="lower right")
    text_lines = []
    for row in geometry[geometry["dataset"].isin(DATASETS)].itertuples():
        short = "GSE174188" if row.dataset == DATASETS[0] else "GSE285773"
        text_lines.append(f"{short}: entropy rank {row.entropy_effective_rank:.1f}; PC1 {row.pc1_variance_fraction:.1%}")
    scree.text(0.02, 0.97, "\n".join(text_lines), transform=scree.transAxes, va="top", fontsize=8, color=MID)
    scree.text(-0.14, 1.04, "B", transform=scree.transAxes, fontsize=14, fontweight="bold")

    donor_geometry = geometry[geometry["dataset"].isin(DATASETS)].set_index("dataset")
    cell_indexed = cell_geometry.set_index("dataset")
    positions = np.arange(2)
    width = 0.34
    donor_values = [donor_geometry.loc[dataset, "entropy_effective_rank"] for dataset in DATASETS]
    cell_values = [cell_indexed.loc[dataset, "entropy_effective_rank_top256"] for dataset in DATASETS]
    rank_axis.bar(
        positions - width / 2,
        cell_values,
        width,
        color=[BLUE, GOLD],
        edgecolor=[BLUE, GOLD],
        label="Cell-level top-256 spectrum",
    )
    rank_axis.bar(
        positions + width / 2,
        donor_values,
        width,
        color="white",
        edgecolor=[BLUE, GOLD],
        linewidth=1.6,
        hatch="//",
        label="Donor means, full spectrum",
    )
    for x, value in zip(positions - width / 2, cell_values):
        rank_axis.text(x, value + 1.2, f"{value:.1f}", ha="center", va="bottom", fontsize=8)
    for x, value in zip(positions + width / 2, donor_values):
        rank_axis.text(x, value + 1.2, f"{value:.1f}", ha="center", va="bottom", fontsize=8)
    rank_axis.set_xticks(positions, ["GSE174188", "GSE285773"])
    rank_axis.set_ylabel("Entropy effective rank")
    rank_axis.set_title("Cell-to-donor dimensional contraction", loc="left")
    rank_axis.grid(axis="y", color="#E8ECEF", linewidth=0.7)
    rank_axis.spines[["top", "right"]].set_visible(False)
    rank_axis.legend(frameon=False, fontsize=8, loc="upper right")
    rank_axis.text(-0.14, 1.04, "C", transform=rank_axis.transAxes, fontsize=14, fontweight="bold")

    fig.suptitle(
        "Cohort shift and effective dimensionality of donor-mean embeddings",
        fontsize=15,
        fontweight="bold",
        x=0.02,
        y=0.965,
        ha="left",
    )
    fig.text(
        0.01,
        0.035,
        "PCA is descriptive and fitted jointly only for visualization. Cell spectra use a fixed 10,000-cell sample and 256 randomized PCs; donor ranks use the full centered spectrum.",
        fontsize=8,
        color=MID,
    )
    for suffix in ("pdf", "png"):
        fig.savefig(output / f"FIG_5090_donor_embedding_geometry.{suffix}", dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    audit, audit_summary = audit_archive(args.archive_root)
    audit.to_csv(args.output / "archive_integrity_audit.tsv", sep="\t", index=False)
    (args.output / "archive_integrity_summary.json").write_text(json.dumps(audit_summary, indent=2) + "\n")

    canonical = args.archive_root / "outputs/distributional_pooling_cap500"
    metrics = pd.read_csv(canonical / "transfer_metrics.tsv", sep="\t")
    comparisons = pd.read_csv(canonical / "paired_comparisons.tsv", sep="\t")
    comparisons["delong_bh_q_14_tests"] = bh_adjust(comparisons["delong_p"])
    comparisons["method_label"] = comparisons["method"].map(METHOD_LABELS)
    comparisons.to_csv(args.output / "paired_comparisons_with_bh.tsv", sep="\t", index=False)

    metrics["method_label"] = metrics["method"].map(METHOD_LABELS)
    metrics.to_csv(args.output / "transfer_metrics_labeled.tsv", sep="\t", index=False)
    prediction_summary = summarize_predictions(canonical, metrics)
    prediction_summary.to_csv(args.output / "prediction_calibration_extremeness.tsv", sep="\t", index=False)

    combined, geometry = load_donor_means(args.archive_root)
    cell_geometry = sample_cell_geometry(args.archive_root)
    combined[["dataset", "donor_id", "case_control", "PC1", "PC2"]].to_csv(
        args.output / "donor_mean_pca_coordinates.tsv", sep="\t", index=False
    )
    geometry.to_csv(args.output / "donor_embedding_geometry.tsv", sep="\t", index=False)
    cell_geometry.to_csv(args.output / "cell_embedding_sample_geometry.tsv", sep="\t", index=False)

    plot_pooling_summary(metrics, comparisons, args.output)
    plot_geometry(combined, geometry, cell_geometry, args.output)

    large = comparisons[comparisons["target_dataset"].eq("SLE_GSE174188_CD4")]
    significant_large = large.loc[large["delong_bh_q_14_tests"] < 0.05, ["method", "delta_auc", "delong_bh_q_14_tests"]]
    report = {
        **audit_summary,
        "local_analysis_status": "complete",
        "canonical_historical_reproduction_passed": bool(
            json.loads((canonical / "manifest.json").read_text())["historical_mean_reproduction_passed"]
        ),
        "methods_evaluated": int(metrics["method"].nunique()),
        "directions_evaluated": int(metrics["direction"].nunique()),
        "methods_promoted_to_cap1000": [],
        "large_target_bh_significant_differences_from_mean": significant_large.to_dict(orient="records"),
        "mean_pooling_small_target_brier": float(
            metrics.loc[metrics["method"].eq("mean") & metrics["target_dataset"].eq("SLE_GSE285773_CD4"), "brier"].iloc[0]
        ),
        "mean_pooling_large_target_brier": float(
            metrics.loc[metrics["method"].eq("mean") & metrics["target_dataset"].eq("SLE_GSE174188_CD4"), "brier"].iloc[0]
        ),
        "cell_to_donor_entropy_rank_ratios": {
            dataset: float(
                cell_geometry.set_index("dataset").loc[dataset, "entropy_effective_rank_top256"]
                / geometry.set_index("dataset").loc[dataset, "entropy_effective_rank"]
            )
            for dataset in DATASETS
        },
        "outputs": sorted(path.name for path in args.output.iterdir()),
    }
    (args.output / "analysis_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
