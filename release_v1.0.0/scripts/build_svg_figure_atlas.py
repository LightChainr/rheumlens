#!/usr/bin/env python3
"""Replot the traceable local evidence as editable, multi-panel SVG drafts."""

from __future__ import annotations

from pathlib import Path
import os
import string

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


RELEASE = Path(__file__).resolve().parents[1]
PROJECT = Path(os.environ.get("SLE_SOURCE_PROJECT", RELEASE / "source_data")).resolve()
RESULTS = PROJECT / "04_重构_投稿定位_20260714/results"
ANALYSIS = PROJECT / "08_Applied_Sciences_special_issue_submission_20260717/analysis"
SVG = RELEASE / "figures/svg"
PNG = RELEASE / "figures/previews"
PDF = RELEASE / "figures/pdf"
PROV = RELEASE / "provenance"

NAVY = "#16324F"
BLUE = "#2878B5"
CYAN = "#4CB3C7"
GOLD = "#E3A018"
CORAL = "#D95F59"
RED = "#A23B3B"
GREEN = "#3A7D44"
PURPLE = "#7566A8"
GREY = "#69727D"
LIGHT = "#EEF2F5"
METHOD_COLORS = {
    "Frozen Geneformer": BLUE,
    "Source-HVG pseudobulk": GOLD,
    "Source-PCA pseudobulk": CORAL,
    "Mean-HVG pseudobulk": GOLD,
    "PCA pseudobulk": CORAL,
    "HVG pseudobulk": GOLD,
    "Geneformer + HVG": GREEN,
    "Geneformer + PCA": PURPLE,
}
COHORT_LABELS = {
    "SLE_GSE135779": "GSE135779",
    "SLE_GSE174188_CD4": "GSE174188 CD4",
    "SLE_GSE285773_CD4": "GSE285773 CD4",
    "GSE135779": "GSE135779",
    "GSE174188_CD4": "GSE174188 CD4",
    "GSE285773": "GSE285773 CD4",
}


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def direction_label(value: str) -> str:
    a, b = value.replace("SLE_", "").split("_to_")
    return f"{COHORT_LABELS.get(a, a)} → {COHORT_LABELS.get(b, b)}"


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.labelsize": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )
    sns.set_theme(style="ticks", rc={"axes.facecolor": "white", "figure.facecolor": "white"})


def panel_labels(fig: plt.Figure, axes, labels=None) -> None:
    labels = labels or list(string.ascii_uppercase)
    for ax, label in zip(np.ravel(axes), labels):
        pos = ax.get_position()
        fig.text(pos.x0 - 0.018, pos.y1 + 0.012, label, fontsize=13, fontweight="bold", color=NAVY)


def finish(fig: plt.Figure, stem: str, sources: list[Path]) -> None:
    for folder in (SVG, PNG, PDF):
        folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(SVG / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(PNG / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(PDF / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    with (PROV / "figure_source_index.tsv").open("a", encoding="utf-8") as handle:
        for source in sources:
            handle.write(f"{stem}\t{source.relative_to(PROJECT)}\n")


def fig01_overview() -> None:
    metric_path = RESULTS / "repeated_strict_fold_benchmark/repeat_level_metrics.tsv"
    metrics = read(metric_path)
    cohorts = metrics.drop_duplicates("dataset_id")[["dataset_id", "n_donors", "n_case", "n_control"]]
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.1), gridspec_kw={"height_ratios": [1, 1.1]})
    ax = axes[0, 0]
    x = np.arange(len(cohorts))
    ax.bar(x, cohorts.n_control, color=CYAN, label="Control")
    ax.bar(x, cohorts.n_case, bottom=cohorts.n_control, color=CORAL, label="SLE")
    for i, row in enumerate(cohorts.itertuples()):
        ax.text(i, row.n_donors + 4, f"n={row.n_donors}", ha="center", fontweight="bold")
    ax.set_xticks(x, [COHORT_LABELS[v] for v in cohorts.dataset_id], rotation=12, ha="right")
    ax.set_ylabel("Donors")
    ax.set_title("Three independent SLE cohorts")
    ax.legend(frameon=False, ncol=2, loc="upper left")

    ax = axes[0, 1]
    ax.axis("off")
    nodes = [(0.03, "Single-cell\ntranscriptomes"), (0.29, "Cell-level\nrepresentation"), (0.55, "Donor-level\naggregation"), (0.81, "Clinical-state\ndiscrimination")]
    for i, (x0, label) in enumerate(nodes):
        box = FancyBboxPatch((x0, 0.33), 0.17, 0.34, boxstyle="round,pad=0.02,rounding_size=0.025", facecolor=[LIGHT, "#E3F0F5", "#F8EAC7", "#F6DAD8"][i], edgecolor=NAVY, linewidth=1.2)
        ax.add_patch(box)
        ax.text(x0 + 0.085, 0.5, label, ha="center", va="center", fontweight="bold")
        if i < len(nodes) - 1:
            ax.add_patch(FancyArrowPatch((x0 + 0.18, 0.5), (nodes[i + 1][0] - 0.01, 0.5), arrowstyle="-|>", mutation_scale=13, color=NAVY, linewidth=1.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Patient-level representation problem")

    ax = axes[1, 0]
    ax.axis("off")
    labels = ["Frozen\nGeneformer", "Mean-HVG\npseudobulk", "PCA\npseudobulk", "scVI / scGPT\nsensitivities"]
    colors = [BLUE, GOLD, CORAL, PURPLE]
    for i, (label, color) in enumerate(zip(labels, colors)):
        x0 = 0.02 + i * 0.245
        ax.add_patch(FancyBboxPatch((x0, 0.36), 0.2, 0.28, boxstyle="round,pad=0.02,rounding_size=0.02", facecolor=color, alpha=0.9, edgecolor="white"))
        ax.text(x0 + 0.1, 0.5, label, color="white", ha="center", va="center", fontweight="bold")
    ax.text(0.5, 0.18, "Same donors • same outcome • donor-disjoint evaluation", ha="center", color=NAVY, fontweight="bold")
    ax.set_title("Matched representation families")

    ax = axes[1, 1]
    labels = ["Repeated\ninternal CV", "Strict source-only\nreciprocal transfer", "Disease-program\nattribution", "Cell and compute\nbudgets"]
    vals = [60, 2, 20, 4]
    colors = [NAVY, CORAL, GREEN, GOLD]
    y = np.arange(4)[::-1]
    ax.barh(y, np.log10(np.array(vals) + 1), color=colors, height=0.62)
    for yi, v in zip(y, vals):
        ax.text(np.log10(v + 1) + 0.02, yi, str(v), va="center", fontweight="bold")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Analysis units (log-scaled display)")
    ax.set_title("Evidence ladder")
    ax.grid(axis="x", color="#D9E0E5", linewidth=0.7)
    fig.suptitle("From single cells to transferable patient-level evidence", fontsize=16, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG01_study_design_and_evidence_ladder", [metric_path])


def fig02_repeated_benchmark() -> None:
    path = RESULTS / "repeated_strict_fold_benchmark/repeat_level_metrics.tsv"
    pair_path = RESULTS / "repeated_strict_fold_benchmark/paired_repeat_auc_differences.tsv"
    df, pairs = read(path), read(pair_path)
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.5))
    ax = axes[0, 0]
    sns.violinplot(data=df, x="dataset_id", y="roc_auc", hue="method", palette=METHOD_COLORS, inner=None, cut=0, linewidth=0.7, ax=ax)
    sns.stripplot(data=df, x="dataset_id", y="roc_auc", hue="method", dodge=True, palette=METHOD_COLORS, size=2.6, alpha=0.6, ax=ax, legend=False)
    ax.set_xticklabels([COHORT_LABELS.get(t.get_text(), t.get_text()) for t in ax.get_xticklabels()], rotation=12, ha="right")
    ax.set_xlabel(""); ax.set_ylabel("ROC AUC"); ax.set_ylim(0.72, 1.01)
    ax.set_title("Repeated donor-level cross-validation")
    ax.legend(frameon=False, ncol=3, loc="lower right")

    ax = axes[0, 1]
    long = pairs.melt(id_vars=["dataset_id", "repeat"], value_vars=["delta_mean_hvg_minus_geneformer", "delta_pca_minus_geneformer"], var_name="comparison", value_name="delta")
    long["comparison"] = long.comparison.map({"delta_mean_hvg_minus_geneformer": "HVG − Geneformer", "delta_pca_minus_geneformer": "PCA − Geneformer"})
    sns.boxplot(data=long, x="dataset_id", y="delta", hue="comparison", palette=[GOLD, CORAL], width=0.58, fliersize=0, ax=ax)
    sns.stripplot(data=long, x="dataset_id", y="delta", hue="comparison", dodge=True, palette=[GOLD, CORAL], size=2.2, alpha=0.65, ax=ax, legend=False)
    ax.axhline(0, color=NAVY, lw=1)
    ax.set_xticklabels([COHORT_LABELS.get(t.get_text(), t.get_text()) for t in ax.get_xticklabels()], rotation=12, ha="right")
    ax.set_xlabel(""); ax.set_ylabel("Paired ΔAUC")
    ax.set_title("Within-split paired differences")
    ax.legend(frameon=False, loc="lower right")

    ax = axes[1, 0]
    mean = df.groupby(["dataset_id", "method"]).roc_auc.mean().unstack()
    rank = mean.rank(axis=1, ascending=False, method="average")
    sns.heatmap(mean, annot=True, fmt=".3f", cmap=LinearSegmentedColormap.from_list("auc", ["#F3F5F6", CYAN, NAVY]), vmin=0.8, vmax=1.0, cbar_kws={"label": "Mean AUC"}, ax=ax)
    ax.set_yticklabels([COHORT_LABELS.get(t.get_text(), t.get_text()) for t in ax.get_yticklabels()], rotation=0)
    ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_title("Mean discrimination across 20 splits")

    ax = axes[1, 1]
    stability = df.groupby(["dataset_id", "method"]).roc_auc.std().reset_index()
    sns.pointplot(data=stability, x="roc_auc", y="dataset_id", hue="method", palette=METHOD_COLORS, dodge=0.35, markers="o", linestyles="none", ax=ax)
    ax.set_yticklabels([COHORT_LABELS.get(t.get_text(), t.get_text()) for t in ax.get_yticklabels()])
    ax.set_xlabel("SD of AUC across splits"); ax.set_ylabel("")
    ax.set_title("Split sensitivity")
    ax.legend(frameon=False, ncol=1, loc="lower right")
    ax.grid(axis="x", color="#D9E0E5")
    fig.suptitle("Internal accuracy is high, but method ordering depends on cohort", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG02_repeated_within_cohort_benchmark", [path, pair_path])


def fig03_transfer() -> None:
    old_path = RESULTS / "strict_source_only_transfer_20260715/strict_source_only_transfer_metrics.tsv"
    old_pred = RESULTS / "strict_source_only_transfer_20260715/strict_source_only_transfer_predictions.tsv"
    new_path = ANALYSIS / "main_transfer_sharedcv_candidate/strict_source_only_transfer_metrics.tsv"
    new_pred = ANALYSIS / "main_transfer_sharedcv_candidate/strict_source_only_transfer_predictions.tsv"
    old, new, pred = read(old_path), read(new_path), read(new_pred)
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.7))
    ax = axes[0, 0]
    work = old.copy(); work["dir"] = work.direction.map(direction_label)
    ypos = np.arange(len(work))[::-1]
    for y, row in zip(ypos, work.itertuples()):
        c = METHOD_COLORS.get(row.method_label, GREY)
        ax.errorbar(row.roc_auc, y, xerr=[[row.roc_auc - row.roc_auc_ci_low], [row.roc_auc_ci_high - row.roc_auc]], fmt="o", color=c, ecolor=c, capsize=3)
        ax.text(row.roc_auc + 0.007, y, f"{row.roc_auc:.3f}", va="center", fontsize=7)
    ax.set_yticks(ypos, [f"{r.dir}\n{r.method_label}" for r in work.itertuples()])
    ax.set_xlim(0.65, 1.02); ax.set_xlabel("Target-cohort ROC AUC (95% bootstrap CI)")
    ax.axvline(0.5, color=GREY, ls="--", lw=0.8)
    ax.set_title("Strict source-only external transfer")
    ax.grid(axis="x", color="#D9E0E5")

    ax = axes[0, 1]
    old2 = old[["direction", "method_id", "roc_auc"]].assign(protocol="Canonical")
    new2 = new.rename(columns={"method_id": "method_id"})[["direction", "method_id", "roc_auc"]].assign(protocol="Shared source-CV folds")
    sens = pd.concat([old2, new2])
    sens["method"] = sens.method_id.map({"frozen_geneformer": "Frozen Geneformer", "source_hvg_pseudobulk": "HVG", "source_pca_pseudobulk": "PCA"})
    sens["dir"] = sens.direction.map(direction_label)
    for (direc, method), g in sens.groupby(["dir", "method"]):
        x = [0, 1] if len(g) == 2 else list(range(len(g)))
        vals = [g[g.protocol == p].roc_auc.iloc[0] for p in ["Canonical", "Shared source-CV folds"]]
        ax.plot(x, vals, marker="o", lw=1.6, color={"Frozen Geneformer": BLUE, "HVG": GOLD, "PCA": CORAL}[method], alpha=0.85)
        ax.text(1.03, vals[-1], f"{method} | {direc.split('→')[-1].strip()}", va="center", fontsize=6.8)
    ax.set_xticks([0, 1], ["Canonical", "Shared source-CV\nfold sensitivity"])
    ax.set_xlim(-0.15, 1.72); ax.set_ylim(0.84, 1.0); ax.set_ylabel("ROC AUC")
    ax.set_title("Protocol sensitivity preserves the ordering")
    ax.grid(axis="y", color="#D9E0E5")

    ax = axes[1, 0]
    direction = "SLE_GSE285773_CD4_to_SLE_GSE174188_CD4"
    wide = pred[pred.direction == direction].pivot(index=["donor_id", "y_true"], columns="method_id", values="prob_case").reset_index()
    ax.scatter(wide.source_pca_pseudobulk, wide.frozen_geneformer, c=np.where(wide.y_true == 1, CORAL, CYAN), s=23, alpha=0.7, edgecolor="white", linewidth=0.3)
    ax.plot([0, 1], [0, 1], color=GREY, ls="--", lw=0.8)
    ax.set_xlabel("PCA pseudobulk probability"); ax.set_ylabel("Geneformer probability")
    ax.set_title("Paired target-donor scores (large target)")
    ax.text(0.03, 0.94, "SLE", color=CORAL, transform=ax.transAxes, fontweight="bold")
    ax.text(0.03, 0.87, "Control", color=CYAN, transform=ax.transAxes, fontweight="bold")

    ax = axes[1, 1]
    calibration = new.copy(); calibration["dir"] = calibration.direction.map(direction_label)
    calibration["method"] = calibration.method_label.str.replace("Source-", "", regex=False)
    sns.scatterplot(data=calibration, x="brier", y="ece_10bin_supplementary", hue="method", style="dir", palette=METHOD_COLORS, s=95, ax=ax)
    for row in calibration.itertuples():
        ax.text(row.brier + 0.004, row.ece_10bin_supplementary + 0.003, f"AUC {row.roc_auc:.2f}", fontsize=6.5)
    ax.set_xlabel("Brier score (lower is better)"); ax.set_ylabel("ECE, 10 bins (lower is better)")
    ax.set_title("Discrimination does not guarantee calibration")
    ax.legend(frameon=False, fontsize=6.4, loc="upper left")
    fig.suptitle("External cohort shift separates ranking performance from probability reliability", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG03_cross_cohort_transfer_and_calibration", [old_path, old_pred, new_path, new_pred])


def fig04_geometry() -> None:
    coord_path = ANALYSIS / "5090_cap500_results/donor_mean_pca_coordinates.tsv"
    donor_path = ANALYSIS / "5090_cap500_results/donor_embedding_geometry.tsv"
    cell_path = ANALYSIS / "5090_cap500_results/cell_embedding_sample_geometry.tsv"
    cap_path = ANALYSIS / "cap_regularization_path/cap500_cap1000_embedding_geometry.tsv"
    coords, donor, cell, cap = map(read, [coord_path, donor_path, cell_path, cap_path])
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.3))
    for ax, dataset in zip(axes[0], coords.dataset.unique()):
        d = coords[coords.dataset == dataset]
        for lab, color in [("control", CYAN), ("case", CORAL)]:
            q = d[d.case_control == lab]
            ax.scatter(q.PC1, q.PC2, s=26, alpha=0.72, c=color, label="SLE" if lab == "case" else "Control", edgecolor="white", linewidth=0.35)
        ax.set_title(f"{COHORT_LABELS.get(dataset, dataset)} donor means")
        ax.set_xlabel("Donor embedding PC1"); ax.set_ylabel("PC2")
        ax.legend(frameon=False)

    ax = axes[1, 0]
    merged = donor.merge(cell, on="dataset", suffixes=("_donor", "_cell"))
    x = np.arange(len(merged)); width = 0.34
    ax.bar(x - width / 2, merged.entropy_effective_rank, width, color=BLUE, label="Donor means")
    ax.bar(x + width / 2, merged.entropy_effective_rank_top256, width, color=GOLD, label="Sampled cells (top 256 PCs)")
    ax.set_xticks(x, [COHORT_LABELS.get(v, v) for v in merged.dataset], rotation=10)
    ax.set_yscale("log"); ax.set_ylabel("Entropy effective rank (log scale)")
    ax.set_title("Mean pooling contracts cell-level geometry")
    ax.legend(frameon=False)
    for i, row in merged.iterrows():
        ax.text(i - width / 2, row.entropy_effective_rank * 1.15, f"{row.entropy_effective_rank:.1f}", ha="center", fontsize=7)
        ax.text(i + width / 2, row.entropy_effective_rank_top256 * 1.08, f"{row.entropy_effective_rank_top256:.1f}", ha="center", fontsize=7)

    ax = axes[1, 1]
    x = np.arange(len(cap)); width = 0.24
    ax.bar(x - width, cap.raw_cosine_median, width, color=CYAN, label="Raw cosine")
    ax.bar(x, cap.centered_cosine_median, width, color=BLUE, label="Centered cosine")
    ax.bar(x + width, cap.donor_distance_matrix_pearson, width, color=GREEN, label="Distance-matrix r")
    ax.set_xticks(x, [COHORT_LABELS.get(v, v) for v in cap.dataset], rotation=10)
    ax.set_ylim(0.985, 1.001); ax.set_ylabel("Geometry preservation")
    ax.set_title("500- and 1000-cell donor spaces nearly coincide")
    ax.legend(frameon=False, fontsize=6.8)
    fig.suptitle("Frozen donor means occupy a low-dimensional, cohort-sensitive subspace", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG04_donor_embedding_geometry", [coord_path, donor_path, cell_path, cap_path])


def fig05_scale_map() -> None:
    curve_path = RESULTS / "source_size_learning_curve_20260715/source_size_learning_curve_summary.tsv"
    budget_path = RESULTS / "current_geneformer_budget_stability/budget_stability.tsv"
    scgpt_path = RESULTS / "archived_scgpt_cell_budget_repeated/scgpt_cell_budget_repeated_summary.tsv"
    reg_path = ANALYSIS / "cap_regularization_path/cap500_cap1000_regularization_path.tsv"
    curve, budget, scgpt, reg = map(read, [curve_path, budget_path, scgpt_path, reg_path])
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.5))
    ax = axes[0, 0]
    curve["method"] = curve.method_id.map({"frozen_geneformer": "Frozen Geneformer", "source_hvg_pseudobulk": "Source-HVG pseudobulk", "source_pca_pseudobulk": "Source-PCA pseudobulk"})
    for method, g in curve.groupby("method"):
        g = g.sort_values("n_source_train")
        ax.plot(g.n_source_train, g.auc_mean, marker="o", color=METHOD_COLORS[method], label=method, lw=1.8)
        ax.fill_between(g.n_source_train, g.auc_q025, g.auc_q975, color=METHOD_COLORS[method], alpha=0.13)
    ax.set_xlabel("Source donors used for training"); ax.set_ylabel("Target ROC AUC")
    ax.set_title("Representation crossover with source sample size")
    ax.legend(frameon=False); ax.grid(color="#D9E0E5")

    ax = axes[0, 1]
    for (direction, capv), g in reg.groupby(["direction", "cell_cap"]):
        color = BLUE if capv == 500 else CORAL
        ls = "-" if "285773_CD4_to" in direction else "--"
        ax.plot(g.C, g.roc_auc, color=color, ls=ls, marker="o", ms=3, label=f"{capv} cells | {direction_label(direction)}")
    ax.set_xscale("log"); ax.set_xlabel("Logistic regularization C"); ax.set_ylabel("Target ROC AUC")
    ax.set_title("Cell-cap effects depend on source-only regularization")
    ax.legend(frameon=False, fontsize=6.5); ax.grid(color="#D9E0E5")

    ax = axes[1, 0]
    x = np.arange(len(budget)); width = 0.35
    ax.bar(x - width / 2, budget.auc_500, width, color=CYAN, label="500 cells")
    ax.bar(x + width / 2, budget.auc_1000, width, color=BLUE, label="1000 cells")
    ax.set_xticks(x, [COHORT_LABELS.get(v, v) for v in budget.dataset_id], rotation=12)
    ax.set_ylim(0.82, 1.0); ax.set_ylabel("Within-cohort ROC AUC")
    ax.set_title("Geneformer performance plateaus by 500 cells")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    d = scgpt[(scgpt.level.astype(str).str.match(r"^\d+$"))].copy()
    d["cells"] = d.level.astype(int)
    for method, g in d.groupby("method"):
        g = g.sort_values("cells")
        label = method.replace("_", " ").title()
        ax.plot(g.cells, g.mean_auc, marker="o", lw=1.5, label=label)
        ax.fill_between(g.cells, g.mean_auc - g.sd_auc, g.mean_auc + g.sd_auc, alpha=0.1)
    ax.set_xscale("log"); ax.set_xticks([25, 50, 100, 200, 500], [25, 50, 100, 200, 500])
    ax.set_xlabel("Cells per donor"); ax.set_ylabel("Mean ROC AUC")
    ax.set_title("Independent scGPT low-cell sensitivity")
    ax.legend(frameon=False, fontsize=6.6); ax.grid(color="#D9E0E5")
    fig.suptitle("A practical decision map spans donor count, cell count, and regularization", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG05_data_scale_and_regularization_map", [curve_path, budget_path, scgpt_path, reg_path])


def fig06_pooling() -> None:
    base = ANALYSIS / "5090_cap500_results/sharedcv_canonical_candidate"
    metric_path, pair_path, gate_path = base / "transfer_metrics.tsv", base / "paired_comparisons.tsv", base / "gate_decisions.tsv"
    cal_path = ANALYSIS / "5090_cap500_results/prediction_calibration_extremeness.tsv"
    metrics, pairs, gate, cal = map(read, [metric_path, pair_path, gate_path, cal_path])
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.0))
    ax = axes[0, 0]
    pivot = metrics.pivot(index="method", columns="direction", values="roc_auc")
    order = pivot.mean(axis=1).sort_values(ascending=False).index
    pivot = pivot.loc[order]
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap=LinearSegmentedColormap.from_list("pool", ["#F7F7F7", GOLD, CORAL, NAVY]), vmin=0.45, vmax=1.0, cbar_kws={"label": "ROC AUC"}, ax=ax)
    ax.set_xticklabels([direction_label(v) for v in pivot.columns], rotation=15, ha="right")
    ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title("Distribution-preserving pooling atlas")

    ax = axes[0, 1]
    q = pairs[pairs.method != "mean"].copy()
    q["dir"] = q.direction.map(direction_label)
    q = q.sort_values("delta_auc")
    y = np.arange(len(q))
    colors = np.where(q.delta_auc >= 0, GREEN, CORAL)
    ax.hlines(y, 0, q.delta_auc, color=colors, lw=2)
    ax.scatter(q.delta_auc, y, c=colors, s=30)
    ax.axvline(0, color=NAVY, lw=1)
    ax.set_yticks(y, [f"{m} | {d.split('→')[-1].strip()}" for m, d in zip(q.method, q.dir)])
    ax.set_xlabel("ΔAUC versus coordinate mean"); ax.set_title("No alternative clears both transfer directions")
    ax.grid(axis="x", color="#D9E0E5")

    ax = axes[1, 0]
    c = cal[cal.method.isin(metrics.method.unique())].copy()
    sns.scatterplot(data=c, x="roc_auc", y="brier", hue="method", style="target_dataset", palette="tab10", s=80, ax=ax)
    ax.axvline(0.5, color=GREY, ls="--", lw=0.8)
    ax.set_xlabel("ROC AUC"); ax.set_ylabel("Brier score")
    ax.set_title("High ranking can coexist with unusable probabilities")
    ax.legend(frameon=False, fontsize=5.8, ncol=2)

    ax = axes[1, 1]
    gate = gate.set_index("method")
    mat = gate[["delta_auc_target_261", "delta_auc_target_26"]].copy()
    mat.columns = ["ΔAUC: target n=261", "ΔAUC: target n=26"]
    sns.heatmap(mat, annot=True, fmt=".3f", center=0, cmap=sns.diverging_palette(240, 10, as_cmap=True), cbar_kws={"label": "ΔAUC"}, ax=ax)
    ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title("Pre-specified bidirectional promotion gate")
    fig.suptitle("Richer aggregation does not rescue frozen embeddings under external shift", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG06_pooling_stress_test", [metric_path, pair_path, gate_path, cal_path])


def fig07_ifn() -> None:
    base = RESULTS / "repeated_ifn_conditioned_geneformer_20260715"
    repeat_path = base / "repeat_level_ifn_metrics.tsv"
    null_path = base / "matched_random_module_repeat_metrics.tsv"
    null_summary_path = base / "matched_random_module_null_summary.tsv"
    score_path = RESULTS / "ifn_conditioned_geneformer/geneformer_ifn_donor_scores.tsv"
    repeat, nulls, null_summary, score = map(read, [repeat_path, null_path, null_summary_path, score_path])
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    ax = axes[0, 0]
    order = ["Geneformer", "ISG only", "IFN-residual Geneformer", "ISG plus residual Geneformer"]
    palette = [BLUE, GOLD, CORAL, GREEN]
    sns.boxplot(data=repeat, x="model", y="roc_auc", order=order, palette=palette, fliersize=0, ax=ax)
    sns.stripplot(data=repeat, x="model", y="roc_auc", order=order, color=NAVY, size=2.2, alpha=0.55, ax=ax)
    ax.set_xticklabels(["Geneformer", "ISG only", "IFN-residual", "ISG + residual"], rotation=12)
    ax.set_xlabel(""); ax.set_ylabel("ROC AUC"); ax.set_title("Repeated disease-program conditioning")

    ax = axes[0, 1]
    module_mean = nulls.groupby("module_id").delta_auc_from_geneformer.mean()
    ax.hist(module_mean, bins=35, color=CYAN, alpha=0.82, edgecolor="white")
    observed = nulls.ifn_delta_auc_same_split.mean()
    ax.axvline(observed, color=RED, lw=2.2, label=f"IFN mean ΔAUC = {observed:.3f}")
    ax.set_xlabel("Mean AUC decrease after module residualization")
    ax.set_ylabel("Matched random modules")
    ax.set_title("IFN decline is not unique among matched modules")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    genes = [c for c in score.columns if c.endswith("_log1p_cpm") and c != "isg_score_log1p_cpm_mean"]
    z = score[genes].apply(lambda x: (x - x.mean()) / (x.std() + 1e-9))
    order_idx = score.sort_values(["case_control", "isg_score_log1p_cpm_mean"]).index
    z = z.loc[order_idx]
    im = ax.imshow(z.T, aspect="auto", cmap=sns.diverging_palette(240, 10, as_cmap=True), vmin=-2.5, vmax=2.5, interpolation="nearest", rasterized=True)
    ax.set_yticks(np.arange(len(genes)), [g.replace("_log1p_cpm", "") for g in genes])
    ax.set_xticks([]); ax.set_xlabel("Donors ordered by phenotype and IFN score")
    ax.set_title("Coordinated 15-gene interferon axis")
    plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="Gene-wise z score")

    ax = axes[1, 1]
    s = null_summary.copy()
    repeats = np.arange(1, len(s) + 1)
    ax.plot(repeats, s.ifn_delta_auc, marker="o", color=RED, label="IFN")
    ax.fill_between(repeats, s.random_delta_mean, s.random_delta_q95, color=CYAN, alpha=0.3, label="Matched modules: mean to q95")
    ax.set_xticks(repeats)
    ax.set_xlabel("Repeated split set"); ax.set_ylabel("AUC decrease")
    ax.set_title("Program attribution is stable but not IFN-specific")
    ax.legend(frameon=False); ax.grid(color="#D9E0E5")
    fig.suptitle("Geneformer discrimination tracks a broad interferon-associated expression axis", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG07_interferon_and_matched_module_attribution", [repeat_path, null_path, null_summary_path, score_path])


def fig08_fusion() -> None:
    base = RESULTS / "source_only_representation_fusion_20260716"
    metric_path, pair_path, pred_path = base / "source_only_fusion_metrics.tsv", base / "source_only_fusion_paired_delong.tsv", base / "source_only_fusion_predictions.tsv"
    metrics, pairs, pred = map(read, [metric_path, pair_path, pred_path])
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.7))
    ax = axes[0, 0]
    metrics["dir"] = metrics.direction.map(direction_label)
    sns.barplot(data=metrics, x="method_label", y="roc_auc", hue="dir", palette=[NAVY, CORAL], ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=18, ha="right")
    ax.set_xlabel(""); ax.set_ylabel("Target ROC AUC"); ax.set_ylim(0.75, 1.0)
    ax.set_title("Early fusion does not exceed the strongest baseline")
    ax.legend(frameon=False, fontsize=6.6)

    ax = axes[0, 1]
    q = pairs.sort_values("delta_hybrid_minus_reference")
    y = np.arange(len(q)); colors = np.where(q.delta_hybrid_minus_reference >= 0, GREEN, CORAL)
    ax.hlines(y, 0, q.delta_hybrid_minus_reference, color=colors, lw=2)
    ax.scatter(q.delta_hybrid_minus_reference, y, c=colors, s=28)
    ax.axvline(0, color=NAVY, lw=1)
    ax.set_yticks(y, [f"{h.replace('_',' ')} vs {r.replace('_',' ')}" for h, r in zip(q.hybrid_method, q.reference_method)], fontsize=6.6)
    ax.set_xlabel("Hybrid − reference ΔAUC"); ax.set_title("Paired external comparisons")
    ax.grid(axis="x", color="#D9E0E5")

    ax = axes[1, 0]
    d = pred[pred.direction.str.contains("285773_CD4_to")]
    wide = d.pivot(index=["donor_id", "y_true"], columns="method_id", values="prob_case")
    corr = wide.drop(columns=[]).corr(method="spearman")
    sns.heatmap(corr, cmap=LinearSegmentedColormap.from_list("corr", ["#F7F7F7", CYAN, NAVY]), vmin=0, vmax=1, annot=True, fmt=".2f", cbar_kws={"label": "Spearman ρ"}, ax=ax)
    ax.set_title("Target-donor score concordance")
    ax.set_xlabel(""); ax.set_ylabel("")

    ax = axes[1, 1]
    sns.scatterplot(data=metrics, x="brier", y="roc_auc", hue="method_label", style="dir", palette=METHOD_COLORS, s=90, ax=ax)
    ax.set_xlabel("Brier score"); ax.set_ylabel("ROC AUC")
    ax.set_title("Fusion changes ranking and calibration together")
    ax.legend(frameon=False, fontsize=6.2, ncol=2)
    fig.suptitle("Simple concatenation adds features without adding transportable signal", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG08_source_only_fusion_sensitivity", [metric_path, pair_path, pred_path])


def fig09_validity_controls() -> None:
    base = RESULTS / "maximal_legacy_integration_20260716"
    cov_path, perm_path = base / "covariate_residualization_summary.tsv", base / "formal_permutation_summary.tsv"
    kernel_path, matrix_path = base / "kernel_fold_level_diagnostics.tsv", base / "legacy_multicohort_auc_matrix.tsv"
    cov, perm, kernel, matrix = map(read, [cov_path, perm_path, kernel_path, matrix_path])
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    ax = axes[0, 0]
    c = cov.dropna(subset=["unadjusted_auc"]).copy()
    x = np.arange(len(c))
    ax.scatter(c.unadjusted_auc, x, color=BLUE, s=60, label="Unadjusted")
    ax.scatter(c.residualized_auc, x, color=CORAL, s=60, label="Residualized")
    for i, row in c.reset_index(drop=True).iterrows():
        ax.plot([row.unadjusted_auc, row.residualized_auc], [i, i], color=GREY, lw=1.5)
    ax.set_yticks(x, c.method); ax.set_xlim(0.75, 1.0); ax.set_xlabel("ROC AUC")
    ax.set_title("Archived fold-contained covariate residualization")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    perm = perm.dropna(subset=["null_q025", "null_q975", "null_mean"]).copy()
    y = np.arange(len(perm))
    ax.hlines(y, perm.null_q025, perm.null_q975, color=CYAN, linewidth=5, alpha=0.65, label="Permutation 95% interval")
    ax.scatter(perm.null_mean, y, marker="x", color=RED, label="Permutation mean")
    ax.scatter(perm.observed_auc, y, marker="o", color=BLUE, s=45, label="Observed AUC")
    ax.set_yticks(y, perm.method_id.str.replace("_", " "))
    ax.set_xlabel("ROC AUC"); ax.set_title("Formal label-permutation checks")
    ax.legend(frameon=False); ax.grid(axis="x", color="#D9E0E5")

    ax = axes[1, 0]
    sns.scatterplot(data=kernel, x="effective_rank", y="feature_norm_cellcount_corr", hue="level", palette="viridis", s=75, ax=ax)
    ax.axhline(0, color=GREY, lw=0.8); ax.set_xlabel("Kernel effective rank"); ax.set_ylabel("Feature norm–cell count correlation")
    ax.set_title("Kernel geometry and cell-yield coupling")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    m = matrix.set_index("label")
    m = m.loc[m.mean(axis=1).sort_values(ascending=False).head(18).index]
    sns.heatmap(m, cmap=LinearSegmentedColormap.from_list("legacy", ["#F7F7F7", GOLD, CORAL, NAVY]), vmin=0.5, vmax=1.0, cbar_kws={"label": "Archived AUC"}, ax=ax)
    ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title("Historical multicohort method landscape")
    fig.suptitle("Validity controls reveal both biological signal and nuisance exposure", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG09_validity_controls_and_historical_landscape", [cov_path, perm_path, kernel_path, matrix_path])


def fig10_archived_resource_atlas() -> None:
    model_path = RESULTS / "archived_scaling_and_aggregation_20260716/geneformer_model_scale_maxlen.tsv"
    agg_path = RESULTS / "archived_scaling_and_aggregation_20260716/scgpt_aggregation_repeat_metrics.tsv"
    atlas_path = RESULTS / "maximal_legacy_integration_20260716/gse174188_legacy_pooling_atlas.tsv"
    full_path = RESULTS / "archived_scaling_and_aggregation_20260716/geneformer_full_cell_probe.tsv"
    model, agg, atlas, full = map(read, [model_path, agg_path, atlas_path, full_path])
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 8.7))
    ax = axes[0, 0]
    model["config"] = model.model.astype(str) + "M / " + model.max_len.astype(str)
    sns.scatterplot(data=model, x="max_len", y="auc", hue="model", size="sample_cap", palette=[GOLD, BLUE], sizes=(70, 150), ax=ax)
    ax.set_xlabel("Maximum sequence length"); ax.set_ylabel("ROC AUC")
    ax.set_title("Observed Geneformer configuration grid")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    sns.boxplot(data=agg, x="method", y="auc", palette="Set2", fliersize=0, ax=ax)
    sns.stripplot(data=agg, x="method", y="auc", color=NAVY, size=2.4, alpha=0.6, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=18, ha="right")
    ax.set_xlabel(""); ax.set_ylabel("ROC AUC"); ax.set_title("Repeated scGPT aggregation sensitivity")

    ax = axes[1, 0]
    a = atlas.sort_values("roc_auc", ascending=False)
    sns.scatterplot(data=a, x="roc_auc", y="label", hue="family", size="pr_auc", palette="tab10", sizes=(30, 110), ax=ax)
    ax.set_xlabel("ROC AUC"); ax.set_ylabel(""); ax.set_title("Twenty-two fixed-split pooling strategies")
    ax.legend(frameon=False, fontsize=5.7, ncol=2)

    ax = axes[1, 1]
    y = np.arange(len(full))
    ax.scatter(full.auc_cap1000, y, color=CYAN, s=70, label="1000 cells")
    ax.scatter(full.auc_all_available, y, color=BLUE, s=70, label="All available cells")
    for i, row in full.iterrows():
        ax.plot([row.auc_cap1000, row.auc_all_available], [i, i], color=GREY)
    ax.set_yticks(y, [COHORT_LABELS.get(v, v) for v in full.dataset])
    ax.set_xlabel("ROC AUC"); ax.set_title("Full-cell probes add little discrimination")
    ax.legend(frameon=False)
    fig.suptitle("Archived engineering sensitivities define a broad performance plateau", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG10_archived_resource_and_aggregation_atlas", [model_path, agg_path, atlas_path, full_path])


def fig11_pseudobulk_biology() -> None:
    de_dir = PROJECT / "04_重构_投稿定位_20260714/delivery/Pseudobulk_vs_Geneformer_SLE_visual_expansion_20260715/de_tables"
    paths = [
        de_dir / "SLE_GSE135779_pseudobulk_de_all_genes_20260709.tsv",
        de_dir / "SLE_GSE174188_CD4_pseudobulk_de_all_genes_20260709.tsv",
        de_dir / "SLE_GSE285773_CD4_pseudobulk_de_all_genes_20260709.tsv",
    ]
    tables = []
    for path in paths:
        d = read(path)
        p = d.welch_p.fillna(1).clip(lower=np.finfo(float).tiny).to_numpy()
        order = np.argsort(p)
        ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        q = np.empty_like(ranked); q[order] = np.clip(ranked, 0, 1)
        d["bh_q"] = q
        d["minus_log10_q"] = -np.log10(np.clip(q, 1e-50, 1))
        tables.append(d)

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.8))
    isg = {"ISG15", "IFI6", "MX1", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT3", "IFI44", "IFI44L", "STAT1", "RSAD2", "IFITM1", "IFITM3", "HERC5", "IFI27"}
    for ax, d, title in zip(np.ravel(axes)[:3], tables, ["GSE135779", "GSE174188 CD4", "GSE285773 CD4"]):
        sig = (d.bh_q < 0.05) & (d.abs_delta >= 0.5)
        colors = np.where(d.gene.isin(isg), RED, np.where(sig, GOLD, "#B7C0C8"))
        ax.scatter(
            d.delta_case_minus_control,
            d.minus_log10_q.clip(upper=50),
            c=colors,
            s=np.where(d.gene.isin(isg), 18, 7),
            alpha=0.72,
            linewidth=0,
            rasterized=True,
        )
        ax.axvline(0, color=GREY, lw=0.8); ax.axhline(-np.log10(0.05), color=GREY, ls="--", lw=0.8)
        labeled = d[d.gene.isin(isg)].nlargest(4, "abs_delta")
        for j, row in enumerate(labeled.itertuples()):
            ax.annotate(
                row.gene,
                (row.delta_case_minus_control, min(row.minus_log10_q, 50)),
                xytext=(5 + (j % 2) * 8, 7 + j * 4),
                textcoords="offset points",
                fontsize=6.5,
                arrowprops={"arrowstyle": "-", "color": GREY, "lw": 0.5},
            )
        ax.set_xlabel("SLE − control mean log1p CPM"); ax.set_ylabel("−log10 BH q")
        ax.set_title(title)

    ax = axes[1, 1]
    sets = [set(d.nlargest(200, "abs_delta").gene) for d in tables]
    combos = []
    for mask in range(1, 8):
        members = [i for i in range(3) if mask & (1 << i)]
        excluded = [i for i in range(3) if not mask & (1 << i)]
        genes = set.intersection(*(sets[i] for i in members))
        for i in excluded:
            genes -= sets[i]
        combos.append((mask, len(genes)))
    combos.sort(key=lambda x: x[1], reverse=True)
    x = np.arange(len(combos)); vals = [v for _, v in combos]
    ax.bar(x, vals, color=[NAVY if bin(m).count("1") > 1 else CYAN for m, _ in combos])
    ax.set_ylabel("Genes in exact intersection"); ax.set_xticks(x, [""] * len(x))
    ax.set_title("Top-200 absolute-effect overlap")
    ymax = max(vals) if vals else 1
    for j, (mask, val) in enumerate(combos):
        for i in range(3):
            ax.scatter(j, -0.14 * ymax - i * 0.09 * ymax, s=26, color=NAVY if mask & (1 << i) else "#D6DCE1", clip_on=False)
        included = [i for i in range(3) if mask & (1 << i)]
        if len(included) > 1:
            ax.plot([j, j], [-0.14 * ymax - min(included) * 0.09 * ymax, -0.14 * ymax - max(included) * 0.09 * ymax], color=NAVY, lw=1.2, clip_on=False)
    for i, label in enumerate(["GSE135779", "GSE174188", "GSE285773"]):
        ax.text(-0.65, -0.14 * ymax - i * 0.09 * ymax, label, ha="right", va="center", fontsize=7, clip_on=False)
    ax.set_ylim(-0.46 * ymax, ymax * 1.08)
    ax.set_yticks(np.arange(0, ymax * 1.01, 25))
    fig.suptitle("Pseudobulk recovers recurrent interferon-associated SLE expression", fontsize=15, fontweight="bold", color=NAVY, y=1.01)
    fig.tight_layout()
    panel_labels(fig, axes)
    finish(fig, "FIG11_pseudobulk_volcano_and_gene_overlap", paths)


def main() -> None:
    style()
    PROV.mkdir(parents=True, exist_ok=True)
    (PROV / "figure_source_index.tsv").write_text("figure\tsource_file\n", encoding="utf-8")
    builders = [
        fig01_overview,
        fig02_repeated_benchmark,
        fig03_transfer,
        fig04_geometry,
        fig05_scale_map,
        fig06_pooling,
        fig07_ifn,
        fig08_fusion,
        fig09_validity_controls,
        fig10_archived_resource_atlas,
        fig11_pseudobulk_biology,
    ]
    for builder in builders:
        print(f"Building {builder.__name__}")
        builder()


if __name__ == "__main__":
    main()
