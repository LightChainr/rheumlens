#!/usr/bin/env python3
"""Rebuild Figure 6 with nonlinear and uncertainty robustness panels."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "main"
RESULTS = ROOT / "results"

BLUE = "#2866b1"
ORANGE = "#e76f3c"
GREEN = "#16866f"
PURPLE = "#7353a6"
GOLD = "#d9a520"
RED = "#bf3f3f"
INK = "#171717"
MUTED = "#666666"
GRID = "#e4e4e4"
COLORS = {
    "geneformer": BLUE,
    "hvg_pseudobulk": ORANGE,
    "pca_pseudobulk": GREEN,
}
LABELS = {
    "geneformer": "Frozen Geneformer",
    "hvg_pseudobulk": "HVG pseudobulk",
    "pca_pseudobulk": "PCA pseudobulk",
}
REPS = list(COLORS)

plt.rcParams.update(
    {
        "font.size": 10.5,
        "axes.titlesize": 12,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 9.2,
        "ytick.labelsize": 9.2,
        "legend.fontsize": 9,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "svg.fonttype": "none",
    }
)


def clean(ax, axis: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=axis, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def panel(ax, letter: str, title: str, subtitle: str | None = None) -> None:
    ax.set_title(
        f"{letter}  {title}",
        loc="left",
        fontweight="bold",
        pad=15 if subtitle else 8,
    )
    if subtitle:
        ax.text(
            0,
            1.01,
            subtitle,
            transform=ax.transAxes,
            color=MUTED,
            fontsize=8.8,
            va="bottom",
        )


summary = pd.read_csv(
    RESULTS / "locked_validity_audit" / "summary.tsv",
    sep="\t",
)
robust = pd.read_csv(
    RESULTS / "plos_robustness" / "residualisation_sensitivity_summary.tsv",
    sep="\t",
)
intervals = pd.read_csv(
    RESULTS / "plos_robustness" / "attenuation_difference_bootstrap.tsv",
    sep="\t",
)

fig, axes = plt.subplots(2, 3, figsize=(18.5, 10.2))


def slope_adjustment(ax, cohort, adjustment, letter, title, subtitle):
    rows = summary[summary["cohort"].eq(cohort)].set_index(["method", "adjustment"])
    for rep in REPS:
        raw = rows.loc[(rep, "unadjusted"), "auc_mean"]
        adjusted = rows.loc[(rep, adjustment), "auc_mean"]
        ax.plot([0, 1], [raw, adjusted], color=COLORS[rep], linewidth=2.6, alpha=0.82)
        ax.scatter(
            [0, 1],
            [raw, adjusted],
            s=75,
            color=COLORS[rep],
            edgecolor="white",
            zorder=3,
        )
        ax.text(
            1.04,
            adjusted,
            f"{LABELS[rep]}  {adjusted - raw:+.3f}",
            va="center",
            color=COLORS[rep],
            fontsize=8.8,
        )
    ax.set_xticks(
        [0, 1],
        ["unadjusted", adjustment.replace("residual_", "residualised: ")],
    )
    ax.set_xlim(-0.12, 1.8)
    ax.set_ylim(0.45, 1.01)
    ax.set_ylabel("donor-level ROC-AUC")
    panel(ax, letter, title, subtitle)
    clean(ax)


slope_adjustment(
    axes[0, 0],
    "GSE174188_CD4",
    "residual_batch",
    "A",
    "Wave residualisation removes 0.23-0.27 AUC",
    "wave-only disease AUC = 0.922",
)
slope_adjustment(
    axes[0, 1],
    "GSE135779",
    "residual_batch",
    "B",
    "Batch exposure is not batch dependence",
    "best representation-to-batch AUC = 0.993; batch-only disease AUC = 0.499",
)

ax = axes[0, 2]
adjustments = [
    ("residual_batch", "batch"),
    ("residual_collection_year", "collection year"),
    ("residual_all", "all design"),
]
x = np.arange(len(adjustments))
width = 0.23
for offset, rep in enumerate(REPS):
    raw = summary[
        summary["cohort"].eq("GSE135779")
        & summary["method"].eq(rep)
        & summary["adjustment"].eq("unadjusted")
    ]["auc_mean"].iloc[0]
    values = []
    for adjustment, _ in adjustments:
        adjusted = summary[
            summary["cohort"].eq("GSE135779")
            & summary["method"].eq(rep)
            & summary["adjustment"].eq(adjustment)
        ]["auc_mean"].iloc[0]
        values.append(adjusted - raw)
    ax.bar(
        x + (offset - 1) * width,
        values,
        width,
        color=COLORS[rep],
        label=LABELS[rep],
        zorder=3,
    )
ax.axhline(0, color=INK, linewidth=0.9)
ax.set_xticks(x, [label for _, label in adjustments])
ax.set_ylabel("change in ROC-AUC")
ax.legend(frameon=False, loc="lower left")
panel(
    ax,
    "C",
    "Damage follows the projected design block",
    "identical folds and classifier pipeline",
)
clean(ax)

ax = axes[1, 0]
design_auc = {
    "GSE174188 wave": 0.921993,
    "GSE135779 batch": 0.499242,
    "GSE135779 year": 0.720730,
    "GSE135779 all": 0.952066,
}
for label, cohort, adjustment, color in [
    ("GSE174188 wave", "GSE174188_CD4", "residual_batch", PURPLE),
    ("GSE135779 batch", "GSE135779", "residual_batch", BLUE),
    ("GSE135779 year", "GSE135779", "residual_collection_year", GOLD),
    ("GSE135779 all", "GSE135779", "residual_all", RED),
]:
    rows = summary[summary["cohort"].eq(cohort)].set_index(["method", "adjustment"])
    loss = np.mean(
        [
            rows.loc[(rep, "unadjusted"), "auc_mean"]
            - rows.loc[(rep, adjustment), "auc_mean"]
            for rep in REPS
        ]
    )
    ax.scatter(
        design_auc[label],
        loss,
        s=120,
        color=color,
        edgecolor="white",
        linewidth=1.1,
        zorder=3,
    )
    ax.text(design_auc[label] + 0.008, loss, label, va="center", fontsize=8.6)
ax.axhline(0, color=MUTED, linewidth=0.9)
ax.axvline(0.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)))
ax.set_xlim(0.45, 1.01)
ax.set_xlabel("disease AUC from projected design block")
ax.set_ylabel("mean representation AUC removed")
panel(
    ax,
    "D",
    "Design-label alignment forecasts attenuation",
    "descriptive relationship across prespecified blocks",
)
clean(ax)

ax = axes[1, 1]
order = [
    ("unadjusted", "none"),
    ("batch_location_adjustment", "batch\nlocation"),
    ("overlap_weighted_full_design", "overlap\nweights"),
    ("random_forest_full_design", "forest\nresiduals"),
    ("ridge_full_design", "ridge\nresiduals"),
]
x = np.arange(len(order))
for rep in REPS:
    selected = robust[robust["representation"].eq(rep)].set_index("adjustment")
    means = np.asarray([selected.loc[key, "mean"] for key, _ in order])
    lower = means - np.asarray([selected.loc[key, "p025"] for key, _ in order])
    upper = np.asarray([selected.loc[key, "p975"] for key, _ in order]) - means
    ax.errorbar(
        x,
        means,
        yerr=np.vstack([lower, upper]),
        color=COLORS[rep],
        marker="o",
        markersize=5.5,
        linewidth=2,
        capsize=2.5,
        label=LABELS[rep],
        zorder=3,
    )
ax.set_xticks(x, [label for _, label in order])
ax.set_ylim(0.45, 1.01)
ax.set_ylabel("donor-level ROC-AUC")
ax.legend(frameon=False, loc="lower left")
panel(
    ax,
    "E",
    "The conclusion survives nonlinear adjustment",
    "GSE135779; bars are 2.5th-97.5th split percentiles",
)
clean(ax)

ax = axes[1, 2]
positions = []
labels = []
colors = []
for rep_index, rep in enumerate(REPS):
    for stratum_index, stratum in enumerate(("dominant_wave_4", "pure_wave_4")):
        positions.append(rep_index * 2.5 + stratum_index * 0.75)
        labels.append(
            f"{LABELS[rep]}\n"
            + ("dominant wave" if stratum_index == 0 else "pure wave")
        )
        colors.append(COLORS[rep])
for position, (_, row), color in zip(
    positions,
    intervals.sort_values(["representation", "restriction_stratum"]).iterrows(),
    colors,
):
    estimate = row["loss_difference"]
    ax.errorbar(
        estimate,
        position,
        xerr=np.asarray(
            [
                [estimate - row["loss_difference_bootstrap_p025"]],
                [row["loss_difference_bootstrap_p975"] - estimate],
            ]
        ),
        fmt="o",
        color=color,
        markersize=7,
        capsize=3,
        zorder=3,
    )
ax.axvline(0, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)))
ax.set_yticks(positions, labels)
ax.set_xlabel("residualisation loss minus matched restriction loss")
ax.set_xlim(-0.02, 0.29)
panel(
    ax,
    "F",
    "The attenuation gap remains positive",
    "paired repeated-split bootstrap; descriptive split-sensitivity interval",
)
clean(ax, "x")

fig.suptitle(
    "Figure 6  Residualisation measures sensitivity, not a causal batch effect",
    fontsize=17,
    fontweight="bold",
    y=1.01,
)
fig.tight_layout()
for suffix in ("svg", "pdf", "png"):
    fig.savefig(
        OUT / f"Figure_6_residualisation_failure_modes.{suffix}",
        dpi=360,
        bbox_inches="tight",
        facecolor="white",
    )
plt.close(fig)
print(OUT / "Figure_6_residualisation_failure_modes.svg")
