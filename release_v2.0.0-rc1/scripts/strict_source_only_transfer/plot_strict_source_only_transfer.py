#!/usr/bin/env python3
"""Render the release-facing strict source-only transfer figure."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
from path_config import FIGURES_ROOT, RESULTS_ROOT

RESULTS = RESULTS_ROOT / "strict_source_only_transfer"
METRICS = RESULTS / "strict_source_only_transfer_metrics.tsv"
COMPARISONS = RESULTS / "strict_source_only_transfer_paired_delong.tsv"
OUT = FIGURES_ROOT / "reproduced" / "strict_source_only_transfer"

METHODS = ("frozen_geneformer", "source_hvg_pseudobulk", "source_pca_pseudobulk")
LABELS = {
    "frozen_geneformer": "Frozen Geneformer",
    "source_hvg_pseudobulk": "Source-HVG\npseudobulk",
    "source_pca_pseudobulk": "Source-PCA\npseudobulk",
}
COLORS = {
    "frozen_geneformer": "#2F78B7",
    "source_hvg_pseudobulk": "#D76A39",
    "source_pca_pseudobulk": "#2FA38B",
}
DIRECTIONS = (
    "SLE_GSE174188_CD4_to_SLE_GSE285773_CD4",
    "SLE_GSE285773_CD4_to_SLE_GSE174188_CD4",
)
DISPLAY = {
    DIRECTIONS[0]: "GSE174188 CD4 -> GSE285773 CD4\nTarget n = 26",
    DIRECTIONS[1]: "GSE285773 CD4 -> GSE174188 CD4\nTarget n = 261",
}


def tidy(axis: plt.Axes) -> None:
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="x", color="#DCE5ED", linewidth=0.75)
    axis.set_axisbelow(True)


def panel_title(axis: plt.Axes, letter: str, title: str) -> None:
    axis.text(-0.18, 1.09, letter, transform=axis.transAxes, fontsize=15, fontweight="bold", va="top")
    axis.set_title(title, loc="left", fontsize=11.5, fontweight="bold", pad=9)


def main() -> None:
    metrics = pd.read_csv(METRICS, sep="\t")
    comparisons = pd.read_csv(COMPARISONS, sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.4), gridspec_kw={"width_ratios": [1.0, 1.0, 1.06]})

    for axis, direction in zip(axes[:2], DIRECTIONS):
        subset = metrics.loc[metrics["direction"].eq(direction)].set_index("method_id").loc[list(METHODS)]
        for position, (method_id, row) in enumerate(subset.iterrows()):
            axis.errorbar(
                row.roc_auc,
                position,
                xerr=[[row.roc_auc - row.roc_auc_ci_low], [row.roc_auc_ci_high - row.roc_auc]],
                fmt="o",
                color=COLORS[method_id],
                markersize=8,
                capsize=3.5,
                linewidth=1.6,
                zorder=3,
            )
            axis.text(0.716, position, f"{row.roc_auc:.3f}", va="center", fontsize=9, color=COLORS[method_id], fontweight="bold")
        axis.set_yticks(range(len(METHODS)), [LABELS[method] for method in METHODS], fontsize=8.6)
        axis.set_xlim(0.70, 1.012)
        axis.set_xlabel("Target-cohort ROC-AUC (95% CI)")
        axis.invert_yaxis()
        tidy(axis)
        panel_title(axis, "A" if direction == DIRECTIONS[0] else "B", DISPLAY[direction])

    axis = axes[2]
    axis.axvline(0, color="#6B7280", linestyle="--", linewidth=1.0, zorder=1)
    y_positions = []
    labels = []
    for group_index, direction in enumerate(DIRECTIONS):
        for offset, baseline in enumerate(("source_hvg_pseudobulk", "source_pca_pseudobulk")):
            row = comparisons.loc[(comparisons.direction.eq(direction)) & (comparisons.baseline_method_id.eq(baseline))].iloc[0]
            y = group_index * 3 + offset
            y_positions.append(y)
            labels.append(f"{DISPLAY[direction].splitlines()[0]}\n{LABELS[baseline].replace(chr(10), ' ')}")
            axis.scatter(row.delta_geneformer_minus_baseline, y, color=COLORS[baseline], s=58, zorder=3)
            axis.text(0.043, y, f"BH p={row.delong_p_bh_across_four_tests:.4g}", va="center", fontsize=8.4, color="#536273")
    axis.set_yticks(y_positions, labels, fontsize=7.7)
    axis.set_xlim(-0.105, 0.072)
    axis.set_xlabel("Paired delta AUC (Geneformer minus baseline)")
    axis.invert_yaxis()
    tidy(axis)
    panel_title(axis, "C", "Paired comparison on identical target donors")

    fig.suptitle("Cross-cohort transfer favors pseudobulk in both evaluated directions", fontsize=14.4, fontweight="bold", y=1.03)
    fig.text(
        0.5,
        -0.025,
        "Shared feature identifiers; source-cohort HVG selection, scaling, PCA, and C selection. "
        "Intervals: 5,000 donor-stratified bootstrap percentile CIs. Paired tests: DeLong with BH adjustment across four prespecified contrasts.",
        ha="center",
        fontsize=8.3,
        color="#52657A",
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(f"{OUT}.png", dpi=320, bbox_inches="tight")
    fig.savefig(f"{OUT}.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
