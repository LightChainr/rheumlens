#!/usr/bin/env python3
"""Generate publication-ready figures for the identifiability upgrade."""

from __future__ import annotations

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from common import WORKSPACE


RESULTS = WORKSPACE / "results"
OUT = WORKSPACE / "figures"
INK = "#18212B"
BLUE = "#2878B5"
BLUE_LIGHT = "#A8CCE2"
ORANGE = "#E07A2D"
ORANGE_LIGHT = "#F3C39F"
GOLD = "#C89B2C"
GREY = "#A8B0B8"
LIGHT_GREY = "#E8ECEF"


def style() -> None:
    sns.set_theme(style="whitegrid")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.edgecolor": INK,
            "axes.linewidth": 0.8,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "axes.labelcolor": INK,
            "grid.color": LIGHT_GREY,
            "grid.linewidth": 0.7,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.11,
        1.075,
        label,
        transform=axis.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
    )


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def heatmap(
    axis: plt.Axes,
    data: pd.DataFrame,
    value: str,
    title: str,
    center: float | None = None,
    cmap: str = "vlag",
) -> None:
    matrix = data.pivot(index="delta", columns="rho", values=value).sort_index(
        ascending=False
    )
    sns.heatmap(
        matrix,
        ax=axis,
        cmap=cmap,
        center=center,
        annot=True,
        fmt=".2f",
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"shrink": 0.72},
        annot_kws={"fontsize": 8},
    )
    axis.set_title(title, loc="left")
    axis.set_xlabel("Design-label association (rho)")
    axis.set_ylabel("Technical effect (delta)")


def figure_simulation() -> None:
    raw = pd.read_csv(
        RESULTS / "simulation" / "simulation_replicates.tsv.gz",
        sep="\t",
    )
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.2))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    info = (
        raw.groupby("rho")["information_fraction"]
        .agg(
            mean="mean",
            low=lambda x: np.quantile(x, 0.025),
            high=lambda x: np.quantile(x, 0.975),
        )
        .reset_index()
    )
    ax_a.plot(info["rho"], info["mean"], color=BLUE, marker="o", lw=2.2)
    ax_a.fill_between(
        info["rho"],
        info["low"],
        info["high"],
        color=BLUE_LIGHT,
        alpha=0.55,
        linewidth=0,
    )
    ax_a.axhline(0, color=INK, lw=0.8, ls="--")
    ax_a.set(
        xlabel="Design-label association (rho)",
        ylabel="Design-adjusted label information",
        ylim=(-0.03, 1.04),
    )
    ax_a.set_title("Label information after design adjustment", loc="left", pad=8)
    panel_label(ax_a, "A")

    selected = (
        raw[raw["beta"].eq(1.0)]
        .assign(
            attenuation=lambda d: d["auc_internal_residualised"]
            - d["auc_internal_unadjusted"]
        )
        .groupby(["delta", "rho"], as_index=False)["attenuation"]
        .mean()
    )
    heatmap(
        ax_b,
        selected,
        "attenuation",
        "Internal AUC change after residualisation\n(biological effect = 1)",
        center=0,
        cmap="vlag",
    )
    panel_label(ax_b, "B")

    selected_external = (
        raw[raw["beta"].eq(1.0)]
        .assign(
            external_gap=lambda d: d["auc_external_source_only"]
            - d["auc_internal_unadjusted"]
        )
        .groupby(["delta", "rho"], as_index=False)["external_gap"]
        .mean()
    )
    heatmap(
        ax_c,
        selected_external,
        "external_gap",
        "External minus internal AUC\n(biological effect = 1)",
        center=0,
        cmap="vlag",
    )
    panel_label(ax_c, "C")

    coefficient = (
        raw[raw["beta"].gt(0)]
        .groupby(["beta", "rho"], as_index=False)["coefficient_l2_error"]
        .median()
    )
    for beta, group in coefficient.groupby("beta"):
        ax_d.plot(
            group["rho"],
            group["coefficient_l2_error"],
            marker="o",
            lw=2,
            label=f"beta = {beta:g}",
        )
    ax_d.set_yscale("log")
    ax_d.set(
        xlabel="Design-label association (rho)",
        ylabel="Median L2 error (log scale)",
    )
    ax_d.set_title("Biological coefficient recovery error", loc="left", pad=8)
    ax_d.legend(title="Biological effect")
    panel_label(ax_d, "D")

    fig.suptitle(
        "Simulated design-label identifiability landscape",
        x=0.06,
        y=0.985,
        ha="left",
        fontsize=17,
        fontweight="bold",
    )
    fig.text(
        0.06,
        0.01,
        (
            "Donor-level additive simulation; 100 replicates per parameter cell. "
            "External targets set rho=0. Heatmaps report mean AUC differences."
        ),
        fontsize=9,
        color="#4E5B66",
    )
    fig.tight_layout(rect=(0.04, 0.05, 0.99, 0.91), h_pad=3.0, w_pad=2.5)
    save(fig, "Figure_ID1_simulated_identifiability_landscape")


def figure_empirical() -> None:
    nulls = pd.read_csv(
        RESULTS / "gse135779_null" / "design_null_distributions.tsv.gz",
        sep="\t",
    )
    observed = pd.read_csv(
        RESULTS / "gse135779_null" / "observed_design_metrics.tsv",
        sep="\t",
    )
    information = pd.read_csv(
        RESULTS / "gse135779_null" / "design_information_fraction.tsv",
        sep="\t",
    )
    composition = pd.read_csv(
        RESULTS / "composition" / "composition_summary.tsv",
        sep="\t",
    )
    restriction = pd.read_csv(
        RESULTS / "composition" / "composition_restriction_metrics.tsv",
        sep="\t",
    )

    order = ["batch", "collection_year", "qc", "demographic", "all"]
    labels = {
        "batch": "Batch",
        "collection_year": "Collection year",
        "qc": "Sequencing/QC",
        "demographic": "Demographics",
        "all": "Complete design",
    }
    fig, axes = plt.subplots(2, 2, figsize=(13.3, 9.6))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    label_nulls = nulls[nulls["null_type"].eq("label_permutation")].copy()
    label_nulls["block"] = pd.Categorical(
        label_nulls["block"], categories=order, ordered=True
    )
    sns.violinplot(
        data=label_nulls,
        y="block",
        x="roc_auc",
        order=order,
        color=GREY,
        inner=None,
        cut=0,
        linewidth=0.6,
        ax=ax_a,
    )
    means = observed.groupby("block")["roc_auc"].mean()
    for index, block in enumerate(order):
        ax_a.scatter(
            means[block],
            index,
            s=72,
            color=BLUE,
            edgecolor="white",
            linewidth=0.9,
            zorder=5,
        )
    ax_a.axvline(0.5, color=INK, ls="--", lw=0.9)
    ax_a.set_yticks(np.arange(len(order)), [labels[x] for x in order])
    ax_a.set(
        xlabel="Cross-validated disease AUC",
        ylabel="",
        xlim=(0.25, 1.01),
    )
    ax_a.set_title("Design-only AUC versus permutation null", loc="left", pad=8)
    panel_label(ax_a, "A")

    info = information.set_index("block").loc[order].reset_index()
    obs = observed.groupby("block", as_index=False)["roc_auc"].mean()
    info = info.merge(obs, on="block")
    y_positions = np.arange(len(order))
    ax_b.hlines(
        y_positions,
        info["information_fraction"],
        info["roc_auc"],
        color=LIGHT_GREY,
        lw=3,
    )
    ax_b.scatter(
        info["information_fraction"],
        y_positions,
        color=ORANGE,
        marker="s",
        s=58,
        label="Information fraction",
        zorder=3,
    )
    ax_b.scatter(
        info["roc_auc"],
        y_positions,
        color=BLUE,
        marker="o",
        s=58,
        label="Design-only AUC",
        zorder=3,
    )
    ax_b.set_yticks(y_positions, [labels[x] for x in order])
    ax_b.invert_yaxis()
    ax_b.set(
        xlabel="Metric value",
        ylabel="",
        xlim=(-0.02, 1.02),
    )
    ax_b.set_title(
        "Design predictability versus remaining label information",
        loc="left",
        pad=8,
    )
    ax_b.legend(loc="lower right")
    panel_label(ax_b, "B")

    variant_order = [
        "raw_proportions",
        "raw_plus_log_cells",
        "clr_composition",
        "clr_plus_log_cells",
    ]
    variant_labels = {
        "raw_proportions": "Raw proportions",
        "raw_plus_log_cells": "Raw + cell yield",
        "clr_composition": "CLR composition",
        "clr_plus_log_cells": "CLR + cell yield",
    }
    for index, variant in enumerate(variant_order):
        rows = composition[composition["variant"].eq(variant)].set_index(
            "adjustment"
        )
        unadjusted = rows.loc["unadjusted"]
        residualised = rows.loc["residual_processing_wave"]
        ax_c.plot(
            [unadjusted["auc_mean"], residualised["auc_mean"]],
            [index, index],
            color=LIGHT_GREY,
            lw=3,
        )
        for row, color, marker in (
            (unadjusted, BLUE, "o"),
            (residualised, ORANGE, "s"),
        ):
            ax_c.errorbar(
                row["auc_mean"],
                index,
                xerr=[
                    [row["auc_mean"] - row["auc_p025"]],
                    [row["auc_p975"] - row["auc_mean"]],
                ],
                fmt=marker,
                color=color,
                capsize=3,
                markersize=6,
                zorder=4,
            )
    ax_c.axvline(0.5, color=INK, ls="--", lw=0.9)
    ax_c.set_yticks(
        np.arange(len(variant_order)),
        [variant_labels[x] for x in variant_order],
    )
    ax_c.invert_yaxis()
    ax_c.set(
        xlabel="Disease AUC",
        ylabel="",
        xlim=(0.45, 0.98),
    )
    ax_c.set_title("CD4 composition after wave adjustment", loc="left", pad=8)
    ax_c.scatter([], [], color=BLUE, marker="o", label="Unadjusted")
    ax_c.scatter([], [], color=ORANGE, marker="s", label="Residualised")
    ax_c.legend(loc="lower right")
    panel_label(ax_c, "C")

    paired = (
        restriction.pivot_table(
            index=["variant", "stratum", "seed"],
            columns="control_type",
            values="roc_auc",
        )
        .dropna()
        .reset_index()
    )
    paired["delta"] = paired["observed"] - paired["size_label_matched"]
    summaries = (
        paired.groupby(["variant", "stratum"])["delta"]
        .agg(
            mean="mean",
            low=lambda x: np.quantile(x, 0.025),
            high=lambda x: np.quantile(x, 0.975),
        )
        .reset_index()
    )
    positions = []
    labels_d = []
    colors = []
    for variant_index, variant in enumerate(variant_order):
        for offset, stratum in ((-0.15, "dominant_wave_4"), (0.15, "pure_wave_4")):
            row = summaries[
                summaries["variant"].eq(variant)
                & summaries["stratum"].eq(stratum)
            ].iloc[0]
            position = variant_index + offset
            positions.append(position)
            labels_d.append(row)
            colors.append(BLUE if stratum == "dominant_wave_4" else GOLD)
    for position, row, color in zip(positions, labels_d, colors):
        ax_d.errorbar(
            row["mean"],
            position,
            xerr=[
                [row["mean"] - row["low"]],
                [row["high"] - row["mean"]],
            ],
            fmt="o",
            color=color,
            capsize=3,
            markersize=6,
        )
    ax_d.axvline(0, color=INK, ls="--", lw=0.9)
    ax_d.set_yticks(
        np.arange(len(variant_order)),
        [variant_labels[x] for x in variant_order],
    )
    ax_d.invert_yaxis()
    ax_d.set(
        xlabel="Observed minus size/label-matched AUC",
        ylabel="",
    )
    ax_d.set_title("CD4 composition within processing wave 4", loc="left", pad=8)
    ax_d.scatter([], [], color=BLUE, label="Dominant wave 4")
    ax_d.scatter([], [], color=GOLD, label="Pure wave 4")
    ax_d.legend(loc="lower right")
    panel_label(ax_d, "D")

    fig.suptitle(
        "Empirical identifiability diagnostics and CD4 composition extension",
        x=0.06,
        y=0.985,
        ha="left",
        fontsize=17,
        fontweight="bold",
    )
    fig.text(
        0.06,
        0.01,
        (
            "GSE135779: 44 donors and full-pipeline nulls. GSE174188: "
            "261 donors; intervals are 2.5th-97.5th percentiles across 20 split seeds."
        ),
        fontsize=9,
        color="#4E5B66",
    )
    fig.tight_layout(rect=(0.04, 0.05, 0.99, 0.91), h_pad=3.0, w_pad=2.5)
    save(fig, "Figure_ID2_empirical_null_and_composition")


def main() -> None:
    style()
    figure_simulation()
    figure_empirical()
    print(f"Wrote identifiability figures to {OUT}")


if __name__ == "__main__":
    main()
