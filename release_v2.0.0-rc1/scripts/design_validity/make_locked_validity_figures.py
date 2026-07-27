#!/usr/bin/env python3
"""Render the locked two-cohort validity figures as SVG, PDF and PNG."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
from path_config import FIGURES_ROOT, RESULTS_ROOT

OUT = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else FIGURES_ROOT / "reproduced_locked"
)
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#2866b1"
ORANGE = "#e76f3c"
GREEN = "#16866f"
PURPLE = "#7353a6"
GOLD = "#d9a520"
RED = "#bf3f3f"
INK = "#171717"
MUTED = "#666666"
GRID = "#e4e4e4"
CASE = "#d95745"
CONTROL = "#4d86c6"
REP_COLORS = {
    "geneformer": BLUE,
    "hvg_pseudobulk": ORANGE,
    "pca_pseudobulk": GREEN,
}
REP_LABELS = {
    "geneformer": "Frozen Geneformer",
    "hvg_pseudobulk": "HVG pseudobulk",
    "pca_pseudobulk": "PCA pseudobulk",
}
COHORT_LABELS = {
    "GSE174188_CD4": "GSE174188 CD4 (n=261)",
    "GSE135779": "GSE135779 (n=44)",
}
REPS = list(REP_COLORS)

plt.rcParams.update(
    {
        "font.size": 10.5,
        "axes.titlesize": 12.5,
        "axes.labelsize": 10.5,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9.5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "svg.fonttype": "none",
    }
)


def clean(ax, axis: str = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=axis, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


def panel(ax, letter: str, title: str, subtitle: str | None = None) -> None:
    ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", pad=16 if subtitle else 8)
    if subtitle:
        ax.text(
            0,
            1.01,
            subtitle,
            transform=ax.transAxes,
            va="bottom",
            color=MUTED,
            fontsize=9,
        )


def save(fig, stem: str) -> None:
    for suffix in ("svg", "pdf", "png"):
        fig.savefig(
            OUT / f"{stem}.{suffix}",
            dpi=360,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)
    print(f"wrote {OUT / (stem + '.svg')}")


summary = pd.read_csv(
    RESULTS_ROOT / "locked_validity_audit" / "summary.tsv", sep="\t"
)
batch_pred = pd.read_csv(
    RESULTS_ROOT
    / "locked_validity_audit"
    / "representation_batch_predictability_summary.tsv",
    sep="\t",
)
subset135 = pd.read_csv(
    RESULTS_ROOT
    / "locked_validity_audit"
    / "gse135779_mixed_batch_matched.tsv",
    sep="\t",
)
meta135 = pd.read_csv(
    RESULTS_ROOT
    / "gse135779_metadata"
    / "gse135779_donor_metadata_restored.tsv",
    sep="\t",
)
ent174 = pd.read_csv(RESULTS_ROOT / "design_entanglement_metrics.tsv", sep="\t")
pure174 = pd.read_csv(RESULTS_ROOT / "pure_wave_stratum.tsv", sep="\t")
comp174 = pd.read_csv(
    RESULTS_ROOT / "raw_h5ad_wave_by_disease_cd4_only.tsv",
    sep="\t",
    index_col=0,
)


def stacked_design(ax, table: pd.DataFrame, labels: list[str], title: str) -> None:
    controls = table["control"].to_numpy()
    cases = table["case"].to_numpy()
    y = np.arange(len(table))
    ax.barh(y, controls, color=CONTROL, label="control", zorder=3)
    ax.barh(y, cases, left=controls, color=CASE, label="SLE", zorder=3)
    for yi, control, case in zip(y, controls, cases):
        if control:
            ax.text(control / 2, yi, str(int(control)), ha="center", va="center", color="white")
        if case:
            ax.text(
                control + case / 2,
                yi,
                str(int(case)),
                ha="center",
                va="center",
                color="white",
                fontweight="bold",
            )
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("donors")
    panel(ax, title[0], title[3:])
    clean(ax)


# Figure 1: two-cohort design channels
fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.5))

t174 = pd.DataFrame(
    {
        "control": comp174["normal"].to_numpy(dtype=int),
        "case": comp174["systemic lupus erythematosus"].to_numpy(dtype=int),
    }
)
stacked_design(
    axes[0, 0],
    t174,
    [f"wave {int(x)}" for x in comp174.index],
    "A  GSE174188 processing waves",
)
axes[0, 0].legend(frameon=False, ncol=2, loc="lower right")

t135 = pd.crosstab(meta135["batch"], meta135["case_control"]).rename_axis(None)
t135 = t135.reindex(columns=["control", "case"], fill_value=0)
stacked_design(
    axes[0, 1],
    t135,
    [str(x) for x in t135.index],
    "B  GSE135779 sequencing batches",
)

tyear = pd.crosstab(meta135["collection_year"], meta135["case_control"]).rename_axis(None)
tyear = tyear.reindex(columns=["control", "case"], fill_value=0)
stacked_design(
    axes[1, 0],
    tyear,
    [str(x) for x in tyear.index],
    "C  GSE135779 collection years",
)

ax = axes[1, 1]
blocks = [
    ("covariates_batch", "batch"),
    ("covariates_collection_year", "collection year"),
    ("covariates_qc", "sequencing/QC"),
    ("covariates_demographic", "demographics"),
    ("covariates_all", "all measured design"),
]
y = np.arange(len(blocks))
offset = {"GSE174188_CD4": -0.12, "GSE135779": 0.12}
colors = {"GSE174188_CD4": PURPLE, "GSE135779": GOLD}
for cohort in ("GSE174188_CD4", "GSE135779"):
    rows = summary[
        summary["cohort"].eq(cohort) & summary["adjustment"].eq("covariates_only")
    ].set_index("method")
    values, low, high = [], [], []
    for method, _ in blocks:
        if method in rows.index:
            values.append(rows.loc[method, "auc_mean"])
            low.append(rows.loc[method, "auc_p025"])
            high.append(rows.loc[method, "auc_p975"])
        else:
            values.append(np.nan)
            low.append(np.nan)
            high.append(np.nan)
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    yy = y[mask] + offset[cohort]
    ax.errorbar(
        values[mask],
        yy,
        xerr=[values[mask] - np.asarray(low)[mask], np.asarray(high)[mask] - values[mask]],
        fmt="o",
        color=colors[cohort],
        markersize=7,
        capsize=3,
        linewidth=1.4,
        label=COHORT_LABELS[cohort],
        zorder=3,
    )
ax.axvline(0.5, color=MUTED, linestyle=(0, (4, 3)), linewidth=1)
ax.set_yticks(y, [label for _, label in blocks])
ax.invert_yaxis()
ax.set_xlim(0.4, 1.01)
ax.set_xlabel("case/control ROC-AUC from design covariates only")
ax.legend(frameon=False, loc="lower right")
panel(ax, "D", "The study design can act as a label channel", "20 repeated donor-level five-fold analyses")
clean(ax)

fig.suptitle(
    "Figure 1  Outcome labels are entangled with study design in two widely reused SLE cohorts",
    fontsize=16,
    fontweight="bold",
    y=1.01,
)
fig.tight_layout()
save(fig, "Figure_1_two_cohort_design_channels")


# Figure 2: disease and design information coexist
fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw={"width_ratios": [1, 1, 1.25]})
for ax, cohort, letter in zip(
    axes[:2],
    ("GSE174188_CD4", "GSE135779"),
    ("A", "B"),
):
    disease = summary[
        summary["cohort"].eq(cohort) & summary["adjustment"].eq("unadjusted")
    ].set_index("method")
    bp = batch_pred[batch_pred["cohort"].eq(cohort)]
    best = bp.groupby("method")["auc_mean"].max()
    for rep in REPS:
        ax.plot(
            [disease.loc[rep, "auc_mean"], best.loc[rep]],
            [rep, rep],
            color=REP_COLORS[rep],
            linewidth=3,
            alpha=0.7,
            zorder=2,
        )
        ax.scatter(
            disease.loc[rep, "auc_mean"],
            rep,
            s=95,
            color=REP_COLORS[rep],
            edgecolor="white",
            linewidth=1.2,
            zorder=3,
        )
        ax.scatter(
            best.loc[rep],
            rep,
            s=95,
            facecolor="white",
            edgecolor=REP_COLORS[rep],
            linewidth=2.2,
            zorder=3,
        )
    ax.set_yticks(REPS, [REP_LABELS[r] for r in REPS])
    ax.set_xlim(0.65 if cohort == "GSE135779" else 0.94, 1.005)
    ax.set_xlabel("ROC-AUC")
    panel(
        ax,
        letter,
        COHORT_LABELS[cohort],
        "filled: disease; open: most predictable recorded batch",
    )
    clean(ax)

ax = axes[2]
matrix_rows = []
row_labels = []
for cohort in ("GSE174188_CD4", "GSE135779"):
    for rep in REPS:
        disease = summary[
            summary["cohort"].eq(cohort)
            & summary["method"].eq(rep)
            & summary["adjustment"].eq("unadjusted")
        ]["auc_mean"].iloc[0]
        best = batch_pred[
            batch_pred["cohort"].eq(cohort) & batch_pred["method"].eq(rep)
        ]["auc_mean"].max()
        matrix_rows.append([disease, best, best - disease])
        row_labels.append(f"{cohort.replace('_CD4', '')} | {REP_LABELS[rep]}")
matrix = np.asarray(matrix_rows)
cmap = LinearSegmentedColormap.from_list("validity", ["#f4efe3", GOLD, RED])
im = ax.imshow(matrix[:, :2], cmap=cmap, vmin=0.70, vmax=1.0, aspect="auto")
for i in range(matrix.shape[0]):
    for j in range(2):
        ax.text(j, i, f"{matrix[i, j]:.3f}", ha="center", va="center", fontweight="bold")
    ax.text(2.0, i, f"{matrix[i, 2]:+.3f}", ha="center", va="center", color=MUTED)
ax.set_xlim(-0.5, 2.5)
ax.set_xticks([0, 1, 2], ["disease", "best batch", "batch - disease"])
ax.set_yticks(np.arange(len(row_labels)), row_labels)
ax.tick_params(length=0)
for spine in ax.spines.values():
    spine.set_visible(False)
panel(ax, "C", "Representations retain design identity", "one-vs-rest batch prediction, five repeated donor-level CV")
fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03, label="ROC-AUC")

fig.suptitle(
    "Figure 2  High disease discrimination coexists with near-perfect design recognition",
    fontsize=16,
    fontweight="bold",
    y=1.02,
)
fig.tight_layout()
save(fig, "Figure_2_disease_and_design_information")


# Figure 3: restriction against matched donor controls
fig, axes = plt.subplots(1, 3, figsize=(17, 5.6), sharex=True)

def restriction_panel(ax, observed, matched, letter, title, subtitle):
    y = np.arange(len(REPS))
    for i, rep in enumerate(REPS):
        o = observed[rep]
        m = matched[rep]
        ax.plot([m["mean"], o["mean"]], [i, i], color=REP_COLORS[rep], linewidth=3, alpha=0.65)
        ax.errorbar(
            m["mean"],
            i,
            xerr=[[m["mean"] - m["low"]], [m["high"] - m["mean"]]],
            fmt="o",
            color=MUTED,
            markerfacecolor="white",
            capsize=3,
            zorder=3,
        )
        ax.scatter(o["mean"], i, s=100, color=REP_COLORS[rep], edgecolor="white", zorder=4)
        ax.text(
            min(1.005, max(o["mean"], m["mean"]) + 0.008),
            i,
            f"{o['mean'] - m['mean']:+.3f}",
            va="center",
            fontsize=9,
            color=REP_COLORS[rep],
            fontweight="bold",
        )
    ax.set_yticks(y, [REP_LABELS[r] for r in REPS])
    ax.set_xlim(0.72, 1.015)
    ax.set_xlabel("ROC-AUC")
    panel(ax, letter, title, subtitle)
    clean(ax)


obs89_rows = ent174[
    ent174["block"].eq("C_within_wave") & ent174["stratum"].eq("wave4")
].set_index("representation")
mat89_rows = ent174[ent174["block"].eq("C_size_matched_random")].set_index("representation")
obs89 = {
    r: {"mean": obs89_rows.loc[r, "auc_mean"]}
    for r in REPS
}
mat89 = {
    r: {
        "mean": mat89_rows.loc[r, "auc_mean"],
        "low": mat89_rows.loc[r, "auc_p025"],
        "high": mat89_rows.loc[r, "auc_p975"],
    }
    for r in REPS
}
restriction_panel(
    axes[0],
    obs89,
    mat89,
    "A",
    "GSE174188 dominant wave 4",
    "n=89; open point/range: size- and label-matched random controls",
)

obs66_rows = pure174[pure174["block"].eq("C2_pure_wave")].set_index("representation")
mat66_rows = pure174[pure174["block"].eq("C2_size_matched_random")].set_index("representation")
obs66 = {r: {"mean": obs66_rows.loc[r, "auc_mean"]} for r in REPS}
mat66 = {
    r: {
        "mean": mat66_rows.loc[r, "auc_mean"],
        "low": mat66_rows.loc[r, "auc_p025"],
        "high": mat66_rows.loc[r, "auc_p975"],
    }
    for r in REPS
}
restriction_panel(
    axes[1],
    obs66,
    mat66,
    "B",
    "GSE174188 pure wave 4",
    "n=66; donors have at least 99% of analysed cells in wave 4",
)

group135 = subset135.groupby(["stratum", "method"])["roc_auc"]
obs135 = {}
mat135 = {}
for rep in REPS:
    observed_values = group135.get_group(("mixed_batches_B2_B6", rep))
    matched_values = group135.get_group(("matched_random_for_mixed_batches_B2_B6", rep))
    obs135[rep] = {"mean": observed_values.mean()}
    mat135[rep] = {
        "mean": matched_values.mean(),
        "low": matched_values.quantile(0.025),
        "high": matched_values.quantile(0.975),
    }
restriction_panel(
    axes[2],
    obs135,
    mat135,
    "C",
    "GSE135779 mixed batches B2-B6",
    "n=36; all-case B1 removed and every retained batch contains both labels",
)

fig.suptitle(
    "Figure 3  Disease discrimination survives design restriction after controlling for lost sample size",
    fontsize=16,
    fontweight="bold",
    y=1.02,
)
fig.tight_layout()
save(fig, "Figure_3_design_restriction_matched_controls")


# Figure 4: residualisation failure modes across cohorts
fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.8))

def slope_adjustment(ax, cohort, adjustment, letter, title, subtitle):
    rows = summary[summary["cohort"].eq(cohort)].set_index(["method", "adjustment"])
    y = np.arange(len(REPS))
    for i, rep in enumerate(REPS):
        raw = rows.loc[(rep, "unadjusted"), "auc_mean"]
        adjusted = rows.loc[(rep, adjustment), "auc_mean"]
        ax.plot([0, 1], [raw, adjusted], color=REP_COLORS[rep], linewidth=2.8, alpha=0.8)
        ax.scatter([0, 1], [raw, adjusted], s=90, color=REP_COLORS[rep], edgecolor="white", zorder=3)
        ax.text(1.04, adjusted, f"{REP_LABELS[rep]}  {adjusted - raw:+.3f}", va="center", color=REP_COLORS[rep])
    ax.set_xticks([0, 1], ["unadjusted", adjustment.replace("residual_", "residualised: ")])
    ax.set_xlim(-0.12, 1.75)
    ax.set_ylim(0.45, 1.01)
    ax.set_ylabel("donor-level ROC-AUC")
    panel(ax, letter, title, subtitle)
    clean(ax, "y")


slope_adjustment(
    axes[0, 0],
    "GSE174188_CD4",
    "residual_batch",
    "A",
    "Batch residualisation removes 0.23-0.27 AUC",
    "batch-only disease AUC = 0.922",
)
slope_adjustment(
    axes[0, 1],
    "GSE135779",
    "residual_batch",
    "B",
    "Batch residualisation changes almost nothing",
    "batch-only disease AUC = 0.499",
)

ax = axes[1, 0]
available_adjustments = [
    ("residual_batch", "batch"),
    ("residual_collection_year", "collection year"),
    ("residual_all", "all measured design"),
]
x = np.arange(len(available_adjustments))
width = 0.23
for j, rep in enumerate(REPS):
    values = []
    raw = summary[
        summary["cohort"].eq("GSE135779")
        & summary["method"].eq(rep)
        & summary["adjustment"].eq("unadjusted")
    ]["auc_mean"].iloc[0]
    for adjustment, _ in available_adjustments:
        row = summary[
            summary["cohort"].eq("GSE135779")
            & summary["method"].eq(rep)
            & summary["adjustment"].eq(adjustment)
        ]
        values.append(np.nan if row.empty else row["auc_mean"].iloc[0] - raw)
    ax.bar(
        x + (j - 1) * width,
        values,
        width,
        color=REP_COLORS[rep],
        label=REP_LABELS[rep],
        zorder=3,
    )
ax.axhline(0, color=INK, linewidth=0.9)
ax.set_xticks(x, [label for _, label in available_adjustments])
ax.set_ylabel("change in ROC-AUC after residualisation")
ax.legend(frameon=False, ncol=1, loc="lower left")
panel(ax, "C", "The answer depends on which design block is projected out", "identical folds and classifier pipeline")
clean(ax, "y")

ax = axes[1, 1]
design_auc = {
    "GSE174188": 0.921993,
    "GSE135779 batch": 0.499242,
    "GSE135779 year": 0.720730,
    "GSE135779 all": 0.952066,
}
loss = {}
for label, cohort, adjustment in [
    ("GSE174188", "GSE174188_CD4", "residual_batch"),
    ("GSE135779 batch", "GSE135779", "residual_batch"),
    ("GSE135779 year", "GSE135779", "residual_collection_year"),
    ("GSE135779 all", "GSE135779", "residual_all"),
]:
    rows = summary[summary["cohort"].eq(cohort)].set_index(["method", "adjustment"])
    deltas = [
        rows.loc[(rep, "unadjusted"), "auc_mean"] - rows.loc[(rep, adjustment), "auc_mean"]
        for rep in REPS
    ]
    loss[label] = np.mean(deltas)
for label, color in zip(design_auc, (PURPLE, BLUE, GOLD, RED)):
    ax.scatter(design_auc[label], loss[label], s=150, color=color, edgecolor="white", linewidth=1.2, zorder=3)
    ax.text(design_auc[label] + 0.008, loss[label], label, va="center", fontsize=9)
ax.axhline(0, color=MUTED, linewidth=0.9)
ax.axvline(0.5, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)))
ax.set_xlim(0.45, 1.01)
ax.set_xlabel("ROC-AUC from the projected-out design block")
ax.set_ylabel("mean representation AUC removed")
panel(ax, "D", "Design predictability forecasts residualisation damage", "descriptive relationship across prespecified blocks")
clean(ax)

fig.suptitle(
    "Figure 4  Residualisation is a sensitivity operation, not a causal batch-effect estimate",
    fontsize=16,
    fontweight="bold",
    y=1.01,
)
fig.tight_layout()
save(fig, "Figure_4_residualisation_failure_modes")
