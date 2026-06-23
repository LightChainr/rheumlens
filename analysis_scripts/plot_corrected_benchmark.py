#!/usr/bin/env python3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "fold_contained"
OUT = ROOT / "manuscript" / "figures_corrected"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {
    "scgpt": "#2F6B9A",
    "expression_pca": "#D9822B",
    "hvg_pseudobulk": "#3A8D5D",
    "scgpt_ifn_residual": "#8D63A8",
}
LABELS = {
    "scgpt": "Frozen scGPT",
    "expression_pca": "Expression PCA",
    "hvg_pseudobulk": "HVG pseudobulk",
}
COHORTS = ["GSE135779", "GSE285773", "GSE174188"]


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


df = pd.read_csv(RESULTS / "all_cohorts_auc_bootstrap_ci.csv")

# Corrected cohort AUC panel.
fig, ax = plt.subplots(figsize=(7.2, 4.4))
x = np.arange(len(COHORTS))
width = 0.22
for j, method in enumerate(["scgpt", "expression_pca", "hvg_pseudobulk"]):
    sub = df[df.comparison == method].set_index("cohort").loc[COHORTS]
    pos = x + (j - 1) * width
    vals = sub.estimate.to_numpy()
    low = vals - sub.ci_low.to_numpy()
    high = sub.ci_high.to_numpy() - vals
    ax.errorbar(pos, vals, yerr=[low, high], fmt="o", ms=7, capsize=3,
                lw=1.5, color=COLORS[method], label=LABELS[method])
ax.axhline(0.5, color="#777777", ls="--", lw=1)
ax.set_xticks(x, ["GSE135779\n(n=33/11)", "GSE285773\n(n=16/10)", "GSE174188\n(n=162/99)"])
ax.set_ylim(0.68, 1.015)
ax.set_ylabel("Donor-level OOF AUC (95% paired-bootstrap CI)")
ax.set_title("Fold-contained donor-level classification")
ax.legend(frameon=False, ncol=3, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "FIG1_corrected_cohort_auc")

# Corrected paired method differences.
fig, ax = plt.subplots(figsize=(7.2, 4.4))
comparisons = [
    ("expression_pca_minus_scgpt", "PCA − scGPT", "#D9822B"),
    ("hvg_pseudobulk_minus_scgpt", "HVG pseudobulk − scGPT", "#3A8D5D"),
]
for j, (key, label, color) in enumerate(comparisons):
    sub = df[df.comparison == key].set_index("cohort").loc[COHORTS]
    pos = x + (j - 0.5) * width
    vals = sub.estimate.to_numpy()
    ax.errorbar(pos, vals, yerr=[vals - sub.ci_low.to_numpy(), sub.ci_high.to_numpy() - vals],
                fmt="o", ms=7, capsize=3, lw=1.5, color=color, label=label)
ax.axhline(0, color="#333333", lw=1)
ax.set_xticks(x, COHORTS)
ax.set_ylabel("Paired AUC difference (95% CI)")
ax.set_title("No demonstrated incremental advantage of frozen scGPT")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "FIG3_corrected_paired_differences")

# Corrected IFN residual result and permutation null.
g = df[(df.cohort == "GSE135779") & df.comparison.isin(["scgpt", "scgpt_ifn_residual"])].set_index("comparison")
null = np.load(RESULTS / "GSE135779_ifn_residual_permutation_null.npy")
perm = pd.read_csv(RESULTS / "GSE135779_ifn_residual_permutation.csv").iloc[0]
fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), gridspec_kw={"width_ratios": [0.8, 1.4]})
methods = ["scgpt", "scgpt_ifn_residual"]
vals = g.loc[methods, "estimate"].to_numpy()
axes[0].errorbar([0, 1], vals,
                 yerr=[vals - g.loc[methods, "ci_low"].to_numpy(), g.loc[methods, "ci_high"].to_numpy() - vals],
                 fmt="o", ms=8, capsize=4, color="#2F6B9A")
axes[0].set_xticks([0, 1], ["Original", "ISG-residual"])
axes[0].set_ylim(0.48, 1.01)
axes[0].set_ylabel("OOF AUC (95% CI)")
axes[0].set_title("Fold-contained IFN removal")
axes[0].spines[["top", "right"]].set_visible(False)
axes[1].hist(null, bins=30, color="#C8D6E5", edgecolor="white")
axes[1].axvline(float(perm.observed_auc), color="#8D2E2E", lw=2,
                label=f"Observed={perm.observed_auc:.3f}\nP={perm.permutation_p:.3f}")
axes[1].set_xlabel("Permuted residual-model AUC")
axes[1].set_ylabel("Frequency")
axes[1].set_title("Donor-label permutation (1,000×)")
axes[1].legend(frameon=False)
axes[1].spines[["top", "right"]].set_visible(False)
fig.tight_layout()
save(fig, "FIG2_corrected_ifn_residual")

print(OUT)
