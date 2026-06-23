#!/usr/bin/env python3
"""Finalize report/figure/manifest for existing P9 structured transfer CSVs."""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path("/Users/lc/Documents/RheumLens")
OUT = ROOT / "manuscript/extension_results_20260623/P9_structured_transfer"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def md_table(df: pd.DataFrame) -> str:
    try:
        return df.to_markdown(index=False)
    except Exception:
        return "```\n" + df.to_string(index=False) + "\n```"


def main() -> None:
    summary = pd.read_csv(OUT / "EXT_P9_structured_transfer_summary.csv")
    failures_path = OUT / "EXT_P9_structured_transfer_failures.csv"
    try:
        failures = pd.read_csv(failures_path)
    except Exception:
        failures = pd.DataFrame()

    # Figure.
    plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False})
    plot = summary.copy()
    plot["label"] = plot["direction"].str.replace("_to_", " → ", regex=False) + "\n" + plot["method_id"]
    plot = plot.sort_values(["direction", "auc"])
    colors = []
    for m in plot["method_id"]:
        if m in {"donor_expression_pca", "donor_mean_hvg"}:
            colors.append("#4c78a8")
        elif m == "scgpt_mean":
            colors.append("#7f8c8d")
        elif m == "focus_lite@scgpt":
            colors.append("#b55242")
        elif m == "moments_mean_var@scgpt":
            colors.append("#8e6bbd")
        else:
            colors.append("#9aa6a6")
    fig, ax = plt.subplots(figsize=(8.0, 5.6), dpi=220)
    ax.barh(np.arange(len(plot)), plot["auc"], color=colors)
    ax.set_yticks(np.arange(len(plot)))
    ax.set_yticklabels(plot["label"])
    ax.set_xlim(0.35, 1.02)
    ax.set_xlabel("Source-only transfer AUC")
    ax.set_title("P9 extension: structured non-neural transfer methods")
    for y, v in enumerate(plot["auc"]):
        ax.text(v + 0.006, y, f"{v:.3f}", va="center", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(OUT / "EXT_FIG_P9_structured_transfer_auc.png")
    fig.savefig(OUT / "EXT_FIG_P9_structured_transfer_auc.pdf")
    plt.close(fig)

    best_by_dir = summary.sort_values("auc", ascending=False).groupby("direction").head(3)
    report = "# P9 structured source-only transfer extension\n\n"
    report += "Non-neural structured methods were fit on source donors only and scored on target donors before metric calculation. "
    report += "This extension does not include neural MIL models because their source-target training protocol requires a separate pre-registered transfer definition.\n\n"
    report += "## Key findings\n\n"
    report += (
        "- GSE285773 → GSE174188 CD4: expression baselines remained strongest; "
        "FOCUS was the best structured scGPT method but did not exceed donor_mean_hvg or donor_expression_pca.\n"
    )
    report += (
        "- GSE174188 CD4 → GSE285773: moments_mean_var@scGPT achieved the highest AUC in this extension, "
        "with focus_lite@scGPT comparable to donor_expression_pca. The pediatric target has only 26 donors, so uncertainty remains wide.\n"
    )
    report += "- KME did not improve transfer performance in either direction.\n\n"
    report += "## Top methods by direction\n\n"
    report += md_table(
        best_by_dir[
            [
                "direction",
                "method_id",
                "auc",
                "auc_ci_low",
                "auc_ci_high",
                "pr_auc",
                "brier",
                "n_donors",
            ]
        ]
    )
    report += "\n\n## Full summary\n\n"
    report += md_table(
        summary[
            [
                "direction",
                "method_id",
                "auc",
                "auc_ci_low",
                "auc_ci_high",
                "pr_auc",
                "brier",
                "n_donors",
                "status",
            ]
        ]
    )
    report += "\n\n## Failures\n\n"
    report += "No failures.\n" if failures.empty else md_table(failures)
    report += "\n\n## Claim boundary\n\n"
    report += (
        "These are retrospective source-only transfer analyses. They support method-prioritization hypotheses but do not establish clinical deployment, "
        "causal biology, or query-specific FOCUS mechanisms. Neural MIL transfer remains a future analysis requiring a separate protocol.\n"
    )
    (OUT / "EXT_P9_STRUCTURED_TRANSFER_REPORT.md").write_text(report)

    files = [p for p in sorted(OUT.rglob("*")) if p.is_file() and p.name != "MANIFEST_SHA256.tsv"]
    pd.DataFrame([{"path": str(p.relative_to(OUT)), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in files]).to_csv(
        OUT / "MANIFEST_SHA256.tsv", sep="\t", index=False
    )
    print(OUT / "EXT_P9_STRUCTURED_TRANSFER_REPORT.md")


if __name__ == "__main__":
    main()
