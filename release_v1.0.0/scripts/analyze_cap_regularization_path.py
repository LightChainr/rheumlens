#!/usr/bin/env python3
"""Compare cap500 and cap1000 Geneformer transfer across a fixed C path."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler


DATASETS = ("SLE_GSE174188_CD4", "SLE_GSE285773_CD4")
CS = np.logspace(-4, 4, 9)
BLUE = "#246B8E"
GOLD = "#C58A1C"
INK = "#20252B"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cap500-reference-root", type=Path, required=True)
    parser.add_argument("--cap1000-model-root", type=Path, required=True)
    parser.add_argument("--method-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    labels = {}
    embeddings = {}
    geometry_rows = []
    for dataset in DATASETS:
        labels[dataset] = (
            pd.read_csv(
                args.cap500_reference_root / dataset / "donor_labels.tsv",
                sep="\t",
                dtype={"donor_id": str},
            )
            .set_index("donor_id")["case_control"]
            .sort_index()
        )
        cap500 = pd.read_parquet(
            args.cap500_reference_root / dataset / "historical_cap500_donor_embedding.parquet"
        )
        cap1000 = pd.read_parquet(
            args.cap1000_model_root / dataset / args.method_id / "donor_embedding.parquet"
        )
        cap500.index = cap500.index.astype(str)
        cap1000.index = cap1000.index.astype(str)
        cap500 = cap500.loc[labels[dataset].index]
        cap1000 = cap1000.loc[labels[dataset].index]
        embeddings[(dataset, 500)] = cap500
        embeddings[(dataset, 1000)] = cap1000
        x = cap500.to_numpy(float)
        y = cap1000.to_numpy(float)
        raw_cosine = np.sum(x * y, axis=1) / (np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1))
        global_mean = np.vstack((x, y)).mean(axis=0)
        x_centered = x - global_mean
        y_centered = y - global_mean
        centered_cosine = np.sum(x_centered * y_centered, axis=1) / (
            np.linalg.norm(x_centered, axis=1) * np.linalg.norm(y_centered, axis=1)
        )
        x_distance = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
        y_distance = np.linalg.norm(y[:, None, :] - y[None, :, :], axis=2)
        upper = np.triu_indices(len(x), 1)
        geometry_rows.append(
            {
                "dataset": dataset,
                "n_donors": len(x),
                "raw_cosine_median": float(np.median(raw_cosine)),
                "raw_cosine_minimum": float(raw_cosine.min()),
                "centered_cosine_median": float(np.median(centered_cosine)),
                "centered_cosine_minimum": float(centered_cosine.min()),
                "coordinate_rmse_median": float(np.median(np.sqrt(np.mean((x - y) ** 2, axis=1)))),
                "donor_distance_matrix_pearson": float(np.corrcoef(x_distance[upper], y_distance[upper])[0, 1]),
            }
        )

    rows = []
    for source, target in ((DATASETS[0], DATASETS[1]), (DATASETS[1], DATASETS[0])):
        source_y = labels[source].eq("case").astype(int).to_numpy()
        target_y = labels[target].eq("case").astype(int).to_numpy()
        for cap in (500, 1000):
            scaler = StandardScaler().fit(embeddings[(source, cap)])
            source_x = scaler.transform(embeddings[(source, cap)])
            target_x = scaler.transform(embeddings[(target, cap)])
            for c in CS:
                classifier = LogisticRegression(
                    C=c,
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=20_000,
                ).fit(source_x, source_y)
                rows.append(
                    {
                        "direction": f"{source}_to_{target}",
                        "source_dataset": source,
                        "target_dataset": target,
                        "n_source": len(source_y),
                        "n_target": len(target_y),
                        "cell_cap": cap,
                        "C": c,
                        "roc_auc": roc_auc_score(target_y, classifier.predict_proba(target_x)[:, 1]),
                    }
                )
    path = pd.DataFrame(rows)
    path.to_csv(args.output / "cap500_cap1000_regularization_path.tsv", sep="\t", index=False)
    pd.DataFrame(geometry_rows).to_csv(
        args.output / "cap500_cap1000_embedding_geometry.tsv", sep="\t", index=False
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.edgecolor": "#AAB2BA",
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "pdf.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.18, top=0.78, wspace=0.22)
    for axis, direction in zip(axes, path["direction"].drop_duplicates()):
        subset = path[path["direction"].eq(direction)]
        for cap, color, marker in ((500, BLUE, "o"), (1000, GOLD, "s")):
            values = subset[subset["cell_cap"].eq(cap)]
            axis.plot(
                values["C"],
                values["roc_auc"],
                color=color,
                marker=marker,
                linewidth=1.6,
                markersize=4,
                label=f"{cap} cells per donor",
            )
        target = subset["target_dataset"].iloc[0].replace("SLE_", "").replace("_CD4", "")
        n_target = int(subset["n_target"].iloc[0])
        axis.set_xscale("log")
        axis.set_xlabel("Logistic-regression C")
        axis.set_ylabel("External ROC-AUC")
        axis.set_title(f"Target: {target} (n={n_target})", loc="left")
        axis.grid(color="#E8ECEF", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False)
    fig.suptitle(
        "Cell-cap comparison across a fixed regularization path",
        x=0.02,
        y=0.96,
        ha="left",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.02,
        0.035,
        "At a fixed C, cap500 and cap1000 transfer similarly; single-split C selection can create larger apparent differences.",
        color="#69727D",
        fontsize=8,
    )
    for suffix in ("pdf", "png"):
        fig.savefig(args.output / f"FIG_cap500_cap1000_regularization_path.{suffix}", dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
