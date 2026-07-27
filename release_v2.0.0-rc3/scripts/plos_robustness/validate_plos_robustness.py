#!/usr/bin/env python3
"""Validate the PLOS robustness extension and its release-facing claims."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "plos_robustness"
errors: list[str] = []


def require(condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


residual = pd.read_csv(OUT / "residualisation_sensitivity_summary.tsv", sep="\t")
design = pd.read_csv(OUT / "nonlinear_design_only_summary.tsv", sep="\t")
negative = pd.read_csv(OUT / "batch_exposure_negative_control.tsv", sep="\t")
intervals = pd.read_csv(OUT / "attenuation_difference_bootstrap.tsv", sep="\t")
positive = pd.read_csv(OUT / "synthetic_positive_control.tsv", sep="\t")
leakage = pd.read_csv(OUT / "confound_leakage_summary.tsv", sep="\t")

require(len(residual) == 15, "expected 15 representation-adjustment summaries")
require(
    set(residual["adjustment"])
    == {
        "unadjusted",
        "ridge_full_design",
        "random_forest_full_design",
        "batch_location_adjustment",
        "overlap_weighted_full_design",
    },
    "residualisation arms are incomplete",
)
for representation in ("geneformer", "hvg_pseudobulk", "pca_pseudobulk"):
    values = residual[residual["representation"].eq(representation)].set_index(
        "adjustment"
    )["mean"]
    require(
        values["random_forest_full_design"] < values["unadjusted"] - 0.07,
        f"nonlinear full-design residualisation did not attenuate {representation}",
    )
    require(
        abs(values["batch_location_adjustment"] - values["unadjusted"]) < 0.04,
        f"batch location adjustment unexpectedly changed {representation}",
    )

require(
    design.set_index("model").loc["logistic", "mean"] > 0.90,
    "locked logistic design-only result was not reproduced",
)
require(
    negative["best_batch_auc"].max() > 0.99,
    "near-perfect representation-to-batch exposure was not reproduced",
)
require(
    abs(negative["batch_only_disease_auc"].iloc[0] - 0.5) < 0.02,
    "batch-only disease negative control is not near chance",
)
require(
    negative["delta_batch_residualised_minus_unadjusted"].abs().max() < 0.01,
    "batch residualisation changed disease AUC by at least 0.01",
)
require(len(intervals) == 6, "expected six attenuation-gap intervals")
require(
    (intervals["loss_difference_bootstrap_p025"] > 0).all(),
    "an attenuation-gap interval crossed zero",
)
require(
    positive["passed"].all(),
    "the low-entanglement synthetic positive control did not pass all checks",
)
require(
    set(leakage["scenario"])
    == {
        "raw_random_forest",
        "global_residualisation_random_forest",
        "fold_contained_residualisation_random_forest",
    },
    "confound-leakage sensitivity scenarios are incomplete",
)

for suffix in ("svg", "pdf", "png"):
    path = ROOT / "figures" / "main" / f"Figure_6_residualisation_failure_modes.{suffix}"
    require(path.is_file() and path.stat().st_size > 0, f"missing {path.name}")

report = {
    "status": "passed" if not errors else "failed",
    "checks": 16,
    "errors": errors,
}
(OUT / "validation.json").write_text(
    json.dumps(report, indent=2) + "\n",
    encoding="utf-8",
)
print(json.dumps(report, indent=2))
raise SystemExit(1 if errors else 0)
