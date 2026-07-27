#!/usr/bin/env python3
"""Scientific and packaging validation for the identifiability upgrade."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import WORKSPACE, sha256, write_json


RESULTS = WORKSPACE / "results"
FIGURES = WORKSPACE / "figures"


def check(condition: bool, message: str, checks: list[dict]) -> None:
    checks.append({"check": message, "passed": bool(condition)})


def main() -> None:
    checks: list[dict] = []

    simulation = pd.read_csv(
        RESULTS / "simulation" / "simulation_replicates.tsv.gz",
        sep="\t",
    )
    manifest = json.loads(
        (RESULTS / "simulation" / "simulation_manifest.json").read_text()
    )
    expected_rows = (
        len(manifest["betas"])
        * len(manifest["deltas"])
        * len(manifest["rhos"])
        * manifest["n_replicates"]
    )
    check(len(simulation) == expected_rows, "simulation row count", checks)
    equivalence = json.loads(
        (RESULTS / "simulation" / "observational_equivalence.json").read_text()
    )
    check(
        equivalence["max_absolute_difference"] == 0,
        "perfect-collinearity observational equivalence",
        checks,
    )
    info_by_rho = simulation.groupby("rho")["information_fraction"].mean()
    check(
        np.all(np.diff(info_by_rho.to_numpy()) <= 1e-10),
        "mean information fraction is non-increasing with rho",
        checks,
    )
    check(
        info_by_rho.loc[1.0] < 1e-12,
        "information fraction vanishes at perfect collinearity",
        checks,
    )
    no_biology_external = simulation[simulation["beta"].eq(0)].groupby(
        "rho"
    )["auc_external_source_only"].mean()
    check(
        float(np.max(np.abs(no_biology_external - 0.5))) < 0.03,
        "external AUC is near chance when biological effect is zero",
        checks,
    )

    observed = pd.read_csv(
        RESULTS / "gse135779_null" / "observed_design_metrics.tsv",
        sep="\t",
    )
    null_summary = pd.read_csv(
        RESULTS / "gse135779_null" / "design_null_summary.tsv",
        sep="\t",
    )
    complete_auc = observed.loc[
        observed["block"].eq("all"), "roc_auc"
    ].mean()
    check(
        abs(complete_auc - 0.952066) < 0.01,
        "GSE135779 complete-design AUC reproduces locked estimate",
        checks,
    )
    complete_label_null = null_summary[
        null_summary["block"].eq("all")
        & null_summary["null_type"].eq("label_permutation")
    ].iloc[0]
    check(
        complete_label_null["empirical_upper_p"] <= 0.01,
        "GSE135779 complete design exceeds label-permutation null",
        checks,
    )

    features = pd.read_csv(
        RESULTS
        / "composition"
        / "gse174188_cd4_composition_features.tsv.gz",
        sep="\t",
    )
    composition = pd.read_csv(
        RESULTS / "composition" / "composition_summary.tsv",
        sep="\t",
    )
    check(len(features) == 261, "composition donor count", checks)
    check(int(features["y_true"].sum()) == 162, "composition case count", checks)
    check(
        set(composition["variant"])
        == {
            "raw_proportions",
            "raw_plus_log_cells",
            "clr_composition",
            "clr_plus_log_cells",
        },
        "all composition variants present",
        checks,
    )
    wide = composition.pivot(
        index="variant",
        columns="adjustment",
        values="auc_mean",
    )
    check(
        bool((wide["unadjusted"] > 0.8).all()),
        "all composition variants have nontrivial internal AUC",
        checks,
    )
    check(
        bool(
            (
                wide["residual_processing_wave"] - wide["unadjusted"]
                < -0.1
            ).all()
        ),
        "wave residualisation attenuates every composition variant",
        checks,
    )

    for stem in (
        "Figure_ID1_simulated_identifiability_landscape",
        "Figure_ID2_empirical_null_and_composition",
    ):
        for suffix in ("svg", "pdf", "png"):
            path = FIGURES / f"{stem}.{suffix}"
            check(path.exists() and path.stat().st_size > 10_000, str(path.name), checks)

    failures = [item for item in checks if not item["passed"]]
    output_hashes = {}
    for path in sorted(RESULTS.rglob("*")) + sorted(FIGURES.rglob("*")):
        if path.is_file():
            output_hashes[str(path.relative_to(WORKSPACE))] = sha256(path)
    report = {
        "status": "passed" if not failures else "failed",
        "n_checks": len(checks),
        "n_failures": len(failures),
        "checks": checks,
        "output_sha256": output_hashes,
    }
    write_json(WORKSPACE / "VALIDATION_REPORT.json", report)
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
