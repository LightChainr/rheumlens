#!/usr/bin/env python3
"""Run the retained main transfer analysis with shared source-CV folds by direction."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import pandas as pd


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("main_transfer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    module = load_module(args.source_script)
    args.output.mkdir(parents=True, exist_ok=True)
    module.OUT = args.output
    original_fit = module.fit_classifier

    def fit_shared(source_x, source_y, target_x, direction_seed):
        direction_index = 0 if direction_seed < module.SEED + 100 else 1
        return original_fit(
            source_x,
            source_y,
            target_x,
            direction_seed=module.SEED + direction_index,
        )

    module.fit_classifier = fit_shared
    module.main()

    metrics_path = args.output / "strict_source_only_transfer_metrics.tsv"
    metrics = pd.read_csv(metrics_path, sep="\t")
    metrics["source_internal_cv"] = (
        "stratified, shuffled, ROC-AUC selection; one fixed shared seed per transfer direction"
    )
    metrics["source_cv_shared_across_methods"] = True
    metrics.to_csv(metrics_path, sep="\t", index=False)

    summary = {
        "status": "complete",
        "source_script": str(args.source_script),
        "source_cv_shared_across_methods": True,
        "direction_seeds": {
            "SLE_GSE174188_CD4_to_SLE_GSE285773_CD4": module.SEED,
            "SLE_GSE285773_CD4_to_SLE_GSE174188_CD4": module.SEED + 1,
        },
        "metrics": metrics[
            ["direction", "method_id", "roc_auc", "roc_auc_ci_low", "roc_auc_ci_high", "selected_c_source_internal_cv"]
        ].to_dict(orient="records"),
    }
    (args.output / "sharedcv_rerun_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
