#!/usr/bin/env python3
"""Add GSE135779 year and full-design residualisation to the locked audit.

The primary locked run residualised both cohorts on the recorded batch variable.
For GSE135779, batch alone is not a useful cross-validated label predictor even
though collection year and the complete measured design are. This targeted
augmentation keeps the same outer folds and classifier pipeline and appends the
two scientifically relevant sensitivity analyses without rerunning GSE174188.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed

import locked_validity_audit as audit


OUT = audit.OUT
N_JOBS = min(int(os.environ.get("N_JOBS", "6")), os.cpu_count() or 1)
ADJUSTMENTS = ("residual_collection_year", "residual_all")


def main() -> None:
    cohort = audit.load_cohort_135779()
    completed = Parallel(
        n_jobs=N_JOBS,
        prefer="processes",
        max_nbytes="10M",
        mmap_mode="r",
        verbose=5,
    )(
        delayed(audit.evaluate_repeat)(
            cohort,
            repeat,
            ("collection_year", "all"),
        )
        for repeat in range(len(audit.SEEDS))
    )

    new_metrics = []
    new_predictions = []
    for metrics, predictions in completed:
        new_metrics.extend(
            row for row in metrics if row["adjustment"] in ADJUSTMENTS
        )
        new_predictions.extend(
            row for row in predictions if row["adjustment"] in ADJUSTMENTS
        )

    metrics_path = OUT / "repeat_metrics.tsv"
    predictions_path = OUT / "oof_predictions.tsv.gz"
    metrics = pd.read_csv(metrics_path, sep="\t")
    predictions = pd.read_csv(predictions_path, sep="\t")

    metrics = metrics[
        ~(
            metrics["cohort"].eq("GSE135779")
            & metrics["adjustment"].isin(ADJUSTMENTS)
        )
    ]
    predictions = predictions[
        ~(
            predictions["cohort"].eq("GSE135779")
            & predictions["adjustment"].isin(ADJUSTMENTS)
        )
    ]
    metrics = pd.concat([metrics, pd.DataFrame(new_metrics)], ignore_index=True)
    predictions = pd.concat(
        [predictions, pd.DataFrame(new_predictions)], ignore_index=True
    )

    metrics.to_csv(metrics_path, sep="\t", index=False)
    predictions.to_csv(predictions_path, sep="\t", index=False, compression="gzip")
    audit.summarise(metrics).to_csv(OUT / "summary.tsv", sep="\t", index=False)

    manifest_path = OUT / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["gse135779_additional_residualisation"] = {
        "blocks": ["collection_year", "all"],
        "implementation": (
            "Same outer folds, feature construction, scaling, inner C selection, "
            "ridge alpha and classifier as the primary locked run."
        ),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )

    selected = audit.summarise(metrics)
    selected = selected[
        selected["cohort"].eq("GSE135779")
        & selected["adjustment"].isin(ADJUSTMENTS)
    ]
    print(selected.to_string(index=False))
    print(f"Wrote augmented locked audit to {OUT}")


if __name__ == "__main__":
    main()
