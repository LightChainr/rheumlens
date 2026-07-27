#!/usr/bin/env python3
"""Fail fast if headline manuscript values diverge from locked result tables."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))
from path_config import RELEASE_ROOT, RESULTS_ROOT

MANUSCRIPT = RELEASE_ROOT / "manuscript" / "manuscript.md"
SUMMARY = RESULTS_ROOT / "locked_validity_audit" / "summary.tsv"
SUBSET = (
    RESULTS_ROOT
    / "locked_validity_audit"
    / "gse135779_mixed_batch_matched.tsv"
)
OUT = RESULTS_ROOT / "locked_validity_audit" / "release_validation.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def value(
    table: pd.DataFrame,
    cohort: str,
    method: str,
    adjustment: str,
) -> float:
    row = table[
        table["cohort"].eq(cohort)
        & table["method"].eq(method)
        & table["adjustment"].eq(adjustment)
    ]
    if len(row) != 1:
        raise AssertionError(
            f"Expected one row for {cohort}/{method}/{adjustment}; found {len(row)}"
        )
    return float(row["auc_mean"].iloc[0])


def require_number(text: str, number: float, decimals: int = 3) -> None:
    token = f"{number:.{decimals}f}"
    if token not in text:
        raise AssertionError(f"Manuscript does not contain locked value {token}")


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    summary = pd.read_csv(SUMMARY, sep="\t")
    subset = pd.read_csv(SUBSET, sep="\t")

    checks: dict[str, float] = {}
    for cohort in ("GSE174188_CD4", "GSE135779"):
        for method in ("geneformer", "hvg_pseudobulk", "pca_pseudobulk"):
            raw = value(summary, cohort, method, "unadjusted")
            checks[f"{cohort}/{method}/unadjusted"] = raw
            require_number(text, raw)

    for method in ("geneformer", "hvg_pseudobulk", "pca_pseudobulk"):
        adjusted = value(summary, "GSE174188_CD4", method, "residual_batch")
        checks[f"GSE174188_CD4/{method}/residual_batch"] = adjusted
        require_number(text, adjusted)

        adjusted_all = value(summary, "GSE135779", method, "residual_all")
        checks[f"GSE135779/{method}/residual_all"] = adjusted_all
        require_number(text, adjusted_all)

    for cohort, method in (
        ("GSE174188_CD4", "covariates_batch"),
        ("GSE174188_CD4", "covariates_all"),
        ("GSE135779", "covariates_batch"),
        ("GSE135779", "covariates_collection_year"),
        ("GSE135779", "covariates_all"),
    ):
        locked = value(summary, cohort, method, "covariates_only")
        checks[f"{cohort}/{method}/covariates_only"] = locked
        require_number(text, locked)

    means = (
        subset.groupby(["stratum", "method"], as_index=False)["roc_auc"]
        .mean()
        .pivot(index="method", columns="stratum", values="roc_auc")
    )
    for method in means.index:
        delta = (
            means.loc[method, "mixed_batches_B2_B6"]
            - means.loc[method, "matched_random_for_mixed_batches_B2_B6"]
        )
        checks[f"GSE135779/{method}/restriction_delta"] = float(delta)
        require_number(text, delta)

    forbidden = [
        r"source-internal cross-validation overestimates independent-target AUC by 0\.10",
        r"pure wave 4 \(n=89\)",
        r"only GSE174188 distributes donor-level design metadata",
        r"5-10 times larger",
        r"5–10 times larger",
    ]
    for pattern in forbidden:
        if re.search(pattern, text, flags=re.IGNORECASE):
            raise AssertionError(f"Forbidden stale claim remains: {pattern}")

    result = {
        "status": "passed",
        "manuscript": str(MANUSCRIPT),
        "manuscript_sha256": sha256(MANUSCRIPT),
        "summary_sha256": sha256(SUMMARY),
        "subset_sha256": sha256(SUBSET),
        "headline_values": checks,
        "forbidden_stale_claims": forbidden,
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
