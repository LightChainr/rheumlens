#!/usr/bin/env python3
"""Validate the self-contained v2.0.0-rc3 research-object package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


VERSION = "2.0.0-rc3"
FIGURES = {
    "Figure_1_two_cohort_design_channels",
    "Figure_2_simulated_identifiability_landscape",
    "Figure_3_empirical_null_and_composition",
    "Figure_4_disease_and_design_information",
    "Figure_5_design_restriction_matched_controls",
    "Figure_6_residualisation_failure_modes",
    "Figure_7_generalisation_ladder",
    "Figure_8_learned_pooling",
}
REQUIRED = {
    "README.md",
    "REPRODUCE.md",
    "RELEASE_NOTES_v2.0.0-rc3.md",
    "CITATION.cff",
    ".zenodo.json",
    "environment-analysis-lock.yml",
    "environment-geneformer-lock.yml",
    "manuscript/manuscript.md",
    "docs/METHODS_LOCK_20260726.md",
    "docs/LITERATURE_AND_NOVELTY_AUDIT_20260726.md",
    "docs/INDEPENDENT_REPRODUCTION_20260727.md",
    "docs/TECHNICAL_AUDIT_20260727.md",
    "docs/TARGET_JOURNAL_STRATEGY_20260727.md",
    "results/locked_validity_audit/release_validation.json",
    "results/locked_validity_audit/summary.tsv",
    "results/locked_validity_audit/repeat_metrics.tsv",
    "results/locked_validity_audit/oof_predictions.tsv.gz",
    "results/gse135779_metadata/gse135779_donor_metadata_restored.tsv",
    "results/strict_source_only_transfer/strict_source_only_transfer_predictions.tsv",
    "results/strict_source_only_transfer/strict_source_only_transfer_metrics.tsv",
    "results/learned_pooling/learned_pooling_metrics.tsv",
    "results/learned_pooling/learned_pooling_paired_tests.tsv",
    "results/plos_robustness/residualisation_sensitivity_metrics.tsv",
    "results/plos_robustness/residualisation_sensitivity_summary.tsv",
    "results/plos_robustness/nonlinear_design_only_metrics.tsv",
    "results/plos_robustness/nonlinear_design_only_summary.tsv",
    "results/plos_robustness/batch_exposure_negative_control.tsv",
    "results/plos_robustness/attenuation_difference_bootstrap.tsv",
    "results/plos_robustness/synthetic_positive_control.tsv",
    "results/plos_robustness/confound_leakage_summary.tsv",
    "results/plos_robustness/manifest.json",
    "results/plos_robustness/validation.json",
    "inputs/design_metadata/gse174188_final_donor_covariates.tsv",
    "inputs/public_metadata/GSE135779_ST1b_donor_clinical.csv",
    "inputs/public_metadata/GSE135779_ST1c_sequencing.csv",
    "inputs/public_metadata/GSE135779_GEO_samples.csv",
    "inputs/public_metadata/GSE135779_study_name_to_donor_id.csv",
    "inputs/donor_level/SLE_GSE135779/donor_labels.tsv",
    "inputs/donor_level/SLE_GSE135779/donor_log1p_cpm.parquet",
    "inputs/donor_level/SLE_GSE135779/donor_embedding.parquet",
    "inputs/donor_level/SLE_GSE174188_CD4/donor_labels.tsv",
    "inputs/donor_level/SLE_GSE174188_CD4/donor_log1p_cpm.parquet",
    "inputs/donor_level/SLE_GSE174188_CD4/donor_embedding.parquet",
    "inputs/donor_level/SLE_GSE285773_CD4/donor_labels.tsv",
    "inputs/donor_level/SLE_GSE285773_CD4/donor_log1p_cpm.parquet",
    "inputs/donor_level/SLE_GSE285773_CD4/donor_embedding.parquet",
    "scripts/path_config.py",
    "scripts/build_plos_submission_pdf.py",
    "scripts/plos_robustness/run_plos_robustness.py",
    "scripts/plos_robustness/validate_plos_robustness.py",
    "scripts/plos_robustness/make_plos_robustness_figure.py",
    "submission/plos_computational_biology/AUTHOR_SUMMARY.md",
    "submission/plos_computational_biology/COVER_LETTER.md",
    "submission/plos_computational_biology/FIGURE_LEGENDS.md",
    "submission/plos_computational_biology/SUPPLEMENTARY_TABLE_INDEX.md",
    "submission/plos_computational_biology/Supplementary_Tables_v2.0.0-rc3.xlsx",
    "submission/plos_computational_biology/Design_identifiability_v2.0.0-rc3_initial_submission.pdf",
    "identifiability_extension/results/simulation/simulation_summary.tsv",
    "identifiability_extension/results/simulation/observational_equivalence.json",
    "identifiability_extension/results/gse135779_null/design_null_summary.tsv",
    "identifiability_extension/results/gse135779_null/design_information_fraction.tsv",
    "identifiability_extension/results/composition/composition_summary.tsv",
    "identifiability_extension/results/composition/composition_restriction_metrics.tsv",
    "identifiability_extension/docs/NOVELTY_AND_RELEASE_DECISION_20260727.md",
}
FORBIDDEN_MANUSCRIPT_PATTERNS = {
    "source-internal cross-validation overestimates independent-target AUC by 0.10",
    "pure wave 4 (n=89)",
    "only GSE174188 distributes donor-level design metadata",
    "5-10 times larger",
    "5–10 times larger",
}
GENERATED = {"RELEASE_MANIFEST.json", "SHA256SUMS"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def row_count(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def fail(errors: list[str], message: str) -> None:
    errors.append(message)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "root",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = parser.parse_args()
    root = args.root.resolve()
    errors: list[str] = []

    for rel in sorted(REQUIRED):
        path = root / rel
        if not path.is_file() or path.stat().st_size == 0:
            fail(errors, f"missing or empty required file: {rel}")

    for stem in sorted(FIGURES):
        for suffix in (".svg", ".pdf", ".png"):
            rel = f"figures/main/{stem}{suffix}"
            path = root / rel
            if not path.is_file() or path.stat().st_size == 0:
                fail(errors, f"missing or empty figure: {rel}")

    manuscript_path = root / "manuscript/manuscript.md"
    if manuscript_path.is_file():
        manuscript = manuscript_path.read_text(encoding="utf-8")
        for pattern in sorted(FORBIDDEN_MANUSCRIPT_PATTERNS):
            if pattern.lower() in manuscript.lower():
                fail(errors, f"forbidden stale manuscript claim: {pattern}")

        locked_path = root / "results/locked_validity_audit/release_validation.json"
        if locked_path.is_file():
            locked = json.loads(locked_path.read_text(encoding="utf-8"))
            if locked.get("status") != "passed":
                fail(errors, "copied scientific validation is not marked passed")
    metadata_path = (
        root
        / "results/gse135779_metadata/gse135779_donor_metadata_restored.tsv"
    )
    if metadata_path.is_file() and row_count(metadata_path) != 44:
        fail(errors, "GSE135779 restored metadata does not contain 44 donor rows")

    zenodo_path = root / ".zenodo.json"
    if zenodo_path.is_file():
        zenodo = json.loads(zenodo_path.read_text(encoding="utf-8"))
        if zenodo.get("version") != VERSION:
            fail(errors, "Zenodo metadata version is not 2.0.0-rc3")
        creators = {item.get("name") for item in zenodo.get("creators", [])}
        expected = {"Ying, Hongyu", "Yun, Dandan", "Liu, Dan"}
        if creators != expected:
            fail(errors, "Zenodo creator metadata does not match the manuscript authors")

    robustness_path = root / "results/plos_robustness/validation.json"
    if robustness_path.is_file():
        robustness = json.loads(robustness_path.read_text(encoding="utf-8"))
        if robustness.get("status") != "passed":
            fail(errors, "PLOS robustness validation is not marked passed")
        if robustness.get("checks") != 16:
            fail(errors, "PLOS robustness validation did not pass all 16 checks")

    author_summary_path = (
        root / "submission/plos_computational_biology/AUTHOR_SUMMARY.md"
    )
    if author_summary_path.is_file():
        author_summary = author_summary_path.read_text(encoding="utf-8")
        words = [
            token
            for token in author_summary.split()
            if token.lower() not in {"#", "author", "summary"}
        ]
        if not 150 <= len(words) <= 200:
            fail(
                errors,
                f"PLOS Author Summary has {len(words)} words; expected 150-200",
            )

    for path in root.rglob("*"):
        if "tmp" in path.relative_to(root).parts:
            continue
        if path.is_symlink():
            fail(errors, f"release candidate contains a symlink: {path.relative_to(root)}")

    forbidden_prefixes = ("/" + "Volumes" + "/", "/" + "Users" + "/")
    for path in (root / "scripts").rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        if any(prefix in text for prefix in forbidden_prefixes):
            fail(
                errors,
                f"machine-specific absolute path remains in {path.relative_to(root)}",
            )

    files = [
        path
        for path in root.rglob("*")
        if path.is_file()
        and path.name not in GENERATED
        and "tmp" not in path.relative_to(root).parts
        and ".pytest_cache" not in path.parts
        and "__pycache__" not in path.parts
        and path.suffix != ".pyc"
    ]
    records = [
        {
            "path": path.relative_to(root).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(files)
    ]
    manifest = {
        "version": VERSION,
        "status": "passed" if not errors else "failed",
        "file_count": len(records),
        "total_bytes": sum(item["bytes"] for item in records),
        "errors": errors,
        "files": records,
    }
    (root / "RELEASE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    (root / "SHA256SUMS").write_text(
        "".join(f"{item['sha256']}  {item['path']}\n" for item in records),
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "status": manifest["status"],
                "file_count": manifest["file_count"],
                "total_bytes": manifest["total_bytes"],
                "errors": errors,
            },
            indent=2,
        )
    )
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
