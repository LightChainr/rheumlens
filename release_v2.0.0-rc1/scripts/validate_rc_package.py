#!/usr/bin/env python3
"""Validate the self-contained v2.0.0-rc1 research-object package."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


VERSION = "2.0.0-rc1"
FIGURES = {
    "Figure_1_two_cohort_design_channels",
    "Figure_2_disease_and_design_information",
    "Figure_3_design_restriction_matched_controls",
    "Figure_4_residualisation_failure_modes",
    "Figure_5_generalisation_ladder",
    "Figure_6_learned_pooling",
}
REQUIRED = {
    "README.md",
    "CITATION.cff",
    ".zenodo.json",
    "manuscript/manuscript.md",
    "docs/METHODS_LOCK_20260726.md",
    "docs/LITERATURE_AND_NOVELTY_AUDIT_20260726.md",
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
    "submission/plos_computational_biology/AUTHOR_SUMMARY.md",
    "submission/plos_computational_biology/COVER_LETTER.md",
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
            expected_hash = locked.get("manuscript_sha256")
            observed_hash = sha256(manuscript_path)
            if expected_hash != observed_hash:
                fail(
                    errors,
                    "candidate manuscript hash differs from the locked scientific manuscript",
                )

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
            fail(errors, "Zenodo metadata version is not 2.0.0-rc1")
        creators = {item.get("name") for item in zenodo.get("creators", [])}
        expected = {"Ying, Hongyu", "Yun, Dandan", "Liu, Dan"}
        if creators != expected:
            fail(errors, "Zenodo creator metadata does not match the manuscript authors")

    for path in root.rglob("*"):
        if path.is_symlink():
            fail(errors, f"release candidate contains a symlink: {path.relative_to(root)}")

    files = [
        path
        for path in root.rglob("*")
        if path.is_file() and path.name not in GENERATED
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
