#!/usr/bin/env python3
"""Validate the self-contained version 1.0.0 public archive."""

from __future__ import annotations

from pathlib import Path
import re
import xml.etree.ElementTree as ET

from PIL import Image
import pandas as pd


RELEASE = Path(__file__).resolve().parents[1]


def main() -> None:
    rows = []
    svgs = sorted((RELEASE / "figures/svg").glob("*.svg"))
    for svg in svgs:
        stem = svg.stem
        png = RELEASE / "figures/previews" / f"{stem}.png"
        root = ET.parse(svg).getroot()
        texts = ["".join(node.itertext()).strip() for node in root.iter() if node.tag.endswith("text")]
        with Image.open(png) as image:
            width, height = image.size
        rows.append(
            {
                "figure": stem,
                "svg_bytes": svg.stat().st_size,
                "svg_xml_valid": True,
                "editable_text_nodes": len(texts),
                "panel_A_to_D_present": all(letter in texts for letter in "ABCD"),
                "png_exists": png.exists(),
                "png_width": width,
                "png_height": height,
            }
        )
    figures = pd.DataFrame(rows)
    figures.to_csv(RELEASE / "provenance/figure_validation.tsv", sep="\t", index=False)

    manuscript = RELEASE / "manuscript/Cross_cohort_SLE_donor_representation_benchmark.md"
    links = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", manuscript.read_text(encoding="utf-8"))
    linked = []
    for link in links:
        target = (manuscript.parent / link).resolve()
        linked.append({"markdown_link": link, "resolved_path": str(target.relative_to(RELEASE)), "exists": target.exists()})
    pd.DataFrame(linked).to_csv(RELEASE / "provenance/manuscript_image_links.tsv", sep="\t", index=False)

    source_index = pd.read_csv(RELEASE / "provenance/figure_source_index.tsv", sep="\t")
    source_index["copied_to_source_data"] = source_index.source_file.map(lambda p: (RELEASE / "source_data" / p).exists())
    source_index.to_csv(RELEASE / "provenance/figure_source_validation.tsv", sep="\t", index=False)

    failures = []
    if len(figures) != 11:
        failures.append(f"expected 11 editable SVG figures, found {len(figures)}")
    if not figures[["svg_xml_valid", "panel_A_to_D_present", "png_exists"]].all().all():
        failures.append("one or more SVG/PNG figure pairs failed validation")
    if not all(row["exists"] for row in linked):
        failures.append("one or more Markdown image links are broken")
    if not source_index.copied_to_source_data.all():
        failures.append("one or more figure-source records are missing")

    report = [
        "# Version 1.0.0 archive validation",
        "",
        f"- Editable SVG figures: {len(figures)}",
        f"- Manuscript image links: {len(linked)} ({sum(row['exists'] for row in linked)} valid)",
        f"- Figure-source records: {len(source_index)} ({int(source_index.copied_to_source_data.sum())} copied)",
        f"- Status: {'PASS' if not failures else 'FAIL'}",
    ]
    if failures:
        report += ["", "## Failures", ""] + [f"- {item}" for item in failures]
    (RELEASE / "provenance/RELEASE_VALIDATION.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    if failures:
        raise SystemExit("; ".join(failures))


if __name__ == "__main__":
    main()
