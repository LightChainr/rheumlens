#!/usr/bin/env python3
"""Extract structured GEO SOFT and GSE135779 supplementary metadata tables."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PUBLIC = ROOT / "results" / "evidence_package" / "public_metadata"
OUT = PUBLIC / "structured"


def parse_soft(path: Path) -> tuple[dict, pd.DataFrame]:
    series: dict[str, list[str]] = {}
    samples = []
    current = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    samples.append(current)
                current = {"geo_accession": line.split(" = ", 1)[1]}
            elif line.startswith("!Series_") and " = " in line:
                key, value = line[1:].split(" = ", 1)
                series.setdefault(key, []).append(value)
            elif current is not None and line.startswith("!Sample_") and " = " in line:
                key, value = line[1:].split(" = ", 1)
                if key == "Sample_characteristics_ch1" and ": " in value:
                    name, val = value.split(": ", 1)
                    current[name.strip().lower().replace(" ", "_")] = val.strip()
                elif key == "Sample_relation":
                    if value.startswith("BioSample:"):
                        current["biosample"] = value.rsplit("/", 1)[-1]
                    elif value.startswith("SRA:"):
                        current["sra_experiment"] = value.split("term=", 1)[-1]
                else:
                    current[key.removeprefix("Sample_").lower()] = value
    if current:
        samples.append(current)
    return series, pd.DataFrame(samples)


def clean_table(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.dropna(axis=1, how="all").dropna(axis=0, how="all")
    frame.columns = [str(c).strip().replace("\n", " ") for c in frame.columns]
    return frame


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for accession in ("GSE135779", "GSE285773", "GSE174188"):
        series, samples = parse_soft(PUBLIC / f"{accession}_family.soft.gz")
        samples.to_csv(OUT / f"{accession}_GEO_samples.csv", index=False)
        pd.DataFrame([{"field": k, "value": v} for k, values in series.items() for v in values]).to_csv(
            OUT / f"{accession}_GEO_series.csv", index=False
        )

    xlsx = PUBLIC / "41590_2020_743_MOESM3_ESM.xlsx"
    sheet_specs = {
        "ST1a-Cilinical table": (3, "GSE135779_ST1a_cohort_summary.csv"),
        "ST1b-Clinical information ": (3, "GSE135779_ST1b_donor_clinical.csv"),
        "ST1c-Sequencing information ": (3, "GSE135779_ST1c_sequencing.csv"),
        "ST4a": (3, "GSE135779_ST4a_official_cluster_counts.csv"),
        "ST4b": (2, "GSE135779_ST4b_official_subcluster_counts.csv"),
        "ST4c": (3, "GSE135779_ST4c_combined_cluster_counts.csv"),
        "ST4d": (2, "GSE135779_ST4d_combined_subcluster_counts.csv"),
    }
    for sheet, (header, filename) in sheet_specs.items():
        clean_table(pd.read_excel(xlsx, sheet_name=sheet, header=header)).to_csv(
            OUT / filename, index=False
        )

    # Join supplement names (cSLE1 etc.) to project donor IDs in GEO titles.
    geo = pd.read_csv(OUT / "GSE135779_GEO_samples.csv")
    extracted = geo["title"].astype(str).str.extract(r"^(?P<study_name>[^ ]+) \[(?P<donor_id>[^]]+)\]")
    mapping = pd.concat([geo[["geo_accession", "title"]], extracted], axis=1)
    mapping.to_csv(OUT / "GSE135779_study_name_to_donor_id.csv", index=False)


if __name__ == "__main__":
    main()
