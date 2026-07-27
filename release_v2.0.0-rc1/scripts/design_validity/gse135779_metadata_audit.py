#!/usr/bin/env python3
"""Restore and audit donor-level design metadata for the GSE135779 analysis cohort.

The 44 childhood donors used in the classification benchmark are linked from the
internal donor IDs to the study sample names, Supplementary Table 1b clinical
metadata, Supplementary Table 1c sequencing metrics, and GEO sample metadata.
The script fails on any incomplete or ambiguous mapping and records source hashes.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve()
WORKSPACE = HERE.parents[1]
PROJECT = HERE.parents[2]
PUBLIC = (
    PROJECT
    / "10_zenodo_release_20260717"
    / "rheumlens"
    / "evidence_package"
    / "public_metadata"
)
EXTRA = (
    PROJECT
    / "03_远程回传"
    / "final_gpu_figure_ready_20260709"
    / "extracted_plus"
    / "RheumLens_GPU_figure_ready_plus_20260709"
    / "extra"
)
LABELS = EXTRA / "pseudobulk" / "SLE_GSE135779" / "donor_labels.tsv"
OUT = WORKSPACE / "results" / "gse135779_metadata"

CLINICAL = PUBLIC / "GSE135779_ST1b_donor_clinical.csv"
SEQUENCING = PUBLIC / "GSE135779_ST1c_sequencing.csv"
GEO = PUBLIC / "GSE135779_GEO_samples.csv"
MAPPING = PUBLIC / "GSE135779_study_name_to_donor_id.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def clean_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace({"ND": np.nan, "": np.nan}), errors="coerce")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    labels = pd.read_csv(LABELS, sep="\t", dtype=str)
    mapping = pd.read_csv(MAPPING, dtype=str)
    clinical = pd.read_csv(CLINICAL, dtype=str)
    sequencing = pd.read_csv(SEQUENCING, dtype=str)
    geo = pd.read_csv(GEO, dtype=str)

    require(labels["donor_id"].is_unique, "Analysis donor IDs are not unique")
    require(mapping["donor_id"].is_unique, "Study-name mapping has duplicate donor IDs")
    require(clinical["Names"].is_unique, "Clinical table has duplicate study names")
    require(sequencing["SampleID"].is_unique, "Sequencing table has duplicate study names")
    require(geo["geo_accession"].is_unique, "GEO table has duplicate accessions")

    restored = (
        labels.merge(
            mapping[["donor_id", "study_name", "geo_accession", "title"]],
            on="donor_id",
            how="left",
            validate="one_to_one",
        )
        .merge(
            clinical,
            left_on="study_name",
            right_on="Names",
            how="left",
            validate="one_to_one",
        )
        .merge(
            sequencing,
            left_on="study_name",
            right_on="SampleID",
            how="left",
            validate="one_to_one",
        )
        .merge(
            geo,
            on="geo_accession",
            how="left",
            validate="one_to_one",
            suffixes=("", "_geo"),
        )
    )

    required = [
        "study_name",
        "geo_accession",
        "Groups",
        "Batch",
        "Age",
        "Gender",
        "Race",
        "Ethnicity",
        "Collection_year",
        "Estimated Number of Cells",
        "Mean Reads per Cell",
        "Median Genes per Cell",
        "Sequencing Saturation",
        "Fraction Reads in Cells",
        "Median UMI Counts per Cell",
    ]
    missing = restored[required].isna().sum()
    require(not missing.any(), f"Incomplete donor mapping:\n{missing[missing.gt(0)]}")
    require(len(restored) == 44, f"Expected 44 analysis donors, found {len(restored)}")
    require(
        set(restored["Groups"]) == {"cSLE", "cHD"},
        f"Unexpected study groups: {sorted(restored['Groups'].unique())}",
    )

    expected_label = restored["Groups"].map({"cSLE": "case", "cHD": "control"})
    require(
        expected_label.equals(restored["case_control"]),
        "Study group and analysis case/control label disagree",
    )

    column_map = {
        "Batch": "batch",
        "Age": "age_years",
        "Gender": "sex",
        "Race": "race",
        "Ethnicity": "ethnicity",
        "Collection_year": "collection_year",
        "Estimated Number of Cells": "estimated_cells",
        "Mean Reads per Cell": "mean_reads_per_cell",
        "Median Genes per Cell": "median_genes_per_cell",
        "Sequencing Saturation": "sequencing_saturation",
        "Fraction Reads in Cells": "fraction_reads_in_cells",
        "Median UMI Counts per Cell": "median_umi_per_cell",
        "instrument_model": "instrument_model",
        "platform_id": "platform_id",
        "library_strategy": "library_strategy",
        "library_source": "library_source",
    }
    keep = [
        "donor_id",
        "study_name",
        "geo_accession",
        "case_control",
        *column_map.keys(),
    ]
    donor = restored[keep].rename(columns=column_map).copy()
    for column in [
        "age_years",
        "collection_year",
        "estimated_cells",
        "mean_reads_per_cell",
        "median_genes_per_cell",
        "sequencing_saturation",
        "fraction_reads_in_cells",
        "median_umi_per_cell",
    ]:
        donor[column] = clean_numeric(donor[column])
    require(
        not donor.isna().any().any(),
        f"Restored analysis table still has missing values:\n{donor.isna().sum()}",
    )

    donor["y_true"] = donor["case_control"].eq("case").astype(int)
    donor["log_estimated_cells"] = np.log1p(donor["estimated_cells"])
    donor["log_mean_reads_per_cell"] = np.log1p(donor["mean_reads_per_cell"])
    donor["log_median_genes_per_cell"] = np.log1p(donor["median_genes_per_cell"])
    donor["log_median_umi_per_cell"] = np.log1p(donor["median_umi_per_cell"])
    donor = donor.sort_values("donor_id").reset_index(drop=True)

    batch_label = pd.crosstab(
        donor["batch"], donor["case_control"], margins=True, margins_name="Total"
    )
    year_label = pd.crosstab(
        donor["collection_year"].astype(int),
        donor["case_control"],
        margins=True,
        margins_name="Total",
    )
    batch_year = pd.crosstab(donor["batch"], donor["collection_year"].astype(int))

    donor.to_csv(OUT / "gse135779_donor_metadata_restored.tsv", sep="\t", index=False)
    batch_label.to_csv(OUT / "gse135779_batch_by_label.tsv", sep="\t")
    year_label.to_csv(OUT / "gse135779_year_by_label.tsv", sep="\t")
    batch_year.to_csv(OUT / "gse135779_batch_by_year.tsv", sep="\t")

    source_files = [LABELS, MAPPING, CLINICAL, SEQUENCING, GEO]
    manifest = {
        "analysis": "GSE135779 childhood-cohort donor metadata restoration",
        "n_donors": int(len(donor)),
        "n_cases": int(donor["y_true"].sum()),
        "n_controls": int((1 - donor["y_true"]).sum()),
        "mapping_complete": True,
        "source_publication": {
            "title": "Mapping systemic lupus erythematosus heterogeneity at the single-cell level",
            "doi": "10.1038/s41590-020-0743-0",
            "metadata_location": "Supplementary Table 1b (clinical/batch) and 1c (sequencing QC)",
        },
        "source_files": {
            str(path.relative_to(PROJECT)): {"sha256": sha256(path), "bytes": path.stat().st_size}
            for path in source_files
        },
        "batch_counts": {
            str(batch): {
                "case": int(group["y_true"].sum()),
                "control": int((1 - group["y_true"]).sum()),
                "total": int(len(group)),
            }
            for batch, group in donor.groupby("batch")
        },
        "claim_boundary": (
            "Batch and collection year are observed design variables. Associations with disease "
            "labels do not identify a purely technical causal effect."
        ),
    }
    (OUT / "gse135779_metadata_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    print(batch_label.to_string())
    print()
    print(year_label.to_string())
    print(f"\nWrote restored metadata and provenance to {OUT}")


if __name__ == "__main__":
    main()
