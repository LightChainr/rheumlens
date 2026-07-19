#!/usr/bin/env python3
"""Build a submission workbook and an upgraded cohort-characteristics table."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill


PROJECT = Path(__file__).resolve().parents[2]
REPO = PROJECT / "10_zenodo_release_20260717" / "rheumlens"
TABLES = REPO / "supplementary_tables"
PUBLIC = REPO / "evidence_package" / "public_metadata"
EVIDENCE = REPO / "evidence_package"
OUTPUT = PROJECT / "11_advanced_visualizations_20260717" / "source_data" / "supplementary_tables"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def numeric_smd(case: pd.Series, control: pd.Series) -> float:
    case = pd.to_numeric(case, errors="coerce").dropna()
    control = pd.to_numeric(control, errors="coerce").dropna()
    if len(case) < 2 or len(control) < 2:
        return np.nan
    pooled = np.sqrt((case.var(ddof=1) + control.var(ddof=1)) / 2)
    return float((case.mean() - control.mean()) / pooled) if pooled > 0 else np.nan


def binary_smd(case: pd.Series, control: pd.Series) -> float:
    case = pd.to_numeric(case, errors="coerce").dropna()
    control = pd.to_numeric(control, errors="coerce").dropna()
    if case.empty or control.empty:
        return np.nan
    p1, p0 = case.mean(), control.mean()
    pooled = np.sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    return float((p1 - p0) / pooled) if pooled > 0 else np.nan


def continuous_display(values: pd.Series) -> str:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if values.empty:
        return "NA"
    q1, median, q3 = values.quantile([0.25, 0.5, 0.75])
    return f"{median:.1f} [{q1:.1f}, {q3:.1f}]"


def binary_display(values: pd.Series) -> str:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if values.empty:
        return "NA"
    return f"{int(values.sum())}/{len(values)} ({100 * values.mean():.1f}%)"


def append_continuous(rows: list[dict], dataset: str, frame: pd.DataFrame, variable: str, column: str) -> None:
    case = frame.loc[frame["label"].eq(1), column]
    control = frame.loc[frame["label"].eq(0), column]
    rows.append(
        {
            "Dataset": dataset,
            "Variable": variable,
            "SLE": continuous_display(case),
            "Control": continuous_display(control),
            "SMD (SLE-control)": numeric_smd(case, control),
            "Missing SLE": int(pd.to_numeric(case, errors="coerce").isna().sum()),
            "Missing control": int(pd.to_numeric(control, errors="coerce").isna().sum()),
        }
    )


def append_binary(rows: list[dict], dataset: str, frame: pd.DataFrame, variable: str, column: str) -> None:
    case = frame.loc[frame["label"].eq(1), column]
    control = frame.loc[frame["label"].eq(0), column]
    rows.append(
        {
            "Dataset": dataset,
            "Variable": variable,
            "SLE": binary_display(case),
            "Control": binary_display(control),
            "SMD (SLE-control)": binary_smd(case, control),
            "Missing SLE": int(pd.to_numeric(case, errors="coerce").isna().sum()),
            "Missing control": int(pd.to_numeric(control, errors="coerce").isna().sum()),
        }
    )


def gse135779() -> pd.DataFrame:
    clinical = pd.read_csv(PUBLIC / "GSE135779_ST1b_donor_clinical.csv", na_values=["ND", "not reported"])
    clinical = clinical.loc[clinical["Groups"].astype(str).str.startswith("c")].copy()
    clinical["label"] = clinical["Groups"].eq("cSLE").astype(int)
    clinical["female"] = clinical["Gender"].map({"F": 1, "M": 0})
    mapping = pd.read_csv(PUBLIC / "GSE135779_study_name_to_donor_id.csv")
    sequencing = pd.read_csv(PUBLIC / "GSE135779_ST1c_sequencing.csv")
    cell_counts = pd.read_csv(
        PROJECT
        / "03_远程回传/final_gpu_figure_ready_20260709/extracted_plus/RheumLens_GPU_figure_ready_plus_20260709/extra/pseudobulk/SLE_GSE135779/donor_cell_counts.tsv",
        sep="\t",
    )
    frame = clinical.merge(mapping[["study_name", "donor_id"]], left_on="Names", right_on="study_name", how="left")
    frame = frame.merge(sequencing, left_on="Names", right_on="SampleID", how="left")
    frame = frame.merge(cell_counts, on="donor_id", how="left")
    return frame


def gse174188() -> pd.DataFrame:
    frame = pd.read_csv(EVIDENCE / "GSE174188_donor_covariates.csv", dtype={"donor_id": str})
    frame["age"] = pd.to_numeric(frame["development_stage"].str.extract(r"(\d+(?:\.\d+)?)")[0], errors="coerce")
    frame["female"] = frame["sex"].map({"female": 1, "male": 0})
    return frame


def gse285773() -> pd.DataFrame:
    return pd.read_csv(EVIDENCE / "GSE285773_donor_covariates.csv", dtype={"donor_id": str})


def cohort_characteristics() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    frames = {"GSE135779": gse135779(), "GSE174188 CD4": gse174188(), "GSE285773 CD4": gse285773()}
    rows: list[dict] = []
    for dataset, frame in frames.items():
        rows.append(
            {
                "Dataset": dataset,
                "Variable": "Donors, n",
                "SLE": str(int(frame["label"].eq(1).sum())),
                "Control": str(int(frame["label"].eq(0).sum())),
                "SMD (SLE-control)": np.nan,
                "Missing SLE": 0,
                "Missing control": 0,
            }
        )
    append_continuous(rows, "GSE135779", frames["GSE135779"], "Age, years", "Age")
    append_binary(rows, "GSE135779", frames["GSE135779"], "Female sex", "female")
    append_continuous(rows, "GSE135779", frames["GSE135779"], "Retained cells per donor", "n_cells_raw")
    append_continuous(rows, "GSE135779", frames["GSE135779"], "Mean reads per cell", "Mean Reads per Cell")
    append_continuous(rows, "GSE135779", frames["GSE135779"], "Median detected genes per cell", "Median Genes per Cell")
    append_continuous(rows, "GSE135779", frames["GSE135779"], "Median UMI counts per cell", "Median UMI Counts per Cell")
    for level in sorted(frames["GSE135779"]["Race"].dropna().astype(str).unique()):
        column = f"race_{level}"
        frames["GSE135779"][column] = frames["GSE135779"]["Race"].eq(level).astype(int)
        append_binary(rows, "GSE135779", frames["GSE135779"], f"Source race category: {level}", column)

    append_continuous(rows, "GSE174188 CD4", frames["GSE174188 CD4"], "Age, years", "age")
    append_binary(rows, "GSE174188 CD4", frames["GSE174188 CD4"], "Female sex", "female")
    for variable, column in [
        ("Retained CD4 cells per donor", "cells_per_donor"),
        ("Mean UMI per cell", "mean_umi_per_cell"),
        ("Mean detected genes per cell", "mean_genes_per_cell"),
        ("Mean mitochondrial fraction, %", "mean_pct_mito"),
        ("Libraries per donor", "n_libraries"),
        ("Suspensions per donor", "n_suspensions"),
    ]:
        append_continuous(rows, "GSE174188 CD4", frames["GSE174188 CD4"], variable, column)
    for level in sorted(frames["GSE174188 CD4"]["ethnicity"].dropna().astype(str).unique()):
        column = f"ethnicity_{level}"
        frames["GSE174188 CD4"][column] = frames["GSE174188 CD4"]["ethnicity"].eq(level).astype(int)
        append_binary(rows, "GSE174188 CD4", frames["GSE174188 CD4"], f"Source ethnicity: {level}", column)

    for variable, column in [
        ("Retained CD4 cells per donor", "cells_per_donor"),
        ("Mean UMI per cell", "mean_umi_per_cell"),
        ("Mean detected genes per cell", "mean_genes_per_cell"),
        ("Mean mitochondrial fraction, %", "mean_pct_mito"),
    ]:
        append_continuous(rows, "GSE285773 CD4", frames["GSE285773 CD4"], variable, column)

    result = pd.DataFrame(rows)
    result["SMD (SLE-control)"] = result["SMD (SLE-control)"].round(3)
    return result, frames


def source_tables() -> list[tuple[str, str, pd.DataFrame, Path | None]]:
    records: list[tuple[str, str, pd.DataFrame, Path | None]] = []
    for number in range(1, 12):
        path = next(TABLES.glob(f"SUPP_TABLE_S{number}_*.csv"))
        records.append((f"S{number}", path.stem, pd.read_csv(path), path))

    path = TABLES / "SUPP_MECH_GSE174188_CD4_donor_level_gene_discrimination.csv"
    records.append(("S12", "Donor-level gene discrimination", pd.read_csv(path), path))

    pieces = []
    s13_paths = [
        TABLES / "SUPP_MECH_TOP100_SLE_HIGH_GENES.csv",
        TABLES / "SUPP_MECH_TOP100_CONTROL_HIGH_GENES.csv",
        TABLES / "SUPP_MECH_TOP100_DISCRIMINATIVE_GENES.csv",
    ]
    for path in s13_paths:
        frame = pd.read_csv(path)
        frame.insert(0, "list_type", path.stem.replace("SUPP_MECH_TOP100_", ""))
        pieces.append(frame)
    records.append(("S13", "Top gene lists", pd.concat(pieces, ignore_index=True), None))

    path = TABLES / "SUPP_MECH_GSE174188_CD4_donor_ISG_scores.csv"
    records.append(("S14", "Donor-level ISG scores", pd.read_csv(path), path))
    for number in range(15, 19):
        path = next(TABLES.glob(f"SUPP_TABLE_S{number}_*.csv"))
        records.append((f"S{number}", path.stem, pd.read_csv(path), path))
    return records


def style_workbook(writer: pd.ExcelWriter) -> None:
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions
        worksheet.row_dimensions[1].height = 28
        for cell in worksheet[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="1F4E78")
            cell.alignment = Alignment(wrap_text=True, vertical="center")
        for column in worksheet.columns:
            letter = column[0].column_letter
            width = min(max(len(str(cell.value or "")) for cell in column) + 2, 48)
            worksheet.column_dimensions[letter].width = max(width, 10)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    characteristics, donor_frames = cohort_characteristics()
    characteristics.to_csv(OUTPUT / "SUPP_TABLE_S19_COHORT_CHARACTERISTICS.csv", index=False)

    scaled_path = (
        PROJECT
        / "11_advanced_visualizations_20260717/source_data/scaled_internal_sensitivity/tuned_c/summary.csv"
    )
    scaled = pd.read_csv(scaled_path)
    scaled.to_csv(OUTPUT / "SUPP_TABLE_S20_SCALED_TUNED_INTERNAL_CV.csv", index=False)

    workbook = OUTPUT / "Applied_Sciences_Supplementary_Tables_S1-S20.xlsx"
    index_rows = []
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        readme = pd.DataFrame(
            {
                "Item": ["Scope", "S19 display", "S19 SMD", "Missing metadata", "S20 purpose"],
                "Description": [
                    "Tables S1-S18 reproduce the retained public repository tables; S19 adds case-control cohort characteristics; S20 reports the fold-scaled, inner-C-tuned internal sensitivity.",
                    "Continuous variables are median [Q1, Q3]; categorical variables are n/N (%).",
                    "Continuous SMD uses the pooled standard deviation; binary SMD uses the pooled Bernoulli variance. Positive values indicate higher values or proportions in SLE.",
                    "Age, sex, treatment, disease activity, ancestry, and batch were not available in the archived GSE285773 donor object and are not imputed.",
                    "S20 addresses feature-scale sensitivity of the archived fixed-C internal benchmark; the strict external transfer already uses source-only scaling and remains the primary model-selection evidence.",
                ],
            }
        )
        readme.to_excel(writer, sheet_name="README", index=False)
        for number, title, frame, path in source_tables():
            sheet = f"Table_{number}"
            frame.to_excel(writer, sheet_name=sheet, index=False)
            index_rows.append(
                {
                    "Table": number,
                    "Worksheet": sheet,
                    "Description": title,
                    "Rows": len(frame),
                    "Source SHA256": sha256(path) if path else "assembled from named source tables",
                }
            )
        characteristics.to_excel(writer, sheet_name="Table_S19", index=False)
        scaled.to_excel(writer, sheet_name="Table_S20", index=False)
        index_rows.extend(
            [
                {"Table": "S19", "Worksheet": "Table_S19", "Description": "Case-control cohort characteristics and available QC metrics", "Rows": len(characteristics), "Source SHA256": "derived by accompanying script"},
                {"Table": "S20", "Worksheet": "Table_S20", "Description": "Fold-scaled, inner-C-tuned repeated internal sensitivity", "Rows": len(scaled), "Source SHA256": sha256(scaled_path)},
            ]
        )
        pd.DataFrame(index_rows).to_excel(writer, sheet_name="INDEX", index=False)
        donor_frames["GSE135779"].to_excel(writer, sheet_name="S19_GSE135779_source", index=False)
        donor_frames["GSE174188 CD4"].to_excel(writer, sheet_name="S19_GSE174188_source", index=False)
        donor_frames["GSE285773 CD4"].to_excel(writer, sheet_name="S19_GSE285773_source", index=False)
        style_workbook(writer)
    print(workbook)


if __name__ == "__main__":
    main()
