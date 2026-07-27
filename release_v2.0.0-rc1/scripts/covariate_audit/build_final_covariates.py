#!/usr/bin/env python3
"""Build final-donor technical and demographic covariates from raw CD4 counts."""

from __future__ import annotations

import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from common import H5AD, WORKSPACE, label_path, load_labels, sha256_file


DATASET = "SLE_GSE174188_CD4"
OUT = WORKSPACE / "results" / "01_inputs"
CELL_TYPE = "CD4-positive, alpha-beta T cell"
CHUNK_SIZE = 5_000


def age_in_years(value: object) -> float:
    match = re.search(r"(\d+(?:\.\d+)?)", str(value))
    return float(match.group(1)) if match else np.nan


def single_value(group: pd.Series) -> object:
    values = group.dropna().astype(str).unique()
    return values[0] if len(values) == 1 else ";".join(sorted(values))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    labels = load_labels(DATASET)
    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs
    selected = obs["cell_type"].astype(str).eq(CELL_TYPE).to_numpy()
    rows = np.flatnonzero(selected)
    cd4 = obs.iloc[rows].copy()
    cd4["donor_id"] = cd4["donor_id"].astype(str)
    if set(cd4["donor_id"]) != set(labels.index):
        raise RuntimeError("Final CD4 donor universe does not match retained labels")

    donors = labels.index.to_list()
    donor_index = {donor: index for index, donor in enumerate(donors)}
    row_donor = cd4["donor_id"].map(donor_index).to_numpy(dtype=int)
    n = len(donors)
    n_cells = np.zeros(n, dtype=np.int64)
    umi_sum = np.zeros(n, dtype=np.float64)
    genes_sum = np.zeros(n, dtype=np.float64)
    mito_sum = np.zeros(n, dtype=np.float64)
    pct_mito_sum = np.zeros(n, dtype=np.float64)

    feature_names = adata.var["feature_name"].astype(str)
    mito_index = np.flatnonzero(feature_names.str.upper().str.startswith("MT-").to_numpy())
    if len(mito_index) == 0:
        raise RuntimeError("No mitochondrial features found in raw matrix")

    raw = adata.raw.X
    for start in range(0, len(rows), CHUNK_SIZE):
        stop = min(start + CHUNK_SIZE, len(rows))
        matrix = raw[rows[start:stop], :]
        cell_umi = np.asarray(matrix.sum(axis=1)).ravel()
        cell_genes = matrix.getnnz(axis=1).astype(float)
        cell_mito = np.asarray(matrix[:, mito_index].sum(axis=1)).ravel()
        cell_pct_mito = np.divide(cell_mito, cell_umi, out=np.zeros_like(cell_mito), where=cell_umi > 0)
        index = row_donor[start:stop]
        np.add.at(n_cells, index, 1)
        np.add.at(umi_sum, index, cell_umi)
        np.add.at(genes_sum, index, cell_genes)
        np.add.at(mito_sum, index, cell_mito)
        np.add.at(pct_mito_sum, index, cell_pct_mito)

    grouped = cd4.groupby("donor_id", observed=True)
    covariates = pd.DataFrame(index=pd.Index(donors, name="donor_id"))
    covariates["case_control"] = labels
    covariates["y_true"] = labels.eq("case").astype(int)
    covariates["cells_per_donor"] = n_cells
    covariates["mean_umi_per_cell"] = umi_sum / n_cells
    covariates["mean_genes_per_cell"] = genes_sum / n_cells
    covariates["aggregate_pct_mito"] = mito_sum / umi_sum
    covariates["mean_cell_pct_mito"] = pct_mito_sum / n_cells
    covariates["n_libraries"] = grouped["library_uuid"].nunique().reindex(donors).to_numpy()
    covariates["n_samples"] = grouped["sample_uuid"].nunique().reindex(donors).to_numpy()
    covariates["n_suspensions"] = grouped["suspension_uuid"].nunique().reindex(donors).to_numpy()
    covariates["n_processing_cohorts"] = grouped["Processing_Cohort"].nunique().reindex(donors).to_numpy()
    covariates["age_years"] = grouped["development_stage"].apply(lambda x: np.nanmedian([age_in_years(v) for v in x])).reindex(donors)
    covariates["sex"] = grouped["sex"].apply(single_value).reindex(donors)
    covariates["ethnicity"] = grouped["self_reported_ethnicity"].apply(single_value).reindex(donors)

    processing = pd.crosstab(cd4["donor_id"], cd4["Processing_Cohort"], normalize="index")
    processing = processing.reindex(donors).fillna(0.0)
    processing.columns = [f"processing_fraction_{str(value).replace('.', '_')}" for value in processing.columns]
    covariates = covariates.join(processing)
    covariates.to_csv(OUT / "gse174188_final_donor_covariates.tsv", sep="\t")

    dictionary_rows = [
        ("cells_per_donor", "technical", "Number of retained CD4 cells"),
        ("mean_umi_per_cell", "technical", "Mean raw UMI count per retained CD4 cell"),
        ("mean_genes_per_cell", "technical", "Mean detected genes per retained CD4 cell"),
        ("aggregate_pct_mito", "technical", "Donor aggregate mitochondrial UMI fraction"),
        ("mean_cell_pct_mito", "technical", "Mean per-cell mitochondrial UMI fraction"),
        ("n_libraries", "technical", "Distinct library UUIDs represented for donor"),
        ("n_samples", "technical", "Distinct sample UUIDs represented for donor"),
        ("n_suspensions", "technical", "Distinct suspension UUIDs represented for donor"),
        ("n_processing_cohorts", "technical", "Distinct Processing_Cohort levels represented for donor"),
        ("age_years", "demographic", "Age parsed from development_stage"),
        ("sex", "demographic", "Registered donor sex"),
        ("ethnicity", "demographic_sensitivity", "Self-reported ethnicity; sparse levels excluded from primary adjustment"),
    ]
    dictionary_rows.extend((column, "technical", "Fraction of retained CD4 cells in processing cohort") for column in processing.columns)
    pd.DataFrame(dictionary_rows, columns=["field", "block", "definition"]).to_csv(
        OUT / "covariate_dictionary.tsv", sep="\t", index=False
    )
    summary = covariates.describe(include="all").transpose()
    summary.to_csv(OUT / "covariate_summary.tsv", sep="\t")
    manifest = {
        "dataset": DATASET,
        "cell_type": CELL_TYPE,
        "n_cells": int(len(rows)),
        "n_donors": int(n),
        "n_mito_features": int(len(mito_index)),
        "h5ad": {"path": str(H5AD), "size_bytes": H5AD.stat().st_size, "mtime_ns": H5AD.stat().st_mtime_ns},
        "labels": {"path": str(label_path(DATASET)), "sha256": sha256_file(label_path(DATASET))},
        "raw_layer": "adata.raw.X",
    }
    (OUT / "input_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    adata.file.close()


if __name__ == "__main__":
    main()
