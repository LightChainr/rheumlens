#!/usr/bin/env python3
"""Donor-level confounding audit for the two high-AUC replication cohorts."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from fold_contained_benchmark import COHORTS, classifier, folds_from_json


ROOT = Path(__file__).resolve().parents[1]
FEATURE_DIR = ROOT / "data" / "donor_features"
OUT = ROOT / "results" / "evidence_package"
GSE174 = Path(COHORTS["GSE174188"]["h5ad"])
GSE285 = Path(COHORTS["GSE285773"]["h5ad"])


def safe_mode(x: pd.Series):
    values = x.dropna().astype(str)
    return values.mode().iloc[0] if len(values) else np.nan


def build_gse285_covariates() -> pd.DataFrame:
    a = ad.read_h5ad(GSE285, backed="r")
    obs = a.obs.copy()
    a.file.close()
    obs["label"] = (obs["disease"].astype(str).str.upper() == "SLE").astype(int)
    numeric = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    donor = obs.groupby(obs["donor_id"].astype(str)).agg(
        label=("label", "first"),
        cells_per_donor=("label", "size"),
        mean_umi_per_cell=(numeric[0], "mean"),
        mean_genes_per_cell=(numeric[1], "mean"),
        mean_pct_mito=(numeric[2], "mean"),
    ).reset_index(names="donor_id")
    donor["metadata_limitations"] = (
        "archived object has no batch, site, ancestry, age, sex, treatment, chemistry, lane, or processing-date fields"
    )
    return donor


def build_gse174_covariates() -> pd.DataFrame:
    cache = OUT / "GSE174188_donor_covariates.csv"
    if cache.exists():
        return pd.read_csv(cache)
    a = ad.read_h5ad(GSE174, backed="r")
    cd4 = a.obs["cell_type"].astype(str).eq("CD4-positive, alpha-beta T cell").to_numpy()
    idx_all = np.flatnonzero(cd4)
    obs = a.obs.iloc[idx_all].copy()
    donor_cell = obs["donor_id"].astype(str).to_numpy()
    donors = np.asarray(sorted(pd.unique(donor_cell)), dtype=object)
    dmap = {d: i for i, d in enumerate(donors)}
    codes = np.asarray([dmap[d] for d in donor_cell], dtype=np.int32)

    raw = a.raw.X
    raw_var = a.raw.var
    symbols = (raw_var["feature_name"].astype(str).to_numpy()
               if "feature_name" in raw_var else raw_var.index.astype(str).to_numpy())
    mt_idx = np.flatnonzero(np.char.startswith(np.char.upper(symbols.astype(str)), "MT-"))
    n = len(donors)
    cells = np.zeros(n, dtype=np.int64)
    umi_sum = np.zeros(n, dtype=np.float64)
    gene_sum = np.zeros(n, dtype=np.float64)
    mt_fraction_sum = np.zeros(n, dtype=np.float64)
    chunk = 4096
    for start in range(0, len(idx_all), chunk):
        stop = min(start + chunk, len(idx_all))
        x = raw[idx_all[start:stop]].tocsr()
        genes = np.asarray(x.getnnz(axis=1)).ravel()
        keep = genes >= 200
        if not keep.any():
            continue
        x = x[keep]
        c = codes[start:stop][keep]
        umi = np.asarray(x.sum(axis=1)).ravel()
        mt = np.asarray(x[:, mt_idx].sum(axis=1)).ravel() / umi if len(mt_idx) else np.zeros(len(umi))
        np.add.at(cells, c, 1)
        np.add.at(umi_sum, c, umi)
        np.add.at(gene_sum, c, genes[keep])
        np.add.at(mt_fraction_sum, c, mt)
        if start == 0 or stop == len(idx_all) or (start // chunk) % 25 == 0:
            print(f"GSE174188 raw-QC aggregation {stop:,}/{len(idx_all):,}", flush=True)

    frame = pd.DataFrame({
        "donor_id": donors,
        "cells_per_donor": cells,
        "mean_umi_per_cell": umi_sum / cells,
        "mean_genes_per_cell": gene_sum / cells,
        "mean_pct_mito": 100 * mt_fraction_sum / cells,
    })
    obs["donor_id"] = obs["donor_id"].astype(str)
    meta_rows = []
    for donor, group in obs.groupby("donor_id"):
        meta_rows.append({
            "donor_id": donor,
            "label": int(group["disease"].astype(str).eq("systemic lupus erythematosus").iloc[0]),
            "sex": safe_mode(group["sex"]),
            "ethnicity": safe_mode(group["self_reported_ethnicity"]),
            "processing_cohort": safe_mode(group["Processing_Cohort"]),
            "development_stage": safe_mode(group["development_stage"]),
            "disease_state": safe_mode(group["disease_state"]),
            "n_samples": group["sample_uuid"].nunique(),
            "n_libraries": group["library_uuid"].nunique(),
            "n_suspensions": group["suspension_uuid"].nunique(),
        })
    a.file.close()
    donor = frame.merge(pd.DataFrame(meta_rows), on="donor_id", validate="one_to_one")
    donor.to_csv(cache, index=False)
    return donor


def cross_tables(cohort: str, donor: pd.DataFrame, categorical: list[str]) -> None:
    rows = []
    for cov in categorical:
        table = pd.crosstab(donor[cov].fillna("NA"), donor["label"])
        for level, values in table.iterrows():
            rows.append({"cohort": cohort, "covariate": cov, "level": level,
                         "HC": int(values.get(0, 0)), "SLE": int(values.get(1, 0))})
    pd.DataFrame(rows).to_csv(OUT / f"{cohort}_disease_covariate_crosstabs.csv", index=False)


def encoded_covariates(train: pd.DataFrame, test: pd.DataFrame,
                       numeric: list[str], categorical: list[str]):
    transformer = ColumnTransformer([
        ("num", make_pipeline(SimpleImputer(strategy="median"), StandardScaler()), numeric),
        ("cat", make_pipeline(SimpleImputer(strategy="most_frequent"),
                              OneHotEncoder(handle_unknown="ignore", sparse_output=False)), categorical),
    ], sparse_threshold=0)
    return transformer.fit_transform(train), transformer.transform(test)


def model_audit(cohort: str, donor: pd.DataFrame, numeric: list[str],
                categorical: list[str]) -> pd.DataFrame:
    z = np.load(FEATURE_DIR / f"{cohort}_donor_features.npz", allow_pickle=True)
    donors = z["donors"].astype(str)
    y = z["y"].astype(int)
    scgpt = z["scgpt"].astype(np.float64)
    meta = donor.set_index("donor_id").loc[donors].reset_index()
    if not np.array_equal(meta["label"].to_numpy(), y):
        raise ValueError(f"{cohort}: metadata labels do not match frozen donor features")
    folds = folds_from_json(Path(COHORTS[cohort]["folds"]), donors, y)
    p_tech = np.full(len(y), np.nan)
    p_resid = np.full(len(y), np.nan)
    for tr, te in folds:
        ctr, cte = encoded_covariates(meta.iloc[tr], meta.iloc[te], numeric, categorical)
        p_tech[te] = classifier().fit(ctr, y[tr]).predict_proba(cte)[:, 1]
        residualizer = LinearRegression().fit(ctr, scgpt[tr])
        xtr = scgpt[tr] - residualizer.predict(ctr)
        xte = scgpt[te] - residualizer.predict(cte)
        scaler = StandardScaler().fit(xtr)
        p_resid[te] = classifier().fit(scaler.transform(xtr), y[tr]).predict_proba(scaler.transform(xte))[:, 1]
    pred = pd.DataFrame({"cohort": cohort, "donor_id": donors, "label": y,
                         "technical_covariates_only": p_tech,
                         "scgpt_covariate_residual": p_resid})
    pred.to_csv(OUT / f"{cohort}_confounding_oof_predictions.csv", index=False)
    frozen = pd.read_csv(ROOT / "results" / "fold_contained" / f"{cohort}_oof_predictions.csv")
    return pd.DataFrame([
        {"cohort": cohort, "model": "technical_covariates_only", "auc": roc_auc_score(y, p_tech)},
        {"cohort": cohort, "model": "scgpt_unadjusted_frozen", "auc": roc_auc_score(y, frozen["scgpt"])},
        {"cohort": cohort, "model": "scgpt_residualized_for_available_covariates", "auc": roc_auc_score(y, p_resid)},
    ])


def pc_associations(cohort: str, donor: pd.DataFrame, numeric: list[str],
                    categorical: list[str]) -> None:
    z = np.load(FEATURE_DIR / f"{cohort}_donor_features.npz", allow_pickle=True)
    donors = z["donors"].astype(str)
    expression = z["expression"].astype(np.float64)
    meta = donor.set_index("donor_id").loc[donors].reset_index()
    hvg = np.argsort(-np.var(expression, axis=0, ddof=1), kind="mergesort")[:2000]
    pcs = PCA(n_components=5, svd_solver="full").fit_transform(StandardScaler().fit_transform(expression[:, hvg]))
    rows = []
    variables = ["label", *numeric, *categorical]
    for pc_i in range(5):
        score = pcs[:, pc_i]
        for var in variables:
            if var in categorical:
                x = pd.get_dummies(meta[var].fillna("NA").astype(str), drop_first=False).to_numpy(float)
                r2 = LinearRegression().fit(x, score).score(x, score)
                p = stats.f_oneway(*(score[meta[var].fillna("NA").astype(str).to_numpy() == level]
                                     for level in meta[var].fillna("NA").astype(str).unique())).pvalue
                metric = "eta_squared_r2"
                value = r2
            else:
                vals = pd.to_numeric(meta[var], errors="coerce").to_numpy()
                keep = np.isfinite(vals)
                value, p = stats.pearsonr(vals[keep], score[keep])
                metric = "pearson_r"
            rows.append({"cohort": cohort, "pc": pc_i + 1, "variable": var,
                         "metric": metric, "value": value, "p_value_unadjusted": p})
    pd.DataFrame(rows).to_csv(OUT / f"{cohort}_PC1_PC5_covariate_associations.csv", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d285 = build_gse285_covariates()
    d285.to_csv(OUT / "GSE285773_donor_covariates.csv", index=False)
    n285 = ["cells_per_donor", "mean_umi_per_cell", "mean_genes_per_cell", "mean_pct_mito"]
    cross_tables("GSE285773", d285, [])
    m285 = model_audit("GSE285773", d285, n285, [])
    pc_associations("GSE285773", d285, n285, [])

    d174 = build_gse174_covariates()
    d174.to_csv(OUT / "GSE174188_donor_covariates.csv", index=False)
    n174 = ["cells_per_donor", "mean_umi_per_cell", "mean_genes_per_cell", "mean_pct_mito",
            "n_samples", "n_libraries", "n_suspensions"]
    c174 = ["sex", "ethnicity", "processing_cohort", "development_stage"]
    cross_tables("GSE174188", d174, [*c174, "disease_state"])
    m174 = model_audit("GSE174188", d174, n174, c174)
    pc_associations("GSE174188", d174, n174, c174)
    pd.concat([m285, m174], ignore_index=True).to_csv(OUT / "confounding_model_summary.csv", index=False)

    limitations = {
        "GSE285773": "Only cell-count and per-cell QC summaries are present in the archived object. Batch/site/ancestry/age/sex/treatment/chemistry/lane/date cannot be tested until source metadata are obtained.",
        "GSE174188": "Available fields include sex, self-reported ethnicity, processing cohort, disease state, sample/library/suspension IDs and raw-count QC. Exact age, sequencing lane, processing date, and site are absent from the CELLxGENE object.",
        "interpretation": "These are within-cohort confounding diagnostics, not proof that discrimination is causal disease biology.",
    }
    (OUT / "confounding_metadata_limitations.json").write_text(json.dumps(limitations, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
