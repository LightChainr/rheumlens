#!/usr/bin/env python3
"""Complete available-covariate and cohort-structure sensitivity analyses."""

from __future__ import annotations

from contextlib import redirect_stdout
import io
from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from confounding_audit import encoded_covariates
from fold_contained_benchmark import (
    COHORTS, classifier, folds_from_json, run_oof, select_hvg,
)


ROOT = Path(__file__).resolve().parents[1]
FEATURES = ROOT / "data" / "donor_features"
EVIDENCE = ROOT / "results" / "evidence_package"
PUBLIC = EVIDENCE / "public_metadata" / "structured"
OUT = EVIDENCE / "confounding_complete"


def gse135_covariates() -> pd.DataFrame:
    h5ad = Path(COHORTS["GSE135779"]["h5ad"])
    a = ad.read_h5ad(h5ad, backed="r")
    obs = a.obs.copy()
    a.file.close()
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["label"] = obs["disease"].astype(str).str.upper().eq("SLE").astype(int)
    qc = obs.groupby("donor_id").agg(
        label=("label", "first"), cells_per_donor=("label", "size"),
        mean_umi_per_cell=("total_counts", "mean"),
        mean_genes_per_cell=("n_genes_by_counts", "mean"),
        mean_pct_mito=("pct_counts_mt", "mean"),
    ).reset_index()
    clinical = pd.read_csv(PUBLIC / "GSE135779_ST1b_donor_clinical.csv")
    mapping = pd.read_csv(PUBLIC / "GSE135779_study_name_to_donor_id.csv")
    clinical = clinical.merge(mapping[["study_name", "donor_id"]], left_on="Names", right_on="study_name")
    keep = ["donor_id", "Batch", "Age", "Gender", "Race", "Ethnicity", "Collection_year",
            "SLEDAI", "MMF", "OS", "MTX", "Plaquenil", "Neph_all"]
    donor = qc.merge(clinical[keep], on="donor_id", how="left", validate="one_to_one")
    return donor


def parse_age(value) -> float:
    match = re.search(r"(\d+)", str(value))
    return float(match.group(1)) if match else np.nan


def available_covariates(cohort: str):
    if cohort == "GSE135779":
        donor = gse135_covariates()
        numeric = ["cells_per_donor", "mean_umi_per_cell", "mean_genes_per_cell", "mean_pct_mito", "Age"]
        categorical = ["Batch", "Gender", "Race", "Ethnicity", "Collection_year"]
    elif cohort == "GSE285773":
        donor = pd.read_csv(EVIDENCE / "GSE285773_donor_covariates.csv")
        numeric = ["cells_per_donor", "mean_umi_per_cell", "mean_genes_per_cell", "mean_pct_mito"]
        categorical = []
    else:
        donor = pd.read_csv(EVIDENCE / "GSE174188_donor_covariates.csv")
        donor["age"] = donor["development_stage"].map(parse_age)
        numeric = ["cells_per_donor", "mean_umi_per_cell", "mean_genes_per_cell", "mean_pct_mito",
                   "age", "n_samples", "n_libraries", "n_suspensions"]
        categorical = ["sex", "ethnicity", "processing_cohort"]
    return donor, numeric, categorical


def ci(y: np.ndarray, p: np.ndarray, n_boot: int = 10000):
    rng = np.random.default_rng(20260619)
    case, control = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    values = np.empty(n_boot)
    for i in range(n_boot):
        idx = np.concatenate([rng.choice(case, len(case), replace=True),
                              rng.choice(control, len(control), replace=True)])
        values[i] = roc_auc_score(y[idx], p[idx])
    return np.quantile(values, .025), np.quantile(values, .975)


def adjusted_oof(cohort: str):
    z = np.load(FEATURES / f"{cohort}_donor_features.npz", allow_pickle=True)
    donors, y = z["donors"].astype(str), z["y"].astype(int)
    expression, scgpt = z["expression"].astype(float), z["scgpt"].astype(float)
    donor, numeric, categorical = available_covariates(cohort)
    meta = donor.set_index("donor_id").loc[donors].reset_index()
    folds = folds_from_json(Path(COHORTS[cohort]["folds"]), donors, y)
    pred = {name: np.full(len(y), np.nan) for name in (
        "covariates_only", "scgpt_covariate_residual", "pca_covariate_residual", "hvg_covariate_residual"
    )}
    for tr, te in folds:
        ctr, cte = encoded_covariates(meta.iloc[tr], meta.iloc[te], numeric, categorical)
        pred["covariates_only"][te] = classifier().fit(ctr, y[tr]).predict_proba(cte)[:, 1]

        def residual_classify(xtr, xte):
            reg = LinearRegression().fit(ctr, xtr)
            rtr, rte = xtr - reg.predict(ctr), xte - reg.predict(cte)
            model = make_pipeline(StandardScaler(), classifier()).fit(rtr, y[tr])
            return model.predict_proba(rte)[:, 1]

        pred["scgpt_covariate_residual"][te] = residual_classify(scgpt[tr], scgpt[te])
        hvg = select_hvg(expression[tr], 2000)
        pca = PCA(n_components=min(25, len(tr) - 1), svd_solver="full").fit(expression[tr][:, hvg])
        pred["pca_covariate_residual"][te] = residual_classify(
            pca.transform(expression[tr][:, hvg]), pca.transform(expression[te][:, hvg])
        )
        pred["hvg_covariate_residual"][te] = residual_classify(
            expression[tr][:, hvg], expression[te][:, hvg]
        )
    frame = pd.DataFrame({"cohort": cohort, "donor_id": donors, "label": y, **pred})
    frame.to_csv(OUT / f"{cohort}_available_covariate_adjusted_oof.csv", index=False)
    rows = []
    for method, p in pred.items():
        lo, hi = ci(y, p)
        rows.append({"cohort": cohort, "analysis": "archived_folds_available_covariates",
                     "method": method, "n": len(y), "auc": roc_auc_score(y, p),
                     "ci_low": lo, "ci_high": hi,
                     "numeric_covariates": ";".join(numeric),
                     "categorical_covariates": ";".join(categorical)})
    return pd.DataFrame(rows)


def subset_cv(name: str, donor_mask: pd.Series):
    cohort = "GSE174188"
    z = np.load(FEATURES / f"{cohort}_donor_features.npz", allow_pickle=True)
    donors, y = z["donors"].astype(str), z["y"].astype(int)
    meta = pd.read_csv(EVIDENCE / "GSE174188_donor_covariates.csv").set_index("donor_id").loc[donors]
    keep = donor_mask.reindex(meta.index).fillna(False).to_numpy(bool)
    expression, scgpt = z["expression"].astype(float)[keep], z["scgpt"].astype(float)[keep]
    ys, ds = y[keep], donors[keep]
    folds = list(StratifiedKFold(5, shuffle=True, random_state=42).split(ds, ys))
    with redirect_stdout(io.StringIO()):
        pred, _ = run_oof(expression, scgpt, ys, folds)
    pd.DataFrame({"analysis": name, "donor_id": ds, "label": ys, **pred}).to_csv(
        OUT / f"GSE174188_{name}_oof.csv", index=False
    )
    rows = []
    for method, p in pred.items():
        lo, hi = ci(ys, p)
        rows.append({"cohort": "GSE174188", "analysis": name, "method": method,
                     "n": len(ys), "n_sle": int(ys.sum()), "n_hc": int((ys == 0).sum()),
                     "auc": roc_auc_score(ys, p), "ci_low": lo, "ci_high": hi})
    return pd.DataFrame(rows)


def leave_processing_cohort_out():
    z = np.load(FEATURES / "GSE174188_donor_features.npz", allow_pickle=True)
    donors, y = z["donors"].astype(str), z["y"].astype(int)
    expression, scgpt = z["expression"].astype(float), z["scgpt"].astype(float)
    meta = pd.read_csv(EVIDENCE / "GSE174188_donor_covariates.csv").set_index("donor_id").loc[donors]
    groups = meta["processing_cohort"].to_numpy()
    rows = []
    for group in sorted(np.unique(groups)):
        te, tr = np.flatnonzero(groups == group), np.flatnonzero(groups != group)
        if np.unique(y[te]).size < 2 or np.unique(y[tr]).size < 2:
            rows.append({"processing_cohort": group, "method": "all", "n_test": len(te),
                         "n_sle": int(y[te].sum()), "n_hc": int((y[te] == 0).sum()),
                         "auc": np.nan, "status": "one_class_in_train_or_test"})
            continue
        scaler = StandardScaler().fit(scgpt[tr])
        p = classifier().fit(scaler.transform(scgpt[tr]), y[tr]).predict_proba(scaler.transform(scgpt[te]))[:, 1]
        rows.append({"processing_cohort": group, "method": "scgpt", "n_test": len(te),
                     "n_sle": int(y[te].sum()), "n_hc": int((y[te] == 0).sum()),
                     "auc": roc_auc_score(y[te], p), "status": "estimated"})
        hvg = select_hvg(expression[tr], 2000)
        pca = PCA(n_components=25, svd_solver="full").fit(expression[tr][:, hvg])
        for method, xtr, xte in (
            ("expression_pca", pca.transform(expression[tr][:, hvg]), pca.transform(expression[te][:, hvg])),
            ("donor_mean_hvg", expression[tr][:, hvg], expression[te][:, hvg]),
        ):
            s = StandardScaler().fit(xtr)
            p = classifier().fit(s.transform(xtr), y[tr]).predict_proba(s.transform(xte))[:, 1]
            rows.append({"processing_cohort": group, "method": method, "n_test": len(te),
                         "n_sle": int(y[te].sum()), "n_hc": int((y[te] == 0).sum()),
                         "auc": roc_auc_score(y[te], p), "status": "estimated"})
    return pd.DataFrame(rows)


def sample_state_map():
    a = ad.read_h5ad(COHORTS["GSE174188"]["h5ad"], backed="r")
    obs = a.obs[a.obs["cell_type"].astype(str).eq("CD4-positive, alpha-beta T cell")].copy()
    a.file.close()
    table = obs.groupby(["donor_id", "sample_uuid", "disease_state"], observed=True).agg(
        n_cd4_cells=("observation_joinid", "size"),
        n_libraries=("library_uuid", "nunique"),
        n_suspensions=("suspension_uuid", "nunique"),
    ).reset_index()
    table.to_csv(OUT / "GSE174188_donor_sample_state_map.csv", index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summaries = [adjusted_oof(c) for c in ("GSE135779", "GSE285773", "GSE174188")]
    d174 = pd.read_csv(EVIDENCE / "GSE174188_donor_covariates.csv").set_index("donor_id")
    summaries.append(subset_cv("managed_or_control_only", d174["disease_state"].isin(["managed", "na"])))
    summaries.append(subset_cv("single_sample_donors_only", d174["n_samples"].eq(1)))
    pd.concat(summaries, ignore_index=True).to_csv(OUT / "confounding_sensitivity_summary.csv", index=False)
    leave_processing_cohort_out().to_csv(OUT / "GSE174188_leave_processing_cohort_out.csv", index=False)
    sample_state_map()
    print(pd.concat(summaries, ignore_index=True).to_string(index=False))


if __name__ == "__main__":
    main()
