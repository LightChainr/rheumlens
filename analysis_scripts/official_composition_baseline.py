#!/usr/bin/env python3
"""GSE135779 composition baseline from author-provided cluster counts."""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from fold_contained_benchmark import COHORTS, classifier, folds_from_json


ROOT = Path(__file__).resolve().parents[1]
PUBLIC = ROOT / "results" / "evidence_package" / "public_metadata" / "structured"
OUT = ROOT / "results" / "evidence_package" / "secondary_foldcontained"


counts = pd.read_csv(PUBLIC / "GSE135779_ST4a_official_cluster_counts.csv")
counts = counts.rename(columns={counts.columns[0]: "study_name"})
mapping = pd.read_csv(PUBLIC / "GSE135779_study_name_to_donor_id.csv")
mapping = mapping.dropna(subset=["study_name", "donor_id"])
counts = counts.merge(mapping[["study_name", "donor_id"]], on="study_name", validate="one_to_one")
features = [c for c in counts.columns if c.startswith("Cluster ")]
counts[features] = counts[features].apply(pd.to_numeric, errors="coerce").fillna(0)
counts[features] = counts[features].div(counts[features].sum(axis=1), axis=0)

z = np.load(ROOT / "data" / "donor_features" / "GSE135779_donor_features.npz", allow_pickle=True)
donors = z["donors"].astype(str)
y = z["y"].astype(int)
x = counts.set_index("donor_id").loc[donors, features].to_numpy(float)
folds = folds_from_json(Path(COHORTS["GSE135779"]["folds"]), donors, y)

prob = np.full(len(y), np.nan)
for tr, te in folds:
    model = make_pipeline(StandardScaler(), classifier()).fit(x[tr], y[tr])
    prob[te] = model.predict_proba(x[te])[:, 1]
observed = roc_auc_score(y, prob)

rng = np.random.default_rng(20260619)
case, control = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
frozen = pd.read_csv(ROOT / "results" / "fold_contained" / "GSE135779_oof_predictions.csv")
scgpt = frozen["scgpt"].to_numpy()
n_boot = 10000
boot_comp, boot_delta = np.empty(n_boot), np.empty(n_boot)
for i in range(n_boot):
    idx = np.concatenate([rng.choice(case, len(case), replace=True),
                          rng.choice(control, len(control), replace=True)])
    boot_comp[i] = roc_auc_score(y[idx], prob[idx])
    boot_delta[i] = roc_auc_score(y[idx], scgpt[idx]) - boot_comp[i]

n_perm = 1000
null = np.empty(n_perm)
for i in range(n_perm):
    yp = rng.permutation(y)
    p = np.full(len(y), np.nan)
    if any(np.unique(yp[tr]).size < 2 for tr, _ in folds):
        null[i] = 0.5
        continue
    for tr, te in folds:
        model = make_pipeline(StandardScaler(), classifier()).fit(x[tr], yp[tr])
        p[te] = model.predict_proba(x[te])[:, 1]
    null[i] = roc_auc_score(yp, p)
exceed = int((null >= observed).sum())

OUT.mkdir(parents=True, exist_ok=True)
pd.DataFrame({"donor_id": donors, "label": y, "official_cluster_composition_probability": prob}).to_csv(
    OUT / "GSE135779_official_composition_oof_predictions.csv", index=False
)
pd.DataFrame([{
    "cohort": "GSE135779", "model": "author_cluster_proportions_20",
    "observed_auc": observed, "ci_low": np.quantile(boot_comp, .025),
    "ci_high": np.quantile(boot_comp, .975), "n_boot": n_boot,
    "permutation_p": (exceed + 1) / (n_perm + 1), "n_perm": n_perm,
    "scgpt_minus_composition": roc_auc_score(y, scgpt) - observed,
    "delta_ci_low": np.quantile(boot_delta, .025),
    "delta_ci_high": np.quantile(boot_delta, .975),
}]).to_csv(OUT / "GSE135779_official_composition_summary.csv", index=False)
np.savez_compressed(OUT / "GSE135779_official_composition_resampling.npz",
                    composition_auc=boot_comp, scgpt_minus_composition=boot_delta,
                    permutation_null=null)
print(pd.read_csv(OUT / "GSE135779_official_composition_summary.csv").to_string(index=False))
