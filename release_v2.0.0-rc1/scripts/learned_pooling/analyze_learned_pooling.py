"""Inference layer for the learned-pooling experiment.

Consumes `results/pred_*.npy` written by `run_learned_pooling.py` and produces, for
both transfer directions:

  * target ROC-AUC with 5,000-resample case/control-stratified donor bootstrap CIs;
  * a PCA32 donor-mean + logistic control fitted with the same locked selector, which
    separates "source-fitted dimension reduction" from "learned pooling";
  * paired DeLong tests against the two Geneformer baselines (1152-d mean pooling and
    the PCA32 donor-mean control) -- PRIMARY family, BH-adjusted across all of them;
  * paired DeLong tests against the archived HVG / PCA pseudobulk transfer predictions
    -- SECONDARY family, BH-adjusted separately.

Usage:  python3 scripts/analyze_learned_pooling.py inputs results
"""
import csv
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import protocol

IN, OUT = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from path_config import RESULTS_ROOT

ARCH = os.environ.get(
    "RHEUMLENS_TRANSFER_PREDICTIONS",
    str(
        RESULTS_ROOT
        / "strict_source_only_transfer"
        / "strict_source_only_transfer_predictions.tsv"
    ),
)
MODELS = ("deepsets", "gatedmil", "pma")
DIRECTIONS = (("GSE174188", "GSE285773"), ("GSE285773", "GSE174188"))
BOOT = 5000


def donors(tag):
    f = IN + ("/GSE174188_donor_index_v2.csv" if tag == "GSE174188" else "/GSE285773_donor_index.csv")
    with open(f) as fh:
        return [r["donor_id"] for r in csv.DictReader(fh)]


arch = pd.read_csv(ARCH, sep="\t")
arch["direction"] = arch["direction"].str.replace("SLE_", "", regex=False).str.replace("_CD4", "", regex=False)

summary = json.load(open(OUT + "/learned_pooling_summary.json"))
cv_lookup = {(r["method"], r["direction"]): r for r in summary["results"]}

rows, tests = [], []
for src, tgt in DIRECTIONS:
    key = "%s_to_%s" % (src, tgt)
    y = np.load(OUT + "/y_target_%s.npy" % key)
    ys = np.load(OUT + "/y_source_%s.npy" % key)
    dn = donors(tgt)
    assert len(dn) == len(y)

    preds = {}
    for m in MODELS:
        preds[m] = np.load(OUT + "/pred_%s_%s.npy" % (m, key))

    # baseline 1: locked-protocol 1152-d mean pooling on the same cap500 embeddings
    preds["mean_pool_1152d"] = np.load(OUT + "/harness_pred_mean_%s.npy" % key)

    # baseline 2: PCA32 donor mean + logistic, same source-fitted projection
    fs = np.load(OUT + "/feat_pca32_mean_source_%s.npy" % key)
    ft = np.load(OUT + "/feat_pca32_mean_target_%s.npy" % key)
    C, _ = protocol.select_C(fs, ys)
    p_pca32, _ = protocol.transfer(fs, ys, ft, C=C)
    preds["pca32_donor_mean"] = p_pca32
    print("%s  PCA32 donor-mean control: C=%g AUC=%.4f" % (key, C, protocol.roc_auc(y, p_pca32)), flush=True)

    # archived pseudobulk comparators, reindexed to the donor order used here
    a = arch[arch.direction == key.replace("GSE", "GSE")]
    for mid, name in (("source_hvg_pseudobulk", "hvg_pseudobulk"), ("source_pca_pseudobulk", "pca_pseudobulk")):
        g = a[a.method_id == mid].set_index("donor_id")
        assert set(g.index) == set(dn), "donor id mismatch for %s / %s" % (key, mid)
        assert np.array_equal(g.loc[dn, "y_true"].to_numpy().astype(float), y), "label mismatch %s" % mid
        preds[name] = g.loc[dn, "prob_case"].to_numpy()

    for name, p in preds.items():
        met = protocol.evaluate(y, p)
        lo, hi = protocol.boot_auc_ci(y, p, n=BOOT)
        cv = cv_lookup.get((name, key))
        rows.append(dict(direction=key, method=name, n_source=len(ys), n_target=len(y),
                         source_cv_auc=(cv or {}).get("source_cv_auc", ""),
                         best_H=(cv or {}).get("best_H", ""), best_wd=(cv or {}).get("best_wd", ""),
                         target_auc=met["roc_auc"], auc_ci_low=lo, auc_ci_high=hi,
                         pr_auc=met["pr_auc"], brier=met["brier"], ece_10bin=met["ece_10bin"]))

    for m in MODELS:
        for base, family in (("mean_pool_1152d", "primary"), ("pca32_donor_mean", "primary"),
                             ("hvg_pseudobulk", "secondary"), ("pca_pseudobulk", "secondary")):
            diff, p = protocol.delong_paired_p(y, preds[m], preds[base])
            lo, hi = protocol.boot_auc_diff_ci(y, preds[m], preds[base], n=BOOT)
            tests.append(dict(direction=key, family=family, method=m, baseline=base,
                              auc_method=protocol.roc_auc(y, preds[m]),
                              auc_baseline=protocol.roc_auc(y, preds[base]),
                              auc_diff=diff, diff_ci_low=lo, diff_ci_high=hi, delong_p=p))

t = pd.DataFrame(tests)
t["delong_bh_q"] = np.nan
for fam in ("primary", "secondary"):
    mask = t.family == fam
    t.loc[mask, "delong_bh_q"] = protocol.bh(t.loc[mask, "delong_p"].to_numpy())

pd.DataFrame(rows).to_csv(OUT + "/learned_pooling_metrics.tsv", sep="\t", index=False)
t.to_csv(OUT + "/learned_pooling_paired_tests.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
print("\n=== target metrics ===")
print(pd.DataFrame(rows)[["direction", "method", "source_cv_auc", "target_auc", "auc_ci_low", "auc_ci_high", "brier"]]
      .to_string(index=False, float_format=lambda v: "%.4f" % v))
print("\n=== paired tests ===")
print(t[["direction", "family", "method", "baseline", "auc_diff", "diff_ci_low", "diff_ci_high", "delong_p", "delong_bh_q"]]
      .to_string(index=False, float_format=lambda v: "%.4f" % v))
