"""Harness validation for the archived fixed pooling operators.

For each of the four fixed operators that can be rebuilt from cell-level embeddings
(mean, coordinate median, 10% coordinate-wise trimmed mean, mean+std), this script

  1. recomputes the donor features locally,
  2. scores both transfer directions under the locked deterministic protocol, and
  3. scans the C grid to check whether the archived AUC is reproducible at SOME C.

Step 3 is the point: if every archived value is hit somewhere on the grid, the donor
features are bit-faithful and any discrepancy comes only from inner-CV regularisation
selection (the archived run used sklearn's LogisticRegressionCV with method- and
direction-specific fold seeds).

Usage:  python3 scripts/validate_harness.py inputs results
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import protocol

IN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

# archived values: 08_.../outputs/distributional_pooling_cap500/transfer_metrics.tsv
ARCHIVED = {
    ("GSE174188_to_GSE285773", "mean"): 0.8624999999999999,
    ("GSE174188_to_GSE285773", "median"): 0.9187500000000002,
    ("GSE174188_to_GSE285773", "trimmed_mean_10"): 0.8562500000000001,
    ("GSE174188_to_GSE285773", "mean_std"): 0.825,
    ("GSE285773_to_GSE174188", "mean"): 0.8842748472378101,
    ("GSE285773_to_GSE174188", "median"): 0.8580246913580246,
    ("GSE285773_to_GSE174188", "trimmed_mean_10"): 0.8918817807706696,
    ("GSE285773_to_GSE174188", "mean_std"): 0.8852724778650705,
}


def donors(tag):
    f = IN + ("/GSE174188_donor_index_v2.csv" if tag == "GSE174188" else "/GSE285773_donor_index.csv")
    with open(f) as fh:
        return [r["donor_id"] for r in csv.DictReader(fh)]


def labels(tag):
    f, sep = ((IN + "/gse174188_final_donor_covariates.tsv", "\t") if tag == "GSE174188"
              else (IN + "/labels_285773.tsv", "\t"))
    with open(f) as fh:
        m = {r["donor_id"]: r["case_control"] for r in csv.DictReader(fh, delimiter=sep)}
    return np.array([1.0 if m[d] == "case" else 0.0 for d in donors(tag)])


def operator(tag, kind):
    mean = np.load(IN + "/%s_pool_mean.npy" % tag)
    if kind == "mean":
        return mean
    if kind == "median":
        return np.load(IN + "/%s_pool_median.npy" % tag)
    if kind == "trimmed_mean_10":
        return np.load(IN + "/%s_pool_trimmed10.npy" % tag)
    if kind == "mean_std":
        return np.hstack([mean, np.load(IN + "/%s_pool_std.npy" % tag)])
    raise ValueError(kind)


Y = {t: labels(t) for t in ("GSE174188", "GSE285773")}
rows = []
for src, tgt in (("GSE174188", "GSE285773"), ("GSE285773", "GSE174188")):
    key = "%s_to_%s" % (src, tgt)
    for kind in ("mean", "median", "trimmed_mean_10", "mean_std"):
        Xs, Xt = operator(src, kind), operator(tgt, kind)
        C, _ = protocol.select_C(Xs, Y[src])
        p, _ = protocol.transfer(Xs, Y[src], Xt, C=C)
        m = protocol.evaluate(Y[tgt], p)
        arch = ARCHIVED[(key, kind)]
        hits = protocol.reproducible_at(Xs, Y[src], Xt, Y[tgt], arch, tol=1e-9)
        rows.append(dict(direction=key, operator=kind, n_features=Xs.shape[1],
                         selected_C=C, protocol_auc=m["roc_auc"], archived_auc=arch,
                         same_as_archived=abs(m["roc_auc"] - arch) < 1e-9,
                         archived_reproducible_at_C=";".join("%g" % c for c, _ in hits) or "none",
                         pr_auc=m["pr_auc"], brier=m["brier"], ece_10bin=m["ece_10bin"]))
        print("%-26s %-16s C=%-8g protocol=%.6f archived=%.6f  archived reproducible at C=[%s]"
              % (key, kind, C, m["roc_auc"], arch, rows[-1]["archived_reproducible_at_C"]), flush=True)
        np.save(OUT + "/harness_pred_%s_%s.npy" % (kind, key), p)

with open(OUT + "/harness_validation.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

n_repro = sum(1 for r in rows if r["archived_reproducible_at_C"] != "none")
print("\n%d/%d archived AUCs reproducible somewhere on the C grid" % (n_repro, len(rows)), flush=True)
