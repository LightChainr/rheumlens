"""Reproduction gate: do the rebuilt donor representations hit the archived AUCs?

Rebuilt cap500 Geneformer donor means (from the 63 GSE174188 shards and the 26
GSE285773 donor shards) are pushed through the locked source-only transfer protocol.
The gate passes only if both directions reproduce the archived cap500 mean-pooling
ROC-AUC exactly:

    GSE174188 -> GSE285773   0.8625000000
    GSE285773 -> GSE174188   0.8842748472

Usage:  python3 scripts/repro_gate.py inputs results
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import protocol

IN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

ARCHIVED = {"GSE174188_to_GSE285773": 0.8625000000000000,
            "GSE285773_to_GSE174188": 0.8842748472378102}


def donors(tag):
    f = IN + ("/GSE174188_donor_index_v2.csv" if tag == "GSE174188" else "/GSE285773_donor_index.csv")
    with open(f) as fh:
        return [r["donor_id"] for r in csv.DictReader(fh)]


def labels(tag):
    if tag == "GSE174188":
        with open(IN + "/gse174188_final_donor_covariates.tsv") as fh:
            m = {r["donor_id"]: r["case_control"] for r in csv.DictReader(fh, delimiter="\t")}
    else:
        with open(IN + "/labels_285773.tsv") as fh:
            m = {r["donor_id"]: r["case_control"] for r in csv.DictReader(fh, delimiter="\t")}
    return np.array([1.0 if m[d] == "case" else 0.0 for d in donors(tag)])


def features(tag):
    f = IN + ("/GSE174188_donor_mean_f64_v2.npy" if tag == "GSE174188" else "/GSE285773_donor_mean_f64.npy")
    return np.load(f)


DATA = {t: dict(X=features(t), y=labels(t), donors=donors(t)) for t in ("GSE174188", "GSE285773")}
for t, d in DATA.items():
    assert len(d["X"]) == len(d["y"]) == len(d["donors"])
    print("%s donors=%d cases=%d dim=%d" % (t, len(d["y"]), int(d["y"].sum()), d["X"].shape[1]), flush=True)

    # the independently rebuilt pooling-operator mean must equal the donor-mean file
    pm = np.load(IN + "/%s_pool_mean.npy" % t)
    with open(IN + "/%s_pool_index.csv" % t) as fh:
        pool_donors = [r["donor_id"] for r in csv.DictReader(fh)]
    assert pool_donors == d["donors"], "pool_index donor order differs from donor_index"
    print("   pool_mean identical to donor_mean: %s (max abs diff %.3g)"
          % (np.allclose(pm, d["X"]), np.abs(pm - d["X"]).max()), flush=True)

rows = []
for src, tgt in (("GSE174188", "GSE285773"), ("GSE285773", "GSE174188")):
    S, T = DATA[src], DATA[tgt]
    C, grid = protocol.select_C(S["X"], S["y"])
    p, _ = protocol.transfer(S["X"], S["y"], T["X"], C=C)
    m = protocol.evaluate(T["y"], p)
    key = "%s_to_%s" % (src, tgt)
    hit = abs(m["roc_auc"] - ARCHIVED[key]) < 1e-9
    rows.append(dict(direction=key, n_source=len(S["y"]), n_target=len(T["y"]), selected_C=C,
                     recomputed_auc=m["roc_auc"], archived_auc=ARCHIVED[key],
                     abs_diff=abs(m["roc_auc"] - ARCHIVED[key]), passed=hit,
                     pr_auc=m["pr_auc"], brier=m["brier"], ece_10bin=m["ece_10bin"]))
    print("%-26s C=%-8g recomputed=%.10f archived=%.10f  %s"
          % (key, C, m["roc_auc"], ARCHIVED[key], "PASS" if hit else "FAIL"), flush=True)
    np.save(OUT + "/repro_gate_pred_%s.npy" % key, p)

with open(OUT + "/repro_gate.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

ok = all(r["passed"] for r in rows)
print("\nGATE", "PASS" if ok else "FAIL", flush=True)
raise SystemExit(0 if ok else 1)
