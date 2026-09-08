#!/usr/bin/env python3
"""Cross-check every claim in the manuscript that names a cohort and a number.

Not a proof of correctness - it checks a list of assertions that were written by
hand and are the ones most likely to drift when the analysis is re-run.
"""
from __future__ import annotations
import re, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TEXT = (ROOT / "manuscript" / "manuscript_v5.md").read_text()


def seed_of(f: Path) -> int:
    for part in reversed(f.parts):
        m = re.fullmatch(r"seed_(\d{8})", part)
        if m:
            return int(m.group(1))
    raise SystemExit(f"cannot read a seed from {f}")


df = pd.concat([pd.read_csv(f, sep="\t").assign(seed=seed_of(f))
                for f in sorted((ROOT / "results/screen").glob("**/design_screen.tsv"))])
wt = pd.read_csv(ROOT / "walkthrough/results/decision_tree_walkthrough.tsv", sep="\t")
cal_a = pd.read_csv(ROOT / "sim/results/calibration_arm_a_summary.tsv", sep="\t")
cal_b = pd.read_csv(ROOT / "sim/results/calibration_arm_b_summary.tsv", sep="\t")
med = pd.read_csv(ROOT / "sim/results/mediator_arm_summary.tsv", sep="\t")
med0 = med[med.rho == 0.0].set_index("admitted")


def r(coh, block, col, fmt="{:.3f}"):
    g = df[(df.cohort == coh) & (df.block == block)][col].dropna()
    lo, hi = g.min(), g.max()
    return fmt.format(lo) if lo == hi else f"{fmt.format(lo)}-{fmt.format(hi)}"


CHECKS: list[tuple[str, str]] = [
    # (description, string that must appear verbatim in the manuscript)
    ("CMV all AUC",            r("CMV_HIHA", "all", "design_auc_linear").replace("-", " to ")),
    ("CMV collection AUC",     r("CMV_HIHA", "batch", "design_auc_linear").replace("-", " to ")),
    ("CMV quality AUC",        r("CMV_HIHA", "qc", "design_auc_linear").replace("-", " to ")),
    ("CMV demographic AUC",    r("CMV_HIHA", "demographic", "design_auc_linear").replace("-", " to ")),
    ("CMV diagnosis AUC",      r("CMV_HIHA", "all", "disease_auc").replace("-", " to ")),
    ("Stephenson all AUC",     r("COVID_STEPHENSON", "all", "design_auc_linear")),
    ("Stephenson quality AUC", r("COVID_STEPHENSON", "qc", "design_auc_linear")),
    ("Stephenson coll AUC",    r("COVID_STEPHENSON", "batch", "design_auc_linear")),
    ("Ren all AUC",            r("COVID_REN", "all", "design_auc_linear")),
    ("Ren collection AUC",     r("COVID_REN", "batch", "design_auc_linear")),
    ("Ren quality AUC",        r("COVID_REN", "qc", "design_auc_linear")),
    ("COMBAT_COVID all AUC",   r("COMBAT_COVID", "all", "design_auc_linear")),
    ("Cross all AUC",          r("COMBAT_CROSS", "all", "design_auc_linear").replace("-", " to ")),
    ("CMV residualisation loss",
     "{:.3f}".format(float(wt[wt.cohort == "CMV_HIHA"].residualisation_loss.iloc[0]))),
    ("CMV walkthrough design AUC",
     "{:.3f}".format(float(wt[wt.cohort == "CMV_HIHA"].design_only_auc.iloc[0]))),
    ("CMV walkthrough disease AUC",
     "{:.3f}".format(float(wt[wt.cohort == "CMV_HIHA"].disease_auc.iloc[0]))),
    ("CMV walkthrough residualised AUC",
     "{:.3f}".format(float(wt[wt.cohort == "CMV_HIHA"].residualised_auc.iloc[0]))),
    ("calibration gamma=0, collection-preserving",
     "{:.1f}%".format(100 * float(cal_a[cal_a.gamma == 0].reject_collection_preserving.iloc[0]))),
    ("calibration gamma=0, free",
     "{:.1f}%".format(100 * float(cal_a[cal_a.gamma == 0].reject_unstratified.iloc[0]))),
    ("arm B loss at 92 columns",
     "{:.3f}".format(float(cal_b[(cal_b.layout == "balanced") &
                                 (cal_b.n_design_col == 92)].loss_mean.iloc[0]))),
    ("mediator arm, covariate excluded",
     "{:.3f}".format(float(med0.loc[False, "design_auc"]))),
    ("mediator arm, covariate admitted",
     "{:.3f}".format(float(med0.loc[True, "design_auc"]))),
    ("mediator arm, flag rate without the covariate",
     "{:.1f}%".format(100 * float(med0.loc[False, "sig"]))),
    ("arm B loss, CMV-shaped",
     "{:.3f}".format(float(cal_b[cal_b.layout == "cmv_ragged"].loss_mean.iloc[0]))),
]

bad = 0
for name, needle in CHECKS:
    ok = needle in TEXT
    if not ok:
        bad += 1
    print(f"{'OK  ' if ok else 'MISS'}  {name:<42} {needle}")
print(f"\n{len(CHECKS) - bad}/{len(CHECKS)} assertions found verbatim in the manuscript.")
sys.exit(1 if bad else 0)
