#!/usr/bin/env python3
"""Stage the simulation and calibration outputs as flat supplementary tables.

The simulation summary is written by pandas with a two-level column index, which
is fine for analysis and unreadable as a released table. Flatten it, and give the
calibration arms their own table rather than leaving them only in sim/results.
"""
from pathlib import Path
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from si_names import rename_columns, check_header

ROOT = Path(__file__).resolve().parents[1]
SUP = ROOT / "supplementary"
RES = ROOT / "sim" / "results"

# ---- Table S2: extended simulation + mediator arm, one flat header ------------
sim = pd.read_csv(RES / "extended_simulation_summary.tsv", sep="\t", header=[0, 1],
                  index_col=[0, 1, 2])
sim.columns = [f"{a}_{b}" for a, b in sim.columns]
sim.index.names = ["arm", "label_type", "rho"]
sim = sim.reset_index()
sim.columns = rename_columns(sim.columns)
med = pd.read_csv(RES / "mediator_arm_summary.tsv", sep="\t")
med.insert(0, "arm", "mediator")
med.insert(1, "label_type", "binary")
med.columns = rename_columns(med.columns)

out = SUP / "Table_S2_simulation_summary.tsv"
with out.open("w") as fh:
    fh.write("# Part 1: extended simulation, 7 arms x 7 rho x 2 label types x 200 reps.\n")
    fh.write("# c = concordance index (= ROC AUC for a binary outcome).\n")
    sim.round(4).to_csv(fh, sep="\t", index=False)
    fh.write("\n# Part 2: mediator arm. 'admitted' = a disease-caused covariate was\n")
    fh.write("# admitted to the design variables. 'sig' = fraction of runs flagged.\n")
    med.round(4).to_csv(fh, sep="\t", index=False)

# ---- Table S5: permutation calibration and residualisation width -------------
a = pd.read_csv(RES / "calibration_arm_a_summary.tsv", sep="\t")
b = pd.read_csv(RES / "calibration_arm_b_summary.tsv", sep="\t")
a.columns = rename_columns(a.columns)
b.columns = rename_columns(b.columns)
out9 = SUP / "Table_S5_calibration.tsv"
with out9.open("w") as fh:
    fh.write("# Arm A: type-I error of both permutation tests when a sample-quality\n")
    fh.write("# variable predicts the diagnosis WITHIN a collection stratum and no\n")
    fh.write("# disease effect exists. 300 replicates per gamma, 200 permutations each.\n")
    a.round(4).to_csv(fh, sep="\t", index=False)
    fh.write("\n# Arm B: residualisation loss when the design is independent of the\n")
    fh.write("# diagnosis by construction. 108 donors, 60 replicates per layout.\n")
    fh.write("# 'cmv_ragged' reuses the observed CMV batch-by-pool level sizes.\n")
    b.round(4).to_csv(fh, sep="\t", index=False)

for f in (out, out9):
    print(f"wrote {f.relative_to(ROOT)}  ({f.stat().st_size} bytes)")

# Nothing released may still carry a retired name.
for f in (out, out9):
    header = [ln for ln in f.read_text().splitlines() if not ln.startswith("#")][0]
    bad = check_header(header.split("\t"))
    if bad:
        raise SystemExit(f"{f.name} still uses retired names: {bad}")
