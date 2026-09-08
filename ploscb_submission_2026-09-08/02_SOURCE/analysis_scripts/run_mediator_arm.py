"""Supplementary: what admitting a disease-caused covariate does to the audit.

Section 4.1 states that covariates caused by the disease (severity, medication,
symptom timing) must be kept out of the design variables, and that admitting them
inflates the result "by construction, not by discovery". This quantifies that claim
on data where the truth is known: the biological effect is held fixed, and the only
thing that changes is whether one disease-caused variable is treated as a design
variable.
"""
from __future__ import annotations
import hashlib
import itertools
import os
from pathlib import Path
import numpy as np, pandas as pd
from joblib import Parallel, delayed
import importlib.util

spec = importlib.util.spec_from_file_location(
    "sim", str(Path(__file__).resolve().parent / "run_extended_simulation.py"))
sim = importlib.util.module_from_spec(spec); spec.loader.exec_module(sim)

OUT = Path(__file__).resolve().parent / "results"; OUT.mkdir(exist_ok=True)
RHOS = [0.0, 0.25, 0.50, 0.75, 0.90]
N_REP = 200


def cell_seed(rho, admit, rep) -> int:
    """Stable across interpreter sessions.

    This used Python's built-in hash(), which is salted per process for str, so the
    seeds - and therefore every number this arm reports - did not reproduce in a
    fresh interpreter. Derived the same way as the main simulation's seeds.
    """
    key = f"{sim.MASTER_SEED}|mediator|{rho:.3f}|{admit}|{rep}".encode()
    return int.from_bytes(hashlib.blake2b(key, digest_size=4).digest(), "big") % (2**31)


def one(rho, admit, rep):
    seed = cell_seed(rho, admit, rep)
    rng = np.random.default_rng(seed)
    X, y, y_lat, D, strata, load = sim.make_cohort(rng, rho, "base", "binary")

    # a covariate caused by the disease: severity is a noisy function of the label
    severity = 1.2 * (y - y.mean()) + rng.standard_normal(len(y))
    severity = (severity - severity.mean()) / severity.std()
    Duse = np.column_stack([D, severity]) if admit else D

    y_c = y.astype(float) - y.mean()
    des = sim.cv_auc(Duse, y, seed)
    i_cv, i_null, p_id = sim.id_contrast(y_c, Duse, seed)
    res = sim.residualised_auc(X, y, Duse, seed)
    dis = sim.cv_auc(X, y, seed)
    return dict(rho=rho, admitted=admit, rep=rep, seed=seed, design_auc=des, I_D_cv=i_cv,
                p_I_D=p_id, disease_auc=dis, residualised_auc=res)


if __name__ == "__main__":
    cells = list(itertools.product(RHOS, [False, True], range(N_REP)))
    print(f"{len(cells)} runs", flush=True)
    n_jobs = int(os.environ.get("RHEUMLENS_JOBS", "4"))
    rows = Parallel(n_jobs=n_jobs, batch_size=32)(delayed(one)(*c) for c in cells)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "mediator_arm_raw.tsv", sep="\t", index=False)
    g = (d.groupby(["rho", "admitted"])
           .agg(design_auc=("design_auc", "mean"), I_D_cv=("I_D_cv", "mean"),
                sig=("p_I_D", lambda s: (s <= .05).mean()),
                disease_auc=("disease_auc", "mean"),
                residualised_auc=("residualised_auc", "mean"))
           .round(3).reset_index())
    g.to_csv(OUT / "mediator_arm_summary.tsv", sep="\t", index=False)
    print(g.to_string(index=False))
