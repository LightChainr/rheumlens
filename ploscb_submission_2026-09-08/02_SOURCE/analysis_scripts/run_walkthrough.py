"""Full decision-tree walkthrough on three real cohorts.

Three cohorts are traced step by step through the decision tree with the evidence for
each branch stated, chosen because they terminate on three different routes: CMV_HIHA
(recorded metadata uninformative), COVID_REN (metadata strongly predictive, mixed strata
available) and COMBAT_INFLUENZA (metadata determines the label, no mixed stratum).

Q4 requires restriction and residualisation on the same donors, and this script runs
both.

The metadata matrix, the strata and the classifier are all IMPORTED rather than
rebuilt here. An earlier version reimplemented the matrix, and the reimplementation
drifted: it admitted every non-label column instead of the screen's whitelist, so CMV
came out with 92 columns where the screen builds 89, and the two halves of the paper
were describing different matrices. Pre-submission review then found the same failure
one level down - this script standardised the whole-cohort matrix before handing it to a
fold-contained ridge, and ridge shrinkage depends on scale, so a step the manuscript
described as inside the fold was partly outside it.

Both are now structural. `pipeline_core.WALKTHROUGH` is `SCREEN_FROZEN` with twenty
repeats instead of one, and `pipeline_core.fold_design` builds the metadata matrix from
the training donors of each split, so the imputation, the encoding and the
standardisation cannot be fitted anywhere else.
"""
from __future__ import annotations
import json, warnings
from dataclasses import replace
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")
REPO = Path(os.environ.get(
    "RHEUMLENS_REPO",
    Path(__file__).resolve().parents[2] / "20_repo_restructure_20260907"))
OUT  = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)

SEED = 20260907
COVMAP = {"COVID_REN": "covid_ren_donor_covariates.tsv",
          "CMV_HIHA": "cmv_hiha_donor_covariates.tsv",
          "COMBAT_INFLUENZA": "combat_influenza_donor_covariates.tsv"}


def load(cohort):
    X = pd.read_parquet(REPO / "inputs/donor_level" / cohort / "donor_log1p_cpm.parquet")
    lab = pd.read_csv(REPO / "inputs/donor_level" / cohort / "donor_labels.tsv", sep="\t")
    cov = pd.read_csv(REPO / "inputs/design_metadata" / COVMAP[cohort], sep="\t")
    X = X.loc[lab.donor_id.values] if X.index.name else X
    cov = cov.set_index("donor_id").loc[lab.donor_id.values].reset_index()
    return X.to_numpy(dtype=float), lab.y_true.to_numpy(int), cov


def _load_screen():
    """Import the screen module so its definitions are used, not copied."""
    import importlib.util
    path = REPO / "scripts" / "cohorts" / "run_design_screen.py"
    spec = importlib.util.spec_from_file_location("rheumlens_screen", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SCREEN = _load_screen()
core = SCREEN.core
SPEC = core.WALKTHROUGH
N_REPEAT = SPEC.n_repeat


def _cov_indexed(cov):
    return cov.set_index("donor_id") if "donor_id" in cov.columns else cov


def design_folds(cov):
    """`get_fold` for the screen's `all` block, built inside the training fold.

    Not standardised here. `fold_contained_auc` standardises it on the training
    donors before the ridge fit, which is the only place a scale can be chosen
    without seeing a held-out donor.
    """
    return core.design_folds(_cov_indexed(cov), "all")


def design_width(cov):
    """How wide the declared matrix is, and how much of the sample it spans.

    Column count is not effective dimension - the CMV collection block has 80
    columns and a centred rank of 45 - so the rank is reported beside it wherever
    the width is used to explain a residualisation loss.
    """
    D = core.design_matrix(_cov_indexed(cov))["all"]
    return int(D.shape[1]), int(np.linalg.matrix_rank(D - D.mean(0)))


def strata_id(cov):
    """Collection stratum, delegated to `collection_strata` in run_design_screen.py.

    Q4 of the decision tree must use the same stratum definition as the
    collection-preserving null at Q2, or the figure traces two different questions.
    """
    idx, label = SCREEN.collection_strata(_cov_indexed(cov))
    return idx.astype(str), label


def cv_auc(X, y, seed, resid_folds=None, n_repeat=None, get_fold=None):
    """Repeated cross-fitted AUC, and the spread over repeats.

    `resid_folds(tr, te)` supplies the adjustment matrix for a split. The adjusted
    and unadjusted arms are the same call with that one argument added, so their
    difference is the cost of adjustment and not a change of representation.
    """
    if len(np.unique(y)) < 2 or min(np.bincount(y)) < 5:
        return np.nan, np.nan
    spec = SPEC if n_repeat is None else replace(SPEC, n_repeat=n_repeat)
    if get_fold is None:
        get_fold = core.array_folds(X)
    per = [core.fold_contained_auc(get_fold, y, replace(spec, n_repeat=1),
                                   seed + r, residualise=resid_folds)
           for r in range(spec.n_repeat)]
    return float(np.mean(per)), float(np.std(per))


def matched_subsets(X, y, n_sub, n_case, seed, n_draw=20):
    """AUC in random subsets of the same size and case/control composition.

    Restriction always costs donors, and a smaller cohort scores differently for
    that reason alone. This is the comparison that separates the two: same size,
    same class balance, drawn without regard to stratum.
    """
    rng = np.random.default_rng(seed)
    ci, co = np.where(y == 1)[0], np.where(y == 0)[0]
    n_ctrl = n_sub - n_case
    if n_case > len(ci) or n_ctrl > len(co) or n_case < 5 or n_ctrl < 5:
        return np.nan, np.nan
    vals = []
    for _ in range(n_draw):
        idx = np.concatenate([rng.choice(ci, n_case, False), rng.choice(co, n_ctrl, False)])
        a, _ = cv_auc(X[idx], y[idx], seed, n_repeat=3)
        if a == a:
            vals.append(a)
    return (float(np.mean(vals)), float(np.std(vals))) if vals else (np.nan, np.nan)


def walk(cohort):
    X, y, cov = load(cohort)
    fold_D = design_folds(cov)
    n_col, rank = design_width(cov)
    strata, scol = strata_id(cov)

    rec = dict(cohort=cohort, n_donor=len(y), n_case=int(y.sum()),
               n_design_col=n_col, n_design_rank=rank, stratum_variable=scol)

    # Q2 evidence
    da, _ = cv_auc(None, y, SEED, get_fold=fold_D)
    rec["design_only_auc"] = round(da, 4)

    # Q3 evidence: strata containing both labels
    tab = pd.crosstab(strata, y)
    mixed = tab[(tab > 0).all(axis=1)]
    rec["n_strata"] = int(tab.shape[0])
    rec["n_mixed_strata"] = int(mixed.shape[0])
    rec["n_donor_in_mixed"] = int(mixed.to_numpy().sum())

    # unadjusted
    ua, us = cv_auc(X, y, SEED)
    rec["disease_auc"] = round(ua, 4); rec["disease_auc_sd"] = round(us, 4)

    # Q4a residualisation on the full metadata block
    ra, rs = cv_auc(X, y, SEED, resid_folds=fold_D)
    rec["residualised_auc"] = round(ra, 4) if ra == ra else None
    rec["residualisation_loss"] = round(ua - ra, 4) if ra == ra else None

    # Q4b restriction to the largest mixed stratum, vs matched random subsets
    # a stratum is usable for restriction only if 5-fold CV can run inside it
    usable = mixed[(mixed >= 5).all(axis=1)]
    rec["n_usable_strata"] = int(usable.shape[0])
    if usable.shape[0] > 0:
        big = usable.sum(axis=1).idxmax()
        m = strata == big
        sa, _ = cv_auc(X[m], y[m], SEED)
        ma, msd = matched_subsets(X, y, int(m.sum()), int(y[m].sum()), SEED)
        rec.update(restriction_stratum=str(big), n_in_stratum=int(m.sum()),
                   n_case_in_stratum=int(y[m].sum()),
                   restricted_auc=round(sa, 4) if sa == sa else None,
                   matched_subset_auc=round(ma, 4) if ma == ma else None,
                   matched_subset_sd=round(msd, 4) if msd == msd else None,
                   restriction_loss=round(sa - ma, 4) if (sa == sa and ma == ma) else None)
    else:
        rec.update(restriction_stratum=None, restricted_auc=None,
                   matched_subset_auc=None, restriction_loss=None)
    return rec


if __name__ == "__main__":
    rows = []
    for c in ["CMV_HIHA", "COVID_REN", "COMBAT_INFLUENZA"]:
        print("running", c, flush=True)
        rows.append(walk(c))
        print(json.dumps(rows[-1], indent=2, ensure_ascii=False), flush=True)
    pd.DataFrame(rows).to_csv(OUT / "decision_tree_walkthrough.tsv", sep="\t", index=False)
    print("wrote", OUT / "decision_tree_walkthrough.tsv")
