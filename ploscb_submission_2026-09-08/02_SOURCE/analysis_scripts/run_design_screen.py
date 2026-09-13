"""Multi-cohort design screen, including the collection-preserving null.

Outputs one row per (cohort, design block) with

    I_D              = 1 - R^2(Y ~ D)      linear label information left after D
    design_auc_lin   cross-fitted logistic  Y ~ D
    design_auc_rf    cross-fitted random forest Y ~ D   (nonlinear design signal)
    disease_auc      cross-fitted logistic  Y ~ pseudobulk
    p_standard       free label permutation p-value
    p_collection_preserving  permutation WITHIN collection strata

These are the analysis-layer column names. tools/si_names.py translates them
to the manuscript vocabulary when a supplementary table is exported, so
`I_D` here is `V_D_insample` in Table S8 and `p_standard` is `p_free`.

Why two nulls
-------------
PLOS reviewer #4 observed that a standard complete-pipeline label permutation
destroys the design-label association along with the biological signal, so it
tests leakage and pipeline validity rather than "can the classifier exploit
acquisition structure". The restricted null permutes the label only within
strata that share a COLLECTION configuration (batch, pool, site, assay), so the
association between the diagnosis and those variables is retained while
disease-specific molecular information is destroyed. A classifier that still
scores above this null is using something beyond those collection strata.

Note the limit, which the manuscript states in the same words: the strata are
built from collection variables only. Sample-quality and demographic columns sit
in the declared design matrix but are NOT conditioned on, so passing this test
does not mean the classifier beat everything recorded.

    p_standard small, p_collection_preserving small  -> signal beyond the
                                                        collection strata
    p_standard small, p_collection_preserving LARGE  -> what looked like signal
                                                        is the collection strata

Usage
    python scripts/cohorts/run_design_screen.py --all
    python scripts/cohorts/run_design_screen.py --cohorts CMV_HIHA COMBAT_CROSS
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from joblib import Parallel, delayed

try:                                     # normal import, when run as a module
    from . import pipeline_core as core
except ImportError:                      # run as a script, or exec'd by path
    # A real import off sys.path, not spec_from_file_location: a joblib worker
    # re-imports by module NAME, and a module that exists only as an exec'd spec
    # in the parent cannot be un-pickled in the child.
    import sys as _sys
    _here = str(Path(__file__).resolve().parent)
    if _here not in _sys.path:
        _sys.path.insert(0, _here)
    import pipeline_core as core

REPO = Path(__file__).resolve().parents[2]
INPUTS = REPO / "inputs"
OUT = REPO / "results" / "multi_cohort_design_screen"

SEED = 20260907
N_SPLITS = 5
N_REPEAT = 5
N_PERM = 1000
CS = np.logspace(-4, 4, 9)
WORKERS = 1        # permutation-stage processes; see permutation_pvalues


def design_matrix(cov: pd.DataFrame) -> dict[str, np.ndarray]:
    """Whole-cohort metadata blocks, delegated to pipeline_core.

    This is the matrix the IN-SAMPLE V_D is defined on, and the one whose width
    is reported as `n_design_feature`. Every cross-fitted quantity is built with
    `core.fold_design` instead, which takes the medians and the categorical level
    set from the training donors alone. Both use `core.block_columns`, so the two
    cannot be assembled from different variables.
    """
    return core.design_matrix(cov)


def label_information_fraction(y: np.ndarray, D: np.ndarray,
                               n_perm: int = 200,
                               seed: int | None = None,
                               fold_D=None) -> dict:
    """Design-adjusted label contrast, reported against its own permutation null.

    The v2 quantity was 1 - R^2(Y ~ D) computed IN SAMPLE. In-sample R^2 inflates
    with the ratio of design features to donors, so its null baseline differs by
    cohort and the raw value cannot be compared across cohorts. Concretely: a
    design block with 89 features on 108 donors has a null I_D near 0.16, while a
    block with 8 features on 182 donors has a null near 0.95. Ranking cohorts by
    the raw number then reports the cohort with the most design columns as the
    most confounded, regardless of any actual association.

    Three quantities are returned:
      i_d_insample : the v2 definition, kept for continuity
      i_d_cv       : cross-fitted, so it is not inflated by feature count
      i_d_null_*   : the permutation null of the in-sample value for THIS design
                     matrix, which is what makes the in-sample number readable
      i_d_cv_null_*: the permutation null of the CROSS-FITTED value, computed by
                     pushing each permuted label through the same folds

    p_i_d is the p-value of the cross-fitted statistic against the cross-fitted
    null - the same quantity on both sides. p_i_d_insample is the in-sample pair,
    kept because the v2 results were reported that way.

    PLOS reviewer #1 raised exactly this: "the boundary 1-R^2(Y~D) depends on the
    parameterisation of the design matrix, which can be problematic."

    `fold_D(tr, te)` supplies the training and held-out metadata matrices for one
    split, imputed and encoded from the training donors only. The in-sample pair
    still uses the whole-cohort matrix `D`, because that is what an in-sample
    statistic is; everything cross-fitted goes through `fold_D`. Passing None
    falls back to slicing `D`, which is the transductive construction the
    pre-submission review flagged and is kept only so the older behaviour can be
    reproduced on demand.
    """
    # Read SEED at call time: a default argument would bind the module-level
    # value at definition time and silently ignore --seed.
    if seed is None:
        seed = SEED

    r2_in = LinearRegression().fit(D, y).score(D, y)
    i_d_in = float(1.0 - max(r2_in, 0.0))

    n = len(y)
    if fold_D is None:
        fold_D = lambda tr, te: (D[tr], D[te])
    oof = np.zeros(n, dtype=float)
    kf = KFold(min(5, n), shuffle=True, random_state=seed)
    folds = list(kf.split(D))
    fold_mats = [fold_D(tr, te) for tr, te in folds]
    for (tr, te), (Dtr, Dte) in zip(folds, fold_mats):
        oof[te] = LinearRegression().fit(Dtr, y[tr]).predict(Dte)
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2_cv = 1.0 - ((y - oof) ** 2).sum() / ss_tot if ss_tot > 0 else 0.0
    i_d_cv = float(1.0 - max(r2_cv, 0.0))

    def _cv(yy: np.ndarray) -> float:
        """Cross-fitted V_D for a label vector, using the folds fixed above."""
        o = np.zeros(n, dtype=float)
        for (tr, te), (Dtr, Dte) in zip(folds, fold_mats):
            o[te] = LinearRegression().fit(Dtr, yy[tr]).predict(Dte)
        sst = ((yy - yy.mean()) ** 2).sum()
        r2 = 1.0 - ((yy - o) ** 2).sum() / sst if sst > 0 else 0.0
        return float(1.0 - max(r2, 0.0))

    rng = np.random.default_rng(seed)
    null_in = np.empty(n_perm, dtype=float)
    null_cv = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        yp = rng.permutation(y)
        null_in[k] = 1.0 - max(LinearRegression().fit(D, yp).score(D, yp), 0.0)
        # The permuted label goes through the SAME cross-fitting, with the same
        # folds, as the observed value. Reporting a cross-fitted statistic against
        # an in-sample null would be comparing two different quantities.
        null_cv[k] = _cv(yp)

    # one-sided: how often does a random label leave AT MOST as much contrast?
    p_in = float((np.sum(null_in <= i_d_in) + 1) / (n_perm + 1))
    p_cv = float((np.sum(null_cv <= i_d_cv) + 1) / (n_perm + 1))

    return {
        "i_d_insample": round(i_d_in, 4),
        "i_d_cv": round(i_d_cv, 4),
        "i_d_null_mean": round(float(null_in.mean()), 4),
        "i_d_null_p025": round(float(np.percentile(null_in, 2.5)), 4),
        "i_d_vs_null": round(i_d_in - float(null_in.mean()), 4),
        "p_i_d_insample": round(p_in, 4),
        "i_d_cv_null_mean": round(float(null_cv.mean()), 4),
        "i_d_cv_null_p025": round(float(np.percentile(null_cv, 2.5)), 4),
        "p_i_d": round(p_cv, 4),
    }


def design_auc_permutation_p(D: np.ndarray, y: np.ndarray, observed: float,
                            n_perm: int = 200, seed: int | None = None,
                            fold_D=None) -> float:
    """Permutation p-value for the cross-fitted design-only AUC.

    Why this exists. V_D is an unpenalised linear R-squared. Cross-fitting it is the
    right thing to do, but it costs power exactly where the design matrix is wide:
    with 26 one-hot columns on 182 donors, cross-fitted OLS can explain nothing at
    all out of fold while a penalised logistic model on the SAME matrix reaches AUC
    0.82. p(V_D) then fails to reject in a comparison where the design plainly does
    predict the diagnosis.

    So the association is also tested with the statistic the paper actually leans
    on - the cross-fitted design-only AUC - against its own permutation null, using
    an identical pipeline on both sides. This is the primary test; p(V_D) is
    reported beside it and reads as a linear-variance summary.
    """
    if seed is None:
        seed = SEED
    if fold_D is None:
        fold_D = core.array_folds(D)
    rng = np.random.default_rng(seed + 991)
    hits = 0
    for _ in range(n_perm):
        yp = rng.permutation(y)
        if len(np.unique(yp)) < 2:
            continue
        # the null refits the encoding and the imputation as well as the
        # classifier, so observed and null are the same quantity throughout
        a = core.fold_contained_auc(fold_D, yp, core.SCREEN_METADATA_FROZEN, seed)
        if a >= observed:
            hits += 1
    return float((hits + 1) / (n_perm + 1))


def cross_fitted_auc(X: np.ndarray, y: np.ndarray, kind: str,
                     n_repeat: int = N_REPEAT,
                     fixed_C: float | None = None,
                     n_top_var: int | None = None,
                     n_pc: int | None = None) -> float:
    """Cross-fitted out-of-fold AUC.

    n_top_var and n_pc put the two unsupervised representation steps - the
    top-variance gene filter and the PCA - inside the training fold. They used to
    be fitted once on all donors before cross-fitting. Neither looks at the label,
    so neither could leak it, but both saw the held-out donors' expression, which
    makes the resulting AUC transductive rather than out-of-fold and makes
    "complete-pipeline permutation" an overstatement of what the null refits.

    fixed_C skips the inner C search. The permutation null uses it because a
    nested search costs ~135 logistic fits per permutation, which makes 1000
    permutations impractical. The observed statistic that the null is compared
    against is recomputed with the SAME frozen setting - comparing a tuned
    observed AUC against an untuned null would bias the p-value downwards.
    """
    aucs = []
    for rep in range(n_repeat):
        skf = StratifiedKFold(N_SPLITS, shuffle=True, random_state=SEED + rep)
        oof = np.zeros(len(y), dtype=float)
        for tr, te in skf.split(X, y):
            if kind == "rf":
                m = RandomForestClassifier(
                    n_estimators=500, min_samples_leaf=2, n_jobs=-1,
                    random_state=SEED + rep, class_weight="balanced")
                m.fit(X[tr], y[tr])
                oof[te] = m.predict_proba(X[te])[:, 1]
                continue

            Xtr, Xte = X[tr], X[te]
            if n_top_var is not None and n_top_var < Xtr.shape[1]:
                keep = np.argsort(Xtr.var(axis=0))[::-1][:n_top_var]
                Xtr, Xte = Xtr[:, keep], Xte[:, keep]
            if n_pc is not None:
                pre = StandardScaler().fit(Xtr)
                pca = PCA(n_components=min(n_pc, Xtr.shape[1], len(tr) - 2),
                          random_state=SEED).fit(pre.transform(Xtr))
                Xtr, Xte = (pca.transform(pre.transform(Xtr)),
                            pca.transform(pre.transform(Xte)))

            sc = StandardScaler().fit(Xtr)
            if fixed_C is not None:
                best = fixed_C
            else:
                best, best_auc = CS[0], -np.inf
                inner = StratifiedKFold(3, shuffle=True, random_state=SEED)
                for C in CS:
                    s_ = []
                    for itr, ite in inner.split(Xtr, y[tr]):
                        mm = LogisticRegression(C=C, solver="liblinear",
                                                class_weight="balanced",
                                                max_iter=5000)
                        mm.fit(sc.transform(Xtr[itr]), y[tr][itr])
                        s_.append(roc_auc_score(
                            y[tr][ite],
                            mm.predict_proba(sc.transform(Xtr[ite]))[:, 1]))
                    if np.mean(s_) > best_auc:
                        best_auc, best = float(np.mean(s_)), C
            m = LogisticRegression(C=best, solver="liblinear",
                                   class_weight="balanced", max_iter=5000)
            m.fit(sc.transform(Xtr), y[tr])
            oof[te] = m.predict_proba(sc.transform(Xte))[:, 1]
        aucs.append(roc_auc_score(y, oof))
    return float(np.mean(aucs))


def collection_strata(cov: pd.DataFrame) -> tuple[np.ndarray, str]:
    """Strata for the restricted permutation, and the variables defining them.

    IMPORTANT, and the reason this function returns its own definition: these
    strata are built from the recorded COLLECTION variables only - batch__* and,
    where it varies, assay. The declared design matrix also holds sample-quality
    and demographic columns, which are NOT conditioned on here. So permuting
    within these strata preserves the association between the diagnosis and the
    collection variables; it does not preserve the association between the
    diagnosis and the full design matrix.

    Everything downstream must be named for that. The test is a
    collection-preserving (restricted) permutation, not a "design-preserving"
    one, and a significant result means the classifier beat what these strata
    can explain - not that it beat everything recorded.

    Where a cohort records no collection variable at all there is nothing to
    condition on, and tertiles of log cells per donor stand in. That is a
    sample-quality variable, so the returned label says so and the test for that
    cohort is a sample-quality-preserving permutation instead.
    """
    cats = [c for c in cov.columns if c.startswith("batch__")]
    if "assay" in cov.columns and cov["assay"].nunique() > 1:
        cats.append("assay")
    if not cats:
        q = pd.qcut(cov["log_cells_per_donor"], 3, labels=False, duplicates="drop")
        return q.to_numpy(), "log_cells_per_donor tertiles (no collection variable recorded)"
    key = cov[cats].astype(str).agg("|".join, axis=1)
    return pd.factorize(key)[0], " x ".join(cats)


N_TOP_VAR = 4000    # top-variance genes, selected within the training fold
PERM_C = 1.0        # frozen regularisation for the permutation pipeline
PERM_PCS = 50       # frozen dimensionality for the permutation pipeline


def _frozen_auc_worker(X: np.ndarray, yy: np.ndarray,
                       n_top_var: int = N_TOP_VAR) -> float:
    """Module-level twin of the frozen_auc closure, so loky can pickle it."""
    if len(np.unique(yy)) < 2:
        return 0.5
    return cross_fitted_auc(X, yy, "linear", n_repeat=1, fixed_C=PERM_C,
                            n_top_var=n_top_var, n_pc=PERM_PCS)


def permutation_pvalues(X: np.ndarray, y: np.ndarray, strata: np.ndarray,
                        n_perm: int = N_PERM) -> dict:
    """Free and collection-preserving permutation p-values.

    Both nulls and the observed statistic they are compared against run the
    identical pipeline, refitted from the raw matrix on every permutation: the
    top-variance gene filter, the standardiser and the PCA to PERM_PCS components
    are all fitted inside the training fold, then logistic regression at PERM_C,
    single cross-fitting repeat.

    Fitting the filter and the PCA once on all donors, as an earlier version did,
    made both the observed statistic and the null transductive. Neither step sees
    the label, so the comparison stayed like-for-like and the p-value was still a
    valid test of that fixed representation - but it was not the complete-pipeline
    test the manuscript described, because the representation was never refitted.
    """
    n_pc = int(min(PERM_PCS, X.shape[1], len(y) - 2))

    def frozen_auc(yy: np.ndarray) -> float:
        if len(np.unique(yy)) < 2:
            return 0.5
        return cross_fitted_auc(X, yy, "linear", n_repeat=1, fixed_C=PERM_C,
                                n_top_var=N_TOP_VAR, n_pc=PERM_PCS)

    # Draw every permuted label first, in exactly the order the sequential
    # version drew them: one free permutation, then one draw per usable stratum,
    # repeated n_perm times. frozen_auc consumes no randomness of its own (each
    # fold split is keyed to SEED), so pulling the draws out of the loop leaves
    # the RNG stream untouched and the result is identical to WORKERS=1. Only
    # then is the expensive part - refitting the whole representation once per
    # permutation - handed to a process pool.
    rng = np.random.default_rng(SEED)
    free, restricted = [], []
    for _ in range(n_perm):
        free.append(rng.permutation(y))

        yp, moved = y.copy(), False
        for s_ in np.unique(strata):
            idx = np.flatnonzero(strata == s_)
            if len(idx) > 1 and len(np.unique(y[idx])) > 1:
                yp[idx] = rng.permutation(y[idx])
                moved = True
        if moved:
            restricted.append(yp)
    usable = len(restricted)

    todo = [y] + free + restricted
    if WORKERS > 1:
        aucs = Parallel(n_jobs=WORKERS, backend="loky", batch_size=4)(
            delayed(_frozen_auc_worker)(X, yy, n_top_var=N_TOP_VAR) for yy in todo)
    else:
        aucs = [frozen_auc(yy) for yy in todo]

    observed = aucs[0]
    std_hits = sum(a >= observed for a in aucs[1:1 + n_perm])
    dpn_hits = sum(a >= observed for a in aucs[1 + n_perm:])

    return {
        "observed_frozen_auc": float(observed),
        "p_standard": float((std_hits + 1) / (n_perm + 1)),
        "p_collection_preserving": (float((dpn_hits + 1) / (usable + 1))
                                if usable else float("nan")),
        "n_usable_stratified_perm": int(usable),
        "perm_n_pc": n_pc,
    }


def screen(cohort: str, n_perm: int) -> list[dict]:
    cov_path = INPUTS / "design_metadata" / f"{cohort.lower()}_donor_covariates.tsv"
    pb_path = INPUTS / "donor_level" / cohort / "donor_log1p_cpm.parquet"
    if not cov_path.exists() or not pb_path.exists():
        print(f"[{cohort}] SKIP - build_donor_level.py has not been run")
        return []

    cov = pd.read_csv(cov_path, sep="\t", dtype={"donor_id": str}).set_index("donor_id")
    pb = pd.read_parquet(pb_path)
    pb.index = pb.index.astype(str)
    donors = [d for d in cov.index if d in pb.index]
    cov, pb = cov.loc[donors], pb.loc[donors]
    y = cov["y_true"].to_numpy(dtype=int)

    # The top-variance filter keeps the donor-level problem well conditioned. It
    # is selected inside each training fold: ranking genes on all donors first is
    # unsupervised, but it still uses the held-out donors' expression.
    Xp = pb.to_numpy(dtype=np.float64)

    print(f"[{cohort}] donors={len(y)} cases={int(y.sum())} controls={int((y==0).sum())}")
    disease_auc = cross_fitted_auc(Xp, y, "linear", n_top_var=N_TOP_VAR)
    print(f"  disease AUC (pseudobulk) = {disease_auc:.3f}")

    strata, strata_def = collection_strata(cov)
    perm = permutation_pvalues(Xp, y, strata, n_perm)
    p_std, p_dpn = perm["p_standard"], perm["p_collection_preserving"]
    print(f"  strata: {strata_def}\n"
          f"  frozen-pipeline AUC={perm['observed_frozen_auc']:.3f} "
          f"(PCA {perm['perm_n_pc']} comps, C={PERM_C})")
    print(f"  p_standard={p_std:.4f}  p_collection_preserving={p_dpn:.4f} "
          f"(usable stratified permutations: {perm['n_usable_stratified_perm']})")

    # Degeneracy gate. A contrast where the design perfectly separates the label
    # is the complete-collinearity case: the design and disease contributions are
    # not separately identifiable, the collection-preserving null cannot be formed
    # (no stratum holds both classes), and a donor-level classifier on thousands
    # of features will separate any labelling. Such a row is not a measurement of
    # confounding strength and must not be ranked alongside the others.
    # NOT ESTIMABLE: the quantity itself has no meaning for this contrast.
    # NOT ESTIMABLE is now reserved for genuine non-identifiability: no stratum
    # holds both labels, so the collection-preserving null cannot be formed at all.
    #
    # The previous version also declared a contrast not estimable when the tuned
    # and frozen classifiers differed by more than 0.20 AUC. That rule is gone.
    # It compared two different algorithms through an arbitrary cut-point, and a
    # tuned model legitimately beating a fixed one is not the same kind of problem
    # as a design that perfectly separates the label. The frozen pipeline is now
    # the inferential statistic everywhere - it is what both permutation nulls are
    # built from - and the tuned AUC is reported beside it as a description. Their
    # gap is released as a column so a reader can see pipeline sensitivity
    # directly instead of through a threshold.
    # ONE criterion, not two. The previous version also declared a contrast not
    # estimable when the metadata-only AUC reached 1.000. Pre-submission review
    # was right that this does not follow: a finite-sample AUC of 1.000 says the
    # fitted model separated these donors completely, which with 21 donors and
    # seven one-hot columns can happen without the metadata determining the label
    # in the population. It is corroborating structure, not a proof of
    # non-identifiability, and it is reported as the metadata-only AUC it is.
    degenerate: list[str] = []
    if np.isnan(p_dpn):
        degenerate.append("no collection stratum holds both labels, so the "
                          "stratified permutation has no non-trivial reference "
                          "distribution and no conditioned estimate can be formed")

    # UNDERPOWERED: the quantity is estimable but the estimate is unstable.
    # Kept in the spectrum, flagged, and to be reported with a seed range rather
    # than a single p-value.
    underpowered: list[str] = []
    minority = min(int(y.sum()), int((y == 0).sum()))
    if len(y) < 40:
        underpowered.append(f"n={len(y)} below 40")
    if minority < 15:
        underpowered.append(f"minority class {minority} below 15")

    rows = []
    for block, D in design_matrix(cov).items():
        # Every cross-fitted metadata quantity is built with fold_D, which takes
        # the medians and the categorical level set from the training donors of
        # that split. The whole-cohort D is still passed in, because the
        # in-sample V_D and the reported column count are defined on it.
        fold_D = core.design_folds(cov, block)
        idm = label_information_fraction(y.astype(float), D, fold_D=fold_D)
        a_lin = core.fold_contained_auc(fold_D, y, core.SCREEN_METADATA_TUNED, SEED)
        a_rf = core.fold_contained_forest_auc(fold_D, y, core.SCREEN_METADATA_FOREST, SEED)
        # observed and null share one frozen pipeline, as everywhere else
        a_lin_frozen = core.fold_contained_auc(
            fold_D, y, core.SCREEN_METADATA_FROZEN, SEED)
        p_auc = design_auc_permutation_p(D, y, a_lin_frozen, fold_D=fold_D)
        flag = "" if idm["i_d_cv"] < idm["i_d_cv_null_p025"] else "  <-- NOT below its own null"
        print(f"  [{block:<12}] V_D_cv={idm['i_d_cv']:.3f} "
              f"(null={idm['i_d_cv_null_mean']:.3f}, p={idm['p_i_d']:.3f}; "
              f"in-sample {idm['i_d_insample']:.3f}, p={idm['p_i_d_insample']:.3f})  "
              f"design AUC lin={a_lin:.3f} rf={a_rf:.3f} "
              f"(frozen {a_lin_frozen:.3f}, p={p_auc:.4f}){flag}")
        rows.append({
            "cohort": cohort, "block": block, "n_donor": len(y),
            "n_case": int(y.sum()), "n_design_feature": int(D.shape[1]),
            "I_D": idm["i_d_insample"],
            "I_D_cv": idm["i_d_cv"],
            "I_D_null_mean": idm["i_d_null_mean"],
            "I_D_vs_null": idm["i_d_vs_null"],
            "p_I_D_insample": idm["p_i_d_insample"],
            "I_D_cv_null_mean": idm["i_d_cv_null_mean"],
            "I_D_cv_null_p025": idm["i_d_cv_null_p025"],
            "p_I_D": idm["p_i_d"],
            "design_auc_linear": round(a_lin, 4),
            "design_auc_rf": round(a_rf, 4),
            "design_auc_frozen": round(a_lin_frozen, 4),
            "p_design_auc": round(p_auc, 4),
            "disease_auc": round(disease_auc, 4),
            "observed_frozen_auc": round(perm["observed_frozen_auc"], 4),
            "frozen_minus_tuned_auc": round(perm["observed_frozen_auc"] - disease_auc, 4),
            "p_standard": round(p_std, 5),
            "p_collection_preserving": None if np.isnan(p_dpn) else round(p_dpn, 5),
            "n_strata": int(len(np.unique(strata))),
            "strata_definition": strata_def,
            "degenerate": bool(degenerate),
            "degenerate_reason": "; ".join(degenerate) if degenerate else "",
            "underpowered": bool(underpowered),
            "underpowered_reason": "; ".join(underpowered) if underpowered else "",
        })
    return rows


DESIGN_COLUMNS = [
    "n_design_feature", "I_D", "I_D_cv", "I_D_null_mean", "I_D_vs_null",
    "p_I_D_insample", "I_D_cv_null_mean", "I_D_cv_null_p025", "p_I_D",
    "design_auc_linear", "design_auc_rf", "design_auc_frozen", "p_design_auc",
]


def refresh_design(cohort: str, prior: pd.DataFrame) -> list[dict]:
    """Recompute the metadata-only columns of an existing run, and nothing else.

    Fitting the imputation and the one-hot encoding inside the training fold
    changes only quantities built from the metadata matrix. `disease_auc`, the
    frozen expression AUC and both permutation p-values are computed from
    expression and the collection strata, and touch no part of `design_matrix`;
    re-running them would burn a thousand permutations per comparison to
    reproduce numbers that cannot move. So they are carried across from `prior`
    unchanged, which also means the released expression results stay bit-identical
    and any difference in the refreshed file is attributable to this one change.

    The degeneracy gate IS recomputed, because one of its two triggers reads the
    metadata-only AUC.
    """
    cov_path = INPUTS / "design_metadata" / f"{cohort.lower()}_donor_covariates.tsv"
    cov = pd.read_csv(cov_path, sep="\t", dtype={"donor_id": str}).set_index("donor_id")
    pb_path = INPUTS / "donor_level" / cohort / "donor_log1p_cpm.parquet"
    donors = pd.read_parquet(pb_path, columns=[]).index.astype(str)
    cov = cov.loc[[d for d in cov.index if d in set(donors)]]
    y = cov["y_true"].to_numpy(dtype=int)

    keep = prior[prior.cohort == cohort]
    if keep.empty:
        raise SystemExit(f"{cohort}: nothing to refresh in the prior run")
    carried = keep.iloc[0]

    reason = carried.get("degenerate_reason", "")
    reason = "" if reason is None or pd.isna(reason) else str(reason)
    # carried across verbatim: the sole criterion is the collection-stratum
    # structure, which this refresh does not touch
    degenerate = [r for r in reason.split("; ") if r]

    rows = []
    for block, D in design_matrix(cov).items():
        fold_D = core.design_folds(cov, block)
        idm = label_information_fraction(y.astype(float), D, fold_D=fold_D)
        a_lin = core.fold_contained_auc(fold_D, y, core.SCREEN_METADATA_TUNED, SEED)
        a_rf = core.fold_contained_forest_auc(fold_D, y, core.SCREEN_METADATA_FOREST, SEED)
        a_lin_frozen = core.fold_contained_auc(
            fold_D, y, core.SCREEN_METADATA_FROZEN, SEED)
        p_auc = design_auc_permutation_p(D, y, a_lin_frozen, fold_D=fold_D)
        print(f"  [{block:<12}] design AUC lin={a_lin:.4f} rf={a_rf:.4f} "
              f"frozen={a_lin_frozen:.4f} p={p_auc:.4f}  V_D_cv={idm['i_d_cv']:.3f}")
        prev = keep[keep.block == block]
        row = dict(prev.iloc[0]) if len(prev) else dict(carried)
        row.update({
            "cohort": cohort, "block": block, "n_donor": len(y),
            "n_case": int(y.sum()), "n_design_feature": int(D.shape[1]),
            "I_D": idm["i_d_insample"], "I_D_cv": idm["i_d_cv"],
            "I_D_null_mean": idm["i_d_null_mean"], "I_D_vs_null": idm["i_d_vs_null"],
            "p_I_D_insample": idm["p_i_d_insample"],
            "I_D_cv_null_mean": idm["i_d_cv_null_mean"],
            "I_D_cv_null_p025": idm["i_d_cv_null_p025"], "p_I_D": idm["p_i_d"],
            "design_auc_linear": round(a_lin, 4), "design_auc_rf": round(a_rf, 4),
            "design_auc_frozen": round(a_lin_frozen, 4),
            "p_design_auc": round(p_auc, 4),
        })
        rows.append(row)
    for r in rows:
        r["degenerate"] = bool(degenerate)
        r["degenerate_reason"] = "; ".join(degenerate) if degenerate else ""
    return rows


def write_outputs(df: pd.DataFrame, out: Path, cohorts: list[str],
                  n_perm: int) -> None:
    """Spectrum, degeneracy gate and config. Shared by a direct run and by
    --merge-from, so a merged sweep cannot drift from a single-process one."""
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "design_screen.tsv", sep="\t", index=False)

    # Sort by cross-fitted design AUC, not by raw I_D: the raw value is not
    # comparable across cohorts with different design-feature counts.
    allb = df[df.block == "all"]
    if "degenerate" in allb.columns and allb["degenerate"].any():
        bad = allb[allb["degenerate"]]
        bad.to_csv(out / "degenerate_contrasts.tsv", sep="\t", index=False)
        print("\n!! excluded from the spectrum as NOT ESTIMABLE:")
        for _, r in bad.iterrows():
            print(f"   {r['cohort']}: {r['degenerate_reason']}")
        allb = allb[~allb["degenerate"]]
    if "underpowered" in allb.columns and allb["underpowered"].any():
        print("\n!! kept but UNDERPOWERED - report a seed range, not a single p:")
        for _, r in allb[allb["underpowered"]].iterrows():
            print(f"   {r['cohort']}: {r['underpowered_reason']}")

    spectrum = (allb
                .sort_values("design_auc_linear")
                [["cohort", "n_donor", "n_design_feature",
                  "design_auc_linear", "design_auc_rf",
                  "I_D", "I_D_cv", "I_D_null_mean", "p_I_D_insample",
                  "I_D_cv_null_mean", "p_I_D", "design_auc_frozen", "p_design_auc",
                  "disease_auc", "p_standard", "p_collection_preserving"]])
    spectrum.to_csv(out / "confounding_spectrum.tsv", sep="\t", index=False)
    (out / "screen_config.json").write_text(json.dumps(
        {"seed": SEED, "n_splits": N_SPLITS, "n_repeat": N_REPEAT,
         "n_perm": n_perm, "cohorts": cohorts,
         "i_d_estimator": "in-sample + cross-fitted + permutation null",
         "spectrum_sorted_by": "design_auc_linear"}, indent=2))

    print("\n=== design-confounding spectrum (block = all) ===")
    print(spectrum.to_string(index=False))
    print(f"\nwrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--seed", type=int, default=None,
                    help="override the global seed (for seed-stability checks)")
    ap.add_argument("--out", default=None,
                    help="override the output directory")
    ap.add_argument("--workers", type=int, default=1,
                    help="processes for the permutation stage; results are "
                         "identical to --workers 1, only the wall time differs")
    ap.add_argument("--refresh-design", default=None, metavar="design_screen.tsv",
                    help="recompute only the metadata-only columns of an existing "
                         "run, carrying the expression columns across unchanged")
    ap.add_argument("--merge-from", nargs="*", default=None,
                    help="concatenate per-cohort design_screen.tsv fragments and "
                         "write the spectrum, instead of running the screen")
    args = ap.parse_args()

    global SEED, OUT, WORKERS
    WORKERS = args.workers
    if args.seed is not None:
        SEED = args.seed
    if args.out:
        OUT = Path(args.out)
    elif args.seed is not None:
        # A seed sweep with no --out wrote every seed to the same
        # design_screen.tsv, so four of five runs were silently destroyed and the
        # seed ranges could not be reproduced from the documented commands. Give
        # each seed its own directory unless the caller names one.
        OUT = OUT.parent / "screen" / f"seed_{SEED}"

    if args.refresh_design:
        prior = pd.read_csv(args.refresh_design, sep="\t")
        cohorts = args.cohorts or sorted(prior.cohort.unique())
        rows: list[dict] = []
        for c in cohorts:
            print(f"[{c}] refreshing metadata-only columns")
            rows.extend(refresh_design(c, prior))
        df = pd.DataFrame(rows)[list(prior.columns)]
        write_outputs(df, OUT, cohorts, args.n_perm)
        return

    if args.merge_from:
        frames = [pd.read_csv(f, sep="\t") for f in sorted(args.merge_from)]
        df = pd.concat(frames, ignore_index=True)
        order = sorted(df.cohort.unique())
        df = (df.assign(_o=df.cohort.map({c: k for k, c in enumerate(order)}))
                .sort_values(["_o"], kind="stable").drop(columns="_o")
                .reset_index(drop=True))
        write_outputs(df, OUT, order, args.n_perm)
        return

    if args.all:
        cohorts = sorted(p.name for p in (INPUTS / "donor_level").iterdir()
                         if p.is_dir())
    elif args.cohorts:
        cohorts = args.cohorts
    else:
        raise SystemExit("pass --cohorts A B or --all")

    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for c in cohorts:
        rows.extend(screen(c, args.n_perm))

    if not rows:
        raise SystemExit("no cohort produced results")

    write_outputs(pd.DataFrame(rows), OUT, cohorts, args.n_perm)


if __name__ == "__main__":
    main()
