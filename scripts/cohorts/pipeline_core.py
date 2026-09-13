"""One definition of every fitted pipeline in the paper.

Four analyses fit classifiers: the multi-cohort screen, the calibration
simulation, the decision-tree walkthrough and the lupus deep dive. Until now
each carried its own copy of "standardise, reduce, fit", and the copies had
drifted - different solvers, different inner-fold counts, and, worse, different
answers to *where* the unsupervised steps were fitted. Two consequences were
found in pre-submission review:

  * `sim/run_calibration.py` arm B compared an unadjusted model whose PCA was
    fitted on all donors against a residualised model whose PCA was fitted
    inside the training fold. The difference between the two arms therefore
    contained the change in representation as well as the residualisation, and
    the whole gap could not be read as the cost of adjustment.
  * The screen's design matrix was median-imputed and one-hot encoded on all
    donors before cross-fitting, while Methods said both happened inside the
    training fold.

So the specifications live here, once, as data; every runner imports them and
`tools/build_hyperparameter_table.py` writes Table S7 straight out of the same
objects. A parameter cannot be described in the manuscript as something other
than what ran, because the description is generated from the object that ran.

Nothing here looks at a held-out donor. `fold_contained_auc` fits imputation,
encoding, feature selection, scaling, PCA, residualisation and the classifier on
training donors only, then applies them.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Callable, Sequence

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

C_GRID = np.logspace(-4, 4, 9)


# --------------------------------------------------------------- specifications
@dataclass(frozen=True)
class PipelineSpec:
    """Everything that decides a number, in one object.

    `fixed_C` set means no inner search: the frozen statistic that both
    permutation nulls are built from. `fixed_C` None means the inner search over
    `C_GRID` runs inside each outer training fold.
    """
    name: str
    purpose: str
    n_splits: int = 5
    n_repeat: int = 1
    fixed_C: float | None = 1.0
    inner_folds: int = 3
    solver: str = "liblinear"
    class_weight: str | None = "balanced"
    max_iter: int = 5000
    n_top_var: int | None = None       # top-variance feature filter, in-fold
    n_pc: int | None = None            # PCA components, in-fold
    score: str = "predict_proba"

    def describe(self) -> dict[str, str]:
        c = "1.0 (frozen)" if self.fixed_C is not None else \
            f"10^-4..10^4 in decade steps, {self.inner_folds}-fold inner search"
        return {
            "pipeline": self.name,
            "used for": self.purpose,
            "outer splits": f"{self.n_splits}-fold stratified x {self.n_repeat} repeat"
                            f"{'s' if self.n_repeat > 1 else ''}",
            "classifier": f"logistic regression, {self.solver}, "
                          f"class_weight={self.class_weight}, max_iter={self.max_iter:,}",
            "regularisation": c,
            "feature filter": "none" if self.n_top_var is None
                              else f"top {self.n_top_var:,} by training-fold variance",
            "dimension reduction": "none" if self.n_pc is None
                                   else f"PCA to {self.n_pc} components, fitted in fold",
            "score": self.score,
        }


@dataclass(frozen=True)
class ForestSpec:
    name: str
    purpose: str
    n_estimators: int
    min_samples_leaf: int
    max_depth: int | None
    n_splits: int = 5
    n_repeat: int = 1
    class_weight: str | None = "balanced"

    def describe(self) -> dict[str, str]:
        return {
            "pipeline": self.name,
            "used for": self.purpose,
            "outer splits": f"{self.n_splits}-fold stratified x {self.n_repeat} repeat"
                            f"{'s' if self.n_repeat > 1 else ''}",
            "classifier": f"random forest, {self.n_estimators} trees, "
                          f"min_samples_leaf={self.min_samples_leaf}, "
                          f"max_depth={self.max_depth if self.max_depth else 'unlimited'}, "
                          f"class_weight={self.class_weight}",
            "regularisation": "fixed, not tuned",
            "feature filter": "none",
            "dimension reduction": "none",
            "score": "predict_proba",
        }


# The screen's linear classifier, in its two roles. `SCREEN_FROZEN` is the
# inferential statistic: both permutation nulls refit exactly this. `SCREEN_TUNED`
# is the description reported beside it.
SCREEN_FROZEN = PipelineSpec(
    name="screen_frozen",
    purpose="the statistic both permutation nulls are built from, on expression "
            "and on the metadata matrix",
    n_repeat=1, fixed_C=1.0, n_top_var=4000, n_pc=50)

SCREEN_TUNED = replace(
    SCREEN_FROZEN, name="screen_tuned",
    purpose="descriptive metadata-only and expression AUC reported beside the "
            "frozen statistic",
    n_repeat=5, fixed_C=None)

# The metadata matrix is small and dense; no filter and no PCA are applied to it.
SCREEN_METADATA_FROZEN = replace(
    SCREEN_FROZEN, name="screen_metadata_frozen",
    purpose="metadata-only AUC and its 200-draw permutation null",
    n_top_var=None, n_pc=None)

SCREEN_METADATA_TUNED = replace(
    SCREEN_TUNED, name="screen_metadata_tuned",
    purpose="descriptive metadata-only AUC (Table 1, `AUC` column source for the "
            "tuned variant reported in Table S4)",
    n_top_var=None, n_pc=None)

SCREEN_METADATA_FOREST = ForestSpec(
    name="screen_metadata_forest",
    purpose="nonlinear metadata-only arm, to show a linear model is not missing "
            "a low-dimensional metadata-label relationship",
    n_estimators=500, min_samples_leaf=2, max_depth=None, n_repeat=5)

# The lupus deep dive predates the screen and runs a wider repeat budget on a
# single dataset. Its settings are NOT the screen's and are no longer described
# as though they were.
SLE_DEEP = PipelineSpec(
    name="sle_deep_dive",
    purpose="GSE174188 residualisation, restriction and cross-cohort arms "
            "(Sections 5.1-5.2)",
    n_repeat=20, fixed_C=None, inner_folds=5, max_iter=20000)

SLE_DEEP_FOREST = ForestSpec(
    name="sle_deep_dive_forest",
    purpose="nonlinear residualiser and nonlinear metadata arm in GSE174188",
    n_estimators=256, min_samples_leaf=3, max_depth=4, n_repeat=20)

# The walkthrough and the calibration simulation are the screen's frozen pipeline
# with a different repeat budget, and nothing else. They used to differ in the
# solver, in class weighting, in whether a variance filter ran at all, and - the
# defect pre-submission review found - in whether the PCA saw the held-out
# donors. Deriving them from SCREEN_FROZEN with `replace` is what makes
# "the same frozen pipeline as the empirical screen" a checkable statement
# rather than a claim in a docstring.
WALKTHROUGH = replace(
    SCREEN_FROZEN, name="walkthrough",
    purpose="decision-tree walkthrough on CMV, Ren and influenza (Section 7); "
            "20 repeats because a single split is too noisy to compare an "
            "adjusted arm against an unadjusted one",
    n_repeat=20)

SIMULATION = replace(
    SCREEN_FROZEN, name="simulation",
    purpose="calibration and extended simulation arms (Section 3)")

# The extended simulation is deliberately NOT derived from SCREEN_FROZEN. Its
# cohorts are 200 donors by 100 features, so the screen's 4,000-gene filter and
# its PCA to 50 components would be a different reduction ratio, not the same
# step. It is listed here so Table S7 describes it too, and because its
# unadjusted, residualised and restricted arms already differ in one operation
# each - the property arm B of the calibration study was missing.
EXTENDED_SIM = PipelineSpec(
    name="extended_simulation",
    purpose="19,600-cohort regime sweep (Sections 3.1-3.3, Figure 2); a check on "
            "estimator behaviour, not a mimic of the empirical pipeline",
    n_repeat=1, fixed_C=1.0, solver="lbfgs", class_weight=None, max_iter=2000,
    n_top_var=None, n_pc=None, score="decision_function")

ALL_SPECS: tuple = (
    SCREEN_FROZEN, SCREEN_TUNED, SCREEN_METADATA_FROZEN, SCREEN_METADATA_TUNED,
    SCREEN_METADATA_FOREST, SLE_DEEP, SLE_DEEP_FOREST, WALKTHROUGH, SIMULATION,
    EXTENDED_SIM)


# ------------------------------------------------------- fold-contained design
QC_COLS = ["log_cells_per_donor", "log_mean_umi_per_cell",
           "log_mean_genes_per_cell", "aggregate_pct_mito"]
DEMO_NUM = ["age_years"]
DEMO_CAT = ["sex", "ethnicity"]


def block_columns(cov: pd.DataFrame) -> dict[str, tuple[list[str], list[str]]]:
    """(numeric, categorical) source columns for each metadata block.

    One place decides what a block contains, so the global matrix used for the
    in-sample statistics and the fold matrices used for everything cross-fitted
    can never be built from different column sets.
    """
    qc = [c for c in QC_COLS if c in cov.columns]
    dnum = [c for c in DEMO_NUM if c in cov.columns]
    dcat = [c for c in DEMO_CAT if c in cov.columns]
    bcat = [c for c in cov.columns if c.startswith("batch__")]
    if "assay" in cov.columns and cov["assay"].nunique() > 1:
        bcat = bcat + ["assay"]
    return {
        "qc": (qc, []),
        "demographic": (dnum, dcat),
        # named "batch" for continuity with the released result files;
        # tools/stage_screen_results.py renames it to "collection" on the way
        # into the supplementary tables.
        "batch": ([], bcat),
        "all": (qc + dnum, dcat + bcat),
    }


def _encode(cov: pd.DataFrame, num: Sequence[str], cat: Sequence[str],
            fit_rows: np.ndarray | None) -> np.ndarray | None:
    """One-hot and impute, with every parameter taken from `fit_rows`.

    `fit_rows` None reproduces the whole-cohort construction, which is still what
    the in-sample V_D is defined on. Otherwise the medians and the categorical
    level sets come from those rows only; a level seen for the first time in a
    held-out donor is encoded as all zeros, that is, as the reference level.
    """
    parts: list[np.ndarray] = []
    ref = cov if fit_rows is None else cov.iloc[fit_rows]

    if num:
        X = np.array(cov[list(num)].to_numpy(dtype=float), copy=True)
        med = np.nanmedian(np.asarray(ref[list(num)].to_numpy(dtype=float)), axis=0)
        keep = []
        for j in range(X.shape[1]):
            if np.isnan(X[:, j]).all() or not np.isfinite(med[j]):
                continue
            X[np.isnan(X[:, j]), j] = med[j]
            keep.append(j)
        if keep:
            parts.append(X[:, keep])

    for c in cat:
        s = cov[c].astype(str)
        levels = sorted(ref[c].astype(str).unique())
        if len(levels) < 2:
            continue
        # drop_first, so the encoding matches pandas.get_dummies(drop_first=True)
        # when fit_rows covers every donor and every level.
        parts.append(np.column_stack([(s == lv).to_numpy(float)
                                      for lv in levels[1:]]))

    if not parts:
        return None
    return np.hstack(parts)


def design_matrix(cov: pd.DataFrame) -> dict[str, np.ndarray]:
    """Whole-cohort metadata blocks. Used for the in-sample V_D and for reporting
    how wide each declared matrix is; every cross-fitted number uses `fold_design`.
    """
    out = {}
    for name, (num, cat) in block_columns(cov).items():
        D = _encode(cov, num, cat, None)
        if D is not None and D.shape[1] > 0:
            out[name] = D
    return out


def fold_design(cov: pd.DataFrame, block: str,
                tr: np.ndarray, te: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Training and held-out metadata matrices, both built from `tr` alone."""
    num, cat = block_columns(cov)[block]
    D = _encode(cov, num, cat, tr)
    if D is None:
        return np.zeros((len(tr), 0)), np.zeros((len(te), 0))
    return D[tr], D[te]


# ------------------------------------------------------------------- estimators
FoldFn = Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]


def _fit_predict(Xtr, Xte, ytr, spec: PipelineSpec, inner_seed: int) -> np.ndarray:
    sc = StandardScaler().fit(Xtr)
    Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
    C = spec.fixed_C
    if C is None:
        best, best_auc = C_GRID[0], -np.inf
        inner = StratifiedKFold(spec.inner_folds, shuffle=True, random_state=inner_seed)
        for cand in C_GRID:
            s = []
            for itr, ite in inner.split(Xtr, ytr):
                m = LogisticRegression(C=cand, solver=spec.solver,
                                       class_weight=spec.class_weight,
                                       max_iter=spec.max_iter)
                m.fit(Xtr[itr], ytr[itr])
                s.append(roc_auc_score(ytr[ite],
                                       m.predict_proba(Xtr[ite])[:, 1]))
            if np.mean(s) > best_auc:
                best_auc, best = float(np.mean(s)), cand
        C = best
    m = LogisticRegression(C=C, solver=spec.solver, class_weight=spec.class_weight,
                           max_iter=spec.max_iter).fit(Xtr, ytr)
    return (m.predict_proba(Xte)[:, 1] if spec.score == "predict_proba"
            else m.decision_function(Xte))


def fold_contained_auc(get_fold: FoldFn, y: np.ndarray, spec: PipelineSpec,
                       seed: int, residualise: FoldFn | None = None) -> float:
    """Cross-fitted AUC where every fitted object comes from the training donors.

    `get_fold(tr, te)` returns the two feature matrices for one split; when the
    features are a fixed array this is a slice, and when they are a metadata
    matrix that has to be imputed and encoded it is `fold_design`. `residualise`
    supplies the adjustment matrix the same way, so the unadjusted and adjusted
    arms of a comparison differ in that argument and in nothing else - the point
    the calibration simulation previously got wrong.
    """
    if len(np.unique(y)) < 2:
        return float("nan")
    aucs = []
    for rep in range(spec.n_repeat):
        skf = StratifiedKFold(spec.n_splits, shuffle=True, random_state=seed + rep)
        oof = np.zeros(len(y), dtype=float)
        for tr, te in skf.split(np.zeros(len(y)), y):
            Xtr, Xte = get_fold(tr, te)
            if Xtr.shape[1] == 0:
                oof[te] = 0.0
                continue
            if residualise is not None:
                Dtr, Dte = residualise(tr, te)
                if Dtr.shape[1] > 0:
                    ds = StandardScaler().fit(Dtr)
                    rg = Ridge(alpha=1.0).fit(ds.transform(Dtr), Xtr)
                    Xtr = Xtr - rg.predict(ds.transform(Dtr))
                    Xte = Xte - rg.predict(ds.transform(Dte))
            if spec.n_top_var is not None and spec.n_top_var < Xtr.shape[1]:
                keep = np.argsort(Xtr.var(axis=0))[::-1][:spec.n_top_var]
                Xtr, Xte = Xtr[:, keep], Xte[:, keep]
            if spec.n_pc is not None:
                pre = StandardScaler().fit(Xtr)
                k = min(spec.n_pc, Xtr.shape[1], len(tr) - 2)
                pca = PCA(n_components=k, random_state=seed).fit(pre.transform(Xtr))
                Xtr = pca.transform(pre.transform(Xtr))
                Xte = pca.transform(pre.transform(Xte))
            oof[te] = _fit_predict(Xtr, Xte, y[tr], spec, seed)
        aucs.append(roc_auc_score(y, oof))
    return float(np.mean(aucs))


def fold_contained_forest_auc(get_fold: FoldFn, y: np.ndarray, spec: ForestSpec,
                              seed: int) -> float:
    if len(np.unique(y)) < 2:
        return float("nan")
    aucs = []
    for rep in range(spec.n_repeat):
        skf = StratifiedKFold(spec.n_splits, shuffle=True, random_state=seed + rep)
        oof = np.zeros(len(y), dtype=float)
        for tr, te in skf.split(np.zeros(len(y)), y):
            Xtr, Xte = get_fold(tr, te)
            if Xtr.shape[1] == 0:
                oof[te] = 0.0
                continue
            m = RandomForestClassifier(
                n_estimators=spec.n_estimators, min_samples_leaf=spec.min_samples_leaf,
                max_depth=spec.max_depth, class_weight=spec.class_weight,
                n_jobs=-1, random_state=seed + rep).fit(Xtr, y[tr])
            oof[te] = m.predict_proba(Xte)[:, 1]
        aucs.append(roc_auc_score(y, oof))
    return float(np.mean(aucs))


def array_folds(X: np.ndarray) -> FoldFn:
    """`get_fold` for a feature matrix that needs no per-fold construction."""
    return lambda tr, te: (X[tr], X[te])


def design_folds(cov: pd.DataFrame, block: str) -> FoldFn:
    """`get_fold` for a metadata block, imputed and encoded inside the fold."""
    return lambda tr, te: fold_design(cov, block, tr, te)
