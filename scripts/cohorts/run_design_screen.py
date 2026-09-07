"""Multi-cohort design screen, including the design-preserving null.

Outputs one row per (cohort, design block) with

    I_D              = 1 - R^2(Y ~ D)      linear label information left after D
    design_auc_lin   cross-fitted logistic  Y ~ D
    design_auc_rf    cross-fitted random forest Y ~ D   (nonlinear design signal)
    disease_auc      cross-fitted logistic  Y ~ pseudobulk
    p_standard       standard label permutation p-value
    p_design_preserving   permutation WITHIN design strata

Why two nulls
-------------
PLOS reviewer #4 observed that a standard complete-pipeline label permutation
destroys the design-label association along with the biological signal, so it
tests leakage and pipeline validity rather than "can the classifier exploit
acquisition structure". The design-preserving null permutes the label only
within strata that share a design configuration: the D-Y association is retained
while disease-specific molecular information is destroyed. A classifier that
still scores above this null is using something beyond the recorded design.

    p_standard small, p_design_preserving small  -> signal beyond design
    p_standard small, p_design_preserving LARGE  -> apparent signal is design

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
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[2]
INPUTS = REPO / "inputs"
OUT = REPO / "results" / "multi_cohort_design_screen"

SEED = 20260907
N_SPLITS = 5
N_REPEAT = 5
N_PERM = 1000
CS = np.logspace(-4, 4, 9)


def design_matrix(cov: pd.DataFrame) -> dict[str, np.ndarray]:
    """Assemble design blocks. Mirrors the v2 GSE174188 block definition."""
    qc_cols = [c for c in ["log_cells_per_donor", "log_mean_umi_per_cell",
                           "log_mean_genes_per_cell", "aggregate_pct_mito"]
               if c in cov.columns]
    demo_num = [c for c in ["age_years"] if c in cov.columns]
    demo_cat = [c for c in ["sex", "ethnicity"] if c in cov.columns]
    batch_cat = [c for c in cov.columns if c.startswith("batch__")]
    if "assay" in cov.columns and cov["assay"].nunique() > 1:
        batch_cat = batch_cat + ["assay"]

    def build(num: list[str], cat: list[str]) -> np.ndarray | None:
        parts = []
        if num:
            # .copy(): pandas >=3.0 can return a read-only view here
            X = np.array(cov[num].to_numpy(dtype=float), copy=True)
            # median-impute; a covariate that is entirely missing is dropped
            for j in range(X.shape[1]):
                col = X[:, j]
                if np.isnan(col).all():
                    continue
                col[np.isnan(col)] = np.nanmedian(col)
            keep = ~np.isnan(X).all(axis=0)
            if keep.any():
                parts.append(X[:, keep])
        for c in cat:
            d = pd.get_dummies(cov[c].astype(str), drop_first=True)
            if d.shape[1]:
                parts.append(d.to_numpy(dtype=float))
        if not parts:
            return None
        return np.hstack(parts)

    blocks = {
        "qc": build(qc_cols, []),
        "demographic": build(demo_num, demo_cat),
        "batch": build([], batch_cat),
        "all": build(qc_cols + demo_num, demo_cat + batch_cat),
    }
    return {k: v for k, v in blocks.items() if v is not None and v.shape[1] > 0}


def label_information_fraction(y: np.ndarray, D: np.ndarray) -> float:
    """I_D = 1 - R^2(Y ~ D), the v2 quantity, computed on centred Y."""
    r2 = LinearRegression().fit(D, y).score(D, y)
    return float(1.0 - max(r2, 0.0))


def cross_fitted_auc(X: np.ndarray, y: np.ndarray, kind: str,
                     n_repeat: int = N_REPEAT,
                     fixed_C: float | None = None) -> float:
    """Cross-fitted out-of-fold AUC.

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

            sc = StandardScaler().fit(X[tr])
            if fixed_C is not None:
                best = fixed_C
            else:
                best, best_auc = CS[0], -np.inf
                inner = StratifiedKFold(3, shuffle=True, random_state=SEED)
                for C in CS:
                    s_ = []
                    for itr, ite in inner.split(X[tr], y[tr]):
                        mm = LogisticRegression(C=C, solver="liblinear",
                                                class_weight="balanced",
                                                max_iter=5000)
                        mm.fit(sc.transform(X[tr][itr]), y[tr][itr])
                        s_.append(roc_auc_score(
                            y[tr][ite],
                            mm.predict_proba(sc.transform(X[tr][ite]))[:, 1]))
                    if np.mean(s_) > best_auc:
                        best_auc, best = float(np.mean(s_)), C
            m = LogisticRegression(C=best, solver="liblinear",
                                   class_weight="balanced", max_iter=5000)
            m.fit(sc.transform(X[tr]), y[tr])
            oof[te] = m.predict_proba(sc.transform(X[te]))[:, 1]
        aucs.append(roc_auc_score(y, oof))
    return float(np.mean(aucs))


def design_strata(cov: pd.DataFrame) -> np.ndarray:
    """Discrete strata used by the design-preserving null."""
    cats = [c for c in cov.columns if c.startswith("batch__")]
    if "assay" in cov.columns and cov["assay"].nunique() > 1:
        cats.append("assay")
    if not cats:                       # fall back to QC tertiles
        q = pd.qcut(cov["log_cells_per_donor"], 3, labels=False, duplicates="drop")
        return q.to_numpy()
    key = cov[cats].astype(str).agg("|".join, axis=1)
    return pd.factorize(key)[0]


PERM_C = 1.0        # frozen regularisation for the permutation pipeline
PERM_PCS = 50       # frozen dimensionality for the permutation pipeline


def permutation_pvalues(X: np.ndarray, y: np.ndarray, strata: np.ndarray,
                        n_perm: int = N_PERM) -> dict:
    """Standard and design-preserving permutation p-values.

    Both nulls and the observed statistic they are compared against run the
    identical frozen pipeline: PCA to PERM_PCS components, logistic regression
    at PERM_C, single cross-fitting repeat.
    """
    from sklearn.decomposition import PCA

    n_pc = int(min(PERM_PCS, X.shape[1], len(y) - 2))
    Xp = PCA(n_components=n_pc, random_state=SEED).fit_transform(
        StandardScaler().fit_transform(X))

    def frozen_auc(yy: np.ndarray) -> float:
        if len(np.unique(yy)) < 2:
            return 0.5
        return cross_fitted_auc(Xp, yy, "linear", n_repeat=1, fixed_C=PERM_C)

    observed = frozen_auc(y)
    rng = np.random.default_rng(SEED)

    std_hits = dpn_hits = usable = 0
    for _ in range(n_perm):
        std_hits += frozen_auc(rng.permutation(y)) >= observed

        yp, moved = y.copy(), False
        for s_ in np.unique(strata):
            idx = np.flatnonzero(strata == s_)
            if len(idx) > 1 and len(np.unique(y[idx])) > 1:
                yp[idx] = rng.permutation(y[idx])
                moved = True
        if moved:
            usable += 1
            dpn_hits += frozen_auc(yp) >= observed

    return {
        "observed_frozen_auc": float(observed),
        "p_standard": float((std_hits + 1) / (n_perm + 1)),
        "p_design_preserving": (float((dpn_hits + 1) / (usable + 1))
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

    # top-variance genes keep the donor-level problem well conditioned
    Xp = pb.to_numpy(dtype=np.float64)
    var = Xp.var(axis=0)
    Xp = Xp[:, np.argsort(var)[::-1][:4000]]

    print(f"[{cohort}] donors={len(y)} cases={int(y.sum())} controls={int((y==0).sum())}")
    disease_auc = cross_fitted_auc(Xp, y, "linear")
    print(f"  disease AUC (pseudobulk) = {disease_auc:.3f}")

    strata = design_strata(cov)
    perm = permutation_pvalues(Xp, y, strata, n_perm)
    p_std, p_dpn = perm["p_standard"], perm["p_design_preserving"]
    print(f"  frozen-pipeline AUC={perm['observed_frozen_auc']:.3f} "
          f"(PCA {perm['perm_n_pc']} comps, C={PERM_C})")
    print(f"  p_standard={p_std:.4f}  p_design_preserving={p_dpn:.4f} "
          f"(usable stratified permutations: {perm['n_usable_stratified_perm']})")

    rows = []
    for block, D in design_matrix(cov).items():
        i_d = label_information_fraction(y.astype(float), D)
        a_lin = cross_fitted_auc(D, y, "linear")
        a_rf = cross_fitted_auc(D, y, "rf")
        print(f"  [{block:<12}] I_D={i_d:.3f}  design AUC lin={a_lin:.3f} rf={a_rf:.3f}")
        rows.append({
            "cohort": cohort, "block": block, "n_donor": len(y),
            "n_case": int(y.sum()), "n_design_feature": int(D.shape[1]),
            "I_D": round(i_d, 4),
            "design_auc_linear": round(a_lin, 4),
            "design_auc_rf": round(a_rf, 4),
            "disease_auc": round(disease_auc, 4),
            "observed_frozen_auc": round(perm["observed_frozen_auc"], 4),
            "p_standard": round(p_std, 5),
            "p_design_preserving": None if np.isnan(p_dpn) else round(p_dpn, 5),
            "n_strata": int(len(np.unique(strata))),
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    args = ap.parse_args()

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

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "design_screen.tsv", sep="\t", index=False)

    spectrum = (df[df.block == "all"]
                .sort_values("I_D")[["cohort", "n_donor", "I_D",
                                     "design_auc_linear", "design_auc_rf",
                                     "disease_auc", "p_standard",
                                     "p_design_preserving"]])
    spectrum.to_csv(OUT / "confounding_spectrum.tsv", sep="\t", index=False)
    (OUT / "screen_config.json").write_text(json.dumps(
        {"seed": SEED, "n_splits": N_SPLITS, "n_repeat": N_REPEAT,
         "n_perm": args.n_perm, "cohorts": cohorts}, indent=2))

    print("\n=== design-confounding spectrum (block = all) ===")
    print(spectrum.to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
