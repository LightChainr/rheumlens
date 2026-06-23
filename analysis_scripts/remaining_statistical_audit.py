#!/usr/bin/env python3
"""Remaining statistical audit for the frozen RheumLens donor benchmark.

Produces:
- permutation P values for all three primary models in all cohorts;
- an auditable PBMC-to-CD4 transfer rerun and provenance record;
- SHA256 manifests for the resulting evidence files.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

from fold_contained_benchmark import COHORTS, classifier, folds_from_json, select_hvg


ROOT = Path(__file__).resolve().parents[1]
FEATURE_DIR = ROOT / "data" / "donor_features"
FROZEN_DIR = ROOT / "results" / "fold_contained"
OUT = ROOT / "results" / "evidence_package"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_array(x: np.ndarray) -> str:
    x = np.ascontiguousarray(x)
    return hashlib.sha256(x.view(np.uint8)).hexdigest()


def prepare_fold_matrices(expression: np.ndarray, scgpt: np.ndarray,
                          folds: list[tuple[np.ndarray, np.ndarray]],
                          n_hvg: int = 2000, n_pcs: int = 25):
    """Fit label-independent transforms once; equivalent under label permutation."""
    prepared: dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]] = {
        "scgpt": [], "expression_pca": [], "donor_mean_hvg": [],
    }
    for tr, te in folds:
        scaler = StandardScaler().fit(scgpt[tr])
        prepared["scgpt"].append((tr, te, scaler.transform(scgpt[tr]), scaler.transform(scgpt[te])))

        hvg = select_hvg(expression[tr], n_hvg)
        pca = PCA(
            n_components=min(n_pcs, len(tr) - 1, len(hvg)),
            svd_solver="full", random_state=42,
        )
        train_pc = pca.fit_transform(expression[tr][:, hvg])
        test_pc = pca.transform(expression[te][:, hvg])
        pc_scaler = StandardScaler().fit(train_pc)
        prepared["expression_pca"].append(
            (tr, te, pc_scaler.transform(train_pc), pc_scaler.transform(test_pc))
        )

        hvg_scaler = StandardScaler().fit(expression[tr][:, hvg])
        prepared["donor_mean_hvg"].append(
            (tr, te, hvg_scaler.transform(expression[tr][:, hvg]),
             hvg_scaler.transform(expression[te][:, hvg]))
        )
    return prepared


def one_permutation(seed: int, y: np.ndarray, prepared) -> dict[str, float] | None:
    rng = np.random.default_rng(seed)
    yp = rng.permutation(y)
    if any(np.unique(yp[tr]).size < 2 for tr, _, _, _ in prepared["scgpt"]):
        return None
    aucs = {}
    for method, fold_data in prepared.items():
        prob = np.full(len(y), np.nan)
        for tr, te, xtr, xte in fold_data:
            model = classifier().fit(xtr, yp[tr])
            prob[te] = model.predict_proba(xte)[:, 1]
        aucs[method] = float(roc_auc_score(yp, prob))
    return aucs


def primary_permutations(cohort: str, n_perm: int, n_jobs: int) -> pd.DataFrame:
    z = np.load(FEATURE_DIR / f"{cohort}_donor_features.npz", allow_pickle=True)
    donors = z["donors"].astype(str)
    y = z["y"].astype(int)
    expression = z["expression"].astype(np.float64)
    scgpt = z["scgpt"].astype(np.float64)
    folds = folds_from_json(Path(COHORTS[cohort]["folds"]), donors, y)
    prepared = prepare_fold_matrices(expression, scgpt, folds)

    frozen = pd.read_csv(FROZEN_DIR / f"{cohort}_oof_predictions.csv")
    observed = {
        "scgpt": roc_auc_score(y, frozen["scgpt"]),
        "expression_pca": roc_auc_score(y, frozen["expression_pca"]),
        "donor_mean_hvg": roc_auc_score(y, frozen["hvg_pseudobulk"]),
    }

    # Generate extra deterministic candidates because a rare permutation can leave
    # a training fold with one class, especially in the 26-donor cohort.
    seeds = np.random.SeedSequence(20260619).generate_state(max(n_perm * 2, n_perm + 100))
    # Threads avoid macOS sandbox semaphore restrictions on loky processes.
    values = Parallel(n_jobs=n_jobs, verbose=5, prefer="threads")(
        delayed(one_permutation)(int(seed), y, prepared) for seed in seeds
    )
    valid = [v for v in values if v is not None][:n_perm]
    if len(valid) < n_perm:
        raise RuntimeError(f"Only {len(valid)} valid permutations for {cohort}")

    null_df = pd.DataFrame(valid)
    null_df.insert(0, "permutation", np.arange(1, n_perm + 1))
    null_path = OUT / f"{cohort}_primary_model_permutation_null.csv.gz"
    null_df.to_csv(null_path, index=False, compression="gzip")

    rows = []
    for method, obs in observed.items():
        null = null_df[method].to_numpy()
        exceed = int(np.sum(null >= obs))
        rows.append({
            "cohort": cohort,
            "method": method,
            "observed_auc": obs,
            "n_perm": n_perm,
            "exceedances": exceed,
            "empirical_p": (exceed + 1) / (n_perm + 1),
            "null_mean": float(null.mean()),
            "null_sd": float(null.std(ddof=1)),
            "permutation_scheme": "global donor-label permutation; fixed archived folds; all fold transforms refit-equivalent because label-independent",
        })
    return pd.DataFrame(rows)


def transfer_provenance() -> None:
    source = np.load(FEATURE_DIR / "GSE135779_donor_features.npz", allow_pickle=True)
    target = np.load(FEATURE_DIR / "GSE285773_donor_features.npz", allow_pickle=True)
    xs = source["scgpt"].astype(np.float64)
    ys = source["y"].astype(int)
    xt = target["scgpt"].astype(np.float64)
    yt = target["y"].astype(int)
    donors = target["donors"].astype(str)

    rows = {"target_donor_id": donors, "label": yt}
    models = []
    for name, scale in (("source_scaled", True), ("unscaled", False)):
        if scale:
            scaler = StandardScaler().fit(xs)
            xtrain, xtest = scaler.transform(xs), scaler.transform(xt)
            scaler_hash = sha256_array(np.concatenate([scaler.mean_, scaler.scale_]))
        else:
            xtrain, xtest = xs, xt
            scaler_hash = "none"
        model = classifier().fit(xtrain, ys)
        prob = model.predict_proba(xtest)[:, 1]
        rows[f"transfer_probability_{name}"] = prob
        models.append({
            "analysis": f"GSE135779_all_PBMC_to_GSE285773_CD4_{name}",
            "source_cell_scope": "all PBMC cells; no cell_type column available",
            "target_cell_scope": "cohort-level CD4 T-cell object; no per-cell subtype column available",
            "n_source_donors": len(ys),
            "n_target_donors": len(yt),
            "auc": roc_auc_score(yt, prob),
            "coefficient_sha256": sha256_array(model.coef_),
            "intercept_sha256": sha256_array(model.intercept_),
            "scaler_sha256": scaler_hash,
            "interpretation": "exploratory unmatched-cell-type transfer; not a matched CD4-to-CD4 validation",
        })

    pred_path = OUT / "GSE135779_to_GSE285773_transfer_predictions.csv"
    pd.DataFrame(rows).to_csv(pred_path, index=False)
    within_path = FROZEN_DIR / "GSE285773_oof_predictions.csv"
    provenance = {
        "frozen_within_cohort": {
            "analysis": "GSE285773 five-fold within-cohort OOF scGPT",
            "auc": float(roc_auc_score(yt, pd.read_csv(within_path)["scgpt"])),
            "prediction_file": str(within_path.relative_to(ROOT)),
            "prediction_file_sha256": sha256_file(within_path),
        },
        "legacy_transfer_claim": {
            "reported_auc": 0.950,
            "status": "unsupported_hardcoded_value",
            "evidence": "legacy p3_transfer.py appends transfer_AUC=0.950 as a literal and no donor prediction file was found",
        },
        "auditable_transfer_rerun": {
            "prediction_file": str(pred_path.relative_to(ROOT)),
            "prediction_file_sha256": sha256_file(pred_path),
            "models": models,
        },
        "decision": "Remove the legacy 0.950 transfer claim from the primary manuscript. The rerun is PBMC-to-CD4, preprocessing-sensitive, and exploratory.",
    }
    (OUT / "transfer_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def write_manifest() -> None:
    files = sorted(p for p in OUT.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
    text = "".join(f"{sha256_file(p)}  {p.name}\n" for p in files)
    (OUT / "SHA256SUMS.txt").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--jobs", type=int, default=-2)
    parser.add_argument("--cohort", action="append", choices=COHORTS)
    parser.add_argument("--skip-permutations", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    transfer_provenance()
    if not args.skip_permutations:
        selected = args.cohort or list(COHORTS)
        for cohort in selected:
            primary_permutations(cohort, args.permutations, args.jobs).to_csv(
                OUT / f"{cohort}_primary_model_permutation_tests.csv", index=False
            )
        existing = [OUT / f"{c}_primary_model_permutation_tests.csv" for c in COHORTS]
        if all(path.exists() for path in existing):
            pd.concat([pd.read_csv(path) for path in existing], ignore_index=True).to_csv(
                OUT / "all_primary_models_permutation_tests.csv", index=False
            )
    write_manifest()


if __name__ == "__main__":
    main()
