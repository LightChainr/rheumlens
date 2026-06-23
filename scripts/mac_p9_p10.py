#!/usr/bin/env python3
"""Run formal P9 transfer and P10 archival integration on the Mac snapshot.

The snapshot is treated as read-only. Every fitted operation in P9 uses source
cohort data only; target labels are accessed only after prediction for metrics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from scipy.stats import norm


SEED = 20260622
METHODS = ("scgpt_mean", "donor_expression_pca", "donor_mean_hvg")
EMBEDDING_METHODS = {"scgpt_mean"}
PAIRS = (("GSE285773", "GSE174188_CD4"), ("GSE174188_CD4", "GSE285773"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def metrics(frame: pd.DataFrame) -> dict[str, float]:
    y = frame.y_true.to_numpy(int)
    raw_score = frame.score.to_numpy(float)
    p = np.clip(raw_score, 1e-6, 1 - 1e-6)
    logit = np.log(p / (1 - p)).reshape(-1, 1)
    calibration = LogisticRegression(C=1e6, solver="lbfgs").fit(logit, y)
    return {
        "auc": float(roc_auc_score(y, raw_score)),
        "pr_auc": float(average_precision_score(y, raw_score)),
        "brier": float(brier_score_loss(y, p)),
        "calibration_intercept": float(calibration.intercept_[0]),
        "calibration_slope": float(calibration.coef_[0, 0]),
        "n_donors": int(len(y)),
        "n_cases": int(y.sum()),
        "n_controls": int((1 - y).sum()),
    }


def auc_delong_ci(frame: pd.DataFrame, alpha: float = 0.95) -> tuple[float, float, float]:
    """Hanley-McNeil/DeLong-style large-sample CI for binary AUC.

    This is used for the main P9 AUC interval because ordinary donor bootstrap
    is unstable for the 26-donor GSE285773 target direction.
    """
    y = frame.y_true.to_numpy(int)
    scores = frame.score.to_numpy(float)
    pos = scores[y == 1]
    neg = scores[y == 0]
    m, n = len(pos), len(neg)
    auc = roc_auc_score(y, scores)
    q1 = auc / (2.0 - auc)
    q2 = 2.0 * auc * auc / (1.0 + auc)
    var = (auc * (1.0 - auc) + (m - 1) * (q1 - auc * auc) + (n - 1) * (q2 - auc * auc)) / (m * n)
    se = float(np.sqrt(max(var, 0.0)))
    z = float(norm.ppf(0.5 + alpha / 2.0))
    return max(0.0, auc - z * se), min(1.0, auc + z * se), se


def bootstrap(frame: pd.DataFrame, reps: int = 10000) -> pd.DataFrame:
    y, p = frame.y_true.to_numpy(int), frame.score.to_numpy(float)
    case, control = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    rng = np.random.default_rng(SEED)
    rows = []
    for rep in range(reps):
        idx = np.r_[rng.choice(case, len(case), True), rng.choice(control, len(control), True)]
        rows.append((rep, roc_auc_score(y[idx], p[idx]), average_precision_score(y[idx], p[idx])))
    return pd.DataFrame(rows, columns=["rep", "auc", "pr_auc"])


def dataset_paths(root: Path, cohort: str) -> dict[str, Path]:
    if cohort == "GSE285773":
        return {
            "scgpt": root / "embeddings/scgpt/GSE285773_v1/GSE285773_scgpt.npz",
            "lognorm": root / "data/processed/GSE285773/lognorm.npz",
        }
    return {
        "scgpt": root / "embeddings/scgpt/GSE174188_CD4_v1/GSE174188_CD4_scgpt.npz",
        "lognorm": root / "data/processed/GSE174188_CD4/lognorm.npz",
    }


def run_p9(root: Path, out: Path, code_root: Path, bootstrap_reps: int) -> None:
    sys.path.insert(0, str(code_root / "src"))
    from rheumlens.data.io import load_npz_dataset
    from rheumlens.evaluation.transfer import run_source_target
    from rheumlens.registry import build_method

    out.mkdir(parents=True, exist_ok=True)
    eligibility = pd.DataFrame([
        {"source": "GSE285773", "target": "GSE174188_CD4", "eligible": True,
         "reason": "Both are CD4-focused; pediatric/adult and study shifts must be reported."},
        {"source": "GSE174188_CD4", "target": "GSE285773", "eligible": True,
         "reason": "Reverse CD4 transfer; small pediatric target requires wide uncertainty intervals."},
        {"source": "GSE135779", "target": "GSE174188_CD4/GSE285773", "eligible": False,
         "reason": "PBMC-to-CD4 compartment mismatch; no validated matched CD4 subset."},
    ])
    eligibility.to_csv(out / "P9_ELIGIBILITY_MATRIX.csv", index=False)
    all_oof, unavailable = [], []
    cache: dict[tuple[str, str], object] = {}
    for source_name, target_name in PAIRS:
        for method_id in METHODS:
            key = "scgpt" if method_id == "scgpt_mean" else "lognorm"
            if method_id not in EMBEDDING_METHODS:
                source_features = np.load(dataset_paths(root, source_name)[key], allow_pickle=True)["feature_names"]
                target_features = np.load(dataset_paths(root, target_name)[key], allow_pickle=True)["feature_names"]
                shared = np.intersect1d(source_features.astype(str), target_features.astype(str))
                if len(shared) == 0:
                    unavailable.append({
                        "direction": f"{source_name}_to_{target_name}",
                        "method_id": method_id,
                        "status": "UNAVAILABLE_FEATURE_ALIGNMENT",
                        "reason": "Critical Mac snapshot has no shared expression feature names for this transfer direction.",
                        "source_n_features": int(len(source_features)),
                        "target_n_features": int(len(target_features)),
                        "shared_n_features": int(len(shared)),
                    })
                    continue
            for cohort in (source_name, target_name):
                ck = (cohort, key)
                if ck not in cache:
                    path = dataset_paths(root, cohort)[key]
                    if not path.exists():
                        raise FileNotFoundError(path)
                    cache[ck] = load_npz_dataset(path)
            registered = build_method(method_id, defaults={"n_hvg": 2000, "donor_pca_dim": 25}, seed=SEED)
            frame = run_source_target(cache[(source_name, key)], cache[(target_name, key)], registered.method,
                                      random_state=SEED)
            frame["direction"] = f"{source_name}_to_{target_name}"
            frame["status"] = "SUCCESS"
            all_oof.append(frame)
    predictions = pd.concat(all_oof, ignore_index=True)
    predictions.to_csv(out / "P9_all_predictions.csv", index=False)
    summaries = []
    for (direction, method_id), frame in predictions.groupby(["direction", "method_id"], sort=False):
            row = {"direction": direction, "method_id": method_id, **metrics(frame)}
            boot = bootstrap(frame, bootstrap_reps)
            boot.to_csv(out / f"bootstrap_{direction}_{method_id}.csv.gz", index=False)
            row["auc_ci_method"] = "Hanley-McNeil"
            row["auc_ci_low"], row["auc_ci_high"], row["auc_se"] = auc_delong_ci(frame)
            row["pr_auc_ci_method"] = "stratified_donor_bootstrap"
            row["pr_auc_ci_low"], row["pr_auc_ci_high"] = np.quantile(boot["pr_auc"], [.025, .975])
            summaries.append(row)
    pd.DataFrame(summaries).to_csv(out / "P9_method_summary.csv", index=False)
    pd.DataFrame(unavailable).to_csv(out / "P9_unavailable_methods.csv", index=False)
    (out / "P9_CLAIM_BOUNDARY.md").write_text(
        "# P9 claim boundary\n\nThese are source-only cross-study CD4 transfer results. "
        "They do not establish clinical deployment, causality, or PBMC-to-CD4 generalization.\n", encoding="utf-8")


def run_p10(root: Path, p9: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    accepted = []
    for path in sorted((root / "results").rglob("*.md")):
        text = path.read_text(errors="ignore")
        if any(x in text for x in ("ACCEPT_RECOMMENDED", "ACCEPTED", "COMPLETED_TECHNICAL")):
            accepted.append({"path": str(path.relative_to(root)), "sha256": sha256(path)})
    pd.DataFrame(accepted).to_csv(out / "P10_ACCEPTED_ASSET_INDEX.csv", index=False)
    p9_summary = pd.read_csv(p9 / "P9_method_summary.csv")
    p9_summary.to_csv(out / "P10_P9_TRANSFER_SUMMARY.csv", index=False)
    unavailable_path = p9 / "P9_unavailable_methods.csv"
    unavailable = pd.read_csv(unavailable_path) if unavailable_path.exists() and unavailable_path.stat().st_size > 1 else pd.DataFrame()
    if unavailable.empty:
        p9_availability_text = "All configured P9 transfer methods completed successfully."
    else:
        p9_availability_text = unavailable.to_markdown(index=False)
    report = ["# RheumLens final integration report", "", "## Formal stage state", "",
              "- P4 Geneformer provenance: accepted.", "- P5 matched-500: accepted (27/27).",
              "- P6 three-cohort benchmarks: accepted with modality-specific boundaries.",
              "- P8.3/P8.4/P8.6/P8.7: accepted with documented caveats.",
              "- P8.5: skipped because no positive KME increment triggered the null audit.",
              "- P8.8: skipped after FOCUS superiority claim retraction; exploratory only.",
              "- P9: source-only CD4 cross-study transfer completed by this run.", "",
              "## P9 results", "", p9_summary.to_markdown(index=False), "",
              "## P9 method availability", "",
              p9_availability_text, "",
              "## Global claim boundary", "",
              "Performance is retrospective and cohort-specific. Covariate sensitivity materially attenuated AUC; "
              "none of the results establishes causal biology or clinical utility."]
    (out / "FINAL_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    files = [p for p in sorted(out.rglob("*")) if p.is_file() and p.name != "MANIFEST_SHA256.tsv"]
    pd.DataFrame([{"path": str(p.relative_to(out)), "sha256": sha256(p), "bytes": p.stat().st_size}
                  for p in files]).to_csv(out / "MANIFEST_SHA256.tsv", sep="\t", index=False)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--snapshot", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--code-root", type=Path, required=True)
    p.add_argument("--bootstrap-reps", type=int, default=10000)
    args = p.parse_args()
    run_p9(args.snapshot, args.output / "P9_transfer", args.code_root, args.bootstrap_reps)
    run_p10(args.snapshot, args.output / "P9_transfer", args.output / "P10_integration")


if __name__ == "__main__":
    main()
