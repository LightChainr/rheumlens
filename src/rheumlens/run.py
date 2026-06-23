from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from rheumlens.analysis.summary import load_result_tree, paired_method_table, summarize_oof
from rheumlens.catalog import catalog_rows
from rheumlens.config import load_yaml
from rheumlens.data.io import load_npz_dataset, save_npz_dataset
from rheumlens.data.validation import validate_embedding, validate_raw_counts
from rheumlens.evaluation.engine import BagMethod, FixedMethod, run_bag_oof, run_fixed_oof
from rheumlens.evaluation.splits import load_folds, make_stratified_folds, save_folds
from rheumlens.evaluation.transfer import run_source_target
from rheumlens.evaluation.workflows import (
    permutation_workflow,
    repeated_cv_workflow,
    successful_auc,
    surrogate_delta_workflow,
)
from rheumlens.registry import RegisteredMethod, build_method
from rheumlens.types import CellDataset
from rheumlens.utils.arrays import ensure_dir
from rheumlens.utils.manifests import write_json


STAGE_METHODS = {
    "reproduce": [
        "donor_mean_hvg",
        "donor_expression_pca",
        "scgpt_mean",
        "geneformer_mean",
        "raw_pseudobulk",
        "isg_scalar",
        "isg_vector",
    ],
    "benchmark_fixed": [
        "composition_clr",
        "celltype_pseudobulk",
        "scfeatures_multiview",
        "moments_mean_var",
        "robust_median_mad",
        "mean_skew",
        "quantiles",
        "tail_fractions",
        "shrinkage_covariance",
        "kme_multiscale",
        "robust_kme",
        "mmd_control",
        "true_bures",
        "sliced_wasserstein",
        "prototype_histogram",
    ],
    "benchmark_learned": ["deep_sets", "attention_mil", "topk_mil", "set_transformer", "mixmil"],
    "originals": ["cckme_u", "cckme_u_weighted", "focus_lite", "red", "gds", "donorclr"],
}


def _dataset_map(cohort_config: dict[str, Any]) -> dict[str, CellDataset]:
    return {key: load_npz_dataset(value) for key, value in cohort_config.get("data", {}).items() if value}


def _load_covariates(cohort_config: dict[str, Any]) -> pd.DataFrame | None:
    path = cohort_config.get("covariates")
    if not path:
        return None
    frame = pd.read_csv(path, dtype={"donor_id": str})
    if "donor_id" not in frame:
        raise ValueError(f"covariate table {path} lacks donor_id")
    return frame


def _donor_table(data: CellDataset) -> tuple[np.ndarray, np.ndarray]:
    mapping = data.donor_label_map()
    donors = np.asarray(list(mapping), dtype=str)
    y = np.asarray([mapping[d] for d in donors], dtype=int)
    return donors, y


def _load_or_make_folds(cohort_config: dict[str, Any], data: CellDataset, output_root: Path) -> list:
    fold_path = cohort_config.get("folds")
    if fold_path and Path(fold_path).exists():
        return load_folds(fold_path)
    if cohort_config.get("require_authoritative_folds", False):
        raise FileNotFoundError(f"authoritative folds are required but missing: {fold_path}")
    donors, y = _donor_table(data)
    folds = make_stratified_folds(
        donors,
        y,
        n_splits=int(cohort_config.get("n_splits", 5)),
        seed=int(cohort_config.get("fold_seed", 42)),
        split_id=str(cohort_config.get("split_id", "generated")),
    )
    path = output_root / "generated_folds.csv"
    save_folds(path, folds)
    return folds


def _defaults(config: dict[str, Any]) -> dict[str, Any]:
    defaults = dict(config.get("method_defaults", {}))
    defaults.setdefault("estimator_C", config.get("evaluation", {}).get("estimator_C", 1.0))
    return defaults


def _query_bank(config: dict[str, Any]) -> dict[str, Any] | None:
    path = config.get("query_bank")
    return load_yaml(path) if path else None


def _run_registered(
    registered: RegisteredMethod,
    data: CellDataset,
    folds: list,
    cohort_name: str,
    estimand: str,
    expression: CellDataset | None,
    covariates: pd.DataFrame | None,
    seed: int,
) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if isinstance(registered.method, FixedMethod):
        return run_fixed_oof(
            data,
            folds,
            registered.method,
            cohort=cohort_name,
            estimand=estimand,
            expression_data=expression,
            donor_covariates=covariates,
            random_state=seed,
        )
    frame = run_bag_oof(
        data,
        folds,
        registered.method,
        cohort=cohort_name,
        estimand=estimand,
        random_state=seed,
    )
    return frame, []


def _run_stage(config: dict[str, Any], stage: str, methods_override: list[str] | None = None) -> None:
    root = Path(config.get("paths", {}).get("root", "."))
    result_root = ensure_dir(config.get("paths", {}).get("results", root / "results"))
    defaults = _defaults(config)
    query_bank = _query_bank(config)
    requested = methods_override or config.get("stages", {}).get(stage) or STAGE_METHODS[stage]
    for cohort_name, cohort_config in config.get("datasets", {}).items():
        datasets = _dataset_map(cohort_config)
        if not datasets:
            print(f"[{cohort_name}] no datasets configured; skipping", file=sys.stderr)
            continue
        anchor = datasets.get("scgpt") or datasets.get("lognorm") or next(iter(datasets.values()))
        folds = _load_or_make_folds(cohort_config, anchor, result_root / cohort_name)
        covariates = _load_covariates(cohort_config)
        estimand = str(cohort_config.get("estimand", "within_cohort"))
        for method_id in requested:
            try:
                registered = build_method(method_id, defaults, query_bank, seed=int(config.get("seed", 0)))
            except Exception as exc:
                print(f"[{cohort_name}] {method_id}: registration failed: {exc}", file=sys.stderr)
                continue
            if registered.data_key not in datasets:
                print(f"[{cohort_name}] {method_id}: missing data key {registered.data_key}; NOT_APPLICABLE")
                continue
            data = datasets[registered.data_key]
            expression = datasets.get(registered.expression_key) if registered.expression_key else None
            output = ensure_dir(result_root / stage / cohort_name / method_id.replace("@", "__"))
            frame, diagnostics = _run_registered(
                registered,
                data,
                folds,
                cohort_name,
                estimand,
                expression,
                covariates,
                int(config.get("seed", 0)),
            )
            frame.to_csv(output / "oof.csv", index=False)
            write_json(output / "diagnostics.json", {"folds": diagnostics})
            print(f"[{cohort_name}] completed {method_id}: {len(frame)} OOF rows; AUC={successful_auc(frame):.4f}")


def _validate(config: dict[str, Any]) -> None:
    reports: dict[str, Any] = {}
    for cohort_name, cohort_config in config.get("datasets", {}).items():
        reports[cohort_name] = {}
        datasets = _dataset_map(cohort_config)
        anchor: CellDataset | None = None
        for key, data in datasets.items():
            report = validate_raw_counts(data) if key == "raw_counts" else validate_embedding(data)
            report["path"] = cohort_config.get("data", {}).get(key)
            reports[cohort_name][key] = report
            anchor = anchor or data
        if anchor is not None:
            folds = _load_or_make_folds(cohort_config, anchor, Path(config.get("paths", {}).get("results", "results")) / cohort_name)
            donor_set = set(anchor.donors.astype(str))
            fold_test = [d for fold in folds for d in fold.test_donors]
            reports[cohort_name]["fold_audit"] = {
                "n_folds": len(folds),
                "all_test_donors_once": len(fold_test) == len(set(fold_test)) == len(donor_set),
                "missing_test_donors": sorted(donor_set - set(fold_test)),
                "unknown_test_donors": sorted(set(fold_test) - donor_set),
            }
    print(json.dumps(reports, indent=2, ensure_ascii=False))


def _final_analysis(config: dict[str, Any]) -> None:
    result_root = Path(config.get("paths", {}).get("results", "results"))
    frame = load_result_tree(result_root)
    if frame.empty:
        raise RuntimeError(f"no oof.csv files under {result_root}")
    output = ensure_dir(result_root / "final_tables")
    summarize_oof(frame).to_csv(output / "method_summary.csv", index=False)
    reference = str(config.get("analysis", {}).get("reference_method", "scgpt_mean"))
    paired_method_table(
        frame,
        reference,
        n_reps=int(config.get("evaluation", {}).get("bootstrap_reps", 10_000)),
        seed=int(config.get("seed", 0)),
    ).to_csv(output / f"paired_vs_{reference}.csv", index=False)
    frame.to_csv(output / "all_oof.csv", index=False)
    pd.DataFrame(catalog_rows()).to_csv(output / "method_catalog.csv", index=False)
    print(f"wrote final tables to {output}")


def _permutation_stage(config: dict[str, Any], methods_override: list[str] | None) -> None:
    result_root = ensure_dir(config.get("paths", {}).get("results", "results"))
    defaults, query_bank = _defaults(config), _query_bank(config)
    settings = config.get("permutation", {})
    methods = methods_override or settings.get("methods", ["scgpt_mean", "geneformer_mean", "donor_expression_pca"])
    n_reps = int(settings.get("n_reps", config.get("evaluation", {}).get("permutation_reps", 1000)))
    for cohort_name, cohort_config in config.get("datasets", {}).items():
        datasets = _dataset_map(cohort_config)
        if not datasets:
            continue
        anchor = datasets.get("scgpt") or datasets.get("lognorm") or next(iter(datasets.values()))
        folds = _load_or_make_folds(cohort_config, anchor, result_root / cohort_name)
        covariates = _load_covariates(cohort_config)
        estimand = str(cohort_config.get("estimand", "within_cohort"))
        for method_id in methods:
            registered = build_method(method_id, defaults, query_bank, seed=int(config.get("seed", 0)))
            if registered.data_key not in datasets:
                continue
            data = datasets[registered.data_key]
            expression = datasets.get(registered.expression_key) if registered.expression_key else None
            observed, _ = _run_registered(
                registered, data, folds, cohort_name, estimand, expression, covariates, int(config.get("seed", 0))
            )
            observed_auc = successful_auc(observed)
            out = result_root / "permutation" / cohort_name / method_id.replace("@", "__")
            result = permutation_workflow(
                data,
                folds,
                registered.method,
                cohort_name,
                estimand,
                observed_auc,
                n_reps,
                int(config.get("seed", 0)),
                expression_data=expression,
                donor_covariates=covariates,
                rebuild_stratified_folds=bool(settings.get("rebuild_stratified_folds", False)),
                output_dir=out,
            )
            print(f"[{cohort_name}] permutation {method_id}: p={result['empirical_p']}")


def _repeated_cv_stage(config: dict[str, Any], methods_override: list[str] | None) -> None:
    result_root = ensure_dir(config.get("paths", {}).get("results", "results"))
    defaults, query_bank = _defaults(config), _query_bank(config)
    settings = config.get("repeated_cv", {})
    methods = methods_override or settings.get("methods", ["scgpt_mean", "geneformer_mean", "donor_expression_pca", "donor_mean_hvg"])
    seeds = list(map(int, settings.get("seeds", list(range(30)))))
    n_splits = int(settings.get("n_splits", 5))
    for cohort_name, cohort_config in config.get("datasets", {}).items():
        datasets = _dataset_map(cohort_config)
        covariates = _load_covariates(cohort_config)
        estimand = str(cohort_config.get("estimand", "within_cohort"))
        for method_id in methods:
            registered = build_method(method_id, defaults, query_bank, seed=int(config.get("seed", 0)))
            if registered.data_key not in datasets:
                continue
            data = datasets[registered.data_key]
            expression = datasets.get(registered.expression_key) if registered.expression_key else None
            out = result_root / "repeated_cv" / cohort_name / method_id.replace("@", "__")
            result = repeated_cv_workflow(
                data,
                registered.method,
                cohort_name,
                estimand,
                seeds,
                n_splits=n_splits,
                expression_data=expression,
                donor_covariates=covariates,
                output_dir=out,
            )
            print(f"[{cohort_name}] repeated CV {method_id}: mean={result.auc.mean():.4f}")


def _surrogate_stage(config: dict[str, Any], methods_override: list[str] | None) -> None:
    result_root = ensure_dir(config.get("paths", {}).get("results", "results"))
    defaults, query_bank = _defaults(config), _query_bank(config)
    settings = config.get("surrogate_nulls", {})
    methods = methods_override or settings.get("methods", ["moments_mean_var", "kme_multiscale", "cckme_u"])
    reference_id = str(settings.get("reference_method", "scgpt_mean"))
    kinds = settings.get("kinds", ["mean_sufficient", "moment2_sufficient"])
    n_reps = int(settings.get("n_reps", 100))
    for cohort_name, cohort_config in config.get("datasets", {}).items():
        datasets = _dataset_map(cohort_config)
        if "scgpt" not in datasets:
            continue
        data = datasets["scgpt"]
        folds = _load_or_make_folds(cohort_config, data, result_root / cohort_name)
        covariates = _load_covariates(cohort_config)
        estimand = str(cohort_config.get("estimand", "within_cohort"))
        reference = build_method(reference_id, defaults, query_bank, seed=int(config.get("seed", 0)))
        if reference.data_key != "scgpt":
            raise ValueError("surrogate reference must use the scgpt representation")
        for method_id in methods:
            registered = build_method(method_id, defaults, query_bank, seed=int(config.get("seed", 0)))
            if registered.data_key != "scgpt" or registered.expression_key is not None:
                print(f"[{cohort_name}] surrogate skip {method_id}: requires different or auxiliary data")
                continue
            for kind in kinds:
                out = result_root / "surrogate_nulls" / cohort_name / method_id.replace("@", "__") / kind
                result = surrogate_delta_workflow(
                    data,
                    folds,
                    registered.method,
                    reference.method,
                    cohort_name,
                    estimand,
                    kind,
                    n_reps,
                    int(config.get("seed", 0)),
                    donor_covariates=covariates,
                    output_dir=out,
                )
                print(f"[{cohort_name}] {kind} {method_id}: q95={result['delta_q95']}")


def _transfer_stage(config: dict[str, Any], methods_override: list[str] | None) -> None:
    result_root = ensure_dir(config.get("paths", {}).get("results", "results"))
    defaults, query_bank = _defaults(config), _query_bank(config)
    all_datasets = {name: _dataset_map(value) for name, value in config.get("datasets", {}).items()}
    pairs = config.get("transfer", {}).get("pairs", [])
    for pair in pairs:
        source_name, target_name = str(pair["source"]), str(pair["target"])
        methods = methods_override or pair.get("methods", ["scgpt_mean", "geneformer_mean", "donor_expression_pca"])
        source_cov = _load_covariates(config["datasets"][source_name])
        target_cov = _load_covariates(config["datasets"][target_name])
        for method_id in methods:
            registered = build_method(method_id, defaults, query_bank, seed=int(config.get("seed", 0)))
            if not isinstance(registered.method, FixedMethod):
                print(f"transfer skip trainable bag method {method_id}")
                continue
            source_map, target_map = all_datasets[source_name], all_datasets[target_name]
            if registered.data_key not in source_map or registered.data_key not in target_map:
                continue
            expression_source = source_map.get(registered.expression_key) if registered.expression_key else None
            expression_target = target_map.get(registered.expression_key) if registered.expression_key else None
            frame = run_source_target(
                source_map[registered.data_key],
                target_map[registered.data_key],
                registered.method,
                expression_source,
                expression_target,
                source_cov,
                target_cov,
                random_state=int(config.get("seed", 0)),
            )
            out = ensure_dir(result_root / "transfer" / f"{source_name}_to_{target_name}" / method_id.replace("@", "__"))
            frame.to_csv(out / "predictions.csv", index=False)
            auc = roc_auc_score(frame.y_true, frame.score) if frame.y_true.nunique() == 2 else np.nan
            write_json(out / "summary.json", {"auc": float(auc), "n_donors": len(frame)})
            print(f"transfer {source_name}->{target_name} {method_id}: AUC={auc:.4f}")


def _synthetic_dataset(seed: int = 0) -> tuple[CellDataset, CellDataset]:
    rng = np.random.default_rng(seed)
    n_donors, cells_per, p = 24, 30, 20
    donor_ids, y, X, expression, cell_types, cell_ids = [], [], [], [], [], []
    genes = np.asarray([f"GENE{i}" for i in range(p)], dtype=str)
    genes[:15] = np.asarray(
        ["ISG15", "IFI6", "MX1", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT3", "IFI44", "IFI44L", "STAT1", "RSAD2", "IFITM1", "IFITM3", "HERC5"]
    )
    projection = rng.normal(0, 0.2, size=(p, p))
    for d in range(n_donors):
        label = d % 2
        base = rng.normal(0, 0.4, p)
        for i in range(cells_per):
            rare = label == 1 and rng.random() < 0.10
            expr = rng.normal(base, 1.0, p)
            expr[:15] += label * 0.25 + rare * 2.5
            emb = np.tanh(expr @ projection)
            X.append(emb)
            expression.append(expr)
            donor_ids.append(f"D{d:03d}")
            y.append(label)
            cell_types.append("CD4_T" if i % 3 else "monocyte")
            cell_ids.append(f"D{d:03d}_C{i:04d}")
    common = dict(
        cell_ids=np.asarray(cell_ids), donor_ids=np.asarray(donor_ids), y=np.asarray(y),
        cell_types=np.asarray(cell_types), feature_names=genes,
    )
    return (
        CellDataset(X=np.asarray(X, dtype=np.float32), name="synthetic_scgpt", **common),
        CellDataset(X=np.asarray(expression, dtype=np.float32), name="synthetic_lognorm", **common),
    )


def _smoke(config: dict[str, Any]) -> None:
    root = Path(config.get("paths", {}).get("root", "/tmp/rheumlens_smoke"))
    data_dir, result_dir = ensure_dir(root / "data"), ensure_dir(root / "results")
    embedding, lognorm = _synthetic_dataset(int(config.get("seed", 0)))
    save_npz_dataset(data_dir / "scgpt.npz", embedding)
    save_npz_dataset(data_dir / "lognorm.npz", lognorm)
    raw = CellDataset(
        X=np.rint(np.exp(np.clip(lognorm.X, -2, 4))).astype(np.int32),
        cell_ids=lognorm.cell_ids, donor_ids=lognorm.donor_ids, y=lognorm.y,
        feature_names=lognorm.feature_names, cell_types=lognorm.cell_types, name="synthetic_raw",
    )
    save_npz_dataset(data_dir / "raw_counts.npz", raw)
    donors, y = _donor_table(embedding)
    fold_path = root / "splits.csv"
    save_folds(fold_path, make_stratified_folds(donors, y, n_splits=3, seed=42, split_id="smoke"))
    smoke_config = {
        **config,
        "paths": {"root": str(root), "results": str(result_dir)},
        "datasets": {"synthetic": {"data": {"scgpt": str(data_dir / "scgpt.npz"), "geneformer": str(data_dir / "scgpt.npz"), "lognorm": str(data_dir / "lognorm.npz"), "raw_counts": str(data_dir / "raw_counts.npz")}, "folds": str(fold_path), "estimand": "smoke"}},
        "method_defaults": {
            "cell_pca_dim": 6, "kme_rff_dim": 12, "kme_scales": [1.0],
            "kme_max_bandwidth_points": 256, "kme_max_diagnostic_points": 64,
            "matched_cells_per_donor": 20, "cckme_subsamples": 3,
        },
        "evaluation": {"bootstrap_reps": 20, "permutation_reps": 5, "estimator_C": 1.0},
    }
    for stage, methods in {
        "reproduce": ["scgpt_mean", "donor_expression_pca", "isg_scalar"],
        "benchmark_fixed": ["moments_mean_var", "quantiles", "kme_multiscale"],
        "originals": ["cckme_u", "focus_lite", "red"],
    }.items():
        _run_stage(smoke_config, stage, methods)
    _final_analysis(smoke_config)
    print(f"smoke test completed under {root}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="RheumLens v4 donor-level benchmark runner")
    parser.add_argument("--config", required=False)
    parser.add_argument(
        "--stage",
        required=True,
        choices=[
            "list_methods", "validate", "smoke", "reproduce", "benchmark_fixed", "benchmark_learned",
            "originals", "permutation", "repeated_cv", "surrogate_nulls", "transfer", "final_analysis", "all",
        ],
    )
    parser.add_argument("--methods", nargs="*", default=None)
    args = parser.parse_args(argv)
    if args.stage == "list_methods":
        print(pd.DataFrame(catalog_rows()).to_string(index=False))
        return
    if not args.config:
        parser.error("--config is required for this stage")
    config = load_yaml(args.config)
    if args.stage == "validate":
        _validate(config)
    elif args.stage == "smoke":
        _smoke(config)
    elif args.stage == "final_analysis":
        _final_analysis(config)
    elif args.stage == "permutation":
        _permutation_stage(config, args.methods)
    elif args.stage == "repeated_cv":
        _repeated_cv_stage(config, args.methods)
    elif args.stage == "surrogate_nulls":
        _surrogate_stage(config, args.methods)
    elif args.stage == "transfer":
        _transfer_stage(config, args.methods)
    elif args.stage == "all":
        _validate(config)
        for stage in ("reproduce", "benchmark_fixed", "benchmark_learned", "originals"):
            _run_stage(config, stage, None)
        if config.get("run_optional", {}).get("repeated_cv", False):
            _repeated_cv_stage(config, None)
        if config.get("run_optional", {}).get("permutation", False):
            _permutation_stage(config, None)
        if config.get("run_optional", {}).get("surrogate_nulls", False):
            _surrogate_stage(config, None)
        if config.get("run_optional", {}).get("transfer", False):
            _transfer_stage(config, None)
        _final_analysis(config)
    else:
        _run_stage(config, args.stage, args.methods)


if __name__ == "__main__":
    main()
