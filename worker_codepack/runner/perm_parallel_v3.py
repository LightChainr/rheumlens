#!/usr/bin/env python3
"""Checkpointed donor-label permutation runner for RheumLens P8.2.

The parent loads one local staged dataset and forks read-only workers. Workers never
read the shared filesystem and never write checkpoints. The parent is the sole
writer and periodically publishes a compact durable checkpoint.
"""

from __future__ import annotations

import os

# Set before importing numerical libraries.
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import argparse
import hashlib
import json
import multiprocessing as mp
import signal
import sys
import tempfile
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml


G: dict[str, Any] = {}
STOP_REQUESTED = False


class GracefulStop(RuntimeError):
    pass


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_json(obj: Any) -> str:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode()).hexdigest()


def atomic_write_bytes(path: Path, payload: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
        dir_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)
    except BaseException:
        try:
            os.unlink(tmp_name)
        except FileNotFoundError:
            pass
        raise


def atomic_write_json(path: Path, obj: Any) -> None:
    atomic_write_bytes(path, (json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n").encode())


def atomic_write_csv(path: Path, frame: pd.DataFrame) -> None:
    atomic_write_bytes(path, frame.to_csv(index=False).encode())


def project_code_hash(project_root: Path) -> str:
    source = Path(project_root) / "src" / "rheumlens"
    files = sorted(source.rglob("*.py"))
    if not files:
        raise RuntimeError(f"no RheumLens Python files below {source}")
    h = hashlib.sha256()
    for path in files:
        rel = path.relative_to(project_root).as_posix().encode()
        h.update(len(rel).to_bytes(8, "big"))
        h.update(rel)
        content = path.read_bytes()
        h.update(len(content).to_bytes(8, "big"))
        h.update(content)
    return h.hexdigest()


def load_seed_table(path: Path, n_reps: int) -> pd.DataFrame:
    raw = pd.read_csv(path)
    aliases = {"rep": "rep_id", "seed": "rep_seed"}
    raw = raw.rename(columns={k: v for k, v in aliases.items() if k in raw.columns})
    required = {"rep_id", "rep_seed"}
    if not required.issubset(raw.columns):
        raise RuntimeError(f"seed table missing {required - set(raw.columns)}")
    frame = raw[["rep_id", "rep_seed"]].copy()
    frame["rep_id"] = pd.to_numeric(frame.rep_id, errors="raise").astype(int)
    frame["rep_seed"] = pd.to_numeric(frame.rep_seed, errors="raise").astype(np.int64)
    frame = frame.sort_values("rep_id").reset_index(drop=True)
    expected = np.arange(n_reps, dtype=int)
    if len(frame) < n_reps or not np.array_equal(frame.rep_id.iloc[:n_reps].to_numpy(), expected):
        raise RuntimeError(f"seed table must contain exact leading rep_ids 0..{n_reps - 1}")
    frame = frame.iloc[:n_reps].copy()
    if frame.rep_seed.nunique() != n_reps:
        raise RuntimeError("seed collision")
    return frame


def load_json(path: Path) -> dict[str, Any]:
    obj = json.loads(Path(path).read_text())
    if not isinstance(obj, dict):
        raise RuntimeError(f"expected JSON object: {path}")
    return obj


def validate_stage_manifest(path: Path, data_key: str) -> tuple[Path, dict[str, Any]]:
    manifest = load_json(path)
    if manifest.get("asset_key") != data_key:
        raise RuntimeError(f"stage manifest asset_key mismatch: {manifest.get('asset_key')} != {data_key}")
    local_path = Path(manifest["local_path"])
    if not local_path.is_file():
        raise RuntimeError(f"staged data missing: {local_path}")
    actual_size = local_path.stat().st_size
    if actual_size != int(manifest["size"]):
        raise RuntimeError(f"staged size mismatch: {actual_size} != {manifest['size']}")
    actual_hash = sha256_file(local_path)
    if actual_hash != manifest["sha256"]:
        raise RuntimeError(f"staged SHA256 mismatch for {local_path}")
    return local_path, manifest


def validate_stage_against_input(input_manifest_path: Path, stage_manifest: dict[str, Any], data_key: str) -> None:
    manifest = load_json(input_manifest_path)
    assets = manifest.get("assets", {})
    if data_key not in assets:
        raise RuntimeError(f"input manifest has no asset {data_key}")
    asset = assets[data_key]
    if int(asset["size"]) != int(stage_manifest["size"]) or asset["sha256"] != stage_manifest["sha256"]:
        raise RuntimeError("stage/input manifest identity mismatch")


def results_frame(rows: dict[int, dict[str, Any]]) -> pd.DataFrame:
    columns = ["rep_id", "rep_seed", "auc", "runtime_sec"]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame([rows[k] for k in sorted(rows)])[columns]


def read_checkpoint(path: Path, seed_table: pd.DataFrame) -> dict[int, dict[str, Any]]:
    if not Path(path).exists():
        return {}
    frame = pd.read_csv(path)
    required = {"rep_id", "rep_seed", "auc", "runtime_sec"}
    if set(frame.columns) != required:
        raise RuntimeError(f"checkpoint columns mismatch at {path}: {list(frame.columns)}")
    if frame.rep_id.duplicated().any():
        raise RuntimeError(f"duplicate rep_id in {path}")
    expected_seeds = seed_table.set_index("rep_id").rep_seed.to_dict()
    rows: dict[int, dict[str, Any]] = {}
    for row in frame.itertuples(index=False):
        rep_id = int(row.rep_id)
        if rep_id not in expected_seeds or int(row.rep_seed) != int(expected_seeds[rep_id]):
            raise RuntimeError(f"seed mismatch for rep {rep_id} in {path}")
        if not np.isfinite(float(row.auc)):
            raise RuntimeError(f"nonfinite AUC for rep {rep_id} in {path}")
        rows[rep_id] = {
            "rep_id": rep_id,
            "rep_seed": int(row.rep_seed),
            "auc": float(row.auc),
            "runtime_sec": float(row.runtime_sec),
        }
    return rows


def merge_checkpoints(*sources: dict[int, dict[str, Any]]) -> dict[int, dict[str, Any]]:
    merged: dict[int, dict[str, Any]] = {}
    for source in sources:
        for rep_id, row in source.items():
            if rep_id in merged:
                old = merged[rep_id]
                if int(old["rep_seed"]) != int(row["rep_seed"]) or not np.isclose(
                    float(old["auc"]), float(row["auc"]), rtol=0, atol=1e-12
                ):
                    raise RuntimeError(f"conflicting checkpoint values for rep {rep_id}")
            merged[rep_id] = row
    return merged


def assert_complete(frame: pd.DataFrame, seed_table: pd.DataFrame, n_reps: int) -> None:
    if len(frame) != n_reps:
        raise AssertionError(f"row count {len(frame)} != {n_reps}")
    if frame.rep_id.duplicated().any():
        raise AssertionError("duplicate rep_id")
    ordered = frame.sort_values("rep_id").reset_index(drop=True)
    if not np.array_equal(ordered.rep_id.to_numpy(dtype=int), np.arange(n_reps)):
        raise AssertionError("rep_id set is not exact 0..n_reps-1")
    if not np.array_equal(
        ordered.rep_seed.to_numpy(dtype=np.int64), seed_table.rep_seed.to_numpy(dtype=np.int64)
    ):
        raise AssertionError("seed table mismatch")
    if not np.isfinite(ordered.auc.to_numpy(dtype=float)).all():
        raise AssertionError("nonfinite AUC")


def make_identity(args: argparse.Namespace, stage_manifest: dict[str, Any]) -> dict[str, Any]:
    runner_hash = sha256_file(Path(__file__).resolve())
    return {
        "schema_version": 1,
        "method": args.method,
        "cohort": args.cohort,
        "estimand": args.estimand,
        "n_reps": int(args.n_reps),
        "base_seed": int(args.base_seed),
        "observed_auc": float(args.observed_auc),
        "test_delay_sec": float(args.test_delay_sec),
        "runner_sha256": runner_hash,
        "project_code_sha256": project_code_hash(Path(args.project_root)),
        "config_sha256": sha256_file(Path(args.config)),
        "folds_sha256": sha256_file(Path(args.folds)),
        "seed_table_sha256": sha256_file(Path(args.seed_table)),
        "input_manifest_sha256": sha256_file(Path(args.input_manifest)),
        "stage_asset_key": stage_manifest["asset_key"],
        "stage_asset_sha256": stage_manifest["sha256"],
    }


def establish_identity(run_dir: Path, args: argparse.Namespace, stage_manifest: dict[str, Any]) -> dict[str, Any]:
    identity = make_identity(args, stage_manifest)
    path = run_dir / "state.json"
    if args.resume:
        if not path.exists():
            raise RuntimeError("cannot resume without state.json")
        old = load_json(path).get("identity")
        if old != identity:
            raise RuntimeError(
                "resume identity mismatch\nold=" + json.dumps(old, sort_keys=True) +
                "\nnew=" + json.dumps(identity, sort_keys=True)
            )
    else:
        if path.exists():
            raise RuntimeError(f"fresh run refuses existing state: {path}")
        atomic_write_json(
            path,
            {
                "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "host": os.uname().nodename,
                "identity": identity,
            },
        )
    return identity


def validate_oof(frame: pd.DataFrame, expected_donors: int) -> None:
    required = {"donor_id", "y_true", "score", "status"}
    if not required.issubset(frame.columns):
        raise RuntimeError(f"OOF missing columns: {required - set(frame.columns)}")
    if len(frame) != expected_donors or frame.donor_id.astype(str).nunique() != expected_donors:
        raise RuntimeError(
            f"OOF donor integrity failure: rows={len(frame)} unique={frame.donor_id.nunique()} expected={expected_donors}"
        )
    if not frame.status.eq("SUCCESS").all():
        counts = frame.status.value_counts(dropna=False).to_dict()
        raise RuntimeError(f"OOF contains failed rows: {counts}")
    if not np.isfinite(pd.to_numeric(frame.score, errors="coerce")).all():
        raise RuntimeError("OOF contains nonfinite scores")
    if frame.y_true.nunique() != 2:
        raise RuntimeError("pooled OOF does not contain both classes")


def worker_task(task: tuple[int, int]) -> dict[str, Any]:
    rep_id, rep_seed = map(int, task)
    started = time.perf_counter()
    try:
        if G.get("test_delay_sec", 0.0):
            time.sleep(float(G["test_delay_sec"]))
        from rheumlens.evaluation.permutation import permute_donor_labels
        from rheumlens.evaluation.engine import run_fixed_oof
        from sklearn.metrics import roc_auc_score

        permuted = permute_donor_labels(G["data"], rep_seed)
        frame, _diagnostics = run_fixed_oof(
            permuted,
            G["folds"],
            G["method"],
            cohort=G["cohort"],
            estimand=G["estimand"] + "_permutation",
            random_state=rep_seed,
        )
        validate_oof(frame, G["expected_donors"])
        auc = float(roc_auc_score(frame.y_true, frame.score))
        if not np.isfinite(auc):
            raise RuntimeError("nonfinite pooled OOF AUC")
        return {
            "ok": True,
            "rep_id": rep_id,
            "rep_seed": rep_seed,
            "auc": auc,
            "runtime_sec": time.perf_counter() - started,
        }
    except BaseException as exc:
        return {
            "ok": False,
            "rep_id": rep_id,
            "rep_seed": rep_seed,
            "error_type": type(exc).__name__,
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "runtime_sec": time.perf_counter() - started,
        }


def worker_init() -> None:
    # The coordinator owns graceful shutdown/checkpointing. Pool termination should
    # kill workers directly rather than running the coordinator signal handler.
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def request_stop(signum: int, _frame: Any) -> None:
    global STOP_REQUESTED
    STOP_REQUESTED = True
    raise GracefulStop(f"received signal {signum}")


def append_error(path: Path, result: dict[str, Any]) -> None:
    # Sole coordinator writer; no cross-process append race.
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(result, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def load_context(args: argparse.Namespace, staged_path: Path) -> None:
    project_root = Path(args.project_root).resolve()
    sys.path.insert(0, str(project_root / "src"))
    from rheumlens.data.io import load_npz_dataset
    from rheumlens.evaluation.splits import load_folds
    from rheumlens.registry import build_method

    config = yaml.safe_load(Path(args.config).read_text())
    defaults = config.get("method_defaults", {})
    registered = build_method(args.method, defaults, None, seed=int(args.base_seed))
    if registered.data_key != args.data_key:
        raise RuntimeError(f"registry data_key={registered.data_key} but job data_key={args.data_key}")
    data = load_npz_dataset(staged_path)
    folds = load_folds(args.folds)
    fold_donors = {str(x) for fold in folds for x in (*fold.train_donors, *fold.test_donors)}
    data_donors = set(map(str, data.donors))
    if fold_donors != data_donors:
        raise RuntimeError(
            f"fold/data donor set mismatch: folds={len(fold_donors)} data={len(data_donors)}"
        )
    G.update(
        data=data,
        folds=folds,
        method=registered.method,
        cohort=args.cohort,
        estimand=args.estimand,
        expected_donors=len(data_donors),
        test_delay_sec=float(args.test_delay_sec),
    )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True, choices=["donor_expression_pca", "kme_multiscale@scgpt"])
    ap.add_argument("--data-key", required=True, choices=["lognorm", "scgpt"])
    ap.add_argument("--cohort", default="GSE174188_CD4")
    ap.add_argument("--estimand", default="matched_500_cells_per_donor")
    ap.add_argument("--n-workers", type=int, required=True)
    ap.add_argument("--n-reps", type=int, required=True)
    ap.add_argument("--base-seed", type=int, default=20260619)
    ap.add_argument("--observed-auc", type=float, required=True)
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--folds", required=True)
    ap.add_argument("--seed-table", required=True)
    ap.add_argument("--input-manifest", required=True)
    ap.add_argument("--stage-manifest", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--local-checkpoint-dir", required=True)
    ap.add_argument("--durable-checkpoint-dir", required=True)
    ap.add_argument("--local-every", type=int, default=5)
    ap.add_argument("--durable-every", type=int, default=10)
    ap.add_argument("--stop-after-reps", type=int)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--test-delay-sec", type=float, default=0.0, help=argparse.SUPPRESS)
    return ap.parse_args()


def run(args: argparse.Namespace) -> int:
    if args.n_workers < 1 or args.n_workers > 64:
        raise RuntimeError("n_workers must be 1..64")
    if args.n_reps < 1 or args.n_reps > 1000:
        raise RuntimeError("n_reps must be 1..1000")
    if args.local_every < 1 or args.durable_every < 1:
        raise RuntimeError("checkpoint intervals must be positive")

    run_dir = Path(args.run_dir).resolve()
    local_dir = Path(args.local_checkpoint_dir).resolve()
    durable_dir = Path(args.durable_checkpoint_dir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    local_dir.mkdir(parents=True, exist_ok=True)
    durable_dir.mkdir(parents=True, exist_ok=True)
    errors_path = run_dir / "errors.jsonl"

    seed_table = load_seed_table(Path(args.seed_table), args.n_reps)
    staged_path, stage_manifest = validate_stage_manifest(Path(args.stage_manifest), args.data_key)
    validate_stage_against_input(Path(args.input_manifest), stage_manifest, args.data_key)
    identity = establish_identity(run_dir, args, stage_manifest)
    load_context(args, staged_path)

    durable_rows = read_checkpoint(durable_dir / "checkpoint.csv", seed_table)
    local_rows = read_checkpoint(local_dir / "checkpoint.csv", seed_table)
    rows = merge_checkpoints(durable_rows, local_rows)
    tasks = [
        (int(row.rep_id), int(row.rep_seed))
        for row in seed_table.itertuples(index=False)
        if int(row.rep_id) not in rows
    ]
    initial_completed = len(rows)

    atomic_write_json(
        run_dir / "runtime_state.json",
        {
            "status": "RUNNING",
            "pid": os.getpid(),
            "n_workers": args.n_workers,
            "n_reps": args.n_reps,
            "initial_completed": initial_completed,
            "identity_sha256": sha256_json(identity),
            "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
    )

    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)
    pool: Any = None
    stopped = False
    failures = 0
    since_start = 0
    try:
        if tasks:
            ctx = mp.get_context("fork")
            pool = ctx.Pool(processes=args.n_workers, initializer=worker_init)
            iterator = pool.imap_unordered(worker_task, tasks, chunksize=1)
            for result in iterator:
                if not result["ok"]:
                    failures += 1
                    append_error(errors_path, result)
                    pool.terminate()
                    pool.join()
                    pool = None
                    raise RuntimeError(f"rep {result['rep_id']} failed; see {errors_path}")
                rep_id = int(result["rep_id"])
                rows[rep_id] = {
                    "rep_id": rep_id,
                    "rep_seed": int(result["rep_seed"]),
                    "auc": float(result["auc"]),
                    "runtime_sec": float(result["runtime_sec"]),
                }
                since_start += 1
                if since_start % args.local_every == 0:
                    atomic_write_csv(local_dir / "checkpoint.csv", results_frame(rows))
                if since_start % args.durable_every == 0:
                    atomic_write_csv(durable_dir / "checkpoint.csv", results_frame(rows))
                    atomic_write_json(
                        run_dir / "progress.json",
                        {
                            "completed": len(rows),
                            "requested": args.n_reps,
                            "failures": failures,
                            "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                        },
                    )
                if args.stop_after_reps is not None and since_start >= args.stop_after_reps:
                    stopped = True
                    break
            if stopped:
                pool.terminate()
            else:
                pool.close()
            pool.join()
            pool = None
    except GracefulStop:
        stopped = True
        if pool is not None:
            pool.terminate()
            pool.join()
            pool = None
    finally:
        atomic_write_csv(local_dir / "checkpoint.csv", results_frame(rows))
        atomic_write_csv(durable_dir / "checkpoint.csv", results_frame(rows))

    if stopped or len(rows) < args.n_reps:
        atomic_write_json(
            run_dir / "runtime_state.json",
            {
                "status": "INTERRUPTED_RECOVERABLE",
                "completed": len(rows),
                "requested": args.n_reps,
                "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        )
        return 75

    frame = results_frame(rows).sort_values("rep_id").reset_index(drop=True)
    assert_complete(frame, seed_table, args.n_reps)
    atomic_write_csv(run_dir / "permutation_results.csv", frame)
    null_auc = frame.auc.to_numpy(dtype=float)
    empirical_p = float((1 + np.sum(null_auc >= args.observed_auc)) / (1 + args.n_reps))
    summary = {
        "schema_version": 1,
        "cohort": args.cohort,
        "method_id": args.method,
        "observed_auc": float(args.observed_auc),
        "n_requested": args.n_reps,
        "n_finite": args.n_reps,
        "empirical_p": empirical_p,
        "null_mean": float(np.mean(null_auc)),
        "null_q025": float(np.quantile(null_auc, 0.025)),
        "null_q975": float(np.quantile(null_auc, 0.975)),
        "identity_sha256": sha256_json(identity),
    }
    atomic_write_json(run_dir / "summary.json", summary)
    atomic_write_json(
        run_dir / "runtime_state.json",
        {
            "status": "COMPLETED_TECHNICAL",
            "completed": args.n_reps,
            "requested": args.n_reps,
            "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
    )
    return 0


def main() -> int:
    args = parse_args()
    try:
        return run(args)
    except BaseException as exc:
        run_dir = Path(getattr(args, "run_dir", "."))
        run_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_json(
            run_dir / "fatal_error.json",
            {
                "error_type": type(exc).__name__,
                "error": repr(exc),
                "traceback": traceback.format_exc(),
                "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        )
        print(traceback.format_exc(), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
