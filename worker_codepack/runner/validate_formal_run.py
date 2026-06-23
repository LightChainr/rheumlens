#!/usr/bin/env python3
"""Strict validator for one completed v3 formal permutation run."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--method", required=True)
    ap.add_argument("--seed-table", required=True)
    ap.add_argument("--expected-observed", type=float, required=True)
    ap.add_argument("--n-reps", type=int, default=1000)
    args = ap.parse_args()
    run = Path(args.run_dir)
    result_path = run / "permutation_results.csv"
    summary_path = run / "summary.json"
    state_path = run / "state.json"
    runtime_path = run / "runtime_state.json"
    exit_path = run / "exit_code.txt"
    for path in [result_path, summary_path, state_path, runtime_path, exit_path]:
        if not path.is_file():
            raise RuntimeError(f"missing required file: {path}")
    if int(exit_path.read_text().strip()) != 0:
        raise RuntimeError("runner exit code is nonzero")
    if (run / "fatal_error.json").exists():
        raise RuntimeError("fatal_error.json exists")
    if (run / "errors.jsonl").exists() and (run / "errors.jsonl").stat().st_size:
        raise RuntimeError("errors.jsonl is nonempty")
    runtime = json.loads(runtime_path.read_text())
    if runtime.get("status") != "COMPLETED_TECHNICAL":
        raise RuntimeError(f"runtime status is {runtime.get('status')}")
    state = json.loads(state_path.read_text())
    identity = state.get("identity", {})
    if identity.get("method") != args.method or int(identity.get("n_reps", -1)) != args.n_reps:
        raise RuntimeError("state identity mismatch")
    frame = pd.read_csv(result_path)
    required = ["rep_id", "rep_seed", "auc", "runtime_sec"]
    if list(frame.columns) != required or len(frame) != args.n_reps or frame.rep_id.duplicated().any():
        raise RuntimeError("result schema/count/duplicate failure")
    frame = frame.sort_values("rep_id").reset_index(drop=True)
    if not np.array_equal(frame.rep_id.to_numpy(dtype=int), np.arange(args.n_reps)):
        raise RuntimeError("rep_id set mismatch")
    if not np.isfinite(frame.auc.to_numpy(dtype=float)).all():
        raise RuntimeError("nonfinite AUC")
    seeds = pd.read_csv(args.seed_table)
    if not np.array_equal(frame.rep_seed.to_numpy(np.int64), seeds.rep_seed.iloc[:args.n_reps].to_numpy(np.int64)):
        raise RuntimeError("seed sequence mismatch")
    summary = json.loads(summary_path.read_text())
    if summary.get("method_id") != args.method or summary.get("n_finite") != args.n_reps:
        raise RuntimeError("summary method/count mismatch")
    observed = float(summary["observed_auc"])
    if not np.isclose(observed, args.expected_observed, atol=1e-12, rtol=0):
        raise RuntimeError("observed AUC mismatch")
    p = float((1 + np.sum(frame.auc.to_numpy() >= observed)) / (1 + args.n_reps))
    if not np.isclose(p, float(summary["empirical_p"]), atol=1e-15, rtol=0):
        raise RuntimeError("empirical P mismatch")

    files = sorted(
        p for p in run.rglob("*")
        if p.is_file() and p.name not in {"MANIFEST_SHA256.txt", "VALIDATION.json"}
    )
    lines = [f"{sha256_file(path)}  {path.relative_to(run).as_posix()}" for path in files]
    manifest_tmp = run / "MANIFEST_SHA256.txt.tmp"
    manifest_tmp.write_text("\n".join(lines) + "\n")
    os.replace(manifest_tmp, run / "MANIFEST_SHA256.txt")
    report = {
        "status": "ACCEPT_RECOMMENDED",
        "method_id": args.method,
        "n_reps": args.n_reps,
        "observed_auc": observed,
        "empirical_p": p,
        "result_sha256": sha256_file(result_path),
        "summary_sha256": sha256_file(summary_path),
    }
    fd, tmp = tempfile.mkstemp(prefix="VALIDATION.json.", suffix=".tmp", dir=run)
    with os.fdopen(fd, "w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, run / "VALIDATION.json")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
