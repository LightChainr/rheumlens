#!/usr/bin/env python3
"""Run serial/parallel, real TERM/resume, bad-hash and resource gates."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd


PACK = Path(__file__).resolve().parents[1]
RUNNER = PACK / "runner/perm_parallel_v3.py"
COMPARE = PACK / "runner/compare_runs.py"


def build(args, run, local, durable, workers, reps, config, resume=False):
    cmd = [
        sys.executable, str(RUNNER), "--method", args.method, "--data-key", args.data_key,
        "--n-workers", str(workers), "--n-reps", str(reps), "--base-seed", "20260619",
        "--observed-auc", str(args.observed_auc), "--project-root", args.project_root,
        "--config", str(config), "--folds", args.folds, "--seed-table", args.seed_table,
        "--input-manifest", args.input_manifest, "--stage-manifest", args.stage_manifest,
        "--run-dir", str(run), "--local-checkpoint-dir", str(local),
        "--durable-checkpoint-dir", str(durable), "--local-every", "1", "--durable-every", "1",
    ]
    if resume:
        cmd.append("--resume")
    if args.test_delay_sec:
        cmd += ["--test-delay-sec", str(args.test_delay_sec)]
    return cmd


def run_logged(cmd, log, expected=(0,)):
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as handle:
        result = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode not in expected:
        raise RuntimeError(f"command failed rc={result.returncode}; see {log}")
    return result.returncode


def compare(a, b, log):
    run_logged([sys.executable, str(COMPARE), str(a), str(b)], log)


def wait_checkpoint(path: Path, minimum: int, timeout: int) -> int:
    deadline = time.time() + timeout
    while time.time() < deadline:
        if path.exists():
            try:
                n = len(pd.read_csv(path))
                if n >= minimum:
                    return n
            except Exception:
                pass
        time.sleep(0.2)
    raise TimeoutError(f"checkpoint did not reach {minimum}: {path}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True)
    ap.add_argument("--data-key", required=True)
    ap.add_argument("--observed-auc", type=float, required=True)
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--folds", required=True)
    ap.add_argument("--seed-table", required=True)
    ap.add_argument("--input-manifest", required=True)
    ap.add_argument("--stage-manifest", required=True)
    ap.add_argument("--gate-root", required=True)
    ap.add_argument("--resource-workers", type=int, default=16)
    ap.add_argument("--resource-reps", type=int, default=50)
    ap.add_argument("--test-delay-sec", type=float, default=0.0, help=argparse.SUPPRESS)
    args = ap.parse_args()
    root = Path(args.gate_root).resolve()
    if root.exists():
        raise RuntimeError(f"gate root already exists: {root}")
    root.mkdir(parents=True)

    serial, parallel = root / "serial8", root / "parallel8"
    run_logged(build(args, serial, root/"local_serial", serial/"checkpoints", 1, 8, args.config), serial/"controller.log")
    run_logged(build(args, parallel, root/"local_parallel", parallel/"checkpoints", 4, 8, args.config), parallel/"controller.log")
    compare(serial, parallel, root/"SERIAL_PARALLEL_COMPARE.log")

    interrupted, reference = root / "term_resume20", root / "reference20"
    cmd = build(args, interrupted, root/"local_term", interrupted/"checkpoints", 8, 20, args.config)
    interrupted.mkdir()
    with (interrupted/"controller.part1.log").open("w") as handle:
        proc = subprocess.Popen(cmd, stdout=handle, stderr=subprocess.STDOUT)
        wait_checkpoint(interrupted/"checkpoints/checkpoint.csv", 3, 1800)
        proc.send_signal(signal.SIGTERM)
        rc = proc.wait(timeout=300)
    if rc != 75:
        raise RuntimeError(f"TERM run expected rc75, got {rc}")
    run_logged(build(args, interrupted, root/"local_term", interrupted/"checkpoints", 8, 20, args.config, True), interrupted/"controller.resume.log")
    run_logged(build(args, reference, root/"local_reference", reference/"checkpoints", 8, 20, args.config), reference/"controller.log")
    compare(interrupted, reference, root/"TERM_RESUME_COMPARE.log")

    bad = root / "bad_hash_resume"
    bad_config = root / "config.bad_hash_test.yaml"
    shutil.copy2(args.config, bad_config)
    cmd = build(args, bad, root/"local_bad", bad/"checkpoints", 4, 8, bad_config)
    cmd += ["--stop-after-reps", "2"]
    run_logged(cmd, bad/"controller.part1.log", expected=(75,))
    with bad_config.open("a") as handle:
        handle.write("\nresume_identity_test_change: true\n")
    rc = run_logged(build(args, bad, root/"local_bad", bad/"checkpoints", 4, 8, bad_config, True), bad/"controller.rejected.log", expected=(1,))
    if rc != 1:
        raise RuntimeError("bad-hash resume was not rejected")

    resource = root / f"resource_{args.resource_workers}w_{args.resource_reps}r"
    run_logged(
        build(args, resource, root/"local_resource", resource/"checkpoints", args.resource_workers, args.resource_reps, args.config),
        resource/"controller.log",
    )
    report = {
        "status": "PASS_ALL_GATES",
        "method": args.method,
        "serial_parallel": "PASS",
        "term_resume": "PASS",
        "bad_hash_resume_rejected": True,
        "resource_workers": args.resource_workers,
        "resource_reps": args.resource_reps,
    }
    (root/"GATE_REPORT.json").write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
