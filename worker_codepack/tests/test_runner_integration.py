from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


PACK = Path(__file__).resolve().parents[1]
RUNNER = PACK / "runner/perm_parallel_v3.py"
PROJECT = PACK / "reference_project"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_fixture(tmp_path: Path) -> dict[str, Path]:
    sys.path.insert(0, str(PROJECT / "src"))
    from rheumlens.data.io import save_npz_dataset
    from rheumlens.evaluation.splits import make_stratified_folds, save_folds
    from rheumlens.types import CellDataset

    rng = np.random.default_rng(7)
    n_donors, cells, features = 20, 5, 8
    donors = np.repeat([f"D{i:02d}" for i in range(n_donors)], cells)
    donor_y = np.asarray([0, 1] * (n_donors // 2), dtype=int)
    y = np.repeat(donor_y, cells)
    X = rng.normal(size=(len(donors), features)).astype(np.float32)
    X[:, :2] += y[:, None] * 0.6
    data = CellDataset(
        X=X,
        cell_ids=np.asarray([f"C{i:04d}" for i in range(len(donors))]),
        donor_ids=donors,
        y=y,
        feature_names=np.asarray([f"G{i}" for i in range(features)]),
        name="toy",
    )
    staged = tmp_path / "lognorm.npz"
    save_npz_dataset(staged, data)
    folds = make_stratified_folds(np.asarray([f"D{i:02d}" for i in range(n_donors)]), donor_y, 5, 11, "toy")
    folds_path = tmp_path / "folds.csv"
    save_folds(folds_path, folds)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump({"method_defaults": {"n_hvg": 8, "donor_pca_dim": 3, "estimator_C": 1.0}}))
    seeds = pd.DataFrame({"rep_id": np.arange(64), "rep_seed": np.arange(101, 165)})
    seed_path = tmp_path / "seeds.csv"
    seeds.to_csv(seed_path, index=False)
    input_manifest = tmp_path / "input.json"
    input_manifest.write_text(json.dumps({"assets": {"lognorm": {"path": str(staged), "size": staged.stat().st_size, "sha256": file_hash(staged)}}}))
    stage_manifest = tmp_path / "stage.json"
    stage_manifest.write_text(json.dumps({"asset_key": "lognorm", "local_path": str(staged), "size": staged.stat().st_size, "sha256": file_hash(staged)}))
    return {"data": staged, "folds": folds_path, "config": config_path, "seeds": seed_path, "input": input_manifest, "stage": stage_manifest}


def command(paths: dict[str, Path], run: Path, local: Path, durable: Path, workers: int, resume=False, stop=None):
    cmd = [
        sys.executable, str(RUNNER), "--method", "donor_expression_pca", "--data-key", "lognorm",
        "--n-workers", str(workers), "--n-reps", "8", "--base-seed", "20260619",
        "--observed-auc", "0.75", "--project-root", str(PROJECT), "--config", str(paths["config"]),
        "--folds", str(paths["folds"]), "--seed-table", str(paths["seeds"]),
        "--input-manifest", str(paths["input"]), "--stage-manifest", str(paths["stage"]),
        "--run-dir", str(run), "--local-checkpoint-dir", str(local),
        "--durable-checkpoint-dir", str(durable), "--local-every", "1", "--durable-every", "1",
    ]
    if resume:
        cmd.append("--resume")
    if stop is not None:
        cmd += ["--stop-after-reps", str(stop)]
    return cmd


def test_serial_parallel_and_resume(tmp_path: Path):
    paths = make_fixture(tmp_path)
    serial = tmp_path / "serial"
    parallel = tmp_path / "parallel"
    resumed = tmp_path / "resumed"
    assert subprocess.run(command(paths, serial, tmp_path / "ls", serial / "checkpoints", 1)).returncode == 0
    assert subprocess.run(command(paths, parallel, tmp_path / "lp", parallel / "checkpoints", 4)).returncode == 0
    a = pd.read_csv(serial / "permutation_results.csv")
    b = pd.read_csv(parallel / "permutation_results.csv")
    pd.testing.assert_frame_equal(a[["rep_id", "rep_seed", "auc"]], b[["rep_id", "rep_seed", "auc"]], check_exact=True)

    first = subprocess.run(command(paths, resumed, tmp_path / "lr", resumed / "checkpoints", 4, stop=3))
    assert first.returncode == 75
    assert subprocess.run(command(paths, resumed, tmp_path / "lr", resumed / "checkpoints", 4, resume=True)).returncode == 0
    c = pd.read_csv(resumed / "permutation_results.csv")
    pd.testing.assert_frame_equal(a[["rep_id", "rep_seed", "auc"]], c[["rep_id", "rep_seed", "auc"]], check_exact=True)


def test_bad_config_resume_rejected(tmp_path: Path):
    paths = make_fixture(tmp_path)
    run = tmp_path / "resume_bad"
    assert subprocess.run(command(paths, run, tmp_path / "local", run / "checkpoints", 2, stop=2)).returncode == 75
    paths["config"].write_text(paths["config"].read_text() + "\nchanged: true\n")
    assert subprocess.run(command(paths, run, tmp_path / "local", run / "checkpoints", 2, resume=True)).returncode != 0
