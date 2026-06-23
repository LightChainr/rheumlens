#!/usr/bin/env python3
"""Create the frozen P8.2 input manifest with streaming hashes."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/autodl-fs/data/rheumlens")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    root = Path(args.root).resolve()
    paths = {
        "raw_counts": root / "data/processed/GSE174188_CD4/raw_counts.npz",
        "lognorm": root / "data/processed/GSE174188_CD4/lognorm.npz",
        "scgpt": root / "embeddings/scgpt/GSE174188_CD4_v1/GSE174188_CD4_scgpt.npz",
        "folds": root / "splits/authoritative_primary/GSE174188_CD4.csv",
        "config": root / "configs/project.p6.gse174188.v1.yaml",
        "method_summary": root / "results/P6_GSE174188_v1/final_tables/method_summary.csv",
        "all_oof": root / "results/P6_GSE174188_v1/final_tables/all_oof.csv",
        "scgpt_permutation": root / "results/P6_GSE174188_v1/permutation/GSE174188_CD4/scgpt_mean/permutation_auc.csv",
        "scgpt_permutation_summary": root / "results/P6_GSE174188_v1/permutation/GSE174188_CD4/scgpt_mean/summary.json",
    }
    assets = {}
    for key, path in paths.items():
        if not path.is_file():
            raise RuntimeError(f"missing asset {key}: {path}")
        assets[key] = {"path": str(path), "size": path.stat().st_size, "sha256": sha256_file(path)}
    obj = {"schema_version": 1, "root": str(root), "assets": assets}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=output.name + ".", suffix=".tmp", dir=output.parent)
    with os.fdopen(fd, "w") as handle:
        json.dump(obj, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, output)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
