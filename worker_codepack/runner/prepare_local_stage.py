#!/usr/bin/env python3
"""Copy one manifest asset to local storage with streaming SHA256 validation."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(obj, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-manifest", required=True)
    ap.add_argument("--asset-key", required=True)
    ap.add_argument("--stage-dir", required=True)
    args = ap.parse_args()

    manifest = json.loads(Path(args.input_manifest).read_text())
    asset = manifest["assets"][args.asset_key]
    source = Path(asset["path"])
    if not source.is_file():
        raise RuntimeError(f"source missing: {source}")
    expected_size = int(asset["size"])
    expected_hash = asset["sha256"]
    stage_dir = Path(args.stage_dir).resolve()
    stage_dir.mkdir(parents=True, exist_ok=True)
    target = stage_dir / source.name
    free = shutil.disk_usage(stage_dir).free
    if free < expected_size + max(5 * 1024**3, expected_size // 5):
        raise RuntimeError(f"insufficient local disk: free={free} required>{expected_size}")

    valid = target.is_file() and target.stat().st_size == expected_size
    if valid:
        valid = sha256_file(target) == expected_hash
    if not valid:
        partial = target.with_suffix(target.suffix + ".partial")
        h = hashlib.sha256()
        copied = 0
        with source.open("rb") as src, partial.open("wb") as dst:
            for chunk in iter(lambda: src.read(8 * 1024 * 1024), b""):
                dst.write(chunk)
                h.update(chunk)
                copied += len(chunk)
            dst.flush()
            os.fsync(dst.fileno())
        if copied != expected_size or h.hexdigest() != expected_hash:
            raise RuntimeError("staged copy size/SHA256 mismatch")
        os.replace(partial, target)

    out = {
        "asset_key": args.asset_key,
        "source_path": str(source),
        "local_path": str(target),
        "size": expected_size,
        "sha256": expected_hash,
    }
    out_path = stage_dir / f"stage_manifest.{args.asset_key}.json"
    atomic_json(out_path, out)
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
