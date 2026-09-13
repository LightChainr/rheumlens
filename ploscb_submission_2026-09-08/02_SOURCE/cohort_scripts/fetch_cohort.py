"""Download registered cohort .h5ad files from CZ CELLxGENE Discover.

Resumable (HTTP Range), checksum-free by design: CELLxGENE does not publish
per-asset digests, so we record size + a local SHA256 after download instead.

    python scripts/cohorts/fetch_cohort.py --all
    python scripts/cohorts/fetch_cohort.py --cohort COVID_REN
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
REGISTRY = REPO / "cohorts" / "registry.yaml"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url: str, dest: Path, expect_gb: float) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    have = dest.stat().st_size if dest.exists() else 0
    req = urllib.request.Request(url)
    if have:
        req.add_header("Range", f"bytes={have}-")
        print(f"  resuming at {have / 1e9:.2f} GB")

    try:
        resp = urllib.request.urlopen(req, timeout=120)
    except urllib.error.HTTPError as exc:
        if exc.code == 416:          # already complete
            print("  already complete")
            return
        raise

    total = have + int(resp.headers.get("Content-Length", 0))
    mode = "ab" if have and resp.status == 206 else "wb"
    if mode == "wb":
        have = 0

    with dest.open(mode) as fh:
        while chunk := resp.read(8 << 20):
            fh.write(chunk)
            have += len(chunk)
            pct = 100 * have / total if total else 0
            print(f"  {have / 1e9:6.2f} / {total / 1e9:.2f} GB ({pct:5.1f}%)",
                  end="\r", flush=True)
    print()
    got = dest.stat().st_size / 1e9
    if abs(got - expect_gb) / max(expect_gb, 1e-9) > 0.10:
        print(f"  ! size {got:.2f} GB differs from registry {expect_gb:.2f} GB")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--include-existing", action="store_true",
                    help="also download cohorts marked status: existing")
    args = ap.parse_args()

    reg = yaml.safe_load(REGISTRY.read_text())
    data_dir = Path(args.data_dir) if args.data_dir else REPO / "data" / "h5ad"

    if args.all:
        targets = [c for c in reg["cohorts"]
                   if c["status"] == "new" or args.include_existing]
    elif args.cohort:
        targets = [c for c in reg["cohorts"] if c["name"] == args.cohort]
        if not targets:
            raise SystemExit(f"unknown cohort {args.cohort}")
    else:
        raise SystemExit("pass --cohort NAME or --all")

    print(f"{len(targets)} cohort(s), "
          f"{sum(c['h5ad_gb'] for c in targets):.1f} GB total -> {data_dir}")

    manifest = []
    for c in targets:
        dest = data_dir / f"{c['name']}.h5ad"
        print(f"\n[{c['name']}] {c['h5ad_gb']} GB  ({c['citation']})")
        fetch(c["h5ad_url"], dest, c["h5ad_gb"])
        digest = sha256(dest)
        print(f"  sha256 {digest}")
        manifest.append({"cohort": c["name"], "dataset_id": c["dataset_id"],
                         "url": c["h5ad_url"], "path": str(dest),
                         "bytes": dest.stat().st_size, "sha256": digest})

    out = data_dir / "DOWNLOAD_MANIFEST.json"
    prev = json.loads(out.read_text()) if out.exists() else []
    keep = [m for m in prev if m["cohort"] not in {m2["cohort"] for m2 in manifest}]
    out.write_text(json.dumps(keep + manifest, indent=2))
    print(f"\nmanifest -> {out}")


if __name__ == "__main__":
    main()
