#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_BIN" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN=python3
  elif command -v python >/dev/null 2>&1; then
    PYTHON_BIN=python
  else
    echo "No python3 or python executable found." >&2
    exit 1
  fi
fi

"$PYTHON_BIN" -m pip install -e .
"$PYTHON_BIN" -m pytest -q
"$PYTHON_BIN" scripts/mac_p9_p10.py --help >/dev/null || true
bash scripts/verify_manifests.sh

"$PYTHON_BIN" - <<'PY'
from pathlib import Path

required = [
    "figures/MANIFEST_SHA256.tsv",
    "supplementary_tables/MANIFEST_SHA256.tsv",
    "docs/DATA_AND_CODE_AVAILABILITY.md",
    "environment.yml",
    "Dockerfile",
    "Singularity.def",
]
missing = [p for p in required if not Path(p).exists()]
if missing:
    raise SystemExit(f"Missing required reproducibility assets: {missing}")
print("Minimal reproduction checks passed.")
PY
