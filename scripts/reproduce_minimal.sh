#!/usr/bin/env bash
set -euo pipefail
python -m pip install -e .
pytest -q
python scripts/mac_p9_p10.py --help >/dev/null || true
