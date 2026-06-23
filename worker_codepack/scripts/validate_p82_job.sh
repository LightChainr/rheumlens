#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 METHOD RUN_ID" >&2
  exit 2
fi
METHOD="$1"
RUN_ID="$2"
ROOT="${RHEUMLENS_ROOT:-/autodl-fs/data/rheumlens}"
PACK="${RHEUMLENS_WORKER_PACK:-$ROOT/supervisor_control/rheumlens_worker_codepack_v2}"
RUN="$ROOT/results/P8_2_formal_permutation_runs/$RUN_ID"
SEEDS="$PACK/frozen/GSE174188_CD4_formal_seed_table_20260619.csv"
case "$METHOD" in
  donor_expression_pca) OBSERVED=0.9856590597331338 ;;
  kme_multiscale@scgpt) OBSERVED=0.977865070457663 ;;
  *) echo "Unsupported method: $METHOD" >&2; exit 2 ;;
esac
PYTHON="${RHEUMLENS_PYTHON:-micromamba run -n rheumlens-core python}"
read -r -a PYTHON_CMD <<< "$PYTHON"
"${PYTHON_CMD[@]}" "$PACK/runner/validate_formal_run.py" \
  --run-dir "$RUN" --method "$METHOD" --seed-table "$SEEDS" \
  --expected-observed "$OBSERVED" --n-reps 1000 | tee "$RUN/validation.log"
