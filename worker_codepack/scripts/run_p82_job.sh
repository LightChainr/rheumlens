#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 4 ]; then
  echo "Usage: $0 METHOD N_WORKERS RUN_ID MODE[fresh|resume]" >&2
  exit 2
fi

METHOD="$1"
N_WORKERS="$2"
RUN_ID="$3"
MODE="$4"
ROOT="${RHEUMLENS_ROOT:-/autodl-fs/data/rheumlens}"
PACK="${RHEUMLENS_WORKER_PACK:-$ROOT/supervisor_control/rheumlens_worker_codepack_v2}"
PROJECT="$ROOT/code/RheumLens_A800_full_code"
CONFIG="$ROOT/configs/project.p6.gse174188.v1.yaml"
FOLDS="$ROOT/splits/authoritative_primary/GSE174188_CD4.csv"
INPUT_MANIFEST="$ROOT/results/P8_2_formal_permutation_runs/manifests/GSE174188_CD4_P8_2_input_manifest.json"
SEEDS="$PACK/frozen/GSE174188_CD4_formal_seed_table_20260619.csv"
RUN="$ROOT/results/P8_2_formal_permutation_runs/$RUN_ID"
LOCAL_BASE="/root/autodl-tmp/rheumlens_p82"
LOCAL_STAGE="$LOCAL_BASE/stage"
LOCAL_CKPT="$LOCAL_BASE/checkpoints/$RUN_ID"
DURABLE_CKPT="$RUN/checkpoints"

case "$METHOD" in
  donor_expression_pca)
    DATA_KEY=lognorm
    OBSERVED_AUC=0.9856590597331338
    ;;
  kme_multiscale@scgpt)
    DATA_KEY=scgpt
    OBSERVED_AUC=0.977865070457663
    ;;
  *) echo "Unsupported method: $METHOD" >&2; exit 2 ;;
esac

test "$MODE" = fresh -o "$MODE" = resume
test -f "$INPUT_MANIFEST"
test -f "$SEEDS"
test -f "$CONFIG"
test -f "$FOLDS"
test -d "$PROJECT/src/rheumlens"
mkdir -p "$RUN" "$LOCAL_CKPT" "$DURABLE_CKPT" "$LOCAL_STAGE"

export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

PYTHON="${RHEUMLENS_PYTHON:-micromamba run -n rheumlens-core python}"
read -r -a PYTHON_CMD <<< "$PYTHON"

STAGE_MANIFEST="$("${PYTHON_CMD[@]}" "$PACK/runner/prepare_local_stage.py" \
  --input-manifest "$INPUT_MANIFEST" --asset-key "$DATA_KEY" --stage-dir "$LOCAL_STAGE")"

{
  date -u
  hostname
  printf 'method=%s\nn_workers=%s\nrun_id=%s\nmode=%s\n' "$METHOD" "$N_WORKERS" "$RUN_ID" "$MODE"
  printf 'project=%s\nconfig=%s\nfolds=%s\ninput_manifest=%s\nseed_table=%s\nstage_manifest=%s\n' \
    "$PROJECT" "$CONFIG" "$FOLDS" "$INPUT_MANIFEST" "$SEEDS" "$STAGE_MANIFEST"
  sha256sum "$PACK/runner/perm_parallel_v3.py" "$CONFIG" "$FOLDS" "$INPUT_MANIFEST" "$SEEDS"
} > "$RUN/command.txt"

{
  "${PYTHON_CMD[@]}" --version
  "${PYTHON_CMD[@]}" - <<'PY'
import numpy, pandas, scipy, sklearn, yaml
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("scipy", scipy.__version__)
print("sklearn", sklearn.__version__)
print("pyyaml", yaml.__version__)
PY
  printf 'CUDA_VISIBLE_DEVICES=%s\nOMP_NUM_THREADS=%s\nMKL_NUM_THREADS=%s\nOPENBLAS_NUM_THREADS=%s\n' \
    "$CUDA_VISIBLE_DEVICES" "$OMP_NUM_THREADS" "$MKL_NUM_THREADS" "$OPENBLAS_NUM_THREADS"
} > "$RUN/environment.txt" 2>&1

RESUME_ARG=()
if [ "$MODE" = resume ]; then RESUME_ARG=(--resume); fi

set +e
"${PYTHON_CMD[@]}" "$PACK/runner/perm_parallel_v3.py" \
  --method "$METHOD" --data-key "$DATA_KEY" \
  --n-workers "$N_WORKERS" --n-reps 1000 --base-seed 20260619 \
  --observed-auc "$OBSERVED_AUC" \
  --project-root "$PROJECT" --config "$CONFIG" --folds "$FOLDS" \
  --seed-table "$SEEDS" --input-manifest "$INPUT_MANIFEST" \
  --stage-manifest "$STAGE_MANIFEST" --run-dir "$RUN" \
  --local-checkpoint-dir "$LOCAL_CKPT" --durable-checkpoint-dir "$DURABLE_CKPT" \
  "${RESUME_ARG[@]}" 2>&1 | tee "$RUN/run.${MODE}.log"
RC=${PIPESTATUS[0]}
set -e
echo "$RC" > "$RUN/exit_code.txt"
exit "$RC"
