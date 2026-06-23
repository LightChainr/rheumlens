#!/usr/bin/env bash
set -euo pipefail

RL_ROOT="${RL_ROOT:-/autodl-fs/data/rheumlens}"
REPO_ROOT="${REPO_ROOT:-$RL_ROOT/code/RheumLens_v3_code}"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$RL_ROOT/tools/micromamba_root}"
export MAMBA_ROOT_PREFIX

mkdir -p "$RL_ROOT"/{code,envs,tools,data/{raw,processed,external,manifests},metadata,immutable,splits,models/{scgpt,geneformer,references},tokenized,embeddings,representations,results,logs,tmp,cache/{huggingface,torch,pip}}

cat > "$RL_ROOT/env.sh" <<EOF
export RL_ROOT=$RL_ROOT
export HF_HOME=$RL_ROOT/cache/huggingface
export HUGGINGFACE_HUB_CACHE=$RL_ROOT/cache/huggingface/hub
export TRANSFORMERS_CACHE=$RL_ROOT/cache/huggingface/transformers
export TORCH_HOME=$RL_ROOT/cache/torch
export PIP_CACHE_DIR=$RL_ROOT/cache/pip
export TMPDIR=$RL_ROOT/tmp
export PYTHONNOUSERSITE=1
export TOKENIZERS_PARALLELISM=false
EOF
# shellcheck disable=SC1090
source "$RL_ROOT/env.sh"

if command -v apt-get >/dev/null 2>&1; then
  apt-get update
  DEBIAN_FRONTEND=noninteractive apt-get install -y \
    tmux git git-lfs rsync curl wget aria2 pigz parallel jq tree htop \
    build-essential cmake pkg-config libhdf5-dev libopenblas-dev liblapack-dev \
    libcurl4-openssl-dev libssl-dev libxml2-dev
  git lfs install
fi

if ! command -v micromamba >/dev/null 2>&1; then
  mkdir -p "$RL_ROOT/tools/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C "$RL_ROOT/tools" bin/micromamba
  ln -sf "$RL_ROOT/tools/bin/micromamba" /usr/local/bin/micromamba
fi

for env_file in \
  "$REPO_ROOT/envs/environment-core.yml" \
  "$REPO_ROOT/envs/environment-geneformer.yml" \
  "$REPO_ROOT/envs/environment-scgpt.yml" \
  "$REPO_ROOT/envs/environment-r.yml"; do
  env_name="$(awk '/^name:/ {print $2}' "$env_file")"
  if micromamba env list | awk '{print $1}' | grep -qx "$env_name"; then
    echo "[skip] environment exists: $env_name"
  else
    micromamba create -y -f "$env_file"
  fi
done

micromamba run -n rheumlens-core pip install -e "$REPO_ROOT[full,dev]"

mkdir -p "$RL_ROOT/results/preflight"
{
  date -Is
  uname -a
  lscpu
  free -h
  df -hT
  nvidia-smi
} > "$RL_ROOT/results/preflight/system_report.txt"

micromamba run -n rheumlens-core python - <<'PY' > "$RL_ROOT/results/preflight/python_gpu_report.txt"
import platform
import torch
print('python:', platform.python_version())
print('torch:', torch.__version__)
print('torch_cuda:', torch.version.cuda)
print('cuda_available:', torch.cuda.is_available())
if torch.cuda.is_available():
    p = torch.cuda.get_device_properties(0)
    print('gpu:', p.name)
    print('memory_gib:', p.total_memory / 1024**3)
PY

echo "Bootstrap complete. Run: source $RL_ROOT/env.sh"
