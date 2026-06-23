#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

check_sha_path_manifest() {
  local manifest="$1"
  local base_dir="$2"
  [[ -f "$manifest" ]] || return 0
  awk 'NF >= 2 {print $1 "\t" $2}' "$manifest" | while IFS=$'\t' read -r sha rel; do
    [[ -n "$rel" && -n "$sha" ]] || continue
    rel="${rel#./}"
    local path="$base_dir/$rel"
    if [[ ! -f "$path" ]]; then
      echo "MISSING $path"
      return 1
    fi
    local got
    got="$(shasum -a 256 "$path" | awk '{print $1}')"
    if [[ "$got" != "$sha" ]]; then
      echo "SHA256_MISMATCH $path expected=$sha got=$got"
      return 1
    fi
  done
}

check_path_bytes_sha_manifest() {
  local manifest="$1"
  local base_dir="$2"
  [[ -f "$manifest" ]] || return 0
  awk -F '\t' 'NR > 1 && NF >= 3 {print $1 "\t" $3}' "$manifest" | while IFS=$'\t' read -r rel sha; do
    [[ -n "$rel" && -n "$sha" ]] || continue
    local path="$base_dir/$rel"
    if [[ ! -f "$path" ]]; then
      echo "MISSING $path"
      return 1
    fi
    local got
    got="$(shasum -a 256 "$path" | awk '{print $1}')"
    if [[ "$got" != "$sha" ]]; then
      echo "SHA256_MISMATCH $path expected=$sha got=$got"
      return 1
    fi
  done
}

if [[ -f figures/MANIFEST_SHA256.tsv ]]; then
  check_sha_path_manifest figures/MANIFEST_SHA256.tsv figures
fi

if [[ -f supplementary_tables/MANIFEST_SHA256.tsv ]]; then
  check_sha_path_manifest supplementary_tables/MANIFEST_SHA256.tsv supplementary_tables
fi

if [[ -f extension_results_20260623/MANIFEST_SHA256.tsv ]]; then
  check_path_bytes_sha_manifest extension_results_20260623/MANIFEST_SHA256.tsv extension_results_20260623
fi

echo "Manifest verification passed."
