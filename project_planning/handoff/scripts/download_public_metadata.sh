#!/usr/bin/env bash
set -euo pipefail

RL_ROOT="${RL_ROOT:-/autodl-fs/data/rheumlens}"
OUT="$RL_ROOT/data/external/geo_metadata"
mkdir -p "$OUT"

accessions=(GSE135779 GSE285773 GSE174188 GSE137029 GSE250024 GSE142016 GSE162577 GSE179633 GSE186476 GSE158055)

for acc in "${accessions[@]}"; do
  mkdir -p "$OUT/$acc"
  prefix="${acc:0:${#acc}-3}nnn"
  base="https://ftp.ncbi.nlm.nih.gov/geo/series/$prefix/$acc"
  aria2c -c -x 4 -s 4 -d "$OUT/$acc" "$base/soft/${acc}_family.soft.gz" || true
  aria2c -c -x 4 -s 4 -d "$OUT/$acc" "$base/miniml/${acc}_family.xml.tgz" || true
done

find "$OUT" -type f -print0 | sort -z | xargs -0 sha256sum > "$OUT/SHA256SUMS.txt"
echo "Metadata download complete: $OUT"
