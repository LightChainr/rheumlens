#!/usr/bin/env bash
# Rebuild every released table, figure and document from the analysis outputs.
#
# It does NOT re-run the screen, the simulation or the calibration - those need a
# 16-vCPU machine and hours, and their raw outputs are in results/ and sim/results/.
# Everything downstream of those outputs is rebuilt here, so a number that drifted
# between the analysis and the manuscript shows up as a failed check rather than
# as a sentence nobody re-read.
#
# Usage:  bash tools/build_all.sh
# Env:    PYTHON=... RSCRIPT=... to override the interpreters.
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON="${PYTHON:-python3}"
RSCRIPT="${RSCRIPT:-Rscript}"

echo "== supplementary tables =="
"$PYTHON" tools/stage_screen_results.py
"$PYTHON" tools/stage_supplementary_sim.py
"$PYTHON" tools/stage_static_tables.py
"$PYTHON" tools/build_registry_table.py
"$PYTHON" audit/build_design_manifest.py

echo
echo "== Table 1 =="
"$PYTHON" tools/build_table1.py | tail -5

echo
echo "== figures =="
for f in figures/src/*.R; do
  printf '%-38s' "$f"
  "$RSCRIPT" "$f" >/dev/null 2>&1 && echo "ok" || { echo "FAILED"; exit 1; }
done

echo
echo "== PLOS TIFFs =="
"$PYTHON" tools/build_tiffs.py | tail -3

echo
echo "== manuscript HTML =="
"$PYTHON" tools/build_html.py

echo
echo "== consistency checks =="
"$PYTHON" tools/check_manuscript_numbers.py | tail -3
"$PYTHON" tools/check_manuscript_structure.py
