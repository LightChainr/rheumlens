#!/usr/bin/env bash
# Re-run the manuscript checks against the files in this research object.
#
# The checkers resolve their inputs relative to the analysis workspace, whose
# layout differs from the submission package. Rather than patching the released
# scripts, assemble that layout from this directory in a temporary tree and run
# them there, so the scripts published here are byte-identical to the ones used
# in the analysis.
#
# Usage:  bash verify.sh
set -euo pipefail
cd "$(dirname "$0")"
HERE=$PWD
PYTHON="${PYTHON:-python3}"

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
mkdir -p "$T"/{manuscript,figures/out,supplementary,tools,results/screen,sim/results} \
         "$T"/walkthrough/results "$T"/submission/01_UPLOAD \
         "$T"/submission/02_SOURCE/analysis_scripts

cp "$HERE/01_UPLOAD/Manuscript.md"            "$T/manuscript/manuscript_v5.md"
cp "$HERE/01_UPLOAD/Cover_Letter.md"          "$T/submission/01_UPLOAD/"
cp "$HERE"/01_UPLOAD/Main_Figures/*       "$T/figures/out/"
cp "$HERE"/01_UPLOAD/Supporting_Figures/* "$T/figures/out/"
cp "$HERE"/01_UPLOAD/Supporting_Tables/*.tsv  "$T/supplementary/"
cp "$HERE"/02_SOURCE/analysis_scripts/*.py    "$T/tools/"
cp "$HERE/02_SOURCE/analysis_scripts/run_design_screen.py" \
   "$T/submission/02_SOURCE/analysis_scripts/"
cp -R "$HERE"/02_SOURCE/per_seed/seed_*       "$T/results/screen/"
for f in extended_simulation_summary mediator_arm_summary \
         calibration_arm_a_summary calibration_arm_b_summary; do
  cp "$HERE/02_SOURCE/result_tables/$f.tsv" "$T/sim/results/"
done
cp "$HERE/02_SOURCE/result_tables/decision_tree_walkthrough.tsv" \
   "$T/walkthrough/results/"

# The figures are released with their manuscript numbers; the checker expects
# the workspace names, which build_html.FIGS maps between.
( cd "$T/figures/out"
  for f in Fig[1-8]_*;  do mv "$f" "Fig_${f#Fig?_}"; done
  for f in FigS[1-7]_*; do n=${f#Fig}; mv "$f" "Fig_${n%%_*}_${n#*_}"; done )

echo "== numbers =="
"$PYTHON" "$T/tools/check_manuscript_numbers.py"
echo
echo "== structure =="
"$PYTHON" "$T/tools/check_manuscript_structure.py"
