#!/usr/bin/env python3
"""Structural checks the journal will apply, run as a test rather than by eye.

Checks: abstract and author-summary length, unstructured abstract, figure and table
citation order against declaration order, every declared item cited at least once,
and every cited figure file present on disk.
"""
from __future__ import annotations
import ast, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S = (ROOT / "manuscript" / "manuscript_v5.md").read_text()
FIGDIR = ROOT / "figures" / "out"
SUPDIR = ROOT / "supplementary"

fails: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(f"{'OK  ' if ok else 'FAIL'}  {msg}")
    if not ok:
        fails.append(msg)


# ---- lengths -----------------------------------------------------------------
abstract = S[S.index("## Abstract"): S.index("**Keywords:**")]
abstract = re.sub(r"^#+.*$", "", abstract, flags=re.M)
n_ab = len(abstract.split())
check(n_ab <= 300, f"abstract is {n_ab} words (PLOS limit 300)")

au = S[S.index("## Author summary"): S.index("## 1.")]
au = re.sub(r"^#+.*$|^-+$", "", au, flags=re.M)
n_au = len(au.split())
check(150 <= n_au <= 200, f"author summary is {n_au} words (PLOS range 150-200)")

check(not re.search(r"^\*\*(Background|Approach|Results|Conclusions)\b", abstract, re.M),
      "abstract is unstructured")
check(S.index("## Abstract") < S.index("## Author summary"),
      "author summary follows the abstract")
check("\n# Results\n" in S, "a top-level Results section exists")
check("\n# Discussion\n" in S and "\n# Methods\n" in S,
      "top-level Discussion and Methods sections exist")

# ---- figures -----------------------------------------------------------------
body = S[: S.index("## Figure legends")]


def cited(pattern: str, text: str) -> list[str]:
    seen: list[str] = []
    for m in re.finditer(pattern, text):
        if m.group(1) not in seen:
            seen.append(m.group(1))
    return seen


fig_cited = cited(r"\bFigures? (\d+)(?![0-9])", body)
fig_decl = re.findall(r"^\*\*Figure (\d+)\.", S, flags=re.M)
check(fig_cited == fig_decl == [str(i) for i in range(1, len(fig_decl) + 1)],
      f"main figures cited and declared in order: {fig_cited} / {fig_decl}")

sfig_cited = cited(r"\bFigures? S(\d+)(?![0-9])", body)
sfig_decl = re.findall(r"^\*\*Figure S(\d+)\.", S, flags=re.M)
check(sfig_cited == sfig_decl == [str(i) for i in range(1, len(sfig_decl) + 1)],
      f"supporting figures cited and declared in order: {sfig_cited} / {sfig_decl}")

si_body = S[: S.index("## Supporting information")]
tab_cited = cited(r"\bTables? S(\d+)(?![0-9])[a-c]?", si_body)
tab_decl = re.findall(r"^\*\*Table S(\d+)\.\*\*", S[S.index("## Supporting information"):],
                      flags=re.M)
check(tab_cited == tab_decl == [str(i) for i in range(1, len(tab_decl) + 1)],
      f"supporting tables cited and declared in order: {tab_cited} / {tab_decl}")

# ---- files exist -------------------------------------------------------------
sys.path.insert(0, str(ROOT / "tools"))
from build_html import FIGS  # noqa: E402

missing = [v for v in FIGS.values() if not (FIGDIR / v).exists()]
check(not missing, f"every figure file present ({len(FIGS)} figures); missing: {missing}")
check(sorted(FIGS) == sorted(fig_decl + ["S" + n for n in sfig_decl]),
      "the HTML builder's figure map covers exactly the declared figures")

want = {f"S{n}" for n in tab_decl}
have = {m.group(1) for f in SUPDIR.glob("Table_S*.tsv")
        if (m := re.match(r"Table_(S\d+)[a-c]?_", f.name))}
check(want <= have, f"every declared table file present; missing: {sorted(want - have)}")

# ---- no leftovers ------------------------------------------------------------
for bad, why in [(r"CALIBRATION_[AB]_RESULT", "unfilled calibration placeholder"),
                 (r"LOSS_(BALANCED|RAGGED)", "unfilled loss placeholder"),
                 (r"\breferees?\b", "reviewer-history language"),
                 (r"design-preserving", "old permutation name"),
                 (r"design-independent component", "withdrawn claim")]:
    hits = re.findall(bad, S, flags=re.I)
    check(not hits, f"no {why} ({len(hits)} hits)")

# ---- reproducibility hazards in the released scripts -------------------------
SCRIPT_DIRS = [ROOT / d for d in ("sim", "walkthrough", "audit", "tools", "figures/src")]
scripts = [f for d in SCRIPT_DIRS if d.exists()
           for f in d.rglob("*") if f.suffix in (".py", ".R", ".sh")]

# Parse rather than grep: the string "hash(" also appears in the comments that
# explain why the built-in is not used, and a grep flags those forever.
hash_seeds = []
for f in scripts:
    if f.suffix != ".py":
        continue
    try:
        tree = ast.parse(f.read_text())
    except SyntaxError:
        continue
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "hash"):
            hash_seeds.append(f"{f.relative_to(ROOT)}:{node.lineno}")
check(not hash_seeds, f"no salted built-in hash() called in any script: {hash_seeds}")

abs_paths = [f"{f.relative_to(ROOT)}:{i}" for f in scripts
             for i, ln in enumerate(f.read_text().splitlines(), 1)
             if re.search(r"[\"'](/Volumes/|/Users/)", ln)]
check(not abs_paths, f"no machine-specific absolute path in a released script: {abs_paths}")

# ---- retired vocabulary anywhere a reader can see it -------------------------
# The manuscript body was cleaned twice while Figure 8 kept the old gate names and
# Table S7 kept "unstratified and stratified". Text-only checks cannot see either,
# so scan the rendered SVGs and the released tables as well.
RETIRED_VISIBLE = [
    (r"design-preserving", "the old permutation name"),
    (r"(?<![a-z-])unstratified", "the old name for the free permutation"),
    (r"(?<![a-z-])stratified permutation", "the old name for the collection-preserving test"),
    (r"design strata", "the old name for collection strata"),
    (r"model-selection artefact", "a withdrawn stop condition"),
    (r"tuned\s*[≈~]\s*frozen", "a gate that was removed from the analysis"),
    (r"identifiability boundary", "a withdrawn framing"),
    (r"\b2,?000 permutation", "a permutation count never run"),
]
def strip_svg(t: str) -> str:
    return re.sub(r"<[^>]+>", " ", t)

# Version-history language belongs in the cover letter, never on a figure of a
# new submission. Checked on figures only: prose may legitimately use these words.
FIGURE_ONLY = [
    (r"\bRevised\b", "version-history language (\"Revised\")"),
    (r"\bv2\b|\bprevious version\b|\bold Figure\b", "a reference to an earlier version"),
    (r"\breferee\b|\breviewer\b", "reviewer-history language"),
]
seen_bad = []
for f in sorted((ROOT / "figures" / "out").glob("*.svg")):
    body = strip_svg(f.read_text(errors="replace"))
    for pat, why in RETIRED_VISIBLE + FIGURE_ONLY:
        if re.search(pat, body, re.I):
            seen_bad.append(f"{f.name}: {why}")
for f in sorted(SUPDIR.glob("Table_S*.tsv")):
    body = f.read_text(errors="replace")
    for pat, why in RETIRED_VISIBLE:
        if re.search(pat, body, re.I):
            seen_bad.append(f"{f.name}: {why}")
check(not seen_bad, f"no retired vocabulary in any figure or released table: {seen_bad}")

# ---- the released tables use the manuscript's column names -------------------
sys.path.insert(0, str(ROOT / "tools"))
from si_names import check_header  # noqa: E402
hdr_bad = []
for f in sorted(SUPDIR.glob("Table_S*.tsv")):
    rows = [ln for ln in f.read_text().splitlines() if not ln.startswith("#")]
    if rows:
        hdr_bad += [f"{f.name}: {b}" for b in check_header(rows[0].split("\t"))]
check(not hdr_bad, f"released table headers use manuscript names: {hdr_bad}")

# ---- counts stated in prose must match the code and the result tables --------
# Every one of these was wrong at some point: the text said six regimes while the
# simulation ran seven, and said 2,000 permutations while the code ran 1,000.
import pandas as pd  # noqa: E402

screen_src = next((q for q in [
    ROOT.parent / "20_repo_restructure_20260907/scripts/cohorts/run_design_screen.py",
    ROOT / "submission/02_SOURCE/analysis_scripts/run_design_screen.py"] if q.exists()), None)
if screen_src is None:
    check(False, "the screen script is available to check permutation counts against")
else:
    n_perm = int(re.search(r"^N_PERM\s*=\s*(\d+)", screen_src.read_text(), re.M).group(1))
    floor = f"1/{n_perm + 1}"
    check(floor in S, f"the complete-pipeline floor quoted in the text is {floor}")
    stale = sorted(set(re.findall(r"(?<![\d.])1/(\d{3,5})\b", S)) - {str(n_perm + 1), "201"})
    check(not stale, f"no permutation floor in the text that the code does not run: {stale}")

sim = pd.read_csv(ROOT / "sim/results/extended_simulation_summary.tsv",
                  sep="\t", header=[0, 1], index_col=[0, 1, 2])
n_arm = sim.index.get_level_values(0).nunique()
n_rho = sim.index.get_level_values(2).nunique()
n_lab = sim.index.get_level_values(1).nunique()
words = {6: "six", 7: "seven", 8: "eight"}
check(f"{words[n_arm]} regimes" in S or f"{words[n_arm]} data-generating regimes" in S,
      f"the text says {words[n_arm]} regimes, matching the simulation output")
wrong = [w for n, w in words.items() if n != n_arm and (f"{w} regimes" in S or
         f"{w} data-generating regimes" in S)]
check(not wrong, f"no other regime count stated anywhere: {wrong}")

# 200 replicates per cell is a fixed setting of the simulation script
total = f"{n_arm * n_rho * n_lab * 200:,}"
check(total in S, f"the total cohort count in the text is {total}")

# ---- short title -------------------------------------------------------------
m = re.search(r"^\*\*Short title:\*\*\s*(.+?)\s*$", S, re.M)
check(m is not None, "a short title is declared")
if m:
    check(len(m.group(1)) <= 70,
          f"short title is {len(m.group(1))} characters (PLOS limit 70)")

# ---- the printed AUC must be the AUC its p-value tests ------------------------
# Table 1 used to print design_auc_linear, the tuned nested-CV model, beside
# p_design_auc, which is computed against design_auc_frozen. The p-value was
# right and the number next to it was a different model - up to 0.057 AUC apart -
# so the table invited the reader to read the p as a test of what they could see.
t1_src = (ROOT / "tools" / "build_table1.py").read_text()
check("design_auc_linear" not in t1_src,
      "Table 1 prints the frozen AUC, the statistic its p-value tests")

# Same trap in the figures: plotting the tuned AUC while Table 1 printed the
# frozen one put two numbers under one name. For COMBAT influenza the diagnosis
# classifier read 0.996 in Figure 3A and 0.475 in Table 1.
INFERENTIAL_FIGS = ["fig_spectrum.R", "fig_supp.R"]
for name in INFERENTIAL_FIGS:
    src = (ROOT / "figures" / "src" / name).read_text()
    stale = ["design_auc_linear"] if "design_auc_linear" in src else []
    # a script may use `disease_auc` as a local name, but only if it assigned the
    # frozen column to it first
    if re.search(r"\bdisease_auc\b", src) and "observed_frozen_auc" not in src:
        stale.append("disease_auc")
    check(not stale, f"{name} plots the frozen pipeline, not the tuned one: {stale}")

# ---- the table in the text is the generated one ------------------------------
gen = (ROOT / "manuscript" / "table1.md").read_text().strip()
i = S.index("| Comparison | Donors |")
check(S[i:S.index(chr(10) * 2, i)].strip() == gen,
      "the Table 1 in the manuscript is the generated one, not a hand-kept copy")

# ---- a restricted permutation is named for what it conditioned on ------------
# Where a comparison records no collection variable the screen falls back to
# tertiles of log cells per donor, which is sample quality. Reporting that column
# as collection-preserving for those rows would misname the test.
import csv
s4 = ROOT / "supplementary" / "Table_S4_full_screen.tsv"
if s4.exists():
    rows = list(csv.DictReader(s4.open(), delimiter="\t"))
    fallback = {r["cohort"] for r in rows
                if r.get("variable_group") == "all"
                and not r.get("strata_definition", "").startswith("batch__")}
    marked = gen.count("\u2020")
    check(marked == len(fallback),
          f"Table 1 marks all {len(fallback)} sample-quality-preserving rows "
          f"(found {marked} daggers)")
    check("dagger" in S,
          "the Table 1 caption explains the dagger")

# ---- cover letter ------------------------------------------------------------
cl = ROOT / "submission" / "01_UPLOAD" / "Cover_Letter.md"
if cl.exists():
    n_cl = len(cl.read_text().split())
    check(n_cl <= 650, f"cover letter is {n_cl} words (aim for one page, <=650)")

print()
if fails:
    print(f"{len(fails)} check(s) failed.")
    sys.exit(1)
print("all structural checks passed.")
