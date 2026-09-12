#!/usr/bin/env python3
"""Render the submission PDF from manuscript_v5.md.

PLOS wants a line-numbered, double-spaced PDF with the main figures embedded.
This replaces the hand-maintained LaTeX body that used to live in
.github/plos_pdf/: that copy was a second source of truth and had already
drifted from the markdown - its abstract still read "both around 0.8" after the
text had been corrected to "both in the mid-0.8s".

The figure map is imported from build_html rather than restated, so the PDF, the
HTML and the TIFFs cannot disagree about which file is Figure 3.

Engine: pdflatex. XeLaTeX would need system fonts that the build container does
not carry; nothing here depends on the engine, and PLOS only needs a legible
line-numbered PDF at this stage.

Usage:  python3 tools/build_pdf.py [--out PATH] [--engine pdflatex]
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_html import FIGS, OUTDIR  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MD = ROOT / "manuscript" / "manuscript_v5.md"
# Rendered into the workspace, not into submission/: build_submission.py
# deletes that tree on every run and copies this file in.
DEFAULT_OUT = ROOT / "manuscript" / "Manuscript.pdf"

# Line numbers, double spacing, and figures that stay where they are put rather
# than floating to the end of a 70-page document.
HEADER = r"""
\usepackage{lineno}
\usepackage{setspace}
\usepackage{float}
\usepackage[margin=1in]{geometry}
\usepackage{pdflscape}
\usepackage{textgreek}
\usepackage{newunicodechar}

% pandoc turns a raw \begin{env} in the markdown into a raw LaTeX BLOCK and stops
% parsing markdown until the matching \end - which fed Table 1 to LaTeX verbatim
% and died on the V_D in its header. A bare macro is passed through instead.
\newcommand{\blandscape}{\begin{landscape}\footnotesize\setlength{\tabcolsep}{3pt}}
\newcommand{\elandscape}{\end{landscape}\normalsize}

% Latin Modern has no glyph for these, and a missing glyph is not an error: the
% first build of this PDF dropped 27 characters in silence, rho 16 times and the
% "approximately equal" in both permutation floors. Map them, and fail the build
% if the engine reports a missing character anyway.
\newunicodechar{ρ}{\textrho}
\newunicodechar{≤}{\ensuremath{\leq}}
\newunicodechar{≥}{\ensuremath{\geq}}
\newunicodechar{≈}{\ensuremath{\approx}}
\newunicodechar{∈}{\ensuremath{\in}}
\newunicodechar{⁴}{\textsuperscript{4}}
\newunicodechar{⁻}{\textsuperscript{\ensuremath{-}}}
\newunicodechar{ȳ}{\ensuremath{\bar{y}}}
\linenumbers
\doublespacing
\makeatletter
\renewcommand{\fps@figure}{H}
\makeatother
\setlength{\emergencystretch}{3em}
"""


WIDE = 8          # columns beyond which a table needs the page turned sideways


def _cells(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def fix_wide_tables(lines: list[str]) -> list[str]:
    """Give wide pipe tables proportional column widths, and lay them sideways.

    pandoc derives a pipe table's relative column widths from the number of
    dashes in its separator row. build_table1.py emits a uniform `|---|---:|`
    separator, so all 15 of Table 1's columns came out the same width: the
    comparison name collided with the donor count and several cells were
    truncated mid-value ("0.982-"). Widths are recomputed here from the widest
    cell in each column, and the table is put in landscape at \footnotesize so
    it fits. The markdown source is not touched - only what pandoc is handed.
    """
    out, i = [], 0
    while i < len(lines):
        sep = i + 1
        if (lines[i].startswith("|") and sep < len(lines)
                and re.fullmatch(r"\|[\s:|-]+", lines[sep] or " ")
                and set(lines[sep]) <= set("|-: ")):
            j = sep + 1
            while j < len(lines) and lines[j].startswith("|"):
                j += 1
            head, body = _cells(lines[i]), [_cells(l) for l in lines[sep + 1:j]]
            n = len(head)
            if n > WIDE:
                width = [6 + max([len(head[k])] + [len(r[k]) for r in body if k < len(r)])
                         for k in range(n)]
                align = _cells(lines[sep])
                bar = "|" + "|".join(
                    ("-" * max(3, w)) + (":" if align[k].endswith(":") else "")
                    for k, w in enumerate(width)) + "|"
                out += ["", r"\blandscape", ""]
                out += [lines[i], bar] + lines[sep + 1:j]
                out += ["", r"\elandscape", ""]
                i = j
                continue
        out.append(lines[i])
        i += 1
    return out


def check_tables(pdf: Path, md_lines: list[str]) -> None:
    """Every cell of every wide table must survive into the PDF.

    The failure this catches is silent: LaTeX overfills the line, the value is
    clipped, and the PDF still builds. Nothing else in the pipeline reads the
    PDF, so without this the manuscript could ship with a truncated Table 1.
    """
    want = set()
    for line in md_lines:
        if line.startswith("|") and len(_cells(line)) > WIDE:
            for c in _cells(line):
                if re.fullmatch(r"[\d.]+(-[\d.]+)?", c):
                    want.update(re.findall(r"[\d.]+", c))
    if not want:
        return
    try:
        txt = subprocess.run(["pdftotext", "-layout", str(pdf), "-"],
                             capture_output=True, text=True, check=True).stdout
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("  (pdftotext unavailable - skipped the table-fidelity check)")
        return
    flat = " ".join(txt.split())
    missing = sorted(v for v in want if v not in flat)
    if missing:
        sys.exit(f"{len(missing)} wide-table value(s) did not survive into the "
                 f"PDF, so a column is being clipped: {missing[:12]}")
    print(f"  table fidelity: all {len(want)} wide-table values present in the PDF")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--engine", default="pdflatex")
    args = ap.parse_args()

    missing = [f for f in FIGS.values() if not (OUTDIR / f).exists()]
    if missing:
        sys.exit("missing figures: " + ", ".join(missing))

    text = re.sub(r"^!\[[^\]]*\]\([^)]*\)\n?", "", MD.read_text(), flags=re.M)

    out, seen = [], []
    for line in text.split("\n"):
        out.append(line)
        m = re.match(r"^\*\*Figure (S?\d+)\.", line)
        if m and m.group(1) in FIGS:
            n = m.group(1)
            out += ["", f"![]({OUTDIR / FIGS[n]}){{width=100%}}"]
            seen.append(n)
    if seen != list(FIGS):
        sys.exit(f"figure order/coverage mismatch: {seen}")

    out = fix_wide_tables(out)

    build = ROOT / "docs" / "_pdf_build"
    build.mkdir(parents=True, exist_ok=True)
    (build / "body.md").write_text("\n".join(out))
    (build / "header.tex").write_text(HEADER)

    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["pandoc", str(build / "body.md"), "-s",
           f"--pdf-engine={args.engine}",
           "-V", "documentclass=article", "-V", "fontsize=11pt",
           "-V", "colorlinks=true", "-V", "linkcolor=black",
           "-V", "urlcolor=blue",
           "-H", str(build / "header.tex"),
           "-o", str(dest)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"pandoc/{args.engine} failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    # A glyph the font lacks is a warning, not an error, so the PDF builds and the
    # character is simply absent. Treat it as fatal.
    dropped = sorted({m for m in re.findall(r"There is no (.) \(U\+[0-9A-F]+\)",
                                            r.stdout + r.stderr)})
    if dropped:
        sys.exit("the engine dropped these characters; add a \\newunicodechar for "
                 f"each and rebuild: {dropped}")
    shutil.rmtree(build, ignore_errors=True)
    print(f"built {dest}  ({dest.stat().st_size/1e6:.2f} MB, {len(seen)} figures)")
    check_tables(dest, out)


if __name__ == "__main__":
    main()
