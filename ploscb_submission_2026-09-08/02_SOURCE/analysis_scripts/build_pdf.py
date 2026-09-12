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
\linenumbers
\doublespacing
\makeatletter
\renewcommand{\fps@figure}{H}
\makeatother
\setlength{\emergencystretch}{3em}
"""


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
    shutil.rmtree(build, ignore_errors=True)
    print(f"built {dest}  ({dest.stat().st_size/1e6:.2f} MB, {len(seen)} figures)")


if __name__ == "__main__":
    main()
