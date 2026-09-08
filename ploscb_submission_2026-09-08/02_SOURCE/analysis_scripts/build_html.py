#!/usr/bin/env python3
"""Render manuscript_v5.md to a single self-contained HTML with figures inline.

Images are inserted immediately after each figure legend header, matched on the
"**Figure N." / "**Figure SN." pattern rather than on any pre-existing image tag,
so a figure that has no tag yet still gets one and the order cannot drift.
"""
import re, subprocess, shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MD   = ROOT / "manuscript" / "manuscript_v5.md"
OUT  = ROOT / "docs" / "manuscript_v5.html"
PKG  = ROOT / "submission" / "01_UPLOAD" / "Manuscript_with_Figures.html"

FIGS = {
    "1":  "Fig_dag.png",
    "2":  "Fig_simulation.png",
    "3":  "Fig_spectrum.png",
    "4":  "Fig_lupus_design.png",
    "5":  "Fig_exposure_vs_use.png",
    "6":  "Fig_adjustment.png",
    "7":  "Fig_transfer.png",
    "8":  "Fig_decision_tree.png",
    "S1": "Fig_S1_permutation_calibration.png",
    "S2": "Fig_S2_residualisation_width.png",
    "S3": "Fig_S3_mediator.png",
    "S4": "Fig_S4_composition.png",
    "S5": "Fig_S5_celltype_composition.png",
    "S6": "Fig_S6_learned_pooling.png",
    "S7": "Fig_S7_seeds.png",
}
OUTDIR = ROOT / "figures" / "out"

css = """body{max-width:44em;margin:2.5em auto;padding:0 1.5em;font:16px/1.62 -apple-system,"Helvetica Neue",Arial,sans-serif;color:#1b1b1b}
h1{font-size:1.55em;line-height:1.3;margin-bottom:.2em}
h2{font-size:1.18em;margin-top:2.2em;border-bottom:1px solid #ddd;padding-bottom:.25em}
h3{font-size:1.02em;margin-top:1.7em;color:#333}
table{border-collapse:collapse;width:100%;font-size:.80em;margin:1em 0}
th,td{border:1px solid #d5d5d5;padding:.35em .55em;text-align:left}
th{background:#f4f4f2}
img{max-width:100%;margin:1.2em 0;border:1px solid #eee}
code{background:#f4f4f2;padding:.1em .3em;border-radius:3px;font-size:.9em}
hr{border:0;border-top:1px solid #ddd;margin:2em 0}
blockquote{color:#555;border-left:3px solid #ddd;margin-left:0;padding-left:1em}
"""
def main() -> None:
    missing = [f for f in FIGS.values() if not (OUTDIR / f).exists()]
    if missing:
        sys.exit("missing figures: " + ", ".join(missing))

    text = MD.read_text()
    # Drop any image tags already present so re-running is idempotent.
    text = re.sub(r"^!\[[^\]]*\]\([^)]*\)\n?", "", text, flags=re.M)

    out, seen = [], []
    for line in text.split("\n"):
        out.append(line)
        m = re.match(r"^\*\*Figure (S?\d+)\.", line)
        if m and m.group(1) in FIGS:
            n = m.group(1)
            out.append("")
            out.append(f"![Figure {n}]({OUTDIR / FIGS[n]})")
            seen.append(n)

    if seen != list(FIGS):
        sys.exit(f"figure order/coverage mismatch: {seen}")

    tmp = ROOT / "docs" / "_with_figures.md"
    tmp.write_text("\n".join(out))
    cssf = ROOT / "docs" / "_style.css"
    cssf.write_text(css)

    subprocess.run(["pandoc", str(tmp), "-s", "--embed-resources", "--standalone",
                    "--metadata", "title=Manuscript v5", "-c", str(cssf),
                    "-o", str(OUT)], check=True)
    tmp.unlink(); cssf.unlink()
    shutil.copy(OUT, PKG)
    print(f"built {OUT} with {len(seen)} figures -> {PKG}")


if __name__ == "__main__":
    main()
