#!/usr/bin/env python3
"""Render the submission figures as PLOS-compliant TIFF.

Editorial Manager accepts TIFF or EPS only, at 300-600 dpi, 789-2250 px wide,
300-2625 px tall and under 10 MB. The PNGs we build for the HTML preview are
400 dpi and 2834 px wide, which is over the width cap, so they cannot be
uploaded as they stand.

Rasterise from the SVG rather than downsampling the PNG: svglite pins every
string with textLength/lengthAdjust, so rsvg reproduces the original glyph
metrics exactly instead of reflowing the text.

Figure numbers come from build_html.FIGS, the same mapping the manuscript and
the submission package use, so this cannot drift from the numbering in the text.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from build_html import FIGS, OUTDIR  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "submission" / "01_UPLOAD" / "TIFF_300dpi"
DPI = 300
# PLOS limits, in pixels.
W_MIN, W_MAX, H_MIN, H_MAX = 789, 2250, 300, 2625
MB_MAX = 10.0


def plos_name(num: str) -> str:
    """PLOS file naming: main figures Fig1..Fig8, supporting S1_Fig..S7_Fig."""
    return f"S{num[1:]}_Fig" if num.startswith("S") else f"Fig{num}"


def main() -> None:
    for tool in ("rsvg-convert", "magick"):
        if not shutil.which(tool):
            sys.exit(f"{tool} not found; install it before building the TIFFs")

    if DEST.exists():
        shutil.rmtree(DEST)
    DEST.mkdir(parents=True)

    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        png = Path(tmp) / "r.png"
        for num, fname in FIGS.items():
            svg = OUTDIR / Path(fname).with_suffix(".svg")
            if not svg.exists():
                sys.exit(f"missing {svg}")
            tif = DEST / f"{plos_name(num)}.tif"
            subprocess.run(["rsvg-convert", "-d", str(DPI), "-p", str(DPI),
                            "-f", "png", "-o", str(png), str(svg)], check=True)
            subprocess.run(["magick", str(png), "-compress", "lzw",
                            "-density", str(DPI), "-units", "PixelsPerInch",
                            str(tif)], check=True)

            w, h, comp = subprocess.run(
                ["magick", "identify", "-format", "%w %h %C", str(tif)],
                capture_output=True, text=True, check=True).stdout.split()
            w, h, mb = int(w), int(h), tif.stat().st_size / 1048576
            ok = (W_MIN <= w <= W_MAX and H_MIN <= h <= H_MAX
                  and mb <= MB_MAX and comp == "LZW")
            print(f"{'OK  ' if ok else 'FAIL'} Figure {num:<3} {tif.name:<12} "
                  f"{w}x{h}  {mb:.2f} MB  {comp}")
            if not ok:
                bad.append(tif.name)

    if bad:
        sys.exit(f"outside PLOS figure limits: {', '.join(bad)}")
    print(f"\nall {len(FIGS)} figures within PLOS limits -> {DEST}")


if __name__ == "__main__":
    main()
