#!/usr/bin/env python3
"""Build a readable initial-submission PDF from the locked Markdown and figures."""

from __future__ import annotations

import html
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from path_config import RELEASE_ROOT


MANUSCRIPT = RELEASE_ROOT / "manuscript" / "manuscript.md"
AUTHOR_SUMMARY = (
    RELEASE_ROOT
    / "submission"
    / "plos_computational_biology"
    / "AUTHOR_SUMMARY.md"
)
FIGURE_LEGENDS = (
    RELEASE_ROOT
    / "submission"
    / "plos_computational_biology"
    / "FIGURE_LEGENDS.md"
)
FIGURE_ROOT = RELEASE_ROOT / "figures" / "main"
OUTPUT = (
    RELEASE_ROOT
    / "submission"
    / "plos_computational_biology"
    / "Design_validity_v2.0.0-rc1_initial_submission.pdf"
)


def clean_inline(text: str) -> str:
    text = text.replace("&lt;sup&gt;", "<super>").replace(
        "&lt;/sup&gt;", "</super>"
    )
    text = text.replace("<sup>", "<super>").replace("</sup>", "</super>")
    text = re.sub(r"`([^`]+)`", r"<font name='Courier'>\1</font>", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"\*([^*]+)\*", r"<i>\1</i>", text)
    return text


def read_blocks(path: Path) -> list[tuple[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    blocks: list[tuple[str, str]] = []
    paragraph: list[str] = []

    def flush() -> None:
        if paragraph:
            blocks.append(("p", " ".join(item.strip() for item in paragraph)))
            paragraph.clear()

    in_table = False
    table_lines: list[str] = []
    for line in lines:
        if line.startswith("|"):
            flush()
            in_table = True
            table_lines.append(line)
            continue
        if in_table:
            blocks.append(("table", "\n".join(table_lines)))
            table_lines.clear()
            in_table = False
        if not line.strip():
            flush()
        elif line.startswith("#### "):
            flush()
            blocks.append(("h4", line[5:].strip()))
        elif line.startswith("### "):
            flush()
            blocks.append(("h3", line[4:].strip()))
        elif line.startswith("## "):
            flush()
            blocks.append(("h2", line[3:].strip()))
        elif line.startswith("# "):
            flush()
            blocks.append(("h1", line[2:].strip()))
        elif re.match(r"^\d+\.\s", line):
            flush()
            blocks.append(("li", re.sub(r"^\d+\.\s+", "", line)))
        else:
            paragraph.append(line)
    flush()
    if table_lines:
        blocks.append(("table", "\n".join(table_lines)))
    return blocks


def parse_markdown_table(raw: str, body_style: ParagraphStyle) -> Table | None:
    rows = []
    for line in raw.splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if cells and all(re.fullmatch(r":?-+:?", cell) for cell in cells):
            continue
        rows.append([Paragraph(clean_inline(cell), body_style) for cell in cells])
    if not rows:
        return None
    available = A4[0] - 4 * cm
    col_width = available / len(rows[0])
    table = Table(rows, colWidths=[col_width] * len(rows[0]), repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#17365D")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#B7C9DC")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F4F7FA")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return table


def page_number(canvas, doc) -> None:
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#D8DEE6"))
    canvas.line(2 * cm, 1.35 * cm, A4[0] - 2 * cm, 1.35 * cm)
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#5A6673"))
    canvas.drawString(
        2 * cm,
        0.9 * cm,
        "Design-validity study v2.0.0-rc1 initial submission",
    )
    canvas.drawRightString(A4[0] - 2 * cm, 0.9 * cm, f"Page {doc.page}")
    canvas.restoreState()


def styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "Title",
            parent=base["Title"],
            fontName="Helvetica-Bold",
            fontSize=18,
            leading=22,
            alignment=TA_LEFT,
            textColor=colors.HexColor("#172230"),
            spaceAfter=12,
        ),
        "h1": ParagraphStyle(
            "Heading1",
            parent=base["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=14,
            leading=18,
            textColor=colors.HexColor("#17365D"),
            spaceBefore=14,
            spaceAfter=6,
        ),
        "h2": ParagraphStyle(
            "Heading2",
            parent=base["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=12,
            leading=15,
            textColor=colors.HexColor("#244F73"),
            spaceBefore=11,
            spaceAfter=5,
        ),
        "h3": ParagraphStyle(
            "Heading3",
            parent=base["Heading3"],
            fontName="Helvetica-Bold",
            fontSize=10.5,
            leading=13,
            textColor=colors.HexColor("#2E5E80"),
            spaceBefore=8,
            spaceAfter=4,
        ),
        "body": ParagraphStyle(
            "Body",
            parent=base["BodyText"],
            fontName="Times-Roman",
            fontSize=10.25,
            leading=14,
            alignment=TA_LEFT,
            textColor=colors.HexColor("#20252B"),
            spaceAfter=7,
        ),
        "small": ParagraphStyle(
            "Small",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=8.7,
            leading=11.5,
            textColor=colors.HexColor("#4C5661"),
            spaceAfter=5,
        ),
        "caption": ParagraphStyle(
            "Caption",
            parent=base["BodyText"],
            fontName="Times-Roman",
            fontSize=9.2,
            leading=12.2,
            textColor=colors.HexColor("#252A30"),
            spaceBefore=6,
            spaceAfter=8,
        ),
    }


def add_blocks(story: list, blocks: list[tuple[str, str]], sty: dict) -> None:
    for kind, raw in blocks:
        text = clean_inline(html.escape(raw, quote=False))
        text = text.replace("&lt;super&gt;", "<super>").replace("&lt;/super&gt;", "</super>")
        text = text.replace("&lt;b&gt;", "<b>").replace("&lt;/b&gt;", "</b>")
        text = text.replace("&lt;i&gt;", "<i>").replace("&lt;/i&gt;", "</i>")
        text = text.replace("&lt;font name='Courier'&gt;", "<font name='Courier'>").replace("&lt;/font&gt;", "</font>")
        if kind == "h1":
            story.append(Paragraph(text, sty["title"]))
        elif kind == "h2":
            story.append(Paragraph(text, sty["h1"]))
        elif kind in {"h3", "h4"}:
            story.append(Paragraph(text, sty["h2"] if kind == "h3" else sty["h3"]))
        elif kind == "li":
            story.append(Paragraph(f"• {text}", sty["body"]))
        elif kind == "table":
            table = parse_markdown_table(raw, sty["small"])
            if table is not None:
                story.extend([table, Spacer(1, 8)])
        else:
            story.append(Paragraph(text, sty["body"]))


def author_summary_blocks() -> list[tuple[str, str]]:
    blocks = read_blocks(AUTHOR_SUMMARY)
    return [
        block
        for block in blocks
        if not (block[0] == "h1" and "author summary" in block[1].lower())
    ]


def legend_map() -> dict[int, str]:
    legends: dict[int, list[str]] = {}
    current: int | None = None
    for kind, text in read_blocks(FIGURE_LEGENDS):
        if kind == "h2":
            match = re.match(r"Figure\s+(\d+)\.\s*(.+)", text)
            if match:
                current = int(match.group(1))
                legends[current] = [f"<b>Figure {current}. {match.group(2)}</b>"]
        elif current is not None and kind == "p":
            legends[current].append(text)
    return {number: " ".join(parts) for number, parts in legends.items()}


def figure_flowables(number: int, legend: str, sty: dict) -> list:
    png = next(FIGURE_ROOT.glob(f"Figure_{number}_*.png"))
    from PIL import Image as PILImage

    with PILImage.open(png) as image:
        width, height = image.size
    max_width = A4[0] - 3.2 * cm
    max_height = A4[1] - 8.2 * cm
    scale = min(max_width / width, max_height / height)
    figure = Image(str(png), width=width * scale, height=height * scale)
    return [PageBreak(), figure, Paragraph(clean_inline(legend), sty["caption"])]


def main() -> None:
    sty = styles()
    story: list = []
    manuscript_blocks = read_blocks(MANUSCRIPT)
    figure_files = next(
        (
            index
            for index, block in enumerate(manuscript_blocks)
            if block[0] == "h2" and block[1] == "Figure files"
        ),
        len(manuscript_blocks),
    )
    manuscript_blocks = manuscript_blocks[:figure_files]
    introduction = next(
        index
        for index, block in enumerate(manuscript_blocks)
        if block[0] == "h2" and block[1].startswith("1. Introduction")
    )
    add_blocks(story, manuscript_blocks[:introduction], sty)
    story.append(Paragraph("Author Summary", sty["h1"]))
    add_blocks(story, author_summary_blocks(), sty)
    add_blocks(story, manuscript_blocks[introduction:], sty)

    legends = legend_map()
    for number in range(1, 7):
        story.extend(figure_flowables(number, legends[number], sty))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc = SimpleDocTemplate(
        str(OUTPUT),
        pagesize=A4,
        rightMargin=2 * cm,
        leftMargin=2 * cm,
        topMargin=1.8 * cm,
        bottomMargin=1.7 * cm,
        title="When design predicts disease",
        author="Hongyu Ying; Dandan Yun; Dan Liu",
        subject="PLOS Computational Biology initial submission",
    )
    doc.build(story, onFirstPage=page_number, onLaterPages=page_number)
    print(OUTPUT)


if __name__ == "__main__":
    main()
