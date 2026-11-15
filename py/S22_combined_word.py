#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build ONE Word file that holds the manuscript and all ten figures.

Why: the submission package keeps text and figures separate because MDPI uploads
them separately, but reading or archiving the paper then means opening eleven
files. This script produces a single reading copy:

    part 1  the full manuscript (A4 portrait, same rendering as the package)
    part 2  one page per figure - four main figures, then supplementary S1-S6

Each figure sits at its native size on a page sized to fit it, so nothing is
shrunk and no label becomes unreadable. Full legends are NOT repeated here: they
already sit in the "Figure Legends" section of part 1, and repeating them would
either shrink the art or split every figure across two pages.

Usage:
    python py/S22_combined_word.py
Output:
    output/Medicina_manuscript_with_figures_<date>.docx
"""

import datetime
import os
import shutil
import sys
import tempfile

import pymupdf
from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Mm, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
FIG = os.path.join(PROJECT, "results", "figures")
MD = os.path.join(PROJECT, "manuscript", "submission", "Medicina",
                  "Medicina_manuscript.md")

FIGURES = [
    ("Figure 1", "Fig1_global_atlas.pdf", 300),
    ("Figure 2", "Fig2_myeloid_depletion.pdf", 300),
    ("Figure 3", "Fig3_function.pdf", 300),
    ("Figure 4", "Fig4_validation.pdf", 300),
    ("Supplementary Figure S1", "Supplementary_Figure_S1.pdf", 200),
    ("Supplementary Figure S2", "Supplementary_Figure_S2.pdf", 200),
    ("Supplementary Figure S3", "Supplementary_Figure_S3.pdf", 200),
    ("Supplementary Figure S4", "Supplementary_Figure_S4.pdf", 200),
    ("Supplementary Figure S5", "Supplementary_Figure_S5.pdf", 250),
    ("Supplementary Figure S6", "Supplementary_Figure_S6.pdf", 300),
]

MARGIN = 8.0          # mm, figure pages only


def rasterise(pdf_path, dpi, out_dir):
    doc = pymupdf.open(pdf_path)
    page = doc[0]
    png = os.path.join(out_dir, os.path.basename(pdf_path) + ".png")
    page.get_pixmap(dpi=dpi).save(png)
    size = (page.rect.width / 72 * 25.4, page.rect.height / 72 * 25.4)
    doc.close()
    return png, size


def build_text_part():
    """Render the manuscript with the same converter the package uses.

    The converted file is opened as the base document - not copied into a blank
    one - so the A4 page setup, Times New Roman 11 pt body and the real
    Title/Heading styles survive.
    """
    sys.path.insert(0, HERE)
    from S14_md_to_docx import convert          # noqa: E402
    tmp = tempfile.mktemp(suffix=".docx")
    convert(MD, tmp)
    return Document(tmp), tmp


def add_figure_page(doc, label, png, size_mm):
    """One figure per page, at native size.

    No label line is printed: every figure PDF already carries its own bold
    title ("Figure 1. ...", "Supplementary Figure S2. ...") as part of the
    artwork, so a label above it would read as a duplicate.
    """
    w, h = size_mm
    sec = doc.add_section(WD_SECTION.NEW_PAGE)
    sec.orientation = WD_ORIENT.LANDSCAPE if w > h else WD_ORIENT.PORTRAIT
    sec.page_width = Mm(w + 2 * MARGIN)
    sec.page_height = Mm(h + 2 * MARGIN)
    for attr in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
        setattr(sec, attr, Mm(MARGIN))

    pic = doc.add_paragraph()
    pic.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pic.paragraph_format.space_before = Pt(0)
    pic.paragraph_format.space_after = Pt(0)
    pic.paragraph_format.line_spacing = 1.0
    pic.add_run().add_picture(png, width=Mm(w))


def main():
    out = os.path.join(
        PROJECT, "output",
        "Medicina_manuscript_with_figures_%s.docx"
        % datetime.date.today().strftime("%Y%m%d"))
    if len(sys.argv) > 1:                  # optional explicit target
        out = sys.argv[1]

    doc, tmp_text = build_text_part()

    work = tempfile.mkdtemp(prefix="figs_")
    try:
        for label, name, dpi in FIGURES:
            png, size = rasterise(os.path.join(FIG, name), dpi, work)
            add_figure_page(doc, label, png, size)
        doc.save(out)
    finally:
        shutil.rmtree(work, ignore_errors=True)
        if os.path.isfile(tmp_text):
            os.remove(tmp_text)

    print("[ok] %s  (%.1f MB)" % (out, os.path.getsize(out) / 1048576))
    return 0


if __name__ == "__main__":
    sys.exit(main())
