#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert a submission markdown file to .docx for MDPI-style upload.

MDPI journals (Genes, Medicina, ...) accept Word or LaTeX, not Markdown, so
every derived package needs a .docx twin. This converter handles only what the
master manuscript actually contains:

  #/##/### headings, paragraphs with soft line wraps, **bold**, *italic*,
  `code`, bullet and numbered lists, pipe tables, and fenced code blocks.

Usage:
    python py/S14_md_to_docx.py <input.md> <output.docx> [--title-lines N]

Figures are NOT embedded: the four main figures and six supplementary figures
are uploaded as PDF files alongside the Word manuscript.
"""

import re
import sys

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Mm, Pt, RGBColor

INLINE = re.compile(r"(\*\*.+?\*\*|\*[^*\s][^*]*?\*|`[^`]+?`)")
TABLE_ROW = re.compile(r"^\s*\|.*\|\s*$")
TABLE_SEP = re.compile(r"^\s*\|[\s:|\-]+\|\s*$")
HEADING = re.compile(r"^(#{1,4})\s+(.*)$")
# Only "-" starts a bullet here. The master writes the correspondence footnote
# as "* Correspondence: ...", and a "* " rule turned it into a bulleted item
# with a visible dot in front of the asterisk.
BULLET = re.compile(r"^\s*-\s+(.*)$")
ORDERED = re.compile(r"^\s*(\d+)\.\s+(.*)$")
# A bullet item whose text already carries its own enumerator, e.g. "- (i) ...".
# Emitting those with a Word bullet gives a doubled marker ("• (i) ...").
LIST_ENUM = re.compile(r"^\((?:[0-9]{1,2}|[a-z]|[ivx]{1,6})\)\s")

# --- gene symbols -----------------------------------------------------------
# MDPI sets gene symbols in italics and protein names in roman. Only symbols
# actually analysed in this manuscript are listed, so that ordinary acronyms
# (GSVA, FDR, AUC, DEG, CPM, UMAP, IQR, LRT), accession numbers (GSE13355,
# GPL570), cell states (LC1, LC2) and package names are never italicised.
# Hyphenated symbols are spelled out because GENE_RE treats '-' as a word
# character on both sides, which keeps "HLA-DQB2" from splitting.
GENES = {
    "CD207", "CD1A", "CD1C", "CD83", "CD274", "LAMP3", "CCR7", "FCER1A",
    "HLA-DRA", "HLA-DQB2", "HLA-A", "B2M", "PSME2", "ACOT7", "FNDC4",
    "HPGDS", "SLC18A2", "DCD", "SCGB2A2", "AQP5", "KRT19", "KRT1", "KRT10",
    "KRT14", "KRT16", "KRT17", "KRT6A", "CPA3", "TPSAB1", "TPSB2", "KIT",
    "MS4A2", "LYZ", "SPI1", "RUNX1", "TCF3", "FOXO1", "FOXO3", "USF1",
    "ETV4", "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B", "IRF1", "RELB",
    "RELA", "NFKB1", "NFKB2", "JUND", "RARA", "NR2F2", "ARNTL", "FOS",
    "PFN1", "SFRP2", "CCL13", "CCL19", "CXCL12", "CCR2", "MARCKSL1",
    "S100A7", "S100A8", "S100A9", "PPA1", "TMEM176B",
}
GENE_RE = re.compile(
    r"(?<![A-Za-z0-9_\-])("
    + "|".join(re.escape(g) for g in sorted(GENES, key=len, reverse=True))
    + r")(?![A-Za-z0-9_\-])")
# "CD207 protein" names the protein, which stays roman.
PROTEIN_TAIL = re.compile(r"^\s*(?:protein|proteins)\b", re.I)
# A table caption must not be left alone at the foot of a page with its table
# overleaf.
CAPTION = re.compile(r"^Table\s+\d+[.:]")
# Superscript digits mark an affiliation line in the author block.
SUPERSCRIPT = set("¹²³⁴⁵⁶⁷⁸⁹⁰")


def add_page_numbers(sec):
    """Centred PAGE field in the footer (Dove: "Page numbers")."""
    par = sec.footer.paragraphs[0]
    par.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = par.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    run._r.append(begin)
    run._r.append(instr)
    run._r.append(end)


def add_line_numbers(sec):
    """Continuous line numbers (Dove: "Line numbers").

    w:lnNumType must sit between w:pgMar and w:cols in the sectPr schema
    order, so it is inserted before w:cols rather than appended.
    """
    sectPr = sec._sectPr
    if sectPr.find(qn("w:lnNumType")) is not None:
        return
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    cols = sectPr.find(qn("w:cols"))
    if cols is not None:
        cols.addprevious(ln)
    else:
        sectPr.append(ln)


def set_base_style(doc):
    # A4 is explicit: python-docx defaults to US Letter (215.9 x 279.4 mm),
    # which reflows every table and page break when the editorial office opens
    # the file in an A4 template. Margins and spacing follow Dove Medical
    # Press ("3-cm margins", "double-spacing", "page numbers", "line numbers");
    # MDPI is happy with the same settings.
    for sec in doc.sections:
        sec.page_width = Mm(210)
        sec.page_height = Mm(297)
        sec.top_margin = Mm(30)
        sec.bottom_margin = Mm(30)
        sec.left_margin = Mm(30)
        sec.right_margin = Mm(30)
        add_page_numbers(sec)
        add_line_numbers(sec)
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style.font.size = Pt(11)
    style.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    pf = style.paragraph_format
    pf.space_before = Pt(0)
    pf.space_after = Pt(0)
    pf.line_spacing = 2.0


def emit_text(par, text, bold=False):
    """Add text as runs, italicising whitelisted gene symbols only.

    Called for every plain span, for the inside of ``backticks`` and for bold
    spans, so gene names are set correctly wherever they appear. A symbol
    immediately followed by the word "protein" stays roman.
    """
    if not text:
        return
    pos = 0
    for m in GENE_RE.finditer(text):
        if PROTEIN_TAIL.match(text[m.end():m.end() + 9]):
            continue
        if m.start() > pos:
            run = par.add_run(text[pos:m.start()])
            run.bold = bold
        run = par.add_run(m.group(1))
        run.italic = True
        run.bold = bold
        pos = m.end()
    if pos < len(text):
        run = par.add_run(text[pos:])
        run.bold = bold


def add_runs(par, text):
    """Emit runs, honouring **bold**, *italic* and `code`."""
    pos = 0
    for m in INLINE.finditer(text):
        if m.start() > pos:
            emit_text(par, text[pos:m.start()])
        tok = m.group(0)
        if tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            emit_text(par, tok[2:-2], bold=True)
        elif tok.startswith("`"):
            # Backticks in the master carry gene symbols (CD207, LAMP3, ...) as
            # often as file and object names. Rendering them all in Courier put
            # 150 monospace islands inside Times body text, which reads as a
            # typesetting error in a Word submission. The span is therefore set
            # in the body face; only fenced blocks stay monospace. Gene symbols
            # inside the span are italicised, paths and package names are not.
            emit_text(par, tok[1:-1])
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            run = par.add_run(tok[1:-1])
            run.italic = True
        else:
            emit_text(par, tok)
        pos = m.end()
    if pos < len(text):
        emit_text(par, text[pos:])


def emit_heading(par, text, size):
    """Heading text with the heading's own formatting plus gene italics.

    The main title names CD207, so headings cannot be emitted with a single
    add_run() or the one gene symbol in the title would stay roman while the
    same symbol is italic everywhere else.
    """
    def styled(chunk, italic=False):
        run = par.add_run(chunk)
        run.bold = True
        run.italic = italic
        run.font.size = Pt(size)
        run.font.color.rgb = RGBColor(0, 0, 0)

    pos = 0
    for m in GENE_RE.finditer(text):
        if m.start() > pos:
            styled(text[pos:m.start()])
        styled(m.group(1), italic=True)
        pos = m.end()
    if pos < len(text):
        styled(text[pos:])


def add_numbered_item(doc, marker, text, hang):
    """One list item with the marker typed as text (never Word auto-numbering).

    Word auto-numbering is shared per numbering instance, so a numbered list in
    the Methods and the reference list end up in ONE sequence: section 2.11
    consumed 1-17 and the references continued at 18. Typing the marker keeps
    every list restarting at its own first number, survives any downstream
    conversion, and lets an existing enumerator such as "(i)" stand alone
    without a bullet in front of it.
    """
    par = doc.add_paragraph()
    pf = par.paragraph_format
    pf.left_indent = Pt(hang)
    pf.first_line_indent = Pt(-hang)
    pf.space_after = Pt(3)
    if marker:
        par.add_run(marker + " ")
    add_runs(par, text)
    return par


def split_cells(line):
    inner = line.strip()
    if inner.startswith("|"):
        inner = inner[1:]
    if inner.endswith("|"):
        inner = inner[:-1]
    return [c.strip() for c in inner.split("|")]


def add_table(doc, rows):
    ncol = max(len(r) for r in rows)
    table = doc.add_table(rows=0, cols=ncol)
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, row in enumerate(rows):
        cells = table.add_row().cells
        for j in range(ncol):
            text = row[j] if j < len(row) else ""
            par = cells[j].paragraphs[0]
            par.paragraph_format.space_after = Pt(0)
            par.paragraph_format.space_before = Pt(0)
            # Tables stay single-spaced: the double-spaced body would otherwise
            # make every cell twice as tall and push the table over several
            # pages (Dove asks whether a table still fits one journal page).
            par.paragraph_format.line_spacing = 1.0
            add_runs(par, text)
            for run in par.runs:
                run.font.size = Pt(9)
                if i == 0:
                    run.bold = True
    # Repeat the header row when a long table breaks across pages, otherwise the
    # reader meets a wall of cells on the next page with no column labels.
    if table.rows:
        trPr = table.rows[0]._tr.get_or_add_trPr()
        trPr.append(OxmlElement("w:tblHeader"))
    # Keep the rows of one table on the same page where they fit, so a table is
    # never split after a single row at the foot of a page, and never split a
    # single row's text across two pages.
    for row in table.rows:
        trPr = row._tr.get_or_add_trPr()
        trPr.insert(0, OxmlElement("w:cantSplit"))
    for row in table.rows[:-1]:
        for cell in row.cells:
            for par in cell.paragraphs:
                par.paragraph_format.keep_with_next = True
    return table


def convert(src, dst):
    with open(src, "r", encoding="utf-8") as fh:
        lines = fh.read().splitlines()

    doc = Document()
    set_base_style(doc)

    i = 0
    n = len(lines)
    first_heading_done = False

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # blank
        if not stripped:
            i += 1
            continue

        # fenced code block
        if stripped.startswith("```"):
            i += 1
            buf = []
            while i < n and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            par = doc.add_paragraph()
            par.paragraph_format.left_indent = Pt(18)
            run = par.add_run("\n".join(buf))
            run.font.name = "Courier New"
            run.font.size = Pt(9)
            continue

        # pipe table
        if TABLE_ROW.match(line):
            block = []
            while i < n and TABLE_ROW.match(lines[i]):
                if not TABLE_SEP.match(lines[i]):
                    block.append(split_cells(lines[i]))
                i += 1
            if block:
                add_table(doc, block)
                doc.add_paragraph()
            continue

        # signature block: after "Sincerely," every non-blank line is its own
        # paragraph. The soft-wrap merger below would otherwise collapse the
        # name, the two affiliations, the e-mail line and the closing sentence
        # into one run-on paragraph, which reads as a formatting error in a
        # cover letter.
        if stripped == "Sincerely,":
            par = doc.add_paragraph("Sincerely,")
            par.paragraph_format.space_after = Pt(2)
            i += 1
            while i < n and lines[i].strip():
                par = doc.add_paragraph(lines[i].strip())
                par.paragraph_format.space_after = Pt(0)
                i += 1
            continue

        # heading
        m = HEADING.match(stripped)
        if m:
            level = len(m.group(1))
            text = m.group(2).strip()
            # Headings carry real Word styles so that the outline/navigation pane
            # works, the editorial system can see the section structure, and the
            # in-text "Section 3.7" references can be located at all.
            if not first_heading_done and level == 1:
                par = doc.add_paragraph(style="Title")
                par.alignment = WD_ALIGN_PARAGRAPH.CENTER
                emit_heading(par, text, 15)
                par.paragraph_format.space_after = Pt(12)
                first_heading_done = True
            else:
                style = {2: "Heading 1", 3: "Heading 2"}.get(level, "Heading 3")
                par = doc.add_paragraph(style=style)
                emit_heading(par, text, 13 if level <= 2 else 11.5)
                par.paragraph_format.space_before = Pt(12 if level <= 2 else 8)
                par.paragraph_format.space_after = Pt(4)
                par.paragraph_format.keep_with_next = True
            i += 1
            continue

        # lists: collected as a block so that one block = one numbering
        mb = BULLET.match(line)
        mo = ORDERED.match(line)
        if mb or mo:
            nums, texts = [], []
            while i < n:
                cur = lines[i]
                mo2 = ORDERED.match(cur)
                mb2 = BULLET.match(cur)
                if mo2:
                    nums.append(int(mo2.group(1)))
                    texts.append(mo2.group(2))
                elif mb2:
                    nums.append(None)
                    texts.append(mb2.group(1))
                elif (cur.strip() and nums and not HEADING.match(cur.strip())
                      and not TABLE_ROW.match(cur)
                      and not cur.strip().startswith("```")):
                    texts[-1] += " " + cur.strip()   # soft-wrapped item
                else:
                    break
                i += 1

            if nums[0] is None:
                marks = [LIST_ENUM.match(t) for t in texts]
                if all(marks):
                    markers = [m.group(0).strip() for m in marks]
                    texts = [t[m.end():] for t, m in zip(texts, marks)]
                else:
                    markers = [""] * len(texts)      # plain bullet block
            else:
                markers = ["%d." % k for k in nums]

            if any(markers):
                hang = min(48.0, max(16.0, 4.6 * max(len(x) for x in markers) + 6.0))
                for mk, tx in zip(markers, texts):
                    add_numbered_item(doc, mk, tx, hang)
            else:
                for tx in texts:
                    par = doc.add_paragraph(style="List Bullet")
                    add_runs(par, tx)
            continue

        # paragraph: gather until blank / heading / table / list
        buf = [stripped]
        i += 1
        while i < n:
            nxt = lines[i]
            if (not nxt.strip() or HEADING.match(nxt.strip())
                    or TABLE_ROW.match(nxt) or BULLET.match(nxt)
                    or ORDERED.match(nxt) or nxt.strip().startswith("```")
                    # an affiliation line starts with a superscript marker and
                    # must stay its own paragraph, not merge into the one above
                    or nxt.strip()[:1] in SUPERSCRIPT):
                break
            buf.append(nxt.strip())
            i += 1
        par = doc.add_paragraph()
        par.alignment = WD_ALIGN_PARAGRAPH.LEFT
        add_runs(par, " ".join(buf))
        if CAPTION.match(par.text.strip()):
            par.paragraph_format.keep_with_next = True

    doc.save(dst)
    return dst


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: python py/S14_md_to_docx.py <input.md> <output.docx>")
        sys.exit(2)
    out = convert(sys.argv[1], sys.argv[2])
    print("[ok] wrote " + out)
