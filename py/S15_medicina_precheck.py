#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pre-submission self-check for the Medicina package.

Reads only files inside manuscript/submission/Medicina/ plus the shared
supplementary archive, and reports PASS/FAIL for every hard requirement that
Medicina (MDPI) states in its Instructions for Authors.

Usage: python py/S15_medicina_precheck.py
Exit 0 only when every check passes.
"""

import hashlib
import os
import re
import sys

from docx import Document

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
SUB = os.path.join(PROJECT, "manuscript", "submission")
PKG = os.path.join(SUB, "Medicina")
FIG = os.path.join(PROJECT, "results", "figures")
ZIP = os.path.join(SUB, "Supplementary_materials_20261112.zip")

MS_MD = os.path.join(PKG, "Medicina_manuscript.md")
MS_DOCX = os.path.join(PKG, "Medicina_manuscript.docx")
CL_MD = os.path.join(PKG, "Medicina_cover_letter.md")
CL_DOCX = os.path.join(PKG, "Medicina_cover_letter.docx")

fails = []
lines = []


def check(ok, msg, detail=""):
    tag = "[PASS] " if ok else "[FAIL] "
    lines.append(tag + msg + (("  -> " + detail) if detail and not ok else ""))
    if not ok:
        fails.append(msg)


def read(p):
    with open(p, encoding="utf-8") as fh:
        return fh.read()


# ---------------------------------------------------------------- 1. contents
EXPECT = ["Medicina_manuscript.md", "Medicina_manuscript.docx",
          "Medicina_cover_letter.md", "Medicina_cover_letter.docx",
          "Figure_1.pdf", "Figure_2.pdf", "Figure_3.pdf", "Figure_4.pdf",
          "Supplementary_Figure_S1.pdf", "Supplementary_Figure_S2.pdf",
          "Supplementary_Figure_S3.pdf", "Supplementary_Figure_S4.pdf",
          "Supplementary_Figure_S5.pdf", "Supplementary_Figure_S6.pdf"]
actual_all = sorted(os.listdir(PKG))
# Word drops a "~$name.docx" lock file while the manuscript is open for editing.
# It is not part of the upload set, but it does mean the .docx may be locked.
actual = [f for f in actual_all if not f.startswith("~$")]
check(sorted(EXPECT) == actual, "package contains exactly the 14 upload files",
      "extra/missing: " + str(set(actual) ^ set(EXPECT)))
lines.append("        " + ", ".join(actual))
if len(actual) != len(actual_all):
    lines.append("        NOTE: Word lock file present -> close the .docx before "
                 "regenerating it")
check(os.path.getmtime(MS_DOCX) >= os.path.getmtime(MS_MD),
      "docx is not older than the markdown it was built from",
      "docx %s < md %s (re-run py/S14_md_to_docx.py)"
      % (os.path.getmtime(MS_DOCX), os.path.getmtime(MS_MD)))

# ------------------------------------------------------------------ 2. abstract
ms = read(MS_MD)
abs_block = ms.split("## Abstract", 1)[1].split("## Keywords", 1)[0]
for label in ("Background and Objectives", "Materials and Methods",
              "Results", "Conclusions"):
    check(label + ":" in abs_block, "abstract has labelled part: " + label)
words = len(re.findall(r"[A-Za-z0-9][A-Za-z0-9\-\+\.]*", abs_block))
check(words <= 300, "abstract <= 300 words", "found %d" % words)
lines.append("        abstract words: %d" % words)

# ---------------------------------------------------------------- 3. keywords
kw = re.split(r"^##\s+(?:\d+\.\s*)?Introduction\s*$",
               ms.split("## Keywords", 1)[1], flags=re.M)[0]
kw_list = [k.strip() for k in kw.replace("\n", " ").split(";") if k.strip()]
check(3 <= len(kw_list) <= 10, "keywords 3-10", "found %d" % len(kw_list))
lines.append("        keywords: " + "; ".join(kw_list))

# ------------------------------------------------------------ 4. declarations
DECL = ["## Supplementary Materials", "## Author Contributions", "## Funding",
        "## Institutional Review Board Statement",
        "## Informed Consent Statement", "## Data Availability Statement",
        "## Acknowledgments", "## Conflicts of Interest"]
for d in DECL:
    check(d in ms, "declaration block present: " + d.strip("# "))
# headings may carry a section number ("5. Conclusions")
check(re.search(r"^##\s+(?:\d+\.\s*)?Conclusions\s*$", ms, flags=re.M) is not None
      and re.search(r"^##\s+(?:\d+\.\s*)?Conclusion\s*$", ms, flags=re.M) is None,
      "conclusion heading is 'Conclusions' (Medicina style)")
check("Case Report" not in
      re.split(r"^##\s+(?:\d+\.\s*)?Introduction\s*$", ms, flags=re.M)[0],
      "submitted as Article, not Case Report")

# --------------------------------------------------------------- 5. references
refs = ms.split("## References", 1)[-1].split("## Figure Legends", 1)[0]
nums = [int(m) for m in re.findall(r"^\s*(\d+)\.\s", refs, flags=re.M)]
check(len(nums) >= 30, "reference list has >= 30 entries", "found %d" % len(nums))
check(nums == list(range(1, len(nums) + 1)),
      "reference numbering is contiguous from 1", "first gaps at %s" % nums[:5])
lines.append("        references: %d" % len(nums))

# ----------------------------------------------------------------- 6. TBD gap
tbd = ms.count("TBD:")
cl_raw = read(CL_MD)
cl_ph = cl_raw.count("[CORRESPONDING AUTHOR")
cl_tbd = cl_raw.count("TBD:")
check(tbd == 0 and cl_ph == 0 and cl_tbd == 0,
      "correspondence details filled (no placeholder in manuscript or cover letter)",
      "manuscript 'TBD:'=%d, cover-letter placeholder=%d, cover-letter 'TBD:'=%d"
      % (tbd, cl_ph, cl_tbd))

# ------------------------------------------------- 6b. section cross-references
# Every "Section 2.13" style reference must point at a heading that exists.
_heads = {m.group(1) for m in
          re.finditer(r"^#{2,3}\s+(\d+(?:\.\d+)?)\.?\s+", ms, flags=re.M)}
_refs = {r for r in re.findall(r"Section (\d+(?:\.\d+)?)", ms)}
_missing = sorted(_refs - _heads, key=lambda x: [int(i) for i in x.split(".")])
check(not _missing, "all 'Section X' cross-references resolve to a heading",
      "missing: %s" % (", ".join(_missing) if _missing else "-"))
lines.append("        section refs: %d cited, %d headings, %d broken"
             % (len(_refs), len(_heads), len(_missing)))

# ------------------------------------------------------------------- 7. docx
doc = Document(MS_DOCX)
texts = [p.text for p in doc.paragraphs]
resid = [t for t in texts if "**" in t or t.strip().startswith("|")]
check(not resid, "docx has no residual markdown markers", str(resid[:2]))
check(len(doc.tables) == 3, "docx carries the three manuscript tables",
      "found %d" % len(doc.tables))
check(len(doc.paragraphs) > 150, "docx body present",
      "paragraphs=%d" % len(doc.paragraphs))

# Word auto-numbering is shared per numbering instance, so the 17 disclosures in
# Section 2.11 and the reference list once ended up in one sequence: 2.11 showed
# 1-17 and the references started at 18. The converter now types every marker as
# literal text, so no paragraph may carry a numbering style any more.
_numbered = [p.text[:40] for p in doc.paragraphs if p.style.name == "List Number"]
check(not _numbered, "no Word auto-numbered list in docx (markers typed as text)",
      "found %d, e.g. %r" % (len(_numbered), _numbered[:1]))
_dbl = [p.text[:40] for p in doc.paragraphs
        if p.style.name == "List Bullet"
        and re.match(r"^\((?:[0-9]{1,2}|[a-z]|[ivx]{1,6})\)", p.text.strip())]
check(not _dbl, "no list item carries a bullet AND its own enumerator",
      "e.g. %r" % (_dbl[:2],))

# The reference list must read 1..N in the rendered file, not continue the
# numbering of an earlier list.
_ref_block, _in_refs = [], False
for p in doc.paragraphs:
    if p.style.name.startswith("Heading"):
        _in_refs = p.text.strip() == "References"
        continue
    if _in_refs and p.text.strip():
        _ref_block.append(p.text.strip())
_rnums = [int(m.group(1)) for m in
          (re.match(r"^(\d+)\.\s", t) for t in _ref_block) if m]
check(_rnums == list(range(1, len(_rnums) + 1)) and len(_rnums) >= 30,
      "docx reference list is numbered 1..N from 1",
      "count=%d, first=%s" % (len(_rnums), _rnums[:5]))
lines.append("        docx: %d paragraphs, %d tables, %d refs"
             % (len(doc.paragraphs), len(doc.tables), len(_rnums)))
cl = Document(CL_DOCX)
cl_text = "\n".join(p.text for p in cl.paragraphs)

# ------------------------------------------------------------- 8. cover letter
check("not under consideration" in cl_text.lower(),
      "cover letter states it is not under consideration elsewhere")
check("All authors have approved" in cl_text,
      "cover letter states all authors approved")
check("Medicina (MDPI)" in cl_text, "cover letter names Medicina (MDPI)")

# -------------------------------------------------------------- 9. figures
for i in range(1, 5):
    a = os.path.join(PKG, "Figure_%d.pdf" % i)
    b = os.path.join(FIG, ["Fig1_global_atlas.pdf", "Fig2_myeloid_depletion.pdf",
                           "Fig3_function.pdf", "Fig4_validation.pdf"][i - 1])
    ha = hashlib.sha256(open(a, "rb").read()).hexdigest()
    hb = hashlib.sha256(open(b, "rb").read()).hexdigest()
    check(ha == hb, "Figure_%d.pdf byte-identical to source" % i)
for i in range(1, 7):
    a = os.path.join(PKG, "Supplementary_Figure_S%d.pdf" % i)
    b = os.path.join(FIG, "Supplementary_Figure_S%d.pdf" % i)
    ha = hashlib.sha256(open(a, "rb").read()).hexdigest()
    hb = hashlib.sha256(open(b, "rb").read()).hexdigest()
    check(ha == hb, "Supplementary_Figure_S%d byte-identical to source" % i)

# ------------------------------------------------------- 10. supplementary ZIP
check(os.path.isfile(ZIP), "shared supplementary archive present")
if os.path.isfile(ZIP):
    lines.append("        ZIP %.2f MB" % (os.path.getsize(ZIP) / 1048576.0))

out = "\n".join(lines)
print(out)
with open(os.path.join(PROJECT, "logs", "Medicina_precheck.log"), "w",
          encoding="utf-8") as fh:
    fh.write("=== Medicina pre-submission self-check ===\n" + out + "\n")
    fh.write("\nRESULT: %s\n" % ("ALL PASS" if not fails
                                 else "BLOCKED (%d)" % len(fails)))
print("\nRESULT: " + ("ALL PASS" if not fails
                      else "BLOCKED (%d) - see FAIL lines above" % len(fails)))
sys.exit(0 if not fails else 1)
