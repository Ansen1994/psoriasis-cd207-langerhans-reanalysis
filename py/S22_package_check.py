# -*- coding: utf-8 -*-
## S22_package_check.py - post-rebuild package verification (replaces the
## ad-hoc 36-assertion run). ASCII output only.
import os, sys, zipfile, datetime
R = "."
md  = "manuscript/submission/CCID/CCID_manuscript.md"
dx  = "manuscript/submission/CCID/CCID_manuscript.docx"
cl  = "manuscript/submission/CCID/CCID_cover_letter.docx"
zipf = "manuscript/submission/Supplementary_materials_20261112.zip"
ok = True
def chk(name, cond, detail=""):
    global ok
    print(("[PASS] " if cond else "[FAIL] ") + name + (" | " + detail if detail else ""))
    if not cond: ok = False

text = open(md, encoding="utf-8").read()
for key in ["### 2.15", "### 3.9", "Supplementary Figure S7:", "Supplementary Figure S8:", "Supplementary Figure S9:",
            "10.5281/zenodo.24369977", "(xvi)", "(xvii)", "GSE226244", "GSE183047",
            "26.4% (NN)", "0.589"]:
    chk("md contains: " + key, key in text)
chk("md legend count S1-S8", all(f"Supplementary Figure S{i}:" in text for i in range(1, 10)))
chk("md no author placeholder", "[AUTHOR NAME 1]" not in text)
chk("md no CHF/APC leftover", "CHF" not in text and "Article Processing Charges" not in text)
chk("md no CN strategy tokens", "Medicina" not in text)

import docx
d = docx.Document(dx)
dtext = "\n".join(p.text for p in d.paragraphs)
chk("docx paragraphs > 200", len(d.paragraphs) > 200, str(len(d.paragraphs)))
for key in ["2.15", "3.9", "Supplementary Figure S8", "GSE226244", "zenodo.24369977"]:
    chk("docx contains: " + key, key in dtext)
m = max(os.path.getmtime(dx), os.path.getmtime(md))
chk("docx newer than md", os.path.getmtime(dx) >= os.path.getmtime(md) - 5)
dcl = docx.Document(cl)
cltext = "\n".join(p.text for p in dcl.paragraphs)
chk("cover letter: submitting author", "Huiyu Zhang" in cltext)
chk("cover letter: corresponding author", "Dongmei Shi" in cltext)
chk("cover letter: treatment recovery highlighted",
    ("IL-17A" in cltext or "secukinumab" in cltext.lower()))

figs = [f for f in os.listdir("manuscript/submission/CCID") if f.endswith(".png")]
chk("CCID folder has 13 figure PNGs", len(figs) == 13, str(len(figs)))
from PIL import Image
for f in figs:
    im = Image.open(os.path.join("manuscript/submission/CCID", f))
    chk("no alpha: " + f, im.mode == "RGB")
    dpix = im.info.get("dpi", (0, 0))
    chk("600dpi: " + f, dpix[0] >= 599.9, str(dpix))

chk("supp zip exists", os.path.exists(zipf))
zf = zipfile.ZipFile(zipf)
names = zf.namelist()
need = ["tables/S18_deconvolution_tests.csv", "tables/S19_treatment_bulk_tests.csv",
        "tables/S20_treatment_sc_samples.csv",
        "tables/S23_recovery_timeline.csv",
        "tables/S23_pasi_correlation.csv", "tables/S17_sub2_trend_genes_source.csv",
        "figures/Supplementary_Figure_S7.pdf", "figures/Supplementary_Figure_S8.pdf",
        "figures/Supplementary_Figure_S9.pdf",
        "R/S18_deconvolution.R", "R/S19_treatment_bulk.R", "R/S20_treatment_sc.R",
        "R/S21_new_figures.R", "R/S23_durability_pasi.R", "logs/S18_deconvolution.log",
        "logs/S23_durability_pasi.log"]
for n in need:
    chk("zip has: " + n, n in names)
bad = [n for n in names if "Ansheng" in n or "cqmu" in n or "CCID" in n or "checklist" in n]
chk("zip no forbidden tokens", not bad, str(bad[:3]))
chk("zip entry count > 125", len(names) > 125, str(len(names)))

print("\nRESULT:", "ALL PASS" if ok else "FAILURES PRESENT")
sys.exit(0 if ok else 1)
