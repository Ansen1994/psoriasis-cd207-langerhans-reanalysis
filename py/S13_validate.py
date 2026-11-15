# -*- coding: utf-8 -*-
"""S13 validation: every claim added to the manuscript is traced back to a table.

Run from the project root:
    python py/S13_validate.py

Exit code 0 and "ALL PASS" only when every check succeeds.
"""
import csv
import hashlib
import os
import re
import sys
import zipfile

TAB = "results/tables"
FIG = "results/figures"
SUB = "manuscript/submission"
ZIP = os.path.join(SUB, "Supplementary_materials_20261112.zip")

ms_files = [os.path.join(SUB, "Genes", "Genes_manuscript.md"),
            os.path.join(SUB, "CCID", "CCID_manuscript.md")]
ms = {f: open(f, encoding="utf-8").read() for f in ms_files}

fails = []
lines = []


def check(ok, msg):
    lines.append(("[PASS] " if ok else "[FAIL] ") + msg)
    if not ok:
        fails.append(msg)


def rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- 1. tables
EXPECT = ["S13_tf_network_summary.csv", "S13_tf_activity_samples.csv",
          "S13_tf_activity_tests.csv", "S13_tf_contrast.csv",
          "S13_tf_depth_confounding.csv", "S13_tf_consistency.csv"]
for name in EXPECT:
    p = os.path.join(TAB, name)
    check(os.path.isfile(p), "table present: " + name)
    if os.path.isfile(p):
        r = rows(p)
        lines.append("        %-40s %5d rows" % (name, len(r)))

tests = rows(os.path.join(TAB, "S13_tf_activity_tests.csv"))
key = {(r["layer"], r["network"], r["tf"]): r for r in tests}

# ------------------------------------------------- 2. numbers used in text
HITS = {"IRF1": 1.91, "STAT1": 1.79, "STAT2": 2.05, "STAT3": 1.48,
        "NFKB1": 2.45, "RELB": 1.32, "JUN": 1.45, "JUND": 1.02, "NFKB2": -0.86}
for tf, val in HITS.items():
    r = key.get(("myeloid_all", "AB", tf))
    ok = r is not None and abs(round(float(r["paired_delta"]), 2) - val) <= 0.005
    check(ok, "paired delta %s = %+.2f matches manuscript" % (tf, val))

robust = [r for r in tests if r["robust"] == "TRUE"]
n_rob_mye = sum(1 for r in robust if r["layer"] == "myeloid_all" and r["network"] == "AB")
n_sig_mye = sum(1 for r in tests if r["layer"] == "myeloid_all" and r["network"] == "AB"
                and float(r["paired_FDR"]) < 0.05)
check(n_rob_mye == 12, "12 robust TFs in pooled A+B layer (found %d)" % n_rob_mye)
check(n_sig_mye == 23, "23 paired-significant TFs in pooled A+B (found %d)" % n_sig_mye)

for tf in ("FOS", "RELA"):
    r = key.get(("myeloid_all", "AB", tf))
    check(r is not None and abs(float(r["rho_delta_depth"])) >= 0.6,
          "%s excluded by |rho| >= 0.6 (rho = %s)" % (tf, r and r["rho_delta_depth"]))
r = key.get(("myeloid_all", "AB", "ARNTL"))
check(r is not None and r["robust"] == "FALSE" and float(r["paired_p"]) > 0.05,
      "ARNTL is not robust despite depth-adjusted significance")

# specification dependence
for lay, net in (("myeloid_all", "ABC"), ("subcluster_0", "AB"), ("subcluster_0", "ABC")):
    n = sum(1 for x in tests if x["layer"] == lay and x["network"] == net
            and x["robust"] == "TRUE")
    check(n == 0, "no robust TF in %s | %s" % (lay, net))

cons = rows(os.path.join(TAB, "S13_tf_consistency.csv"))
cmap = {(r["network"], r["check"]): r for r in cons}
for net, chk, rho in (("AB", "sample_level_myeloid_all_vs_subcluster0", 0.82),
                      ("AB", "sample_level_vs_contrast_level", 0.54)):
    r = cmap.get((net, chk))
    check(r is not None and abs(round(float(r["rho"]), 2) - rho) <= 0.005,
          "consistency %s rho = %.2f" % (chk, rho))

# ---------------------------------------------------------- 3. manuscript
for f, txt in ms.items():
    tag = os.path.basename(f)
    for probe in ("Post hoc transcription-factor activity inference (S13)",
                  "Post hoc transcription-factor activity: a coherent",
                  "(xv) The TF-activity layer (Section 3.8)",
                  "**Supplementary Figure S6. Transcription-factor activity",
                  "**Supplementary Figure S5. Cluster-level atlas"):
        check(probe in txt, tag + " contains: " + probe[:46])
    check("Supplementary Figures S1–S6" in txt, tag + " availability lists S1-S6")
    check("S13" in txt, tag + " mentions stage S13")
    # numbers must agree with the tables
    for probe in ("IRF1 +1.91", "STAT1 +1.79", "STAT2 +2.05", "STAT3 +1.48",
                  "NFKB1 +2.45", "RELB +1.32", "JUND +1.02", "NFKB2 moving in the opposite direction (−0.86)"):
        check(probe in txt, tag + " states: " + probe)
    check("12 of the 23" in txt, tag + " states 12 of 23 robust")
    check("23 of 111" in txt, tag + " states 23 of 111")
    check("none** meeting all three criteria" in txt or "**none** meeting all three criteria" in txt,
          tag + " states A+B+C yields none robust")

# ------------------------------------------------------------- 4. figures
for ext in (".pdf", ".png"):
    p = os.path.join(FIG, "Supplementary_Figure_S6" + ext)
    check(os.path.isfile(p), "figure present: Supplementary_Figure_S6" + ext)
    if os.path.isfile(p):
        lines.append("        %-40s %8d bytes" % (os.path.basename(p), os.path.getsize(p)))
check(os.path.getsize(os.path.join(FIG, "Supplementary_Figure_S6.png")) > 150000,
      "S6 PNG is 300 dpi raster (>150 kB)")

# ------------------------------------------------------------ 5. packages
for journal in ("Genes", "CCID"):
    d = os.path.join(SUB, journal)
    files = sorted(os.listdir(d))
    check(len(files) == 12, "%s package has 12 files (found %d)" % (journal, len(files)))
    check("Supplementary_Figure_S6.pdf" in files, journal + " ships S6")
    # main figures byte-identical to source
    for i in range(1, 5):
        a = os.path.join(d, "Figure_%d.pdf" % i)
        b = os.path.join(FIG, ["Fig1_global_atlas.pdf", "Fig2_myeloid_depletion.pdf",
                               "Fig3_function.pdf", "Fig4_validation.pdf"][i - 1])
        ha = hashlib.sha256(open(a, "rb").read()).hexdigest()
        hb = hashlib.sha256(open(b, "rb").read()).hexdigest()
        check(ha == hb, "%s Figure_%d.pdf hash matches source" % (journal, i))
    for i in range(1, 7):
        a = os.path.join(d, "Supplementary_Figure_S%d.pdf" % i)
        b = os.path.join(FIG, "Supplementary_Figure_S%d.pdf" % i)
        ha = hashlib.sha256(open(a, "rb").read()).hexdigest()
        hb = hashlib.sha256(open(b, "rb").read()).hexdigest()
        check(ha == hb, "%s Supplementary_Figure_S%d hash matches source" % (journal, i))

check(os.path.isfile(ZIP), "supplementary ZIP present")
if os.path.isfile(ZIP):
    zf = zipfile.ZipFile(ZIP)
    names = zf.namelist()
    lines.append("        ZIP %d entries, %.2f MB" %
                 (len(names), os.path.getsize(ZIP) / 1048576.0))
    check(len(names) >= 117, "ZIP has at least 117 entries (found %d)" % len(names))
    for probe in ("tables/S13_tf_activity_tests.csv", "figures/Supplementary_Figure_S6.pdf",
                  "R/S13_tf_activity.R", "R/S13b_tf_figures.R", "logs/S13_tf_activity.log"):
        check(probe in names, "ZIP contains " + probe)
    BAD = ("S7_author_info.log", "author_info_v1.md", "Ansheng Xu", "Cheng Chen",
           "Guobing Yin", "Jing Luo", "checklist_Genes.md")
    for n in names:
        check(not any(b in n for b in BAD), "ZIP entry clean: " + n)

# ------------------------------------------------------------ 6. word count
def body_words(path):
    t = open(path, encoding="utf-8").read()
    start = t.find("## Introduction")
    if start < 0:
        start = t.find("## 1. Introduction")
    end = t.find("## References")
    if end < 0:
        end = t.find("## 5. References")
    seg = t[start:end]
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9\-\+\.]*", seg))

for f in ms_files:
    n = body_words(f)
    lines.append("        body words %-30s %6d" % (os.path.basename(f), n))

out = "\n".join(lines)
print(out)
with open(os.path.join("logs", "S13_validation.log"), "w", encoding="utf-8") as fh:
    fh.write("=== S13 validation ===\n" + out + "\n")
    fh.write("\nRESULT: %s\n" % ("ALL PASS" if not fails else "FAILED (%d)" % len(fails)))
print("\nRESULT: " + ("ALL PASS" if not fails else "FAILED (%d)" % len(fails)))
sys.exit(0 if not fails else 1)
