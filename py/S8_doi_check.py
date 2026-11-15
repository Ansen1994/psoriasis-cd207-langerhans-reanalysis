#!/usr/bin/env python3
"""S8 DOI verification for the dermatology-sci manuscript (v1.1).

Queries the CrossRef REST API by bibliographic string and records the best
match with a token-overlap score. Nothing is invented: when confidence is
below the threshold the row is marked REVIEW and no DOI is transcribed into
the manuscript. Software / book / R citations are marked NO_DOI_EXPECTED.

Usage (from the project root or anywhere):
    python py/S8_doi_check.py

Outputs:
    results/tables/S8_doi_check.csv
    logs/S8_doi_check.log
"""

import csv
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_CSV = os.path.join(ROOT, "results", "tables", "S8_doi_check.csv")
OUT_LOG = os.path.join(ROOT, "logs", "S8_doi_check.log")
API = "https://api.crossref.org/works"

# (ref_id, exact title, expected first author surname, expected year, note)
QUERIES = [
    ("R01", "Psoriasis", "Griffiths", "2021", "Lancet review; short title, verify manually"),
    ("R03", "Increase in TNF-alpha and inducible nitric oxide synthase-expressing dendritic cells in psoriasis and reduction with efalizumab (anti-CD11a)", "Lowes", "2005", ""),
    ("R07", "NCBI GEO: archive for functional genomics data sets--update", "Barrett", "2013", ""),
    ("R08", "org.Hs.eg.db: Genome wide annotation for Human", "Carlson", "2024", "R package, no DOI expected"),
    ("R09", "Doublet identification in single-cell sequencing data using scDblFinder", "Germain", "2021", ""),
    ("R10", "Integrated analysis of multimodal single-cell data", "Hao", "2021", ""),
    ("R11", "Fast, sensitive and accurate integration of single-cell data with Harmony", "Korsunsky", "2019", ""),
    ("R12", "Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage", "Aran", "2019", ""),
    ("R13", "A distribution-free k-sample test against ordered alternatives", "Jonckheere", "1954", "old article; may lack DOI"),
    ("R14", "Controlling the false discovery rate: a practical and powerful approach to multiple testing", "Benjamini", "1995", ""),
    ("R15", "Fitting linear mixed-effects models using lme4", "Bates", "2015", ""),
    ("R16", "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2", "Love", "2014", ""),
    ("R17", "clusterProfiler: an R package for comparing biological themes among gene clusters", "Yu", "2012", ""),
    ("R18", "clusterProfiler 4.0: A universal enrichment tool for interpreting omics data", "Wu", "2021", ""),
    ("R19", "Fast gene set enrichment analysis", "Korotkevich", "2021", "bioRxiv preprint"),
    ("R20", "The Molecular Signatures Database Hallmark gene set collection", "Liberzon", "2015", ""),
    ("R21", "GSVA: gene set variation analysis for microarray and RNA-seq data", "Haenzelmann", "2013", ""),
    ("R22", "Comprehensive Integration of Single-Cell Data", "Stuart", "2019", ""),
    ("R23", "pROC: an open-source package for R and S+ to analyze and compare ROC curves", "Robin", "2011", ""),
    ("R24", "ggplot2: Elegant Graphics for Data Analysis", "Wickham", "2016", "book, DOI optional"),
    ("R25", "R: A Language and Environment for Statistical Computing", "R Core Team", "2024", "software, no DOI expected"),
    ("R26", "The epidermal immune microenvironment plays a dominant role in psoriasis development, as revealed by mass cytometry", "Zhou", "2022", ""),
    ("R27", "Signaling pathways and targeted therapies for psoriasis", "Guo", "2023", ""),
    ("N01", "Specific roles for dendritic cell subsets during initiation and progression of psoriasis", "Glitzner", "2014", "NEW: strongest prior experimental evidence"),
    ("N02", "huSA: a comprehensive database for multi-dimensional resolution of bulk, single cell and spatial transcription profiles in skin diseases", "Zheng", "2026", "NEW: same GSE173706 direction"),
    ("N03", "Distinct human Langerhans cell subsets orchestrate reciprocal functions and require different developmental regulation", "Liu", "2021", "NEW: LC1/LC2/aLC/migLC subsets"),
    ("N04", "Distinct migration behaviors of monocyte-derived versus resident Langerhans cells drive the pathogenesis of imiquimod-induced skin inflammation", "Cai", "2026", "NEW: migration mechanism, in press"),
    ("N05", "A transcriptional continuum from clinically normal to lesional skin: Single-cell trajectory analysis of psoriasis", "", "2026", "NEW: bioRxiv preprint; CrossRef may not index it"),
    ("N06", "Single-cell and bulk RNAseq unveils the immune infiltration landscape and targeted therapeutic biomarkers of psoriasis", "", "2024", "NEW: sc + GSE13355 precedent"),
]

NO_DOI_EXPECTED = {"R08", "R25"}


def norm(s):
    return set(re.findall(r"[a-z0-9]+", (s or "").lower()))


def overlap(a, b):
    A, B = norm(a), norm(b)
    if not A or not B:
        return 0.0
    return len(A & B) / float(max(len(A), len(B)))


def fetch(query):
    qs = urllib.parse.urlencode({
        "query.bibliographic": query,
        "rows": 3,
        "select": "DOI,title,container-title,volume,issue,page,published,author,type",
    })
    req = urllib.request.Request(API + "?" + qs,
                                 headers={"User-Agent": "dermatology-sci-doi-check/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.load(resp)["message"]["items"]


def year_of(item):
    dp = item.get("published", {}).get("date-parts", [[None]])
    try:
        return str(dp[0][0])
    except Exception:
        return ""


def main():
    rows = []
    log = []
    for rid, title, author, year, note in QUERIES:
        if rid in NO_DOI_EXPECTED:
            rows.append(dict(ref_id=rid, title=title, expected_first_author=author,
                             expected_year=year, crossref_title="", doi="",
                             journal="", volume="", issue="", page="", year="",
                             first_author="", score="", status="NO_DOI_EXPECTED",
                             note=note))
            log.append("%s NO_DOI_EXPECTED (%s)" % (rid, note))
            continue
        query = title + " " + author + " " + year
        try:
            items = fetch(query)
        except Exception as exc:
            rows.append(dict(ref_id=rid, title=title, expected_first_author=author,
                             expected_year=year, crossref_title="", doi="",
                             journal="", volume="", issue="", page="", year="",
                             first_author="", score="", status="FETCH_FAILED",
                             note=note + " | " + str(exc)[:120]))
            log.append("%s FETCH_FAILED %s" % (rid, exc))
            time.sleep(0.3)
            continue
        best, best_s = None, -1.0
        for it in items:
            ct = (it.get("title") or [""])[0]
            s = overlap(title + " " + author, ct)
            fa = ""
            au = it.get("author") or []
            if au:
                fa = au[0].get("family", "")
            if author and fa and author.lower() not in fa.lower() and fa.lower() not in author.lower():
                s -= 0.20
            if year and year_of(it) and year_of(it) != year:
                s -= 0.10
            if s > best_s:
                best, best_s = it, s
        if best is None:
            rows.append(dict(ref_id=rid, title=title, expected_first_author=author,
                             expected_year=year, crossref_title="", doi="",
                             journal="", volume="", issue="", page="", year="",
                             first_author="", score="", status="NO_MATCH", note=note))
            log.append("%s NO_MATCH" % rid)
        else:
            au = best.get("author") or []
            fa = au[0].get("family", "") if au else ""
            status = "MATCH" if best_s >= 0.72 else "REVIEW"
            rows.append(dict(
                ref_id=rid, title=title, expected_first_author=author,
                expected_year=year,
                crossref_title=(best.get("title") or [""])[0],
                doi=best.get("DOI", ""),
                journal=(best.get("container-title") or [""])[0],
                volume=best.get("volume", ""), issue=best.get("issue", ""),
                page=best.get("page", ""), year=year_of(best),
                first_author=fa, score=round(best_s, 3), status=status,
                note=note))
            log.append("%s %-14s score=%.3f doi=%s | %s" %
                       (rid, status, best_s, best.get("DOI", ""),
                        (best.get("title") or [""])[0][:90]))
        time.sleep(0.25)

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(OUT_LOG, "w", encoding="utf-8") as fh:
        fh.write("S8 DOI verification - CrossRef REST API\n")
        fh.write("rows=%d  MATCH=%d  REVIEW=%d  NO_MATCH=%d  FETCH_FAILED=%d\n\n" % (
            len(rows),
            sum(r["status"] == "MATCH" for r in rows),
            sum(r["status"] == "REVIEW" for r in rows),
            sum(r["status"] == "NO_MATCH" for r in rows),
            sum(r["status"] == "FETCH_FAILED" for r in rows)))
        fh.write("\n".join(log) + "\n")
    print("WROTE", OUT_CSV)
    print("\n".join(log))
    return 0


if __name__ == "__main__":
    sys.exit(main())
