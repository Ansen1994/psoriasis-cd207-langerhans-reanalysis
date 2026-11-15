"""S4_direct_markers.py - label-transfer-free external direction check.

Motivation (recorded in logs/gates.md): the S4 Seurat label transfer used a
myeloid-only reference against a WHOLE-SKIN query. Every query cell therefore
received a myeloid label, 98% were called macrophage/DC, and the
prediction.score.max filter (0.3/0.5/0.7) removed only 0.1-4% of cells, so it
cannot distinguish a genuine myeloid cell from a keratinocyte. The internal
5-fold CV AUC (0.9988) measures reproducibility of the same clustering, not
external validation.

This script provides an independent, mapping-free readout: direct count of
CD207+ cells and of myeloid-marker-positive cells in GSE162183, by sample and
group. The rule is fixed here, applied identically to all six samples, and
does not use any disease label to define it. No parameter is tuned.

Read-only on data/raw/. Output: results/tables/S4_query_direct_markers.csv
"""
import gzip
import csv
import os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "raw", "PSOR_SC_GSE162183",
                   "GSE162183_Raw_gene_counts_matrix.tab.gz")
OUT = os.path.join(ROOT, "results", "tables", "S4_query_direct_markers.csv")

TARGET = ["CD207", "FCER1A", "LAMP3", "LYZ", "C1QA", "CD14", "AIF1",
          "CD68", "HLA-DQB2", "HPGDS", "KRT14", "DCD", "KRT19"]
MYELOID = ["LYZ", "C1QA", "CD14", "AIF1", "CD68", "FCER1A", "CD207"]


def main():
    vals = {}
    with gzip.open(SRC, "rt", errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        barcodes = hdr[1:]
        for line in fh:
            p = line.rstrip("\n").split("\t")
            g = p[0]
            if g in TARGET and g not in vals:
                row = []
                for x in p[1:]:
                    try:
                        row.append(float(x))
                    except ValueError:
                        row.append(0.0)
                vals[g] = row

    n = len(barcodes)
    present = [g for g in TARGET if g in vals]
    missing = [g for g in TARGET if g not in vals]
    print("genes present:", ", ".join(present))
    print("genes missing:", ", ".join(missing))
    for g in present:
        if len(vals[g]) != n:
            raise RuntimeError("length mismatch for " + g)
    if "CD207" not in vals:
        raise RuntimeError("CD207 absent from validation matrix; cannot run this check")

    cd207 = vals["CD207"]
    myeloid_rows = [vals[g] for g in MYELOID if g in vals]
    fcer1a = vals.get("FCER1A")
    lamp3 = vals.get("LAMP3")

    sample = [b.split("_")[0] for b in barcodes]
    group = ["NS" if s.startswith("Ctrl") else "PP" for s in sample]

    cand = [any(row[i] > 0 for row in myeloid_rows) for i in range(n)]
    cd_pos = [v > 0 for v in cd207]
    lc_pos = [cd_pos[i] or (fcer1a is not None and lamp3 is not None
                            and fcer1a[i] > 0 and lamp3[i] > 0)
              for i in range(n)]

    agg = defaultdict(lambda: defaultdict(int))
    for i in range(n):
        a = agg[sample[i]]
        a["group"] = group[i]
        a["n_cells"] += 1
        a["n_CD207_pos"] += int(cd_pos[i])
        a["n_LC_marker_pos"] += int(lc_pos[i])
        a["n_myeloid_candidate"] += int(cand[i])
        a["n_CD207_in_myeloid"] += int(cd_pos[i] and cand[i])
        a["n_LC_in_myeloid"] += int(lc_pos[i] and cand[i])

    rows = []
    for sid in sorted(agg):
        a = agg[sid]
        rows.append({
            "sample_id": sid,
            "group": a["group"],
            "n_cells": a["n_cells"],
            "n_CD207_pos": a["n_CD207_pos"],
            "pct_CD207_all": round(100 * a["n_CD207_pos"] / a["n_cells"], 4),
            "n_LC_marker_pos": a["n_LC_marker_pos"],
            "pct_LC_marker_all": round(100 * a["n_LC_marker_pos"] / a["n_cells"], 4),
            "n_myeloid_candidate": a["n_myeloid_candidate"],
            "n_CD207_in_myeloid": a["n_CD207_in_myeloid"],
            "pct_CD207_in_myeloid": round(
                100 * a["n_CD207_in_myeloid"] / a["n_myeloid_candidate"], 4)
                if a["n_myeloid_candidate"] else "",
            "n_LC_in_myeloid": a["n_LC_in_myeloid"],
        })

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\nwritten:", OUT)
    for r in rows:
        print("{sample_id:8s} {group:3s} CD207+ {n_CD207_pos:5d}/{n_cells:6d} "
              "({pct_CD207_all:7.4f}%)  myeloid_cand {n_myeloid_candidate:6d} "
              "CD207_in_myeloid {n_CD207_in_myeloid:5d}".format(**r))

    for key in ("pct_CD207_all", "pct_LC_marker_all", "pct_CD207_in_myeloid"):
        for g in ("NS", "PP"):
            v = [r[key] for r in rows if r["group"] == g and r[key] != ""]
            print("group {} {}: n={} values={} median={}".format(
                g, key, len(v), v, sorted(v)[len(v) // 2] if v else "NA"))


if __name__ == "__main__":
    main()
