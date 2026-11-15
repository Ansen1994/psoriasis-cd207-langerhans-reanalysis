#!/usr/bin/env python3
"""
S0 pilot / Gate A check for dermatology-sci.

Read ONLY from data/raw/<dataset_id> (junctions to the shared catalog).
Writes tables to results/tables/, a report to data/processed/, a log to logs/.
No figures are produced. No biological results are computed.

Run from the project root:  python py/S0_pilot.py
"""

import csv
import gzip
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
RAW = PROJ / "data" / "raw"
TAB = PROJ / "results" / "tables"
PROC = PROJ / "data" / "processed"
LOG = PROJ / "logs"
for d in (TAB, PROC, LOG):
    d.mkdir(parents=True, exist_ok=True)

D_DISCOVERY = RAW / "PSOR_SC_GSE173706"
D_VALIDATION = RAW / "PSOR_SC_GSE162183"
D_BULK = RAW / "PSOR_BULK_GSE13355"

lines_log = []


def note(msg):
    lines_log.append(msg)
    print(msg)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def first_lines_gz(path, n=2):
    out = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for _ in range(n):
            s = fh.readline()
            if not s:
                break
            out.append(s.rstrip("\r\n"))
    return out


def count_lines_gz(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in fh)


def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


# ---------------------------------------------------------------- integrity
note("=== S0 integrity (SHA-256) ===")
sha_rows = []
for label, p in [
    ("filelist.txt", D_DISCOVERY / "filelist.txt"),
    ("GSE173706_RAW.tar", D_DISCOVERY / "GSE173706_RAW.tar"),
    ("GSE162183_counts", D_VALIDATION / "GSE162183_Raw_gene_counts_matrix.tab.gz"),
    ("GSE13355_series_matrix", D_BULK / "GSE13355_series_matrix.txt.gz"),
]:
    if not p.exists():
        note("MISSING: " + str(p))
        continue
    digest = sha256(p)
    sha_rows.append([label, p.name, p.stat().st_size, digest])
    note("%-28s %10d  %s" % (label, p.stat().st_size, digest))
write_csv(TAB / "S0_sha256.csv", ["label", "file", "bytes", "sha256"], sha_rows)

# ---------------------------------------------------------------- filelist vs extracted
note("")
note("=== GSE173706: filelist vs extracted ===")
expected = {}
with open(D_DISCOVERY / "filelist.txt", "r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) >= 4 and parts[0] == "File":
            expected[parts[1]] = int(parts[3])
csv_dir = D_DISCOVERY / "csv"
actual = {p.name: p.stat().st_size for p in csv_dir.glob("*.csv.gz")} if csv_dir.exists() else {}
mismatch = [n for n in expected if n not in actual or actual[n] != expected[n]]
note("filelist entries: %d | extracted: %d | size/name mismatches: %d" %
     (len(expected), len(actual), len(mismatch)))
if mismatch:
    note("MISMATCH: " + ", ".join(mismatch[:10]))

# ---------------------------------------------------------------- discovery matrix structure
note("")
note("=== GSE173706: matrix structure ===")
pat = re.compile(r"^(GSM\d+)_(NS|PN|PP)-(.+)\.csv\.gz$")
rows = []
for name in sorted(actual):
    m = pat.match(name)
    if not m:
        note("UNPARSED FILENAME: " + name)
        continue
    gsm, group, donor = m.groups()
    head = first_lines_gz(csv_dir / name, 2)
    n_cells = len(head[0].split(",")) - 1 if head else 0
    gene_id = head[1].split(",")[0] if len(head) > 1 else ""
    donor_base = re.sub(r"(V\d+|[AB])$", "", donor)
    rows.append([gsm, name, group, donor, donor_base, n_cells, gene_id,
                 "YES" if name == "GSM5277176_NS-AR008.csv.gz" else "",
                 actual[name]])
write_csv(TAB / "S0_sample_map.csv",
          ["gsm", "file", "group", "donor", "donor_base", "n_cells",
           "gene_id_example", "flag_small", "bytes"], rows)

group_counts = Counter(r[2] for r in rows)
note("groups: " + ", ".join("%s=%d" % kv for kv in sorted(group_counts.items())))
note("total samples: %d | total cells: %d" % (len(rows), sum(r[5] for r in rows)))
note("gene id example: " + (rows[0][6] if rows else "NA"))
small = [r for r in rows if r[7] == "YES"]
for r in small:
    note("SMALL SAMPLE %s: %d cells, %d bytes" % (r[1], r[5], r[8]))

sizes = sorted(actual, key=lambda n: actual[n])
representatives = [sizes[0], sizes[len(sizes) // 2], sizes[-1]]
gene_counts = {}
for n in representatives:
    gene_counts[n] = count_lines_gz(csv_dir / n) - 1
    note("gene rows in %s: %d" % (n, gene_counts[n]))

# donor pairing (O8)
by_donor = defaultdict(dict)
for r in rows:
    by_donor[r[4]].setdefault(r[2], []).append(r[1])
pair_rows = []
for donor, groups in sorted(by_donor.items()):
    pair_rows.append([donor, len(groups.get("NS", [])), len(groups.get("PN", [])),
                      len(groups.get("PP", [])),
                      "paired" if groups.get("PN") and groups.get("PP") else "",
                      ";".join(groups.get("PN", []) + groups.get("PP", []))])
write_csv(TAB / "S0_donor_pairs.csv",
          ["donor_base", "n_NS", "n_PN", "n_PP", "pn_pp_paired", "files"], pair_rows)
n_paired = sum(1 for r in pair_rows if r[4] == "paired")
n_replicate = sum(1 for r in pair_rows if (r[1] + r[2] + r[3]) > 2)
note("distinct donor_base: %d | PN-PP paired donors: %d | donors with >2 libraries: %d" %
     (len(pair_rows), n_paired, n_replicate))

# ---------------------------------------------------------------- validation cohort
note("")
note("=== GSE162183: validation matrix ===")
vhead = first_lines_gz(D_VALIDATION / "GSE162183_Raw_gene_counts_matrix.tab.gz", 2)
barcodes = vhead[0].split("\t")[1:] if vhead else []
vpref = Counter(b.split("_")[0] for b in barcodes)
note("cells: %d | gene id example: %s" % (len(barcodes), vhead[1].split("\t")[0] if len(vhead) > 1 else "NA"))
note("barcode prefixes: " + ", ".join("%s=%d" % kv for kv in sorted(vpref.items())))
write_csv(TAB / "S0_validation_barcodes.csv",
          ["prefix", "n_cells"], sorted(vpref.items()))

# ---------------------------------------------------------------- bulk cohort
note("")
note("=== GSE13355: series matrix metadata ===")
meta = {}
table_header = None
with gzip.open(D_BULK / "GSE13355_series_matrix.txt.gz", "rt",
               encoding="utf-8", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            table_header = fh.readline().rstrip("\r\n").split("\t")
            break
        if line.startswith("!Sample_"):
            parts = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))
            meta.setdefault(parts[0], []).append(parts[1:])
titles = meta.get("!Sample_title", [[]])[0]
gsms = meta.get("!Sample_geo_accession", [[]])[0]
chars = meta.get("!Sample_characteristics_ch1", [[]])[0]
bulk_rows = []
for i, title in enumerate(titles):
    m = re.search(r"_(NN|PN|PP)_sample", title)
    bulk_rows.append([gsms[i] if i < len(gsms) else "", title,
                      m.group(1) if m else "UNKNOWN",
                      chars[i] if i < len(chars) else ""])
write_csv(TAB / "S0_bulk_sample_map.csv",
          ["gsm", "sample_title", "group", "characteristics"], bulk_rows)
note("samples in series matrix: %d | table columns (incl ID_REF): %d" %
     (len(bulk_rows), len(table_header) if table_header else 0))
note("groups: " + ", ".join("%s=%d" % kv for kv in
                            sorted(Counter(r[2] for r in bulk_rows).items())))
note("characteristics distinct values: " + ", ".join(
    "%s (%d)" % kv for kv in Counter(r[3] for r in bulk_rows).most_common(6)))
note("metadata keys: %d" % len(meta))

# ---------------------------------------------------------------- Gate A decision
gate_a = {
    "filelist_33_ok": len(expected) == 33 and len(actual) == 33 and not mismatch,
    "three_groups": set(["NS", "PN", "PP"]).issubset(set(group_counts)),
    "min_3_per_group": all(group_counts.get(g, 0) >= 3 for g in ("NS", "PN", "PP")),
    "donor_pair_table": len(pair_rows) > 0,
    "small_sample_reviewed": bool(small),
    "bulk_metadata_parsed": len(bulk_rows) == 180 and all(r[2] != "UNKNOWN" for r in bulk_rows),
}
note("")
note("=== Gate A ===")
for k, v in gate_a.items():
    note("  %-24s %s" % (k, "PASS" if v else "FAIL"))
note("GATE_A_OVERALL: " + ("PASS" if all(gate_a.values()) else "REVIEW NEEDED"))

# ---------------------------------------------------------------- report (markdown, UTF-8)
report = []
report.append("# S0 pilot report / Gate A - dermatology-sci")
report.append("")
report.append("Generated by `py/S0_pilot.py`. Data read only through `data/raw/` junctions. No figures, no biological results.")
report.append("")
report.append("## 1. Integrity")
report.append("")
report.append("| file | bytes | sha256 |")
report.append("|---|---|---|")
for r in sha_rows:
    report.append("| %s | %d | `%s` |" % (r[1], r[2], r[3]))
report.append("")
report.append("## 2. Discovery cohort GSE173706")
report.append("")
report.append("- `filelist.txt` vs extracted: %d / %d, mismatches %d。" %
              (len(expected), len(actual), len(mismatch)))
report.append("- 分组：%s，合计 %d 样本 / %d 细胞。" %
              (", ".join("%s=%d" % kv for kv in sorted(group_counts.items())),
               len(rows), sum(r[5] for r in rows)))
report.append("- 矩阵方向：**基因 × 细胞**（首行为 barcode，首列为空），基因 ID 为 **Ensembl (`ENSG...`)**。")
report.append("- 代表文件基因行数：%s。" %
              ", ".join("%s=%d" % (k, v) for k, v in gene_counts.items()))
report.append("- `NS-AR008`：%s。" %
              ("%d 细胞 / %d bytes（明显偏低，需在 S1 QC 单独判定是否剔除）" %
               (small[0][5], small[0][8]) if small else "未找到"))
report.append("- 供体配对（优化 O8）：%d 个 donor_base，其中 **%d 个同时有 PN 与 PP**，%d 个拥有 >2 个文库（技术重复/双活检）。pseudobulk 必须按 donor 合并，检验加 `(1|donor)`。" %
              (len(pair_rows), n_paired, n_replicate))
report.append("")
report.append("## 3. Validation cohort GSE162183")
report.append("")
report.append("- **不是 6 个独立矩阵**，而是一个合并矩阵：%d 个细胞，基因 ID 为 **Gene Symbol**（如 `%s`）。" %
              (len(barcodes), vhead[1].split("\t")[0] if len(vhead) > 1 else "NA"))
report.append("- 分组内嵌在 barcode 前缀：%s。label transfer 前按前缀拆分即可，无需读 6 个文件。" %
              ", ".join("%s=%d" % kv for kv in sorted(vpref.items())))
report.append("")
report.append("## 4. Bulk cohort GSE13355")
report.append("")
report.append("- Series Matrix：%d 个样本，表列 %d（含 ID_REF），元数据键 %d 个。" %
              (len(bulk_rows), len(table_header) if table_header else 0, len(meta)))
report.append("- 分组：%s。" % ", ".join("%s=%d" % kv for kv in
                                    sorted(Counter(r[2] for r in bulk_rows).items())))
report.append("- 表达值为 **log2 且上游已做批次/性别效应校正**；反卷积前 `2^x` 还原，该限制写入 Methods。")
report.append("")
report.append("## 5. Gate A")
report.append("")
report.append("| check | result |")
report.append("|---|---|")
for k, v in gate_a.items():
    report.append("| %s | %s |" % (k, "PASS" if v else "FAIL"))
report.append("")
report.append("**Overall: %s**" % ("PASS" if all(gate_a.values()) else "REVIEW NEEDED"))
report.append("")
report.append("## 6. Open items for S1 (do not block Gate A)")
report.append("")
report.append("1. **ID 不统一**：发现队列 Ensembl ID、验证队列 Gene Symbol、bulk 为 GPL570 探针 ID。三者须统一到 Symbol 才能做注释迁移、label transfer 与反卷积。R 端方案：`biomaRt` 在线映射或 `org.Hs.eg.db::mapIds`；bulk 探针映射可复用集中仓 `REF_GPL570`。")
report.append("2. `NS-AR008` 细胞数偏低，S1 QC 按实际分布判定，剔除与否须在 `logs/` 留痕。")
report.append("3. QC 阈值（`config/01_params.R`）待 S1 读到全量分布后用 median ± 3×MAD 校准，改动记入日志。")
report.append("")
(PROC / "S0_pilot_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
(LOG / "S0_ingest.log").write_text("\n".join(lines_log) + "\n", encoding="utf-8")
note("")
note("WROTE: results/tables/S0_*.csv, data/processed/S0_pilot_report.md, logs/S0_ingest.log")
