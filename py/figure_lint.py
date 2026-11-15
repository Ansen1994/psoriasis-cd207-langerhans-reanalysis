#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
figure_lint.py - 图件版式断言检查（把"看图猜"换成"跑脚本列全部违规项"）

用法:
    python py/figure_lint.py                  # 检查全部主图 + 补充图
    python py/figure_lint.py --only Fig3      # 只查文件名含 Fig3 的
    python py/figure_lint.py --main-only      # 只查主图
    python py/figure_lint.py --save-baseline  # 把当前结果存为基线
    python py/figure_lint.py --baseline       # 与基线对比，新增违规即 FAIL

检查项:
    [1] 最小字号 >= MIN_PT (6.0 pt, 印刷可读下限)
    [2] 文字不出画布 (bbox 不得超出 page rect)
    [3] 文字不互相重叠 (两 span 交叠面积 > 较小者的 OVERLAP_FRAC)
    [4] 主图宽度在 175-185 mm (目标 180 mm; MDPI 单栏上限 170 mm)
    [5] 每张 PDF 都有同名 PNG 预览, 且有效分辨率 >= MIN_DPI

退出码: 0 = 全 PASS, 1 = 有 FAIL
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    import pymupdf as fitz
except ImportError:  # pragma: no cover
    try:
        import fitz
    except ImportError:
        sys.exit("需要 pymupdf:  pip install pymupdf")

MIN_PT = 6.0
OVERLAP_FRAC = 0.50   # 交叠占较小"墨迹框"的比例；真碰撞接近 100%
INK_SHRINK = 0.40     # em 框上下各裁 20%，近似 cap-height 墨迹带
EDGE_TOL_PT = 1.5
MIN_DPI = 200.0      # PNG 预览的有效分辨率下限（目标 300）
WIDTH_RANGE_MM = (175.0, 185.0)
PT_PER_MM = 72.0 / 25.4

ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "results" / "figures"
LOG_DIR = ROOT / "logs"
BASELINE = LOG_DIR / "figure_lint_baseline.json"


def collect(path: Path):
    """返回 (页面尺寸(pt), span 列表)。span = dict(text,size,bbox)。"""
    doc = fitz.open(str(path))
    page = doc[0]
    rect = page.rect
    spans = []
    for blk in page.get_text("dict")["blocks"]:
        for line in blk.get("lines", []):
            dx, dy = line.get("dir", (1.0, 0.0))
            for sp in line.get("spans", []):
                txt = sp["text"].strip()
                if not txt:
                    continue
                x0, y0, x1, y1 = sp["bbox"]
                spans.append(
                    {"text": txt, "size": round(sp["size"], 2),
                     "bbox": (x0, y0, x1, y1), "dir": (round(dx, 3), round(dy, 3))}
                )
    n_pages = doc.page_count
    doc.close()
    return rect, spans, n_pages


def rect_overlap(a, b):
    """交叠面积占较小 bbox 面积的比例。"""
    ax0, ay0, ax1, ay1 = a
    bx0, by0, bx1, by1 = b
    ix0, iy0 = max(ax0, bx0), max(ay0, by0)
    ix1, iy1 = min(ax1, bx1), min(ay1, by1)
    if ix1 <= ix0 or iy1 <= iy0:
        return 0.0
    inter = (ix1 - ix0) * (iy1 - iy0)
    area_a = (ax1 - ax0) * (ay1 - ay0)
    area_b = (bx1 - bx0) * (by1 - by0)
    smaller = min(area_a, area_b)
    return inter / smaller if smaller > 0 else 0.0


def lint_one(path: Path, check_width: bool):
    rect, spans, n_pages = collect(path)
    w_mm = rect.width / PT_PER_MM
    h_mm = rect.height / PT_PER_MM
    fails = []

    # [1] 字号
    small = [s for s in spans if s["size"] < MIN_PT - 1e-6]
    if small:
        worst = sorted(small, key=lambda s: s["size"])[:5]
        fails.append(
            "[1] 字号 < %.1f pt: %d 处; 最小 %.1f pt (%s)"
            % (MIN_PT, len(small), worst[0]["size"],
               " | ".join("%r@%.1f" % (s["text"][:24], s["size"]) for s in worst[:3]))
        )

    # [2] 出画布
    out = []
    for s in spans:
        x0, y0, x1, y1 = s["bbox"]
        if (x0 < -EDGE_TOL_PT or y0 < -EDGE_TOL_PT
                or x1 > rect.width + EDGE_TOL_PT
                or y1 > rect.height + EDGE_TOL_PT):
            out.append(s)
    if out:
        fails.append(
            "[2] 文字出画布: %d 处; 例 %s"
            % (len(out), " | ".join("%r" % s["text"][:28] for s in out[:3]))
        )

    # [3] 重叠：只用"墨迹框"比，旋转文字跳过（轴对齐 bbox 对旋转文字无意义）
    dups, skipped_rot = [], 0
    horiz = []
    for s in spans:
        dx, dy = s["dir"]
        if abs(dx) > 0.99 and abs(dy) < 0.02:
            x0, y0, x1, y1 = s["bbox"]
            pad = (y1 - y0) * INK_SHRINK / 2.0
            horiz.append((s["text"], (x0, y0 + pad, x1, y1 - pad)))
        else:
            skipped_rot += 1
    for i in range(len(horiz)):
        for j in range(i + 1, len(horiz)):
            frac = rect_overlap(horiz[i][1], horiz[j][1])
            if frac > OVERLAP_FRAC:
                dups.append((horiz[i][0], horiz[j][0], frac))
    if dups:
        dups.sort(key=lambda t: -t[2])
        fails.append(
            "[3] 文字互相重叠: %d 对; 最重 %r x %r (%.0f%%)"
            % (len(dups), dups[0][0][:24], dups[0][1][:24], dups[0][2] * 100)
        )

    # [4] 画布宽度（只查主图）
    if check_width and not (WIDTH_RANGE_MM[0] <= w_mm <= WIDTH_RANGE_MM[1]):
        fails.append(
            "[4] 画布宽 %.1f mm 超出 %.0f-%.0f mm"
            % (w_mm, WIDTH_RANGE_MM[0], WIDTH_RANGE_MM[1])
        )

    # [5] 每张 PDF 必须同时有同名 PNG 预览（300 dpi），避免只有矢量图没有预览
    png = path.with_suffix(".png")
    dpi = None
    if not png.exists():
        fails.append("[5] 缺同名 PNG 预览: %s" % png.name)
    else:
        try:
            from PIL import Image
            px_w = Image.open(png).width
            dpi = px_w / (w_mm / 25.4)
            if dpi < MIN_DPI:
                fails.append("[5] PNG 预览只有 %.0f dpi (< %.0f)" % (dpi, MIN_DPI))
        except Exception as e:      # pragma: no cover
            fails.append("[5] PNG 无法读取: %s" % type(e).__name__)

    min_size = min((s["size"] for s in spans), default=float("nan"))
    return {
        "file": path.name,
        "pages": n_pages,
        "size_mm": [round(w_mm, 1), round(h_mm, 1)],
        "n_spans": len(spans),
        "n_rotated_skipped": skipped_rot,
        "min_pt": round(min_size, 2),
        "png_dpi": round(dpi) if dpi else None,
        "fails": fails,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="文件名包含该子串")
    ap.add_argument("--main-only", action="store_true")
    ap.add_argument("--save-baseline", action="store_true")
    ap.add_argument("--baseline", action="store_true")
    args = ap.parse_args()

    pdfs = sorted(FIG_DIR.glob("Fig*_*.pdf"))
    if not args.main_only:
        pdfs += sorted(FIG_DIR.glob("Supplementary_Figure_*.pdf"))
    if args.only:
        pdfs = [p for p in pdfs if args.only in p.name]
    if not pdfs:
        sys.exit("没找到图件 PDF: %s" % FIG_DIR)

    results = []
    print("=" * 78)
    print("figure_lint  |  min %.1f pt | overlap > %.0f%% | width %.0f-%.0f mm | png >= %.0f dpi"
          % (MIN_PT, OVERLAP_FRAC * 100, *WIDTH_RANGE_MM, MIN_DPI))
    print("=" * 78)
    for p in pdfs:
        is_main = p.name.startswith("Fig")
        r = lint_one(p, check_width=is_main)
        results.append(r)
        flag = "PASS" if not r["fails"] else "FAIL"
        print("%-4s %-34s %5.1f x %5.1f mm  min %4.1f pt  png %3s dpi  spans %3d"
              % (flag, r["file"][:34], r["size_mm"][0], r["size_mm"][1],
                 r["min_pt"], r["png_dpi"] if r["png_dpi"] else "-", r["n_spans"]))
        for f in r["fails"]:
            print("       " + f)

    total_fail = sum(1 for r in results if r["fails"])
    print("-" * 78)
    if total_fail == 0:
        print("RESULT: ALL PASS (%d figures)" % len(results))
    else:
        print("RESULT: %d / %d FIGURES FAIL" % (total_fail, len(results)))

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    snapshot = {r["file"]: r for r in results}

    if args.save_baseline:
        BASELINE.write_text(json.dumps(snapshot, indent=1, ensure_ascii=False),
                            encoding="utf-8")
        print("baseline saved -> %s" % BASELINE.relative_to(ROOT))

    if args.baseline and BASELINE.exists():
        old = json.loads(BASELINE.read_text(encoding="utf-8"))
        print("\n--- 与基线对比 ---")
        regress = 0
        for name, cur in snapshot.items():
            prev = old.get(name)
            if prev is None:
                print("  NEW    %s" % name)
                continue
            if cur["min_pt"] < prev["min_pt"] - 0.05:
                print("  WORSE  %s  min_pt %.1f -> %.1f"
                      % (name, prev["min_pt"], cur["min_pt"]))
                regress += 1
            if len(cur["fails"]) > len(prev["fails"]):
                print("  WORSE  %s  违规 %d -> %d"
                      % (name, len(prev["fails"]), len(cur["fails"])))
                regress += 1
        print("  新增倒退: %d 项" % regress)
        if regress:
            total_fail += regress

    return 1 if total_fail else 0


if __name__ == "__main__":
    sys.exit(main())
