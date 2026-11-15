#!/usr/bin/env python3
"""Build the LOCAL supplementary-materials ZIP for dermatology-sci (S8 / P1-4).

Creates manuscript/submission/Supplementary_materials_20261112.zip and a
manifest outside the journal folders. Nothing is uploaded anywhere: a public
repository DOI requires the user's account and is left as an explicit TBD.

Usage (from anywhere):  python py/S8_build_supplementary.py
"""

import csv
import glob
import hashlib
import io
import os
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUB = os.path.join(ROOT, "manuscript", "submission")
ZIP = os.path.join(SUB, "Supplementary_materials_20261112.zip")
MANIFEST = os.path.join(SUB, "Supplementary_manifest.csv")

# (glob relative to project root, destination folder inside the ZIP)
# REVISION POLICY (after internal review):
#   - no internal governance documents (Chinese planning files, gate logs that
#     name fallback journals, submission-strategy notes) are shipped;
#   - no manuscript-generator script is shipped (it embeds an earlier title);
#   - the analysis log shipped here is the sanitised English summary below.
GROUPS = [
    ("results/tables/S0_*.csv", "tables"),
    ("results/tables/S1_*.csv", "tables"),
    ("results/tables/S2_*.csv", "tables"),
    ("results/tables/S3_*.csv", "tables"),
    ("results/tables/S4_*.csv", "tables"),
    ("results/tables/S8_*.csv", "tables"),
    ("results/tables/S9_*.csv", "tables"),
    ("results/tables/S10_*.csv", "tables"),
    ("results/tables/S11_*.csv", "tables"),
    ("results/tables/S13_*.csv", "tables"),
    ("results/tables/S16_*.csv", "tables"),
    ("config/00_paths.R", "config"),
    ("config/01_params.R", "config"),
    ("config/02_s8_params.R", "config"),
    ("R/*.R", "R"),
    ("py/S4_direct_markers.py", "py"),
    ("py/S8_doi_check.py", "py"),
    ("scripts/00_fetch_data.sh", "scripts"),
    ("logs/ANALYSIS_SUMMARY.md", "logs"),
    ("logs/S9_sessionInfo_*.txt", "logs"),
    ("logs/S9_prop_zero.log", "logs"),
    ("logs/S10_signature_confounding.log", "logs"),
    ("logs/S11_lc_maturation.log", "logs"),
    ("logs/S12_revision_figures.log", "logs"),
    ("logs/S13_tf_activity.log", "logs"),
    ("logs/S16_cd207_tf_correlation.log", "logs"),
    ("results/figures/Supplementary_Figure_S1.pdf", "figures"),
    ("results/figures/Supplementary_Figure_S2.pdf", "figures"),
    ("results/figures/Supplementary_Figure_S3.pdf", "figures"),
    ("results/figures/Supplementary_Figure_S4.pdf", "figures"),
    ("results/figures/Supplementary_Figure_S5.pdf", "figures"),
    ("results/figures/Supplementary_Figure_S6.pdf", "figures"),
]

# Files that must never reach a journal / public repository.
EXCLUDE_BASENAMES = {
    "S7_author_info.log", "author_info_v1.md",
    "checklist_Genes.md", "checklist_CCID.md",
    # Local build utilities: they list forbidden tokens or embed an earlier
    # title, and they are not analysis artefacts.
    "S8_build_supplementary.py", "S8_build_master_v1_1.py",
    "apply_author_info.py",
}

# The FOUR authors below were supplied on 2026-11-12 and then explicitly
# retracted by the user at 21:34 the same day. Their contact details are
# obsolete personal data. A hard failure here prevents accidentally shipping
# them inside a supplementary ZIP. Bare digit strings are NOT used as tokens
# because expression tables legitimately contain numbers such as 190448.
PRIVACY_TOKENS = (
    "Ansheng Xu", "Cheng Chen", "Guobing Yin", "Jing Luo",
    "@hospital.cqmu.edu.cn", "0009-0007-5465-8122",
    "Chongqing Medical University",
    # submission-strategy leakage
    "CCID", "Arch Dermatol Res", "BMC Med Genomics", "BMC Medical Genomics",
)

README = """# Supplementary materials - dermatology-sci

Re-analysis of public psoriasis skin transcriptomes. This archive contains the
derived result tables, the analysis and figure scripts, the frozen parameter
files, the environment snapshot and Supplementary Figures S1-S6.

## Contents

| Folder | What it holds |
|---|---|
| `tables/` | All result tables: `S0_*` to `S4_*` (frozen analysis), `S8_*` (post hoc supplementary layer), `S9_*` (structural-zero correction), `S10_*` (composition controls and empirical null), `S11_*` (maturation axis), `S13_*` (transcription-factor activity) and `S16_*` (TF-CD207 correlation) |
| `config/` | `00_paths.R`, the frozen `01_params.R`, and the additive `02_s8_params.R` |
| `R/` | All analysis and figure scripts, including `S9_prop_zero.R`, `S10_signature_confounding.R`, `S11_lc_maturation.R`, `S12_revision_figures.R`, `S13_tf_activity.R`, `S13b_tf_figures.R`, `S16_cd207_tf_correlation.R` and `S24_jt_crosscheck.R` (the last is a read-only verification of the self-implemented Jonckheere–Terpstra test against an independent implementation) |
| `py/` | Data-pilot, mapping-free marker counting and DOI verification scripts |
| `scripts/` | `00_fetch_data.sh`: downloads the five GEO series and lays them out exactly as the scripts expect, with SHA-256 checks |
| `logs/` | `ANALYSIS_SUMMARY.md` (stage order, gates, corrections), the per-stage environment snapshot (`S9_sessionInfo_*.txt`) and the logs of the structural-zero, composition, maturation and TF-activity analyses |
| `figures/` | Supplementary Figures S1-S6 (vector PDF) |
| `manifest.csv` | Per-file size and SHA-256 |

## How to reproduce

1. `bash scripts/00_fetch_data.sh` - downloads GSE173706, GSE162183, GSE13355,
   GSE14905 and GSE30999 and verifies their SHA-256 checksums. This also creates
   the GPL570 annotation path that `R/S4_bulk.R` and `R/S8_bulk_orthogonal.R`
   read. The annotation itself must be obtained from GEO (see the script).
2. Run the R scripts in the order listed in `logs/ANALYSIS_SUMMARY.md`. Each one
   sources `config/00_paths.R` first and uses seed 20261112.
3. `R/06_figures.R` renders Figures 1-4 from the frozen objects and result
   tables. It requires the two `.rds` objects produced by stage S1 and S2
   (about 1.4 GB and 32 MB); these are not included in this archive because of
   their size, and are regenerated by steps S1-S2 of the pipeline.

## Scope notes

- S8 is post hoc supplementary work added after the frozen analysis and does not
  enter any pre-registered gate.
- S9 corrects a data-handling error in the proportion table: pseudo-samples with
  zero cells in a subcluster had been omitted instead of contributing 0%. The
  affected statistics are reported under four specifications in the manuscript.
- Single-gene log2 levels cannot separate fewer CD207+ cells from lower per-cell
  expression. All three bulk cohorts are upstream-normalised log2 relative
  values, and GSE30999 has no healthy controls. LC-subset scores are descriptive
  only, with no p values.
"""


def main():
    os.makedirs(SUB, exist_ok=True)
    entries = []
    seen = set()
    for pattern, folder in GROUPS:
        for src in sorted(glob.glob(os.path.join(ROOT, pattern))):
            if not os.path.isfile(src):
                continue
            name = os.path.basename(src)
            if name in EXCLUDE_BASENAMES:
                continue
            arc = folder + "/" + name
            if arc in seen:
                continue
            seen.add(arc)
            entries.append((src, arc))
    if not entries:
        raise SystemExit("No supplementary files found")

    # Privacy gate: run over EVERY payload, including the README and manifest,
    # before a single byte is written to the ZIP.
    payloads = [("README_supplementary.md", README.encode("utf-8"))]
    for src, arc in entries:
        payloads.append((arc, open(src, "rb").read()))
    for arc, data in payloads:
        text = data.decode("utf-8", errors="ignore")
        leaked = [tok for tok in PRIVACY_TOKENS if tok in text]
        if leaked:
            raise SystemExit("PRIVACY GATE FAILED in %s: %s" % (arc, ", ".join(leaked)))

    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("README_supplementary.md", README)
        rows = []
        for src, arc in entries:
            data = open(src, "rb").read()
            zf.writestr(arc, data)
            rows.append(dict(zip_path=arc,
                             source_path=os.path.relpath(src, ROOT).replace("\\", "/"),
                             bytes=len(data),
                             sha256=hashlib.sha256(data).hexdigest()))
        rows.sort(key=lambda r: r["zip_path"])
        bio = io.StringIO(newline="")
        w = csv.DictWriter(bio, fieldnames=["zip_path", "source_path", "bytes", "sha256"])
        w.writeheader()
        w.writerows(rows)
        zf.writestr("manifest.csv", bio.getvalue())

    with open(MANIFEST, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["zip_path", "source_path", "bytes", "sha256"])
        w.writeheader()
        w.writerows(rows)
    print("WROTE", ZIP, os.path.getsize(ZIP), "bytes,", len(entries), "files")
    print("WROTE", MANIFEST)
    for folder in sorted(set(a.split("/")[0] for _, a in entries)):
        print("  %-10s %d files" % (folder, sum(1 for _, a in entries
                                                if a.startswith(folder + "/"))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
