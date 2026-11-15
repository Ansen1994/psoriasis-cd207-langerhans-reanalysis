#!/usr/bin/env bash
# 00_fetch_data.sh - rebuild the exact input layout expected by the R/py scripts.
#
# WHY THIS FILE EXISTS
#   In the original working copy every dataset lived in a shared catalogue and was
#   linked into data/raw/ as a directory junction. That is invisible to a reader of
#   the supplementary package, so the pipeline could not be started from scratch.
#   This script downloads the five GEO series and the GPL570 annotation, verifies
#   SHA-256, and lays them out exactly as config/00_paths.R and the analysis scripts
#   expect.
#
# USAGE
#   cd <project root>          # the folder containing config/, R/, results/
#   bash scripts/00_fetch_data.sh
#
# REQUIREMENTS
#   bash, curl, tar, sha256sum (Git Bash on Windows provides all four).
#   Disk: ~0.5 GB of downloads + ~0.5 GB extracted; the GPL570 annotation is 76 MB.
#
# NOTE ON PATHS
#   R/S4_bulk.R and R/S8_bulk_orthogonal.R resolve the GPL570 annotation as
#   file.path(RAW_DIR, "..", "..", "..", "data", "ref", "gpl570", "GPL570-55999.txt").
#   With RAW_DIR = "data/raw" that resolves to <parent of project root>/data/ref/gpl570/.
#   This script creates that path. It looks odd but it is what the frozen scripts read,
#   and the frozen scripts are shipped unchanged so that the code matches the results.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RAW="data/raw"
REF="$(cd "$ROOT/.." && pwd)/data/ref/gpl570"
mkdir -p "$RAW"/PSOR_SC_GSE173706/csv "$RAW"/PSOR_SC_GSE162183 \
         "$RAW"/PSOR_BULK_GSE13355 "$RAW"/PSOR_BULK_GSE14905 "$RAW"/PSOR_BULK_GSE30999 \
         "$REF"

die() { echo "FATAL: $*" >&2; exit 1; }

fetch() {  # url dest expected_bytes expected_sha256
  local url="$1" dest="$2" bytes="$3" sha="$4"
  if [ -f "$dest" ] && [ "$(sha256sum "$dest" | cut -d' ' -f1)" = "$sha" ]; then
    echo "  [skip] $(basename "$dest") already present and verified"
    return 0
  fi
  echo "  [get ] $(basename "$dest")"
  curl -L --fail --retry 3 -C - -o "$dest" "$url" || die "download failed: $url"
  local got_bytes got_sha
  got_bytes=$(wc -c < "$dest" | tr -d ' ')
  got_sha=$(sha256sum "$dest" | cut -d' ' -f1)
  [ "$got_bytes" = "$bytes" ] || die "size mismatch for $dest: $got_bytes != $bytes"
  [ "$got_sha" = "$sha" ]     || die "sha256 mismatch for $dest"
}

echo "== 1/6 GSE13355 bulk series matrix (34 MB) =="
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE13nnn/GSE13355/matrix/GSE13355_series_matrix.txt.gz" \
      "$RAW/PSOR_BULK_GSE13355/GSE13355_series_matrix.txt.gz" \
      35921043 087043c4a620828d65f3055840be9de25bd31306d5caf018ae28ce5d18de6194

echo "== 2/6 GSE14905 bulk series matrix (24 MB) =="
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14905/matrix/GSE14905_series_matrix.txt.gz" \
      "$RAW/PSOR_BULK_GSE14905/GSE14905_series_matrix.txt.gz" \
      24387317 72a61fb743b5bb5cebd7951dac379a50bad0aa9e4d9c0544aea4a72b7547dcb5

echo "== 3/6 GSE30999 bulk series matrix (12 MB) =="
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE30nnn/GSE30999/matrix/GSE30999_series_matrix.txt.gz" \
      "$RAW/PSOR_BULK_GSE30999/GSE30999_series_matrix.txt.gz" \
      12251771 988d6ca437028a6e958da379908c62c8b9bb6db614108c182321488ae934574b

echo "== 4/6 GSE162183 merged count matrix (80 MB) =="
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162183/suppl/GSE162183_Raw_gene_counts_matrix.tab.gz" \
      "$RAW/PSOR_SC_GSE162183/GSE162183_Raw_gene_counts_matrix.tab.gz" \
      80418115 9755f4f3219e1fa3877e84d5cb55f6fa4c34a9bea77796c8d7b4c929891f06b6

echo "== 5/6 GSE173706 raw tar of 33 CSV libraries (264 MB) =="
if [ ! -f "$RAW/PSOR_SC_GSE173706/GSE173706_RAW.tar" ]; then
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173706/suppl/GSE173706_RAW.tar" \
        "$RAW/PSOR_SC_GSE173706/GSE173706_RAW.tar" \
        264120320 2d7d02681cadbdc2e15e38876e6c2569dd5c9779debbea397020520dece339f1
fi
if [ "$(ls -1 "$RAW/PSOR_SC_GSE173706/csv" 2>/dev/null | wc -l)" -ne 33 ]; then
  echo "  [tar ] extracting 33 CSV.gz into csv/"
  tar --force-local -xf "$RAW/PSOR_SC_GSE173706/GSE173706_RAW.tar" \
      -C "$RAW/PSOR_SC_GSE173706/csv"
fi
N_CSV=$(ls -1 "$RAW/PSOR_SC_GSE173706/csv" | wc -l | tr -d ' ')
[ "$N_CSV" = "33" ] || die "expected 33 CSV files, found $N_CSV"

echo "== 6/6 GPL570 probe annotation (76 MB) =="
if [ ! -f "$REF/GPL570-55999.txt" ]; then
  echo "  NOTE: the GPL570 annotation is distributed inside the GPL570 SOFT family,"
  echo "        not as a single table. Obtain it with:"
  echo "          curl -L -o /tmp/GPL570_family.soft.gz \\"
  echo "            https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL570nnn/GPL570/soft/GPL570_family.soft.gz"
  echo "        then export the platform table (22,880 symbols x 54,675 probes) as"
  echo "        a tab-delimited file with a 'Gene Symbol' column."
  die "place GPL570-55999.txt at $REF/GPL570-55999.txt (sha256 ecbfa3a4ac145d06b2b42508a05d0c7123fbb4ab681d2b94d6e6330a9b69aa7f)"
fi
[ "$(sha256sum "$REF/GPL570-55999.txt" | cut -d' ' -f1)" = \
  "ecbfa3a4ac145d06b2b42508a05d0c7123fbb4ab681d2b94d6e6330a9b69aa7f" ] \
  || die "GPL570-55999.txt checksum mismatch"

echo
echo "Layout ready:"
echo "  $RAW/PSOR_SC_GSE173706/csv/            33 CSV.gz"
echo "  $RAW/PSOR_SC_GSE162183/                merged count matrix"
echo "  $RAW/PSOR_BULK_GSE13355/               series matrix"
echo "  $RAW/PSOR_BULK_GSE14905/               series matrix"
echo "  $RAW/PSOR_BULK_GSE30999/               series matrix"
echo "  $REF/GPL570-55999.txt"
echo
echo "Next: Rscript R/01_id_map.R  (then 02_s1_qc.R, 03_s1_integrate.R, 04_gateB.R,"
echo "      05_s2_subcluster.R, 05b_s2_stats.R, S3_function.R, S4_transfer.R,"
echo "      S4_bulk.R, S5_gateD.R, 06_figures.R, S8_*.R, S9_prop_zero.R)"
