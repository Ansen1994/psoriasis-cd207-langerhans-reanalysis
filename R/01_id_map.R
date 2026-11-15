## 01_id_map.R - S1 prerequisite: unify three gene ID systems to Symbol.
## ASCII ONLY. Run with working directory = project root.
## Offline: uses org.Hs.eg.db and the local GPL570 annotation in the shared
## data catalog. biomaRt is deliberately NOT used (requires network).
source("config/00_paths.R")
source("config/01_params.R")

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(data.table)
})

log_file <- file.path(LOG_DIR, "S1_id_map.log")
con <- file(log_file, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

cat("=== S1 ID mapping ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
cat("R: ", R.version.string, "\n", sep = "")
cat("org.Hs.eg.db: ", as.character(packageVersion("org.Hs.eg.db")), "\n", sep = "")

## ------------------------------------------------------------------ 1. discovery
disco_dir <- file.path(RAW_DIR, "PSOR_SC_GSE173706", "csv")
csv_files <- sort(list.files(disco_dir, pattern = "\\.csv\\.gz$", full.names = TRUE))
stopifnot(length(csv_files) == 33L)

## All 33 matrices share the same 33,538 Ensembl rows (verified in S0 pilot).
f0 <- csv_files[1]
ens_raw <- data.table::fread(f0, header = TRUE, select = 1L,
                             colClasses = "character", showProgress = FALSE)[[1L]]
cat("discovery reference file: ", basename(f0), "\n", sep = "")
cat("Ensembl rows read: ", length(ens_raw), "\n", sep = "")
cat("example IDs: ", paste(head(ens_raw, 3), collapse = ", "), "\n", sep = "")

## Strip Ensembl version suffix if present (GEO here has none, keep robust).
ens_clean <- sub("\\.[0-9]+$", "", ens_raw)
cat("IDs with version suffix stripped: ", sum(ens_clean != ens_raw), "\n", sep = "")

sym <- AnnotationDbi::mapIds(org.Hs.eg.db,
                             keys = ens_clean,
                             column = "SYMBOL",
                             keytype = "ENSEMBL",
                             multiVals = "first")
sym[is.na(sym) | sym == ""] <- NA_character_

## org.Hs.eg.db names mitochondrial genes ND1/COX2/ATP8 (NCBI style), not
## MT-ND1 (HGNC style). "^MT-" returns ZERO hits in this annotation package.
## Use the frozen whitelist from config/01_params.R.
mt_hgnc_style <- sum(grepl("^MT-", sym), na.rm = TRUE)
is_mt <- !is.na(sym) & sym %in% MT_SYMBOL_WHITELIST

ens_map <- data.table::data.table(
  ensembl_id = ens_raw,
  ensembl_clean = ens_clean,
  symbol = sym,
  is_mt = is_mt,
  row_index = seq_along(ens_raw)
)

n_mapped <- sum(!is.na(ens_map$symbol))
n_mt <- sum(ens_map$is_mt)
cat("Ensembl -> Symbol mapped: ", n_mapped, " / ", nrow(ens_map),
    " (", round(100 * n_mapped / nrow(ens_map), 1), "%)\n", sep = "")
cat("HGNC '^MT-' pattern hits (diagnostic): ", mt_hgnc_style, "\n", sep = "")
cat("mitochondrial genes via NCBI whitelist: ", n_mt, " / ",
    length(MT_SYMBOL_WHITELIST), " whitelisted\n", sep = "")
cat("MT symbols found: ",
    paste(ens_map[ens_map$is_mt == TRUE, sort(unique(symbol))], collapse = ","),
    "\n", sep = "")
if (n_mt < 10L) warning("Too few MT genes; percent.mt may be unreliable.")

## Duplicate Symbols must be handled before collapsing a count matrix:
## keep a flag and record the collision count; summation happens at matrix build.
dups <- ens_map[!is.na(symbol), .N, by = symbol][N > 1L]
cat("Symbols hit by >1 Ensembl ID: ", nrow(dups),
    " (rows involved: ", sum(dups$N), ")\n", sep = "")

## ------------------------------------------------------------------ 2. validation
val_file <- file.path(RAW_DIR, "PSOR_SC_GSE162183",
                      "GSE162183_Raw_gene_counts_matrix.tab.gz")
val_head <- readLines(gzcon(file(val_file, "rb")), n = 4L)
val_genes_example <- strsplit(val_head[2], "\t")[[1]][1]
cat("validation first gene token: ", val_genes_example,
    " | looks like Ensembl: ", grepl("^ENSG", val_genes_example), "\n", sep = "")

## ------------------------------------------------------------------ 3. bulk GPL570
## Shared catalog: <projects>/data/ref/gpl570/GPL570-55999.txt
ann_file <- file.path(dirname(PROJ_ROOT), "data", "ref", "gpl570",
                      "GPL570-55999.txt")
cat("GPL570 annotation exists: ", file.exists(ann_file), "\n", sep = "")
stopifnot(file.exists(ann_file))

## The table header is line 17; lines 1-16 are '#' column descriptions.
ann <- data.table::fread(ann_file, skip = 16L, header = TRUE, sep = "\t",
                         select = c("ID", "Gene Symbol", "ENTREZ_GENE_ID"),
                         colClasses = "character", fill = TRUE,
                         showProgress = FALSE)
data.table::setnames(ann, c("probe_id", "symbol", "entrez_id"))
ann <- ann[!is.na(probe_id) & probe_id != ""]
ann[, symbol := trimws(symbol)]
ann[symbol %in% c("", "---", "NA"), symbol := NA_character_]
## Multi-mapping probes ("A /// B") are dropped, not silently assigned.
multi <- grepl("///", ann$symbol, fixed = TRUE)
cat("GPL570 probes total: ", nrow(ann), "\n", sep = "")
cat("multi-mapping probes excluded: ", sum(multi), "\n", sep = "")
ann <- ann[!multi]
probe_map <- ann[!is.na(symbol)]
cat("GPL570 probes -> Symbol: ", nrow(probe_map), "\n", sep = "")

## ------------------------------------------------------------------ 4. outputs
out <- list(
  version = PROJ_VERSION,
  created = format(Sys.time()),
  ensembl_map = ens_map,
  probe_map = probe_map,
  mt_symbols = ens_map[ens_map$is_mt == TRUE, sort(unique(symbol))],
  mt_detection = "org.Hs.eg.db NCBI symbols via MT_SYMBOL_WHITELIST (not ^MT-)",
  duplicate_symbols = dups,
  discovery_gene_count = nrow(ens_map),
  validation_id_type = if (grepl("^ENSG", val_genes_example)) "ENSEMBL" else "SYMBOL",
  bulk_id_type = "GPL570_PROBE"
)
saveRDS(out, file.path(PROC_DIR, "S1_id_map.rds"))
data.table::fwrite(ens_map, file.path(TAB_DIR, "S1_id_map_ensembl.csv"))
data.table::fwrite(probe_map, file.path(TAB_DIR, "S1_id_map_gpl570.csv"))

cat("validation ID type recorded: ", out$validation_id_type, "\n", sep = "")
cat("OUTPUT: data/processed/S1_id_map.rds\n")
cat("OUTPUT: results/tables/S1_id_map_ensembl.csv\n")
cat("OUTPUT: results/tables/S1_id_map_gpl570.csv\n")
cat("=== S1 ID mapping DONE ===\n")
