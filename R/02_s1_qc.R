## 02_s1_qc.R - single pass over the 33 discovery matrices.
## Produces per-cell QC metrics, a Symbol-collapsed sparse count matrix,
## and the observed distributions used by Gate A QC recalibration.
## ASCII ONLY. Working directory = project root.
## No figures are produced by this script.
source("config/00_paths.R")
source("config/01_params.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

log_file <- file.path(LOG_DIR, "S1_qc.log")
con <- file(log_file, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

cat("=== S1 QC scan ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
set.seed(SEED)

idm <- readRDS(file.path(PROC_DIR, "S1_id_map.rds"))
emap <- idm$ensembl_map
mt_syms <- idm$mt_symbols

## Keep only rows with a Symbol; collapse duplicate Symbols by summation.
keep <- which(!is.na(emap$symbol))
sym_keep <- emap$symbol[keep]
u_sym <- unique(sym_keep)
grp <- match(sym_keep, u_sym)
G <- sparseMatrix(i = grp, j = seq_along(grp), x = 1,
                  dims = c(length(u_sym), length(grp)))
cat("rows kept (mapped to Symbol): ", length(keep), "\n", sep = "")
cat("unique Symbols after collapsing duplicates: ", length(u_sym), "\n", sep = "")
cat("MT genes used for percent.mt: ", length(mt_syms), "\n", sep = "")

## Sample metadata from S0 (read-only junctions were used there).
smap <- data.table::fread(file.path(TAB_DIR, "S0_sample_map.csv"))
data.table::setorder(smap, file)
cat("samples in S0 map: ", nrow(smap), "\n", sep = "")

mats <- vector("list", nrow(smap))
qc_list <- vector("list", nrow(smap))
t_start <- Sys.time()

for (i in seq_len(nrow(smap))) {
  f <- file.path(RAW_DIR, "PSOR_SC_GSE173706", "csv", smap$file[i])
  sid <- sub("\\.csv\\.gz$", "", smap$file[i])
  t0 <- Sys.time()

  d <- data.table::fread(f, header = TRUE, showProgress = FALSE)
  raw_genes <- d[[1L]]
  m_dense <- as.matrix(d[, -1L, with = FALSE])
  bcs <- colnames(d)[-1L]
  rm(d)

  ## Guard: every matrix must use the identical gene ordering as the map.
  if (length(raw_genes) != nrow(emap) ||
      !identical(as.character(raw_genes), as.character(emap$ensembl_id))) {
    stop("Gene row mismatch in ", smap$file[i])
  }

  M <- Matrix(m_dense[keep, , drop = FALSE], sparse = TRUE)
  rm(m_dense)
  M <- drop0(G %*% M)
  rownames(M) <- u_sym
  colnames(M) <- paste0(sid, "_", bcs)

  nFeature <- diff(M@p)
  nCount <- colSums(M)
  mt_idx <- which(u_sym %in% mt_syms)
  mt_sum <- colSums(M[mt_idx, , drop = FALSE])
  pct_mt <- ifelse(nCount > 0, 100 * mt_sum / nCount, NA_real_)

  qc_list[[i]] <- data.table::data.table(
    cell_id = colnames(M),
    sample_id = sid,
    group = smap$group[i],
    donor = smap$donor[i],
    donor_base = smap$donor_base[i],
    nFeature = as.integer(nFeature),
    nCount = as.numeric(nCount),
    pct_mt = as.numeric(pct_mt)
  )
  mats[[i]] <- M

  cat(sprintf("[%02d/%02d] %-42s cells=%5d nnz=%9d med_nFeature=%5.0f med_pctMT=%5.2f  %.1fs\n",
              i, nrow(smap), sid, ncol(M), length(M@x),
              stats::median(nFeature), stats::median(pct_mt, na.rm = TRUE),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  flush.console()
  gc(FALSE)
}

cat("--- merging ", length(mats), " sparse matrices ---\n", sep = "")
t1 <- Sys.time()
counts <- do.call(cbind, mats)
rm(mats); gc(FALSE)
cat("merged dim: ", nrow(counts), " genes x ", ncol(counts), " cells\n", sep = "")
cat("merged nnz: ", length(counts@x), "\n", sep = "")
cat("merged object MB: ", round(as.numeric(object.size(counts)) / 1048576, 1), "\n", sep = "")
cat("merge_sec: ", round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1), "\n", sep = "")

qc <- data.table::rbindlist(qc_list)
rm(qc_list)

## ------------------------------------------------------------ observed distributions
med_f <- stats::median(qc$nFeature)
mad_f <- stats::mad(qc$nFeature)
med_c <- stats::median(qc$nCount)
mad_c <- stats::mad(qc$nCount)
med_m <- stats::median(qc$pct_mt, na.rm = TRUE)
mad_m <- stats::mad(qc$pct_mt, na.rm = TRUE)

cat("\n=== Gate A QC recalibration basis (median +/- K*MAD, K=", QC_MAD_K, ") ===\n", sep = "")
cat("nFeature: median=", round(med_f, 1), " MAD=", round(mad_f, 1),
    "  -> [", round(med_f - QC_MAD_K * mad_f, 1), ", ",
    round(med_f + QC_MAD_K * mad_f, 1), "]\n", sep = "")
cat("nCount  : median=", round(med_c, 1), " MAD=", round(mad_c, 1),
    "  -> [", round(med_c - QC_MAD_K * mad_c, 1), ", ",
    round(med_c + QC_MAD_K * mad_c, 1), "]\n", sep = "")
cat("pct_mt  : median=", round(med_m, 2), " MAD=", round(mad_m, 2),
    "  -> upper ", round(med_m + QC_MAD_K * mad_m, 2), "\n", sep = "")
cat("quantiles nFeature 1/5/25/50/75/95/99%: ",
    paste(round(as.numeric(stats::quantile(qc$nFeature, c(.01,.05,.25,.5,.75,.95,.99))), 0), collapse = "/"), "\n", sep = "")
cat("quantiles pct_mt  50/75/90/95/99%: ",
    paste(round(as.numeric(stats::quantile(qc$pct_mt, c(.5,.75,.9,.95,.99), na.rm = TRUE)), 2), collapse = "/"), "\n", sep = "")

## ------------------------------------------------------------ flagged small sample
flag_sid <- sub("\\.csv\\.gz$", "", FLAG_SMALL_SAMPLE)
sq <- qc[sample_id == flag_sid]
cat("\n=== Flagged sample ", flag_sid, " ===\n", sep = "")
cat("cells=", nrow(sq), " med_nFeature=", round(stats::median(sq$nFeature), 1),
    " med_nCount=", round(stats::median(sq$nCount), 1),
    " med_pctMT=", round(stats::median(sq$pct_mt, na.rm = TRUE), 2), "\n", sep = "")
other <- qc[sample_id != flag_sid]
cat("others: cells=", nrow(other), " med_nFeature=", round(stats::median(other$nFeature), 1),
    " med_nCount=", round(stats::median(other$nCount), 1), "\n", sep = "")

## ------------------------------------------------------------ outputs
data.table::fwrite(qc, file.path(TAB_DIR, "S1_qc_per_cell.csv.gz"))
summary <- qc[, .(
  n_cells = .N,
  med_nFeature = round(stats::median(nFeature), 1),
  med_nCount = round(stats::median(nCount), 1),
  med_pct_mt = round(stats::median(pct_mt, na.rm = TRUE), 2),
  pct_cells_mt_gt_20 = round(100 * mean(pct_mt > 20, na.rm = TRUE), 2),
  pct_cells_nFeature_lt_200 = round(100 * mean(nFeature < 200), 2),
  pct_cells_nFeature_gt_6000 = round(100 * mean(nFeature > 6000), 2)
), by = .(sample_id, group, donor_base)]
data.table::setorder(summary, sample_id)
data.table::fwrite(summary, file.path(TAB_DIR, "S1_qc_summary.csv"))

meta <- qc[, .(cell_id, sample_id, group, donor, donor_base)]
saveRDS(list(counts = counts, meta = meta, symbols = u_sym,
             version = PROJ_VERSION, created = format(Sys.time())),
        file.path(PROC_DIR, "S1_counts.rds"))

cat("\nOUTPUT: results/tables/S1_qc_per_cell.csv.gz\n")
cat("OUTPUT: results/tables/S1_qc_summary.csv\n")
cat("OUTPUT: data/processed/S1_counts.rds\n")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "\n", sep = "")
cat("=== S1 QC scan DONE ===\n")
