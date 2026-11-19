## S18_deconvolution.R - post hoc NNLS deconvolution of the three bulk cohorts
## against the FROZEN myeloid-subcluster signatures (S2_subcluster.rds).
##
## WHY THIS EXISTS (post hoc supplementary layer, mirrors the S8 precedent):
##   Reviewer-facing gap: limitation (vii)/(xii)/(xiii) state that bulk log2
##   data cannot separate "fewer CD207+ cells" from "lower per-cell CD207".
##   This layer quantifies the CD207+ (subcluster 2) share of the myeloid
##   transcriptional program in each bulk sample by NNLS against the four
##   frozen subcluster mean-expression signatures, then re-runs the SAME
##   group/paired tests used in S4/S8 on the estimated shares.
##
## SCOPE AND GUARDS (additive layer, nothing frozen is touched):
##   - Reads the frozen S2_subcluster.rds and the frozen bulk series matrices;
##     writes only NEW S18_ tables. No frozen table, figure or parameter is
##     modified; no Gate criterion uses these outputs.
##   - The share is RELATIVE to discrimination among the four myeloid
##     signatures, NOT an absolute tissue fraction: whole-skin bulk contains
##     keratinocytes/fibroblasts etc. that are not represented here, so all
##     weights absorb background roughly proportionally and the sub2 share
##     remains an internal-ratio readout. Reported with this caveat.
##   - Subcluster 1 (sweat-gland/secretory epithelium) is INCLUDED as a
##     competing signature so that eccrine signal is absorbed there rather
##     than inflating sub0/sub2 shares.
##   - GSE13355/GSE14905/GSE30999 values are log2 RELATIVE (negative values
##     present); 2^x is a monotone rescaling only, as in frozen S4/S8.
##   - Signatures are linear-scale per-cell means (expm1 of log-normalized
##     sc values); NNLS shares are invariant to a global unit rescaling of
##     either the mixture or all signatures, so the sc/array unit mismatch
##     does not bias the ratios.
##
## ASCII only (LC_CTYPE=C). All outputs stay inside the project.
source("config/00_paths.R")
source("config/01_params.R")
source("config/02_s8_params.R")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(nnls)
})

logf <- file.path(LOG_DIR, "S18_deconvolution.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S18 post hoc NNLS deconvolution vs frozen myeloid subclusters ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
set.seed(SEED)
t0 <- Sys.time()

SUBS <- c("0", "1", "2", "3")

## ---------------------------------------------------------------- 1. signatures
cat("\n[1] building linear mean signatures from frozen S2_subcluster.rds\n")
s2 <- readRDS("data/processed/S2_subcluster.rds")
lab <- as.character(s2$subcluster)
stopifnot(setequal(unique(lab), SUBS))
dta <- GetAssayData(s2, assay = "RNA", slot = "data")   # log-normalized
## expm1 in place on the sparse x slot keeps sparsity (expm1(0) = 0), so the
## linear matrix is never materialized densely (commit-memory safe).
stopifnot(inherits(dta, "dgCMatrix"))
dta@x <- expm1(dta@x)                                   # linear scale
lin <- dta
gc(FALSE)
SIG <- matrix(NA_real_, nrow(lin), length(SUBS),
              dimnames = list(rownames(lin), paste0("sub", SUBS)))
for (s in SUBS) {
  idx <- which(lab == s)
  SIG[, paste0("sub", s)] <- Matrix::rowMeans(lin[, idx, drop = FALSE])
  cat("  sub", s, ": ", length(idx), " cells\n", sep = "")
}
rm(lin); gc(FALSE)
drop <- Matrix::rowSums(SIG) <= 0
SIG <- SIG[!drop, ]
cat("  signature matrix: ", nrow(SIG), " genes x ", ncol(SIG), "\n", sep = "")

## ---------------------------------------------------------------- 2. helpers
meta_line <- function(ln, key) {
  i <- grep(paste0("^!", key), ln)[1]
  if (is.na(i)) return(NULL)
  x <- strsplit(ln[i], "\t", fixed = TRUE)[[1]]
  trimws(gsub('"', '', x[-1]))
}

## Streaming reader: scans the gz in chunks and writes only the expression
## table block to a temp file, so the ~250 MB decompressed text is never held
## in R memory (commit-memory safe on this machine).
read_geo_matrix <- function(f) {
  cat("\n[read] ", f, "\n", sep = "")
  con <- gzfile(f, open = "rt")
  on.exit(close(con), add = TRUE)
  tf <- tempfile(fileext = ".txt")
  con2 <- file(tf, open = "wt")
  on.exit(close(con2), add = TRUE)
  ttl <- NULL; gsm <- NULL
  in_tbl <- FALSE
  repeat {
    ch <- readLines(con, n = 50000, warn = FALSE)
    if (!length(ch)) break
    if (!in_tbl) {
      i <- grep("^!Sample_title", ch)[1]
      if (!is.na(i)) ttl <- trimws(gsub('"', '', strsplit(ch[i], "\t", fixed = TRUE)[[1]][-1]))
      i <- grep("^!Sample_geo_accession", ch)[1]
      if (!is.na(i)) gsm <- trimws(gsub('"', '', strsplit(ch[i], "\t", fixed = TRUE)[[1]][-1]))
      b <- grep("^!series_matrix_table_begin", ch)[1]
      if (!is.na(b)) {
        in_tbl <- TRUE
        rest <- ch[-seq_len(b)]
        e <- grep("^!series_matrix_table_end", rest)[1]
        if (!is.na(e)) {
          if (e > 1) cat(rest[seq_len(e - 1)], file = con2, sep = "\n")
          break
        }
        if (length(rest)) cat(rest, file = con2, sep = "\n")
      }
    } else {
      e <- grep("^!series_matrix_table_end", ch)[1]
      if (!is.na(e)) {
        if (e > 1) cat(ch[seq_len(e - 1)], file = con2, sep = "\n")
        break
      }
      cat(ch, file = con2, sep = "\n")
    }
  }
  stopifnot(in_tbl, !is.null(ttl))
  flush(con2)
  cat("  titles: ", length(ttl), " | table lines written to temp\n", sep = "")
  ## read.delim parses incrementally (no mmap) — commit-memory safe here.
  ex <- read.delim(tf, header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "\"", stringsAsFactors = FALSE)
  unlink(tf)
  probe <- gsub('"', '', ex[[1]])
  gsms <- gsub('"', '', colnames(ex)[-1])
  mat <- as.matrix(ex[, -1, drop = FALSE])
  rm(ex); gc(FALSE)
  storage.mode(mat) <- "double"
  rownames(mat) <- probe
  colnames(mat) <- gsms
  cat("  matrix: ", nrow(mat), " probes x ", ncol(mat), " samples\n", sep = "")
  stopifnot(length(ttl) == ncol(mat), identical(gsm, gsms))
  list(mat = mat, title = ttl, gsm = gsms)
}

annf <- normalizePath(file.path(RAW_DIR, "..", "..", "..", "data", "ref",
                               "gpl570", "GPL570-55999.txt"), mustWork = TRUE)
ann <- fread(annf, skip = 16, sep = "\t", header = TRUE, check.names = FALSE,
             quote = "", fill = TRUE, select = c("ID", "Gene Symbol"),
             showProgress = FALSE)

collapse_probes <- function(mat) {
  probe <- rownames(mat)
  sym <- ann[["Gene Symbol"]][match(probe, ann$ID)]
  sym <- sub(" ///.*$", "", sym)
  sym[is.na(sym) | sym %in% c("", "---", "NA")] <- NA_character_
  idx <- !is.na(sym)
  cnt <- table(sym[idx])
  agg <- rowsum(mat[idx, , drop = FALSE], group = sym[idx], reorder = TRUE)
  agg <- agg / as.numeric(cnt[rownames(agg)])
  cat("  probes with Symbol: ", sum(idx), "/", length(sym),
      " -> ", nrow(agg), " gene-level rows\n", sep = "")
  list(expr = agg, nprobe = cnt)
}

smap_13355 <- function(gsms, titles) {
  m0 <- fread(file.path(TAB_DIR, "S0_bulk_sample_map.csv"))
  setnames(m0, "gsm", "sample_id")
  m0[, donor_id := sub("^Individual_([0-9]+)_.*$", "\\1", sample_title)]
  m0[, sample_title := NULL][, characteristics := NULL]
  m0
}

smap_external <- function(cohort, gsms, titles) {
  if (cohort == "GSE14905") {
    grp <- ifelse(grepl("^Normal Skin", titles), "NN",
           ifelse(grepl("^Uninvolved Skin", titles), "PN",
           ifelse(grepl("^Lesion Skin", titles), "PP", NA_character_)))
    donor <- ifelse(grp == "NN",
                    paste0("NN_", sub("^Normal Skin, Normal-", "", titles)),
             ifelse(grp == "PN",
                    paste0("P_", sub("^Uninvolved Skin, NS-", "", titles)),
                    paste0("P_", sub("^Lesion Skin, LS-", "", titles))))
  } else {
    grp <- ifelse(grepl("_NL$", titles), "PN",
           ifelse(grepl("_LS$", titles), "PP", NA_character_))
    donor <- sub("_[^_]+$", "", titles)
  }
  data.table(sample_id = gsms, group = grp, donor_id = donor)
}

deconvolve <- function(expr_lin, sig) {
  g <- intersect(rownames(sig), rownames(expr_lin))
  cat("  NNLS genes: ", length(g), "\n", sep = "")
  Sm <- sig[g, , drop = FALSE]
  Bm <- expr_lin[g, , drop = FALSE]
  W <- matrix(NA_real_, ncol(Sm), ncol(Bm),
              dimnames = list(colnames(Sm), colnames(Bm)))
  ok <- rep(TRUE, ncol(Bm))
  for (j in seq_len(ncol(Bm))) {
    b <- Bm[, j]
    if (any(!is.finite(b)) || max(b, na.rm = TRUE) <= 0) { ok[j] <- FALSE; next }
    fit <- tryCatch(nnls::nnls(Sm, b), error = function(e) NULL)
    if (is.null(fit)) { ok[j] <- FALSE; next }
    W[, j] <- as.numeric(fit$x)
  }
  Wtot <- colSums(W[, ok, drop = FALSE]); Wtot[Wtot <= 0] <- NA_real_
  SH <- sweep(W[, ok, drop = FALSE], 2, Wtot, "/")
  list(W = W, shares = SH, ok = ok)
}

grp_test <- function(v, grp, a = "PP", b = "NN") {
  x <- v[grp == a]; y <- v[grp == b]
  if (length(x) >= 3 && length(y) >= 3 &&
      sum(is.finite(x)) >= 3 && sum(is.finite(y)) >= 3) {
    wilcox.test(x, y, exact = FALSE)$p.value
  } else NA_real_
}

pair_test <- function(v, grp, donor, a = "PP", b = "PN") {
  d <- data.table(v = v, grp = grp, donor = donor)
  d <- d[grp %in% c(a, b) & is.finite(v)]
  w <- dcast(d, donor ~ grp, value.var = "v")
  if (!a %in% names(w)) w[, (a) := NA_real_]
  if (!b %in% names(w)) w[, (b) := NA_real_]
  ok <- is.finite(w[[a]]) & is.finite(w[[b]])
  n <- sum(ok)
  list(n_pairs = n,
       median_a = if (n) round(median(w[[a]][ok]), 4) else NA_real_,
       median_b = if (n) round(median(w[[b]][ok]), 4) else NA_real_,
       delta = if (n) round(median(w[[a]][ok] - w[[b]][ok]), 4) else NA_real_,
       p = if (n >= 6) suppressWarnings(wilcox.test(
             w[[a]][ok], w[[b]][ok], paired = TRUE, exact = FALSE)$p.value)
           else NA_real_)
}

## ---------------------------------------------------------------- 3. run
COHORTS <- list(
  GSE13355 = list(dir = "PSOR_BULK_GSE13355", file = BULK_SOURCE_FILE,
                  smap = "frozen"),
  GSE14905 = S8_EXTERNAL$GSE14905,
  GSE30999 = S8_EXTERNAL$GSE30999
)
READOUTS <- c("share_sub0", "share_sub1", "share_sub2", "share_sub3")
DIAG_GENES <- c("CD207", "CD1A", "DCD", "LAMP3")

all_frac <- list(); all_tests <- list(); all_diag <- list()
for (co in names(COHORTS)) {
  cfg <- COHORTS[[co]]
  cat("\n=== ", co, " ===\n", sep = "")
  rr <- read_geo_matrix(file.path(RAW_DIR, cfg$dir, cfg$file))
  cc <- collapse_probes(rr$mat)
  expr_lin <- 2^cc$expr                  # monotone rescaling, as frozen S4/S8
  sm <- if (identical(cfg$smap, "frozen")) smap_13355(rr$gsm, rr$title)
        else smap_external(co, rr$gsm, rr$title)
  cat("  group counts: "); print(sm[, .N, by = group])

  dv <- deconvolve(expr_lin, SIG)
  fw <- as.data.table(t(dv$W), keep.rownames = "sample_id")
  setnames(fw, c("sample_id", paste0("w_sub", SUBS)))
  fs <- as.data.table(t(dv$shares), keep.rownames = "sample_id")
  setnames(fs, c("sample_id", paste0("share_sub", SUBS)))
  frac <- merge(fw, fs, by = "sample_id", all.x = TRUE)
  
  
  frac[, w_sub2_of_myeloid := fifelse(
        is.finite(w_sub0 + w_sub2 + w_sub3) & (w_sub0 + w_sub2 + w_sub3) > 0,
        w_sub2 / (w_sub0 + w_sub2 + w_sub3), NA_real_)]
  frac <- merge(frac, sm, by = "sample_id", all.x = TRUE, sort = FALSE)
  frac[, cohort := co]

  ## diagnostics: CD207/CD1A vs share_sub2; DCD vs share_sub1
  guse <- intersect(DIAG_GENES, rownames(cc$expr))
  for (g in guse) {
    v <- as.numeric(cc$expr[g, colnames(cc$expr)])
    names(v) <- colnames(cc$expr)
    i <- match(frac$sample_id, names(v))
    frac[, paste0("expr_", g) := v[i]]
    for (sh in c("share_sub2", "share_sub1")) {
      okd <- is.finite(frac[[sh]]) & is.finite(frac[[paste0("expr_", g)]])
      rho <- if (sum(okd) >= 10) suppressWarnings(
        cor.test(frac[[sh]][okd], frac[[paste0("expr_", g)]][okd],
                 method = "spearman", exact = FALSE)) else NULL
      all_diag[[length(all_diag) + 1]] <- data.table(
        cohort = co, share = sh, gene = g, n = sum(okd),
        spearman_rho = if (!is.null(rho)) round(unname(rho$estimate), 4) else NA_real_,
        p = if (!is.null(rho)) rho$p.value else NA_real_)
    }
  }

  for (ro in c(READOUTS, "w_sub2_of_myeloid")) {
    v <- frac[[ro]]
    pp <- pair_test(v, frac$group, frac$donor_id)
    all_tests[[length(all_tests) + 1]] <- data.table(
      cohort = co, readout = ro,
      n_NN = sum(frac$group == "NN"), n_PN = sum(frac$group == "PN"),
      n_PP = sum(frac$group == "PP"),
      median_NN = if (sum(frac$group == "NN")) round(median(v[frac$group == "NN"], na.rm = TRUE), 4) else NA_real_,
      median_PN = if (sum(frac$group == "PN")) round(median(v[frac$group == "PN"], na.rm = TRUE), 4) else NA_real_,
      median_PP = if (sum(frac$group == "PP")) round(median(v[frac$group == "PP"], na.rm = TRUE), 4) else NA_real_,
      p_wilcox_PP_vs_NN = grp_test(v, frac$group),
      n_pairs = pp$n_pairs, median_PN_paired = pp$median_b,
      median_PP_paired = pp$median_a,
      delta_paired_PP_minus_PN = pp$delta,
      p_signed_rank_paired = pp$p)
  }
  all_frac[[co]] <- frac
  rm(rr, cc, expr_lin, dv); gc(FALSE)
}

fr <- rbindlist(all_frac, fill = TRUE)
fwrite(fr[order(cohort, group, sample_id)],
       file.path(TAB_DIR, "S18_deconvolution_fractions.csv"))
ts <- rbindlist(all_tests, fill = TRUE)
ts[, `:=`(FDR_BH_wilcox = p.adjust(p_wilcox_PP_vs_NN, method = S8_FDR),
          FDR_BH_paired = p.adjust(p_signed_rank_paired, method = S8_FDR)),
   by = cohort]
setorder(ts, cohort, readout)
fwrite(ts, file.path(TAB_DIR, "S18_deconvolution_tests.csv"))
fwrite(rbindlist(all_diag, fill = TRUE),
       file.path(TAB_DIR, "S18_deconvolution_diagnostics.csv"))

cat("\n=== S18 tests ===\n")
print(ts[, .(cohort, readout, median_NN, median_PN, median_PP,
             p_wilcox_PP_vs_NN, FDR_BH_wilcox, n_pairs,
             delta_paired_PP_minus_PN, p_signed_rank_paired, FDR_BH_paired)])
cat("\n=== S18 diagnostics (Spearman) ===\n")
print(rbindlist(all_diag, fill = TRUE))
cat("\nelapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1),
    " min\n", sep = "")
cat("=== S18 done ===\n")
