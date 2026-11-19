## S19_treatment_bulk.R - post hoc analysis of the secukinumab treatment
## cohort GSE226244 (GPL570, same platform as GSE13355; Kim & Krueger,
## Front Immunol 2023, trial NCT03131570).
##
## QUESTION (post hoc supplementary layer, S8/S18 precedent):
##   does the CD207+ Langerhans/DC signature recover after effective
##   IL-17A blockade? Readouts, all mirroring the frozen pipeline:
##   (a) single-gene CD207/CD1A/FCER1A/LAMP3/CD1C (frozen S8_SINGLE_GENES);
##   (b) the frozen LC_CD207_signature / Mac_DC_signature / Cycling_signature
##       ssGSEA scores (same settings as S8: alpha = 0.25, normalize = TRUE);
##   (c) S18-style NNLS deconvolution against the four frozen myeloid
##       subcluster signatures.
##
## GROUPS (from series metadata; NO patient identifier exists, so the PRE vs
##   POST comparison is UNPAIRED - stated in the manuscript):
##   NN      : disease state = Control (8 arrays, incl. re-analyzed GSE78097)
##   PRE     : disease Psoriasis + treatment Pretreatment (33 arrays, Baseline)
##   POST    : disease Psoriasis + treatment Posttreatment (28 arrays; 18 on12,
##             9 on24, 1 off20 - the timeline is kept per sample in the output)
##
## CAVEATS carried into the manuscript: values are log2 RELATIVE; 2^x is a
##   monotone rescaling only; PRE/POST arms are unpaired and the POST arm
##   mixes on-treatment time points; deconvolution shares are relative to
##   the four myeloid signatures (see S18 header for the collinearity note).
##
## ASCII only (LC_CTYPE=C). All outputs stay inside the project.
source("config/00_paths.R")
source("config/01_params.R")
source("config/02_s8_params.R")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(nnls); library(GSVA)
})

logf <- file.path(LOG_DIR, "S19_treatment_bulk.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S19 GSE226244 treatment bulk (post hoc) ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
set.seed(SEED)
t0 <- Sys.time()

## ---------------------------------------------------------------- 1. helpers
read_geo_matrix <- function(f) {
  cat("\n[read] ", f, "\n", sep = "")
  con <- gzfile(f, open = "rt")
  on.exit(close(con), add = TRUE)
  tf <- tempfile(fileext = ".txt")
  con2 <- file(tf, open = "wt")
  on.exit(close(con2), add = TRUE)
  ttl <- NULL; gsm <- NULL; in_tbl <- FALSE; chlines <- character(0)
  repeat {
    ch <- readLines(con, n = 50000, warn = FALSE)
    if (!length(ch)) break
    if (!in_tbl) {
      i <- grep("^!Sample_title", ch)[1]
      if (!is.na(i)) ttl <- trimws(gsub('"', '', strsplit(ch[i], "\t", fixed = TRUE)[[1]][-1]))
      i <- grep("^!Sample_geo_accession", ch)[1]
      if (!is.na(i)) gsm <- trimws(gsub('"', '', strsplit(ch[i], "\t", fixed = TRUE)[[1]][-1]))
      ci <- grep("^!Sample_characteristics_ch1", ch)
      if (length(ci)) chlines <- c(chlines, ch[ci])
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
  cat("  titles: ", length(ttl), "\n", sep = "")
  ## characteristics: split each !Sample_characteristics_ch1 line by \t
  chars <- lapply(chlines, function(x)
    trimws(gsub('"', '', strsplit(x, "\t", fixed = TRUE)[[1]][-1])))
  names(chars) <- vapply(chlines, function(x)
    trimws(gsub('"', '', strsplit(x, "\t", fixed = TRUE)[[1]][1])), "")
  stopifnot(all(lengths(chars) == length(ttl)))
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
  list(mat = mat, title = ttl, gsm = gsms, chars = chars)
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

## ---------------------------------------------------------------- 2. signatures
cat("\n[2] building linear mean signatures from frozen S2_subcluster.rds\n")
s2 <- readRDS("data/processed/S2_subcluster.rds")
lab <- as.character(s2$subcluster)
SUBS <- c("0", "1", "2", "3")
stopifnot(setequal(unique(lab), SUBS))
dta <- GetAssayData(s2, assay = "RNA", slot = "data")
stopifnot(inherits(dta, "dgCMatrix"))
dta@x <- expm1(dta@x)                       # linear scale, sparse-safe
SIG <- matrix(0, nrow(dta), length(SUBS),
              dimnames = list(rownames(dta), paste0("sub", SUBS)))
for (s in SUBS) SIG[, paste0("sub", s)] <-
  Matrix::rowMeans(dta[, which(lab == s), drop = FALSE])
rm(dta); gc(FALSE)
SIG <- SIG[Matrix::rowSums(SIG) > 0, ]
cat("  signature matrix: ", nrow(SIG), " genes x ", ncol(SIG), "\n", sep = "")

sigsets <- fread(file.path(TAB_DIR, "S3_signatures.csv"))
sets <- split(sigsets$gene, sigsets$set_name)
cat("  frozen signature sets: ", paste(names(sets), collapse = ", "), "\n", sep = "")

deconvolve <- function(expr_lin, sig) {
  g <- intersect(rownames(sig), rownames(expr_lin))
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

## ---------------------------------------------------------------- 3. data
rr <- read_geo_matrix(file.path(RAW_DIR, "PSOR_BULK_GSE226244",
                                "GSE226244_series_matrix.txt.gz"))
cc <- collapse_probes(rr$mat)
expr_lin <- 2^cc$expr                        # monotone rescaling, as frozen S4/S8
cat("  value range after 2^x: ", round(min(expr_lin, na.rm = TRUE), 4), " to ",
    round(max(expr_lin, na.rm = TRUE), 1), "\n", sep = "")

## metadata parsed from the series matrix characteristics lines (no SOFT file
## is loaded). Values are matched by content prefix, not by column position.
keyval <- function(chars, gsm, pfx) {
  hit <- which(vapply(chars, function(v) grepl(paste0("^", pfx), v[1]), TRUE))
  if (!length(hit)) stop("characteristics prefix not found: ", pfx)
  v <- unlist(chars[[hit[1]]], use.names = FALSE)
  names(v) <- gsm
  v
}
disease_v  <- keyval(rr$chars, rr$gsm, "disease state")
lesional_v <- keyval(rr$chars, rr$gsm, "lesional")
treat_v    <- keyval(rr$chars, rr$gsm, "treatment:")
timeline_v <- keyval(rr$chars, rr$gsm, "treatment_timeline")
smap <- data.table(sample_id = rr$gsm, title = rr$title,
                   disease = disease_v, lesional = lesional_v,
                   treatment = treat_v, timeline = timeline_v)
smap[, group := fifelse(grepl("Control", disease), "NN",
        fifelse(grepl("Pretreatment", treatment) & grepl("Baseline", timeline),
                "PRE", "POST"))]
cat("  group counts: "); print(smap[, .N, by = group])
cat("  POST timeline split: "); print(smap[group == "POST", .N, by = timeline])

## ---------------------------------------------------------------- 4. readouts
res_samples <- data.table(sample_id = rr$gsm, title = rr$title,
                          group = smap$group, timeline = smap$timeline)

## (a) single genes
guse <- intersect(S8_SINGLE_GENES, rownames(cc$expr))
cat("  single genes present: ", paste(guse, collapse = ", "), "\n", sep = "")
for (g in guse) {
  v <- as.numeric(cc$expr[g, colnames(cc$expr)]); names(v) <- colnames(cc$expr)
  res_samples[, paste0("expr_", g) := v[sample_id]]
}

## (b) ssGSEA frozen signature scores
sets_use <- lapply(sets, function(g) intersect(unique(g), rownames(expr_lin)))
keep <- sets_use[vapply(sets_use, length, 1L) >= 5]
cat("  ssGSEA sets: ", paste(names(keep), collapse = ", "), "\n", sep = "")
par <- GSVA::ssgseaParam(exprData = expr_lin, geneSets = keep,
                         minSize = 5, maxSize = 500, alpha = 0.25, normalize = TRUE)
sc <- as.data.table(as.data.frame(GSVA::gsva(par)), keep.rownames = "gene_set")
long <- melt(sc, id.vars = "gene_set", variable.name = "sample_id",
             value.name = "score")
long <- dcast(long, sample_id ~ gene_set, value.var = "score")
res_samples <- merge(res_samples, long, by = "sample_id", all.x = TRUE, sort = FALSE)

## (c) NNLS deconvolution shares
dv <- deconvolve(expr_lin, SIG)
fw <- as.data.table(t(dv$W), keep.rownames = "sample_id")
setnames(fw, c("sample_id", paste0("w_sub", SUBS)))
fs <- as.data.table(t(dv$shares), keep.rownames = "sample_id")
setnames(fs, c("sample_id", paste0("share_sub", SUBS)))
res_samples <- merge(res_samples, merge(fw, fs, by = "sample_id", all.x = TRUE),
                     by = "sample_id", all.x = TRUE, sort = FALSE)
res_samples[, w_sub2_of_myeloid := fifelse(
  is.finite(w_sub0 + w_sub2 + w_sub3) & (w_sub0 + w_sub2 + w_sub3) > 0,
  w_sub2 / (w_sub0 + w_sub2 + w_sub3), NA_real_)]

setorder(res_samples, group, title)
fwrite(res_samples, file.path(TAB_DIR, "S19_treatment_bulk_samples.csv"))

## ---------------------------------------------------------------- 5. tests
READOUTS <- c(paste0("expr_", guse), names(keep),
              paste0("share_sub", SUBS), "w_sub2_of_myeloid")
ts <- rbindlist(lapply(READOUTS, function(ro) {
  v <- res_samples[[ro]]
  post <- v[res_samples$group == "POST"]; pre <- v[res_samples$group == "PRE"]
  nn   <- v[res_samples$group == "NN"]
  data.table(readout = ro, n_NN = sum(!is.na(nn)), n_PRE = sum(!is.na(pre)),
             n_POST = sum(!is.na(post)),
             median_NN = if (sum(!is.na(nn))) round(median(nn, na.rm = TRUE), 4) else NA_real_,
             median_PRE = round(median(pre, na.rm = TRUE), 4),
             median_POST = round(median(post, na.rm = TRUE), 4),
             p_POST_vs_PRE = if (sum(is.finite(post)) >= 3 && sum(is.finite(pre)) >= 3)
               wilcox.test(post, pre, exact = FALSE)$p.value else NA_real_,
             p_PRE_vs_NN = if (sum(is.finite(nn)) >= 3 && sum(is.finite(pre)) >= 3)
               wilcox.test(pre, nn, exact = FALSE)$p.value else NA_real_)
}))
ts[, `:=`(FDR_POST_vs_PRE = p.adjust(p_POST_vs_PRE, method = S8_FDR),
          FDR_PRE_vs_NN = p.adjust(p_PRE_vs_NN, method = S8_FDR))]
setorder(ts, readout)
fwrite(ts, file.path(TAB_DIR, "S19_treatment_bulk_tests.csv"))
cat("\n=== S19 tests ===\n")
print(ts)

## ---------------------------------------------------------------- 6. diagnostics
diag <- rbindlist(lapply(c("share_sub2", "share_sub1", "w_sub2_of_myeloid"),
  function(sh) lapply(intersect(c("CD207", "CD1A", "DCD", "LAMP3"),
                                paste0("expr_", guse)), function(gl) {
    okd <- is.finite(res_samples[[sh]]) & is.finite(res_samples[[gl]])
    if (sum(okd) < 10) return(NULL)
    rho <- suppressWarnings(cor.test(res_samples[[sh]][okd],
                 res_samples[[gl]][okd], method = "spearman", exact = FALSE))
    data.table(share = sh, gene = gsub("^expr_", "", gl), n = sum(okd),
               spearman_rho = round(unname(rho$estimate), 4), p = rho$p.value)
  })))
fwrite(diag, file.path(TAB_DIR, "S19_treatment_bulk_diagnostics.csv"))
cat("\n=== S19 diagnostics (Spearman) ===\n")
print(diag)
cat("\nelapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min\n", sep = "")
cat("=== S19 done ===\n")
