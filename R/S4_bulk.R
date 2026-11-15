## S4_bulk.R - bulk cohort (GSE13355) core-subcluster signature scoring.
## O6 downgrade chain: BayesPrism NOT installed, CIBERSORTx needs a web token,
## so the frozen fallback DECONV_FALLBACK = "ssGSEA_module_score" is used.
##
## Gate D direction criterion (pre-specified here, before seeing results):
##   the LC_CD207_signature score in PP must be LOWER than in NN and the
##   Wilcoxon BH-FDR must be < 0.05. PN is descriptive only.
##   Mac_DC_signature is a specificity control (no direction requirement).
##
## Data caveat (must be in Methods): the GSE13355 series matrix contains
## batch/sex-adjusted log2 RELATIVE expression and includes NEGATIVE values
## (e.g. -0.253). BULK_LOG2 = TRUE is honoured by applying 2^x, but this is a
## monotone rescaling, NOT recovery of absolute linear expression. ssGSEA is
## rank based, so scores are invariant to that transformation. The upstream
## adjustment is irreversible and is reported as a limitation.
##
## ASCII only (LC_CTYPE=C). All outputs stay inside the project.
source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({library(data.table); library(GSVA)})

logf <- file.path(LOG_DIR, "S4_bulk.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S4 bulk signature scoring (offline fallback) ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
cat("fallback: ", DECONV_FALLBACK, " | BULK_LOG2 = ", BULK_LOG2, "\n", sep = "")
cat("Gate D direction: LC_CD207_signature PP < NN with BH-FDR < 0.05.\n")
set.seed(SEED)
t0 <- Sys.time()

## ---------------------------------------------------------------- 1. signatures
sig <- fread(file.path(TAB_DIR, "S3_signatures.csv"))
sets <- split(sig$gene, sig$set_name)
cat("signature sets: ", paste(names(sets), collapse = ", "), "\n", sep = "")

## ---------------------------------------------------------------- 2. expression
f <- file.path(RAW_DIR, "PSOR_BULK_GSE13355", BULK_SOURCE_FILE)
ln <- readLines(gzfile(f, open = "rt"), warn = FALSE)
b <- grep("^!series_matrix_table_begin", ln)[1]
e <- grep("^!series_matrix_table_end", ln)[1]
stopifnot(!is.na(b), !is.na(e), e > b)
cat("series matrix table lines: ", b + 1, " to ", e - 1,
    " (", e - b - 1, " rows)\n", sep = "")
txt <- paste(ln[(b + 1):(e - 1)], collapse = "\n")
rm(ln); gc(FALSE)
ex <- fread(text = txt, sep = "\t", header = TRUE, check.names = FALSE,
            showProgress = FALSE)
rm(txt); gc(FALSE)
probe <- gsub('"', '', ex[[1]])
gsms <- gsub('"', '', colnames(ex)[-1])
mat <- as.matrix(ex[, -1, with = FALSE])
rm(ex); gc(FALSE)
storage.mode(mat) <- "double"
colnames(mat) <- gsms
cat("expression matrix: ", nrow(mat), " probes x ", ncol(mat), " samples\n", sep = "")
cat("NA values: ", sum(is.na(mat)), " | range: ",
    round(min(mat, na.rm = TRUE), 3), " to ", round(max(mat, na.rm = TRUE), 3),
    " | negative values: ", sum(mat < 0, na.rm = TRUE), "\n", sep = "")

## ---------------------------------------------------------------- 3. probe -> Symbol
annf <- file.path(RAW_DIR, "..", "..", "..", "data", "ref", "gpl570", "GPL570-55999.txt")
annf <- normalizePath(annf, mustWork = TRUE)
cat("GPL570 annotation: ", annf, "\n", sep = "")
ann <- fread(annf, skip = 16, sep = "\t", header = TRUE, check.names = FALSE,
             quote = "", fill = TRUE, select = c("ID", "Gene Symbol"),
             showProgress = FALSE)
sym <- ann[["Gene Symbol"]][match(probe, ann$ID)]
sym <- sub(" ///.*$", "", sym)
sym[is.na(sym) | sym %in% c("", "---", "NA")] <- NA_character_
idx <- !is.na(sym)
cat("probes with a Symbol: ", sum(idx), " / ", length(sym), "\n", sep = "")
cnt <- table(sym[idx])
agg <- rowsum(mat[idx, , drop = FALSE], group = sym[idx], reorder = TRUE)
agg <- agg / as.numeric(cnt[rownames(agg)])
cat("gene-level matrix: ", nrow(agg), " Symbols x ", ncol(agg),
    " samples (multi-probe Symbols averaged on log2 scale)\n", sep = "")

## Monotone 2^x rescaling (see header caveat).
if (isTRUE(BULK_LOG2)) expr <- 2^agg else expr <- agg
rm(mat, agg); gc(FALSE)

## ---------------------------------------------------------------- 4. coverage
sets_use <- lapply(sets, function(g) intersect(unique(g), rownames(expr)))
cov <- data.table(set_name = names(sets_use),
                  n_input = vapply(sets, length, 1L),
                  n_mapped = vapply(sets_use, length, 1L))
cov[, fraction_mapped := round(n_mapped / n_input, 3)]
cov[, used_in_ssGSEA := n_mapped >= 5]
fwrite(cov, file.path(TAB_DIR, "S4_bulk_signature_coverage.csv"))
cat("\n--- signature coverage on GPL570 ---\n"); print(cov)
for (n in names(sets_use)) cat(" ", n, ": ",
  paste(intersect(sets_use[[n]], rownames(expr)), collapse = ", "), "\n", sep = "")
sets_keep <- sets_use[cov$used_in_ssGSEA]

## ---------------------------------------------------------------- 5. ssGSEA
score <- NULL
if (length(sets_keep) > 0) {
  par <- GSVA::ssgseaParam(exprData = expr, geneSets = sets_keep,
                           minSize = 5, maxSize = 500,
                           alpha = 0.25, normalize = TRUE)
  score <- tryCatch(as.matrix(GSVA::gsva(par)),
                    error = function(err) {
                      cat("ssGSEA FAILED: ", conditionMessage(err), "\n", sep = ""); NULL })
}
if (is.null(score)) {
  cat("WARNING: no ssGSEA scores produced; Gate D bulk criterion = FAIL.\n")
  fwrite(data.table(gene_set = character(), sample_id = character(),
                    group = character(), score = numeric()),
         file.path(TAB_DIR, T_S4_SCORES))
  fwrite(data.table(gene_set = character(), n_NN = integer(), n_PP = integer(),
                    median_NN = numeric(), median_PP = numeric(),
                    delta_PP_minus_NN = numeric(), p_wilcox = numeric(),
                    FDR_BH = numeric(), direction_expected = logical()),
         file.path(TAB_DIR, "S4_bulk_group_tests.csv"))
} else {
  cat("\nssGSEA score matrix: ", nrow(score), " sets x ", ncol(score), " samples\n", sep = "")
  smap <- fread(file.path(TAB_DIR, "S0_bulk_sample_map.csv"))
  sc <- as.data.table(as.data.frame(score), keep.rownames = "gene_set")
  long <- melt(sc, id.vars = "gene_set", variable.name = "gsm", value.name = "score")
  long <- merge(long, smap[, .(gsm, group, sample_title)], by = "gsm", all.x = TRUE)
  long <- long[!is.na(group)]
  long[, donor_id := sub("^Individual_([0-9]+)_.*$", "\\1", sample_title)]
  setnames(long, "gsm", "sample_id")
  fwrite(long[order(gene_set, group, sample_id)],
         file.path(TAB_DIR, T_S4_SCORES))
  cat("group counts: "); print(long[, .(n_samples = uniqueN(sample_id)), by = group])

  tests <- long[group %in% c("NN", "PP"), .(
    n_NN = sum(group == "NN"), n_PP = sum(group == "PP"),
    median_NN = round(median(score[group == "NN"]), 4),
    median_PP = round(median(score[group == "PP"]), 4),
    delta_PP_minus_NN = round(median(score[group == "PP"]) -
                              median(score[group == "NN"]), 4),
    p_wilcox = wilcox.test(x = score[group == "PP"],
                           y = score[group == "NN"], exact = FALSE)$p.value),
    by = gene_set]
  tests[, FDR_BH := p.adjust(p_wilcox, method = "BH")]
  pn <- long[group == "PN", .(median_PN = round(median(score), 4)), by = gene_set]
  tests <- merge(tests, pn, by = "gene_set", all.x = TRUE)
  tests[, direction_expected := FALSE]
  tests[gene_set == "LC_CD207_signature" & delta_PP_minus_NN < 0 & FDR_BH < 0.05,
        direction_expected := TRUE]
  tests[, is_core_signature := gene_set == "LC_CD207_signature"]
  setorder(tests, -is_core_signature, p_wilcox)
  fwrite(tests, file.path(TAB_DIR, "S4_bulk_group_tests.csv"))
  cat("\n--- bulk group tests (PP vs NN) ---\n"); print(tests)
  BULK_DIR_OK <- any(tests$direction_expected, na.rm = TRUE)
  cat("\nBULK_DIRECTION_CONSISTENT = ", BULK_DIR_OK, "\n", sep = "")

  ## Supportive paired analysis required by frozen route S4: 58 PN-PP pairs
  ## come from the same donors. This is descriptive; it does NOT enter Gate D.
  pd <- long[group %in% c("PN", "PP") & !is.na(donor_id),
             .(score = mean(score)), by = .(gene_set, donor_id, group)]
  w <- dcast(pd, gene_set + donor_id ~ group, value.var = "score")
  if (!"PN" %in% colnames(w)) w[, PN := NA_real_]
  if (!"PP" %in% colnames(w)) w[, PP := NA_real_]
  pair_tab <- w[, {
    ok <- !is.na(PN) & !is.na(PP)
    .(n_pairs = sum(ok),
      median_PN = round(median(PN[ok]), 4),
      median_PP = round(median(PP[ok]), 4),
      delta_PP_minus_PN = round(median(PP[ok]) - median(PN[ok]), 4),
      p_signed_rank = if (sum(ok) >= 6)
        suppressWarnings(wilcox.test(x = PP[ok], y = PN[ok],
                                     paired = TRUE, exact = FALSE)$p.value)
        else NA_real_)
  }, by = gene_set]
  pair_tab[, FDR_BH := p.adjust(p_signed_rank, method = "BH")]
  fwrite(pair_tab, file.path(TAB_DIR, "S4_bulk_paired_PN_PP.csv"))
  cat("\n--- supportive paired PN-PP signed-rank (descriptive, not Gate D) ---\n")
  print(pair_tab)
}

cat("\n=== S4 bulk DONE ===\n")
cat("OUTPUT: results/tables/S4_module_scores.csv\n")
cat("OUTPUT: results/tables/S4_bulk_group_tests.csv\n")
cat("OUTPUT: results/tables/S4_bulk_signature_coverage.csv\n")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")
