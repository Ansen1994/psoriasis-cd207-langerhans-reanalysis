## S8_bulk_orthogonal.R - S8 / P1-1 + P1-2 supplementary analyses.
##
## NEW analysis authorised by the user on 2026-11-12 (route B: P0 + P1).
## This is an ADDITIVE layer. It does NOT:
##   - overwrite any frozen S1-S4 result table (all outputs use the S8_ prefix),
##   - change any frozen parameter in config/01_params.R,
##   - redraw any of the four main figures,
##   - flip the one-sided "greater" alternative, change resolution, or delete samples.
##
## Purpose. The strongest frozen evidence is a rank-based 10-gene ssGSEA score
## (S4). Such a score can fall because there are fewer CD207+ cells OR because
## per-cell expression falls OR because composition shifts. Two orthogonal,
## annotation-independent readouts are added here:
##   P1-1  single-gene log2 levels for CD207/CD1A/FCER1A/LAMP3/CD1C in GSE13355;
##   P1-2  the SAME 10-gene signatures plus the same single genes in two
##         additional independent GPL570 bulk cohorts (GSE14905, GSE30999).
## A single-gene result still cannot distinguish cell number from per-cell
## expression; that limitation is stated in the log and in the manuscript.
##
## ASCII only (LC_CTYPE=C). All outputs stay inside the project.
source("config/00_paths.R")
source("config/01_params.R")
source("config/02_s8_params.R")
suppressPackageStartupMessages({library(data.table); library(GSVA)})

logf <- file.path(LOG_DIR, "S8_bulk_orthogonal.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S8 supplementary bulk analyses (P1-1 + P1-2) ===\n")
cat("timestamp: ", format(Sys.time()), " | S8 version ", S8_VERSION, "\n", sep = "")
cat("FROZEN params unchanged; one-sided 'greater' NOT flipped; no new main figure.\n")
cat("Caveat: GEO series matrices contain upstream-normalised log2 RELATIVE\n")
cat("values. 2^x is a monotone rescaling only; rank tests are invariant to it.\n")
set.seed(SEED)
t0 <- Sys.time()

## ---------------------------------------------------------------- helpers ---
meta_line <- function(ln, key) {
  i <- grep(paste0("^!", key, "\t"), ln)[1]
  if (is.na(i)) return(NULL)
  x <- strsplit(ln[i], "\t", fixed = TRUE)[[1]]
  trimws(gsub('"', '', x[-1]))
}

read_geo_matrix <- function(f) {
  cat("\n[read] ", f, "\n", sep = "")
  gz <- gzfile(f, open = "rt")
  ln <- readLines(gz, warn = FALSE)
  close(gz)
  b <- grep("^!series_matrix_table_begin", ln)[1]
  e <- grep("^!series_matrix_table_end", ln)[1]
  stopifnot(!is.na(b), !is.na(e), e > b)
  ttl <- meta_line(ln, "Sample_title")
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
  rownames(mat) <- probe
  colnames(mat) <- gsms
  cat("  matrix: ", nrow(mat), " probes x ", ncol(mat),
      " samples | titles: ", length(ttl), "\n", sep = "")
  stopifnot(length(ttl) == ncol(mat))
  list(mat = mat, title = ttl, gsm = gsms)
}

annf <- normalizePath(file.path(RAW_DIR, "..", "..", "..", "data", "ref",
                               "gpl570", "GPL570-55999.txt"), mustWork = TRUE)
cat("GPL570 annotation: ", annf, "\n", sep = "")
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

grp_test <- function(d, vcol = "value") {
  a <- d[group == "PP", get(vcol)]
  b <- d[group == "NN", get(vcol)]
  if (length(a) >= 3 && length(b) >= 3 && sum(!is.na(a)) >= 3 &&
      sum(!is.na(b)) >= 3) {
    wilcox.test(a, b, exact = FALSE)$p.value
  } else NA_real_
}

pair_test <- function(d, vcol = "value") {
  w <- d[group %in% c("PN", "PP"),
         .(v = mean(get(vcol), na.rm = TRUE)), by = .(donor_id, group)]
  ww <- dcast(w, donor_id ~ group, value.var = "v")
  if (!"PN" %in% names(ww)) ww[, PN := NA_real_]
  if (!"PP" %in% names(ww)) ww[, PP := NA_real_]
  ok <- !is.na(ww$PN) & !is.na(ww$PP)
  n <- sum(ok)
  list(n_pairs = n,
       median_PN = if (n) round(median(ww$PN[ok]), 4) else NA_real_,
       median_PP = if (n) round(median(ww$PP[ok]), 4) else NA_real_,
       delta = if (n) round(median(ww$PP[ok] - ww$PN[ok]), 4) else NA_real_,
       p = if (n >= 6) suppressWarnings(wilcox.test(
             ww$PP[ok], ww$PN[ok], paired = TRUE, exact = FALSE)$p.value)
           else NA_real_)
}

single_gene_table <- function(cohort, expr, nprobe, smap) {
  genes <- intersect(S8_SINGLE_GENES, rownames(expr))
  cat("  single genes present: ", paste(genes, collapse = ", "), "\n", sep = "")
  out <- list()
  for (g in genes) {
    v <- as.numeric(expr[g, ])
    d <- data.table(sample_id = colnames(expr), value = v)
    d <- merge(d, smap[, .(sample_id, group, donor_id)],
               by = "sample_id", all.x = TRUE, sort = FALSE)
    d <- d[!is.na(group)]
    tab <- d[, .(n = .N, median = round(median(value), 4),
                 q25 = round(quantile(value, 0.25), 4),
                 q75 = round(quantile(value, 0.75), 4)), by = group]
    tab[, `:=`(cohort = cohort, gene = g)]
    out[[length(out) + 1]] <- tab
    nN <- d[group == "NN", .N]; nP <- d[group == "PP", .N]
    pp <- pair_test(d)
    tests <- data.table(
      cohort = cohort, gene = g,
      n_probes = as.integer(nprobe[g]),
      n_NN = nN, n_PP = nP,
      median_NN = if (nN) round(d[group == "NN", median(value)], 4) else NA_real_,
      median_PP = if (nP) round(d[group == "PP", median(value)], 4) else NA_real_,
      delta_PP_minus_NN = if (nN && nP) round(
        d[group == "PP", median(value)] - d[group == "NN", median(value)], 4)
        else NA_real_,
      p_wilcox_PP_vs_NN = grp_test(d),
      n_pairs = pp$n_pairs, median_PN = pp$median_PN,
      median_PP_paired = pp$median_PP,
      delta_paired_PP_minus_PN = pp$delta,
      p_signed_rank_paired = pp$p)
    out[[length(out) + 1]] <- tests
  }
  groups <- rbindlist(out[seq(1, length(out), 2)], fill = TRUE)
  tests  <- rbindlist(out[seq(2, length(out), 2)], fill = TRUE)
  tests[, FDR_BH_wilcox := p.adjust(p_wilcox_PP_vs_NN, method = S8_FDR)]
  tests[, FDR_BH_paired := p.adjust(p_signed_rank_paired, method = S8_FDR)]
  list(groups = groups, tests = tests)
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
  data.table(sample_id = gsms, sample_title = titles,
             group = grp, donor_id = donor)
}

## ------------------------------------------------------- P1-1: GSE13355 ---
cat("\n=== P1-1 single-gene orthogonal check, GSE13355 ===\n")
f13355 <- file.path(RAW_DIR, "PSOR_BULK_GSE13355", BULK_SOURCE_FILE)
r <- read_geo_matrix(f13355)
cc <- collapse_probes(r$mat)
r$mat <- NULL; gc(FALSE)
m0 <- fread(file.path(TAB_DIR, "S0_bulk_sample_map.csv"))
m0[, donor_id := sub("^Individual_([0-9]+)_.*$", "\\1", sample_title)]
setnames(m0, "gsm", "sample_id")
sm0 <- m0[, .(sample_id, group, donor_id)]
sg0 <- single_gene_table("GSE13355", cc$expr, cc$nprobe, sm0)
fwrite(sg0$groups, file.path(TAB_DIR, "S8_bulk_single_gene_groups.csv"))
fwrite(sg0$tests,  file.path(TAB_DIR, "S8_bulk_single_gene_tests.csv"))
cat("\n--- GSE13355 single-gene tests ---\n")
print(sg0$tests[, .(gene, n_probes, n_NN, n_PP, delta_PP_minus_NN,
                    p_wilcox_PP_vs_NN, FDR_BH_wilcox, n_pairs,
                    delta_paired_PP_minus_PN, p_signed_rank_paired,
                    FDR_BH_paired)])
rm(cc); gc(FALSE)

## ------------------------------------------- P1-2: two independent cohorts ---
sig <- fread(file.path(TAB_DIR, "S3_signatures.csv"))
sets <- split(sig$gene, sig$set_name)
cat("\n=== P1-2 independent bulk cohorts ===\n")
cat("signature sets: ", paste(names(sets), collapse = ", "), "\n", sep = "")

all_scores <- list(); all_sig_tests <- list(); all_gene_tests <- list()
for (co in names(S8_EXTERNAL)) {
  cfg <- S8_EXTERNAL[[co]]
  cat("\n--- ", co, ": ", cfg$note, " ---\n", sep = "")
  rr <- read_geo_matrix(file.path(RAW_DIR, cfg$dir, cfg$file))
  cc <- collapse_probes(rr$mat)
  expr2 <- 2^cc$expr                     # monotone rescaling, as in frozen S4
  sm <- smap_external(co, rr$gsm, rr$title)
  cat("  group counts: "); print(sm[, .N, by = group])
  cat("  unclassified titles: ", sum(is.na(sm$group)), "\n", sep = "")

  sets_use <- lapply(sets, function(g) intersect(unique(g), rownames(expr2)))
  cat("  signature coverage: ",
      paste(paste0(names(sets_use), "=", vapply(sets_use, length, 1L)),
            collapse = " "), "\n", sep = "")
  keep <- sets_use[vapply(sets_use, length, 1L) >= 5]
  par <- GSVA::ssgseaParam(exprData = expr2, geneSets = keep,
                           minSize = 5, maxSize = 500,
                           alpha = 0.25, normalize = TRUE)
  score <- as.matrix(GSVA::gsva(par))
  sc <- as.data.table(as.data.frame(score), keep.rownames = "gene_set")
  long <- melt(sc, id.vars = "gene_set", variable.name = "sample_id",
               value.name = "score")
  long <- merge(long, sm[, .(sample_id, group, donor_id)],
                by = "sample_id", all.x = TRUE, sort = FALSE)
  long <- long[!is.na(group)]
  long[, cohort := co]
  all_scores[[co]] <- long

  for (gs in unique(long$gene_set)) {
    d <- long[gene_set == gs, .(sample_id, group, donor_id, value = score)]
    nN <- d[group == "NN", .N]; nP <- d[group == "PP", .N]
    pp <- pair_test(d)
    all_sig_tests[[length(all_sig_tests) + 1]] <- data.table(
      cohort = co, gene_set = gs,
      n_NN = nN, n_PN = d[group == "PN", .N], n_PP = nP,
      median_NN = if (nN) round(d[group == "NN", median(value)], 4) else NA_real_,
      median_PP = if (nP) round(d[group == "PP", median(value)], 4) else NA_real_,
      delta_PP_minus_NN = if (nN && nP) round(
        d[group == "PP", median(value)] - d[group == "NN", median(value)], 4)
        else NA_real_,
      p_wilcox_PP_vs_NN = grp_test(d),
      n_pairs = pp$n_pairs, median_PN = pp$median_PN,
      median_PP_paired = pp$median_PP,
      delta_paired_PP_minus_PN = pp$delta,
      p_signed_rank_paired = pp$p,
      is_core_signature = gs == "LC_CD207_signature")
  }

  sg <- single_gene_table(co, cc$expr, cc$nprobe, sm)
  all_gene_tests[[co]] <- sg$tests
  cat("\n--- ", co, " single-gene tests ---\n", sep = "")
  print(sg$tests[, .(gene, n_probes, delta_PP_minus_NN, p_wilcox_PP_vs_NN,
                     n_pairs, delta_paired_PP_minus_PN, p_signed_rank_paired)])
  rm(rr, cc, expr2, score, sc, long); gc(FALSE)
}

scores <- rbindlist(all_scores, fill = TRUE)
fwrite(scores[order(cohort, gene_set, group, sample_id)],
       file.path(TAB_DIR, "S8_external_bulk_scores.csv"))
sig_tests <- rbindlist(all_sig_tests, fill = TRUE)
sig_tests[, `:=`(FDR_BH_wilcox = p.adjust(p_wilcox_PP_vs_NN, method = S8_FDR),
                 FDR_BH_paired = p.adjust(p_signed_rank_paired, method = S8_FDR)),
          by = cohort]
setorder(sig_tests, cohort, -is_core_signature, gene_set)
fwrite(sig_tests, file.path(TAB_DIR, "S8_external_bulk_tests.csv"))
gene_tests <- rbindlist(all_gene_tests, fill = TRUE)
gene_tests[, `:=`(FDR_BH_wilcox = p.adjust(p_wilcox_PP_vs_NN, method = S8_FDR),
                  FDR_BH_paired = p.adjust(p_signed_rank_paired, method = S8_FDR)),
           by = cohort]
fwrite(gene_tests, file.path(TAB_DIR, "S8_external_single_gene_tests.csv"))

cat("\n=== ALL external signature tests ===\n")
print(sig_tests[, .(cohort, gene_set, n_NN, n_PP, delta_PP_minus_NN,
                    p_wilcox_PP_vs_NN, FDR_BH_wilcox, n_pairs,
                    delta_paired_PP_minus_PN, p_signed_rank_paired,
                    FDR_BH_paired)])

cat("\n=== S8 bulk DONE ===\n")
cat("OUTPUT: S8_bulk_single_gene_groups.csv / S8_bulk_single_gene_tests.csv\n")
cat("OUTPUT: S8_external_bulk_scores.csv / S8_external_bulk_tests.csv\n")
cat("OUTPUT: S8_external_single_gene_tests.csv\n")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    "\n", sep = "")
