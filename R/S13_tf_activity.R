## S13_tf_activity.R
## Transcription-factor activity inference added 2026-11-13 in response to the
## nine-dimension internal review (dimension 5: mechanistic depth).
##
## Motivation
##   The manuscript reports a CD207+ / CD1A- / LAMP3+ shift in the myeloid
##   compartment but stops at marker genes. A reviewer reasonably asks whether
##   the shift is accompanied by a coherent regulatory programme. This script
##   asks that question with an off-the-shelf, reproducible method.
##
## Method
##   decoupleR::run_ulm (univariate linear model) over the dorothea human
##   regulon, at two levels:
##     (a) SAMPLE level  - one TF-activity vector per donor x state pseudobulk,
##         then donor-level tests (KW across three states; paired signed rank
##         PN vs PP). Donor is the unit of inference, as everywhere else.
##     (b) CONTRAST level - one TF-activity vector per DEG contrast, using the
##         DESeq2 Wald statistic. Used only as an internal consistency check,
##         because S3 supports a single two-group contrast (subcluster 0).
##
## Because single-sample pseudobulk depth varies several fold across groups, a
## raw activity difference can be a depth artefact. Two depth diagnostics are
## therefore reported for every TF, and a paired depth-adjusted model is fitted:
##   d_score = a + b * d_log10(n_cells);  H0: a = 0.
##
## HARD LIMITS, all disclosed in the manuscript
##   1. TF activity is inferred from mRNA of regulon targets. It is NOT measured
##      protein activity, NOT DNA binding, NOT a perturbation.
##   2. dorothea A+B covers 118 human TFs. A TF outside the regulon is not
##      "unchanged"; it is not assessable.
##   3. Scores re-weight the same pseudobulk matrix that produced the DEG table,
##      so they are not independent evidence for the DEGs. They summarise them.
##   4. Pseudo-samples hold few cells (median 17 / 34 / 73 for NS / PN / PP), so
##      single-sample activity is noisy and depth-confounded; see point 5.
##   5. Everything here is EXPLORATORY and post hoc. It does not enter Gate C/D.
##
## Inputs (read-only)
##   data/processed/S3_pseudobulk.rds
##   results/tables/S3_two_group_DEG.csv
## Outputs
##   results/tables/S13_tf_network_summary.csv
##   results/tables/S13_tf_activity_samples.csv
##   results/tables/S13_tf_activity_tests.csv   (incl. depth diagnostics)
##   results/tables/S13_tf_contrast.csv
##   results/tables/S13_tf_depth_confounding.csv
##   results/tables/S13_tf_consistency.csv
##   logs/S13_tf_activity.log

source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({
  library(data.table); library(decoupleR); library(dorothea); library(Matrix)
})

logf <- file.path(LOG_DIR, "S13_tf_activity.log")
cat("=== S13 TF activity inference (exploratory, post hoc) ===\n", file = logf)
say <- function(...) {
  m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
says <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))

say("run at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("decoupleR ", as.character(packageVersion("decoupleR")),
    " | dorothea ", as.character(packageVersion("dorothea")))

## ------------------------------------------------------------------ 1. regulon
data(dorothea_hs, package = "dorothea")
doro <- as.data.table(dorothea_hs)
NETS <- list(
  AB  = doro[confidence %in% c("A", "B"), .(tf, target, mor)],
  ABC = doro[confidence %in% c("A", "B", "C"), .(tf, target, mor)]
)
for (nm in names(NETS)) {
  say(sprintf("regulon %-3s : %d edges, %d TFs", nm, nrow(NETS[[nm]]),
              uniqueN(NETS[[nm]]$tf)))
}

## ------------------------------------------------- 2. donor x state pseudobulk
## S3_pseudobulk.rds holds counts already summed per donor x state x subcluster.
## Subcluster 1 is the sweat-gland contamination cluster removed in S2 and is
## absent from this object.
x <- readRDS(file.path(PROC_DIR, "S3_pseudobulk.rds"))
pb <- x$pb
info <- as.data.table(x$info)
## S3 stored info sorted by (subcluster, combo) while matrix columns follow
## sort(unique(combo)); align info to the matrix before anything else.
stopifnot(setequal(colnames(pb), info$combo), !anyDuplicated(info$combo))
info <- info[match(colnames(pb), combo)]
info[, c("pid", "sub") := tstrsplit(combo, "|", fixed = TRUE)]
stopifnot(!any(info$sub == "1"))
say("\npseudobulk: ", nrow(pb), " genes x ", ncol(pb),
    " donor x state x subcluster units")

cpm_log2 <- function(m) log2(t(t(m) / colSums(m)) * 1e6 + 1)

## (a) all QC-eligible myeloid cells pooled over subclusters
idx <- split(seq_len(ncol(pb)), info$pid)
MYE <- vapply(idx, function(i) rowSums(pb[, i, drop = FALSE]), numeric(nrow(pb)))
colnames(MYE) <- names(idx)
CPM <- cpm_log2(MYE)
meta_pid <- info[, .(donor = donor_base[1], group = group[1],
                     n_cells = sum(n_cells)), by = .(pid)]
setorder(meta_pid, group, pid)

## (b) dominant subcluster 0 only, least affected by composition change
SUB0 <- as.matrix(pb[, info$sub == "0", drop = FALSE])
colnames(SUB0) <- info$pid[info$sub == "0"]
CPM0 <- cpm_log2(SUB0)
meta0 <- info[sub == "0", .(donor = donor_base[1], group = group[1],
                            n_cells = n_cells), by = .(pid)]
setorder(meta0, group, pid)

say("\n--- pseudo-samples per layer and group ---")
says(rbind(meta_pid[, .(layer = "myeloid_all", n_pseudo = .N, cells_min = as.numeric(min(n_cells)),
                        cells_median = as.numeric(median(n_cells)), cells_max = as.numeric(max(n_cells))), by = group],
           meta0[, .(layer = "subcluster_0", n_pseudo = .N, cells_min = as.numeric(min(n_cells)),
                     cells_median = as.numeric(median(n_cells)), cells_max = as.numeric(max(n_cells))), by = group])[
         order(layer, group)])

## ------------------------------------------------------- 3. network x matrix
filter_net <- function(net, genes, minsize = 5) {
  n <- net[target %in% genes]
  n[, n_tgt := .N, by = tf]
  n <- n[n_tgt >= minsize][, n_tgt := NULL]
  n
}

net_summary <- list()
run_layer <- function(mat, label) {
  genes <- rownames(mat)[rowSums(mat > 0) >= 5]
  out <- list()
  for (nm in names(NETS)) {
    net <- filter_net(NETS[[nm]], genes)
    net_summary[[length(net_summary) + 1]] <<- data.table(
      network = nm, layer = label,
      n_tf_total = uniqueN(NETS[[nm]]$tf), n_edge_total = nrow(NETS[[nm]]),
      n_genes_in_matrix = length(genes),
      n_tf_used = uniqueN(net$tf), n_edge_used = nrow(net),
      median_targets = median(net[, .N, by = tf]$N))
    act <- as.data.table(run_ulm(mat = mat[genes, , drop = FALSE], net = net,
                                 .source = "tf", .target = "target", .mor = "mor",
                                 minsize = 5))
    act <- act[statistic == "ulm"]
    act[, network := nm]
    out[[nm]] <- act[, .(network, tf = source, condition, score, p_value)]
  }
  rbindlist(out)
}

say("\n--- TF activity, sample level (regulon coverage after minsize = 5) ---")
act_mye <- run_layer(CPM, "myeloid_all")
act_sub0 <- run_layer(CPM0, "subcluster_0")
ns <- rbindlist(net_summary)
says(ns)

## ---------------------------------------- 4. donor-level tests + depth control
test_layer <- function(act, meta, lay) {
  w <- dcast(act, tf + network ~ condition, value.var = "score")
  long <- melt(w, id.vars = c("tf", "network"), variable.name = "pid",
               value.name = "score")
  long <- merge(long, meta[, .(pid, donor, group, n_cells)], by = "pid")
  res <- long[, {
    d <- .SD
    kwp <- tryCatch(
      kruskal.test(score ~ factor(group, levels = c("NS", "PN", "PP")))$p.value,
      error = function(e) NA_real_)
    sc <- dcast(d, donor ~ group, value.var = "score")
    nc <- dcast(d, donor ~ group, value.var = "n_cells")
    sc <- sc[!is.na(PN) & !is.na(PP)]
    nc <- nc[!is.na(PN) & !is.na(PP)]
    setorder(sc, donor); setorder(nc, donor)
    npr <- nrow(sc)
    p_pair <- if (npr >= 5) suppressWarnings(
      wilcox.test(sc$PP, sc$PN, paired = TRUE)$p.value) else NA_real_
    ## paired depth-adjusted model on the difference scale
    ds <- sc$PP - sc$PN
    dl <- log10(nc$PP) - log10(nc$PN)
    fit <- if (npr >= 6) lm(ds ~ dl) else NULL
    int_p <- if (!is.null(fit)) summary(fit)$coefficients[1, 4] else NA_real_
    slope_p <- if (!is.null(fit)) summary(fit)$coefficients[2, 4] else NA_real_
    ct <- if (npr >= 6) suppressWarnings(
      cor.test(ds, dl, method = "spearman", exact = FALSE)) else NULL
    data.table(
      n_pseudo = .N,
      median_NS = if (any(d$group == "NS")) median(d[group == "NS"]$score) else NA_real_,
      median_PN = if (any(d$group == "PN")) median(d[group == "PN"]$score) else NA_real_,
      median_PP = if (any(d$group == "PP")) median(d[group == "PP"]$score) else NA_real_,
      n_pairs_PN_PP = npr,
      paired_delta = if (npr >= 5) mean(ds) else NA_real_,
      paired_median_delta = if (npr >= 5) median(ds) else NA_real_,
      paired_p = p_pair,
      depth_adjusted_intercept_p = int_p,
      depth_slope_p = slope_p,
      rho_delta_depth = if (!is.null(ct)) unname(ct$estimate) else NA_real_,
      p_rho_delta_depth = if (!is.null(ct)) ct$p.value else NA_real_,
      kw_p = kwp)
  }, by = .(network, tf)]
  res[, layer := lay]
  res[, kw_FDR := p.adjust(kw_p, "BH"), by = network]
  res[, paired_FDR := p.adjust(paired_p, "BH"), by = network]
  res[, depth_adj_FDR := p.adjust(depth_adjusted_intercept_p, "BH"), by = network]
  res[]
}

t_mye <- test_layer(act_mye, meta_pid, "myeloid_all")
t_sub0 <- test_layer(act_sub0, meta0, "subcluster_0")
tests <- rbind(t_mye, t_sub0)
setorder(tests, layer, network, paired_p)

## "robust" requires all three: the paired group test survives BH, the paired
## depth-adjusted intercept survives BH, and the activity difference is not
## itself strongly correlated with the depth difference. Without the third arm a
## depth-driven TF can pass the first two (e.g. ARNTL below, paired p = 0.11 but
## depth-adjusted p = 0.004 and rho = 0.75).
tests[, robust := !is.na(paired_FDR) & paired_FDR < 0.05 &
                !is.na(depth_adj_FDR) & depth_adj_FDR < 0.05 &
                !is.na(rho_delta_depth) & abs(rho_delta_depth) < 0.6]
say("\n--- 'robust' TF count (paired FDR<0.05 AND depth-adj FDR<0.05 AND |rho|<0.6) ---")
says(tests[, .(n_robust = sum(robust), n_paired_FDR05 = sum(paired_FDR < 0.05, na.rm = TRUE),
               n_tf = .N), by = .(layer, network)])
say("\n--- TFs that pass paired + depth-adjusted BH but are excluded by depth ",
    "correlation (|rho| >= 0.6): depth artefacts, must not be reported as findings ---")
says(tests[paired_FDR < 0.05 & depth_adj_FDR < 0.05 & abs(rho_delta_depth) >= 0.6,
           .(layer, network, tf, paired_delta = round(paired_delta, 2),
             paired_p = signif(paired_p, 2), rho_delta_depth = round(rho_delta_depth, 2))])

## per-TF depth confound (pooled correlation, complementary to the paired model)
depth <- rbindlist(lapply(list(list("myeloid_all", act_mye, meta_pid),
                               list("subcluster_0", act_sub0, meta0)),
  function(z) {
    layer <- z[[1]]; act <- z[[2]]; meta <- z[[3]]
    m <- merge(act, meta[, .(pid, n_cells)], by.x = "condition", by.y = "pid")
    m[, log_n := log10(n_cells)]
    m[, {
      ct <- suppressWarnings(cor.test(score, log_n, method = "spearman", exact = FALSE))
      .(layer = layer, n = .N, rho_depth = unname(ct$estimate), p_depth = ct$p.value)
    }, by = .(network, tf)]
  }))
say("\n--- depth confound, pooled: |rho| between TF score and log10(cells) ---")
says(depth[, .(median_abs_rho = round(median(abs(rho_depth)), 3),
               max_abs_rho = round(max(abs(rho_depth)), 3),
               n_tf_abs_rho_gt_0.6 = sum(abs(rho_depth) > 0.6), n_tf = .N),
           by = .(layer, network)])
say("\n--- depth confound, paired: |rho| between d_score and d_log10(cells) ---")
says(tests[, .(median_abs_rho = round(median(abs(rho_delta_depth), na.rm = TRUE), 3),
               max_abs_rho = round(max(abs(rho_delta_depth), na.rm = TRUE), 3),
               n_tf_abs_rho_gt_0.6 = sum(abs(rho_delta_depth) > 0.6, na.rm = TRUE),
               n_tf = .N), by = .(layer, network)])

## ------------------------------------------------------- 5. contrast level
deg <- fread(file.path(TAB_DIR, "S3_two_group_DEG.csv"))
say("\ncontrast-level input: ", nrow(deg), " genes; subclusters tested: ",
    paste(sort(unique(deg$subcluster)), collapse = ", "),
    "; contrasts: ", paste(unique(deg$contrast), collapse = ", "))
say("NOTE: S3 only supports a two-group contrast in subcluster 0 (subclusters 2",
    " and 3 were skipped for insufficient pseudo-samples), so the contrast-level",
    " readout is a single vector and is used only for internal consistency.")
sm <- dcast(deg, gene ~ subcluster, value.var = "stat")
smat <- as.matrix(sm[, -1]); rownames(smat) <- sm$gene
smat <- smat[complete.cases(smat), , drop = FALSE]
colnames(smat) <- paste0("subcluster_", colnames(smat))

act_con <- list()
for (nm in names(NETS)) {
  net <- filter_net(NETS[[nm]], rownames(smat))
  a <- as.data.table(run_ulm(mat = smat, net = net, .source = "tf",
                             .target = "target", .mor = "mor", minsize = 5))
  a <- a[statistic == "ulm"]
  a[, `:=`(network = nm, FDR = p.adjust(p_value, "BH"))]
  act_con[[nm]] <- a[, .(network, contrast = condition, tf = source,
                         score, p_value, FDR)]
}
act_con <- rbindlist(act_con)
say("\n--- contrast-level TF activity, hits at FDR < 0.05 ---")
says(act_con[FDR < 0.05][order(network, p_value),
      .(network, tf, score = round(score, 2), p = signif(p_value, 3),
        FDR = signif(FDR, 3))])

## --------------------------------------------------------- 6. consistency
cons <- rbindlist(list(
  ## sample level, myeloid_all vs subcluster 0
  merge(tests[layer == "myeloid_all", .(network, tf, d_all = paired_delta)],
        tests[layer == "subcluster_0", .(network, tf, d_sub0 = paired_delta)],
        by = c("network", "tf"))[, {
          ct <- suppressWarnings(cor.test(d_all, d_sub0, method = "spearman",
                                          exact = FALSE))
          .(check = "sample_level_myeloid_all_vs_subcluster0", n = .N,
            rho = unname(ct$estimate), p = ct$p.value)
        }, by = network],
  ## sample level (subcluster 0) vs contrast level
  merge(tests[layer == "subcluster_0", .(network, tf, d = paired_delta)],
        act_con[, .(network, tf, sc = score)], by = c("network", "tf"))[, {
          ct <- suppressWarnings(cor.test(d, sc, method = "spearman", exact = FALSE))
          .(check = "sample_level_vs_contrast_level", n = .N,
            rho = unname(ct$estimate), p = ct$p.value)
        }, by = network]))
say("\n--- internal consistency (Spearman across TFs) ---")
says(cons[, .(network, check, n, rho = round(rho, 3), p = signif(p, 3))])

## ------------------------------------------------------------- 7. headline
hdr <- function(lay, net, fdr = 0.05) {
  d <- tests[layer == lay & network == net & paired_FDR < fdr][order(paired_p)]
  say(sprintf("\n--- %s | %s | paired PN-PP, FDR < %.2f : %d of %d TFs ---",
              lay, net, fdr, nrow(d), tests[layer == lay & network == net, .N]))
  says(d[, .(tf, median_NS = round(median_NS, 2), median_PN = round(median_PN, 2),
             median_PP = round(median_PP, 2), paired_delta = round(paired_delta, 2),
             paired_p = signif(paired_p, 2), paired_FDR = signif(paired_FDR, 2),
             depth_adj_FDR = signif(depth_adj_FDR, 2))])
  invisible(d)
}
h_mye <- hdr("myeloid_all", "AB")
h_sub0 <- hdr("subcluster_0", "AB")

say("\n--- hits surviving depth adjustment (depth_adj_FDR < 0.05) ---")
says(tests[depth_adj_FDR < 0.05 & !is.na(depth_adj_FDR),
           .(n = .N), by = .(layer, network)])
say("\n--- robust hits (all three criteria) ---")
says(tests[robust == TRUE][order(layer, network, paired_p),
      .(layer, network, tf, paired_delta = round(paired_delta, 2),
        paired_p = signif(paired_p, 2), paired_FDR = signif(paired_FDR, 2),
        depth_adj_FDR = signif(depth_adj_FDR, 2),
        rho_delta_depth = round(rho_delta_depth, 2))])

FOCUS <- c("IRF1", "IRF7", "IRF8", "STAT1", "STAT2", "STAT3", "STAT5A", "STAT5B",
           "RELA", "NFKB1", "NFKB2", "REL", "JUN", "JUNB", "JUND", "FOS", "FOSB",
           "CEBPA", "CEBPB", "SPI1", "MITF", "PPARG", "NR4A1", "BACH1", "NFE2L2",
           "HIF1A", "EGR1", "KLF4", "MYC", "CTNNB1", "TCF7L2", "LEF1", "RUNX1",
           "RUNX3", "ATF4", "XBP1")
say("\n--- focus TFs, reported regardless of FDR (network AB) ---")
says(tests[network == "AB" & tf %in% FOCUS,
           .(layer, tf, median_NS = round(median_NS, 2), median_PN = round(median_PN, 2),
             median_PP = round(median_PP, 2), paired_delta = round(paired_delta, 2),
             paired_p = signif(paired_p, 2), paired_FDR = signif(paired_FDR, 2))][
       order(layer, paired_p)])
say("focus TFs present in regulon: ", sum(FOCUS %in% unique(tests$tf)),
    " of ", length(FOCUS), " (absent ones are not assessable, not negative)")

## ---------------------------------------------------------------- 8. write
samples_out <- rbind(
  merge(act_mye, meta_pid[, .(condition = pid, donor, group, n_cells)],
        by = "condition")[, layer := "myeloid_all"],
  merge(act_sub0, meta0[, .(condition = pid, donor, group, n_cells)],
        by = "condition")[, layer := "subcluster_0"]
)[, .(layer, network, pseudo_id = condition, donor, group, n_cells,
      tf, score, run_ulm_p = p_value)]

fwrite(ns,          file.path(TAB_DIR, "S13_tf_network_summary.csv"))
fwrite(samples_out, file.path(TAB_DIR, "S13_tf_activity_samples.csv"))
fwrite(tests,       file.path(TAB_DIR, "S13_tf_activity_tests.csv"))
fwrite(act_con,     file.path(TAB_DIR, "S13_tf_contrast.csv"))
fwrite(depth,       file.path(TAB_DIR, "S13_tf_depth_confounding.csv"))
fwrite(cons,        file.path(TAB_DIR, "S13_tf_consistency.csv"))
say("\nWROTE S13_tf_{network_summary,activity_samples,activity_tests,contrast,",
    "depth_confounding,consistency}.csv")
say("DONE (exploratory; TF activity = mRNA-weighted regulon score, not measured ",
    "protein activity; depth-adjusted sensitivity reported alongside)")
