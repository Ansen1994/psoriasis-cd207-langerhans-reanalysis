## S11_lc_maturation.R
## Bridge analysis added 2026-11-13 in response to internal review.
##
## Purpose: the manuscript asserts "compartment-selective change" on the basis of
## bulk data only (CD207 down, LAMP3 up). This script tests the same axis inside
## the single-cell data at the donor x state level, so the claim is supported by
## the discovery cohort rather than inferred from bulk alone.
##
## Two readouts per pseudo-sample, because subcluster 2 holds only 352 cells and a
## per-sample pseudobulk of ~12 cells is depth-sensitive:
##   (a) pseudobulk CPM  (sum of counts / total counts * 1e6, log2)
##   (b) detection fraction (share of cells with count > 0)
##
## Two compartments: all QC-eligible myeloid cells (n = 1,666) and the CD207+
## subcluster 2 only (n = 352). Everything here is EXPLORATORY: no alternative was
## pre-registered for these genes, so no one-sided claim is made.
##
## Outputs
##   results/tables/S11_lc_maturation_pseudobulk.csv
##   results/tables/S11_lc_maturation_tests.csv
##   logs/S11_lc_maturation.log

source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({library(Seurat); library(data.table)})

logf <- file.path(LOG_DIR, "S11_lc_maturation.log")
say <- function(...) {
  m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE)
}
cat("=== S11 LC maturation / migration axis (exploratory) ===\n", file = logf)

obj <- readRDS(file.path(PROC_DIR, F_S2_OBJECT))
meta <- obj@meta.data
cnt <- LayerData(obj, assay = "RNA", layer = "counts")
GENES <- c("CD207", "CD1A", "LAMP3", "CD83", "CCR7", "CD1C", "HLA-DRA", "CD274")
GENES <- intersect(GENES, rownames(cnt))
say("genes quantified: ", paste(GENES, collapse = ", "))

meta$pseudo_id <- paste(meta$donor_base, meta$group, sep = "_")

build <- function(cells, label) {
  m <- meta[cells, , drop = FALSE]
  x <- cnt[GENES, cells, drop = FALSE]
  out <- rbindlist(lapply(split(rownames(m), m$pseudo_id), function(ids) {
    tot <- sum(x[, ids, drop = FALSE])
    pb <- log2((rowSums(x[, ids, drop = FALSE]) / tot) * 1e6 + 1)
    det <- rowMeans(x[, ids, drop = FALSE] > 0)
    data.table(compartment = label, pseudo_id = ids[1],
               donor = m[ids[1], "donor_base"], group = m[ids[1], "group"],
               state_ord = as.integer(m[ids[1], "state_ord"]),
               n_cells = length(ids),
               gene = GENES,
               pseudobulk_log2cpm = as.numeric(pb),
               detection_fraction = as.numeric(det))
  }))
  out
}

all_mye <- build(which(meta$subcluster != 1), "myeloid_qc_eligible")
sub2    <- build(which(meta$subcluster == 2), "subcluster2_CD207")
pb <- rbind(all_mye, sub2)
say("\n--- pseudo-samples per compartment / group ---")
print(pb[, .(n_pseudo = uniqueN(pseudo_id), median_cells = median(n_cells)),
         by = .(compartment, group)][order(compartment, group)])

## tests: Kruskal-Wallis across three states, paired PN-PP, and Spearman CD207 vs LAMP3
tests <- list()
for (comp in unique(pb$compartment)) {
  for (metric in c("pseudobulk_log2cpm", "detection_fraction")) {
    d <- dcast(pb[compartment == comp], pseudo_id + donor + group ~ gene,
               value.var = metric)
    for (g in GENES) {
      if (!g %in% names(d)) next
      kw <- kruskal.test(d[[g]] ~ factor(d$group, levels = c("NS", "PN", "PP")))
      w <- dcast(d, donor ~ group, value.var = g)
      w <- w[!is.na(PN) & !is.na(PP)]
      p_pair <- if (nrow(w) >= 5) suppressWarnings(wilcox.test(w$PP, w$PN, paired = TRUE)$p.value) else NA_real_
      tests[[length(tests) + 1]] <- data.table(
        compartment = comp, metric = metric, gene = g,
        n_pseudo = nrow(d),
        median_NS = if (any(d$group == "NS")) median(d[group == "NS"][[g]]) else NA_real_,
        median_PN = if (any(d$group == "PN")) median(d[group == "PN"][[g]]) else NA_real_,
        median_PP = if (any(d$group == "PP")) median(d[group == "PP"][[g]]) else NA_real_,
        kw_p = kw$p.value,
        n_pairs_PN_PP = nrow(w),
        paired_p_PN_PP = p_pair)
    }
  }
}
tt <- rbindlist(tests)
tt[, kw_FDR := p.adjust(kw_p, "BH"), by = .(compartment, metric)]
tt[, paired_FDR := p.adjust(paired_p_PN_PP, "BH"), by = .(compartment, metric)]
setorder(tt, compartment, metric, -kw_p)
say("\n--- tests (exploratory; BH within compartment x metric) ---")
print(tt[, .(compartment, metric, gene, median_NS, median_PN, median_PP,
             kw_p = signif(kw_p, 3), kw_FDR = signif(kw_FDR, 3),
             n_pairs_PN_PP, paired_p = signif(paired_p_PN_PP, 3))])

## the axis question: does LAMP3 rise while CD207 falls?
say("\n--- CD207 vs LAMP3 within pseudo-sample (Spearman, per compartment x metric) ---")
cor_out <- list()
for (comp in unique(pb$compartment)) {
  for (metric in c("pseudobulk_log2cpm", "detection_fraction")) {
    d <- dcast(pb[compartment == comp], pseudo_id + group ~ gene, value.var = metric)
    if (!all(c("CD207", "LAMP3") %in% names(d))) next
    ct <- suppressWarnings(cor.test(d$CD207, d$LAMP3, method = "spearman", exact = FALSE))
    ratio <- d[, .(median_ratio = median(LAMP3 - CD207)), by = group]
    cor_out[[length(cor_out) + 1]] <- data.table(
      compartment = comp, metric = metric, n = nrow(d),
      rho_CD207_LAMP3 = unname(ct$estimate), p = ct$p.value)
    say(sprintf("  %-20s %-18s n=%2d  rho = %+.3f  p = %.3g",
                comp, metric, nrow(d), ct$estimate, ct$p.value))
    say(paste0("      median(LAMP3-CD207) by group: ",
               paste(sprintf("%s %+.2f", ratio$group, ratio$median_ratio), collapse = " | ")))
  }
}
cort <- rbindlist(cor_out)

fwrite(pb, file.path(TAB_DIR, "S11_lc_maturation_pseudobulk.csv"))
fwrite(tt, file.path(TAB_DIR, "S11_lc_maturation_tests.csv"))
fwrite(cort, file.path(TAB_DIR, "S11_lc_maturation_correlation.csv"))
say("\nWROTE S11_lc_maturation_{pseudobulk,tests,correlation}.csv")
say("DONE (exploratory; no pre-registered alternative exists for these genes)")
