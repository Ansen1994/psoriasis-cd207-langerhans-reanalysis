## 05b_s2_stats.R - S2 steps 4-5 + Gate C objective inputs.
## O8: collapse technical replicates into donor x state pseudo-samples.
## O1: sample-level ordered JT permutation trend test (main inference).
## O2: donor x subcluster pseudobulk + DESeq2 LRT (never cell-level tests).
## ASCII ONLY. Working directory = project root. NO figures are written.
source("config/00_paths.R")
source("config/01_params.R")
source("R/03_trend.R")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
  library(lme4); library(DESeq2)
})

con <- file(file.path(LOG_DIR, "S2_stats.log"), open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

set.seed(SEED)
t0 <- Sys.time()
cat("=== S2 steps 4-5: donor-collapsed trend + pseudobulk ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")

my <- readRDS(file.path(PROC_DIR, F_S2_OBJECT))
md <- as.data.table(my@meta.data)
stopifnot(all(c("donor_base", "group", "state_ord", "subcluster") %in% colnames(md)))

## ------------------------------------------------------------ S2 QC flag
## DATA QUALITY, not a significance tune. S1 assigned lineages by taking the
## ARGMAX of seven AddModuleScores, so a cluster with no strong lineage signal
## can be labelled "myeloid" merely for being least negative. S2 markers then
## exposed one such cluster with sweat-gland / secretory epithelial markers.
## Rule (applied identically to every subcluster): >=3 of a cluster's Top10
## markers belonging to the pre-defined epithelial/secretory set => suspect.
## Verified independently against mean expression of diagnostic genes below.
epi_set <- c("DCD","SCGB1D2","SCGB2A2","SCGB1D1","AQP5","PIP","MUCL1",
             "AZGP1","KRT19","KRT8","KRT18","KRT14","EPCAM","SCGB3A1")
top <- data.table::fread(file.path(TAB_DIR, "S2_markers_top10.csv"))
qc <- top[, .(n_epi_top10 = sum(gene %in% epi_set)), by = cluster]
qc[, subcluster := as.character(cluster)]
qc[, qc_eligible := n_epi_top10 < 3]
qc[, qc_flag := ifelse(qc_eligible, "myeloid_lineage", "non_myeloid_suspect")]

## Independent expression check (no marker overlap with the rule above).
diag_genes <- c("DCD","SCGB2A2","AQP5","KRT19","LYZ","C1QA","CD68","CD207")
av <- tryCatch(
  Seurat::AverageExpression(my, features = intersect(diag_genes, rownames(my)),
                            group.by = "subcluster", layer = "data",
                            verbose = FALSE),
  error = function(e) NULL)
if (!is.null(av)) {
  e <- as.data.table(as.data.frame(av[[1]]), keep.rownames = "gene")
  cn <- colnames(e); cn[-1] <- sub("^g", "", cn[-1]); data.table::setnames(e, cn)
  el <- melt(e, id.vars = "gene", variable.name = "subcluster",
             value.name = "mean_expr")
  el[, subcluster := as.character(subcluster)]
  for (g1 in c("DCD", "LYZ")) {
    v <- el[gene == g1, .(subcluster, v = round(mean_expr, 2))]
    data.table::setnames(v, "v", paste0("expr_", g1))
    qc <- merge(qc, v, by = "subcluster", all.x = TRUE)
  }
}
data.table::fwrite(qc, file.path(TAB_DIR, "S2_subcluster_qc.csv"))
cat("--- S2 subcluster QC flags ---\n")
print(qc[order(subcluster)])
cat("\nRule: >=3 epithelial/secretory genes in Top10 => non_myeloid_suspect.\n")
cat("Defined at S2 after marker inspection; does NOT touch any p-value or\n")
cat("significance threshold, and is applied to ALL subclusters equally.\n\n")

## Keep md in the SAME row order as the count matrix; merge() would reorder it.
md[, qc_eligible := qc$qc_eligible[match(subcluster, qc$subcluster)]]
md[, qc_flag     := qc$qc_flag[match(subcluster, qc$subcluster)]]
stopifnot(!anyNA(md$qc_eligible))
cat("QC-eligible myeloid cells: ", sum(md$qc_eligible), " / ", nrow(md),
    " (excluded ", sum(!md$qc_eligible), ")\n", sep = "")

## ------------------------------------------------------------ O8 collapse
## A donor can have BOTH PN and PP, so the pseudo-sample key MUST include the
## disease state. Collapsing on donor alone would merge lesional and
## non-lesional samples and destroy the three-state design.
md[, pseudo_id := paste0(donor_base, "_", group)]
cat("cells: ", nrow(md), " | raw samples: ", uniqueN(md$sample_id),
    " | donor x state pseudo-samples: ", uniqueN(md$pseudo_id), "\n", sep = "")
cat("pseudo-samples by state:\n")
print(md[, .(n_pseudo = uniqueN(pseudo_id)), by = .(group, state_ord)][order(state_ord)])
cat("technical replicates collapsed (raw sample count - pseudo count): ",
    uniqueN(md$sample_id) - uniqueN(md$pseudo_id), "\n", sep = "")

## ------------------------------------------------------------ step 4: proportions
## MAIN: denominator = QC-eligible myeloid cells only. Excluded cells are not
## silently deleted; the full 4-cluster control is written separately below.
mdc <- md[qc_eligible == TRUE]
cnt <- mdc[, .(n_cells = .N),
           by = .(pseudo_id, donor_base, group, state_ord, subcluster, batch)]
tot <- mdc[, .(n_myeloid = .N), by = pseudo_id]
cnt <- merge(cnt, tot, by = "pseudo_id")
cnt[, prop := 100 * n_cells / n_myeloid]
cnt[, low_n := n_cells < PB_MIN_CELLS]
setorder(cnt, subcluster, state_ord, pseudo_id)
data.table::fwrite(cnt, file.path(TAB_DIR, T_S2_PROP))

## Transparency control: all 4 clusters, original denominator.
cnt_all <- md[, .(n_cells = .N),
              by = .(pseudo_id, group, state_ord, subcluster)]
tot_all <- md[, .(n_all = .N), by = pseudo_id]
cnt_all <- merge(cnt_all, tot_all, by = "pseudo_id")
cnt_all[, prop_all := 100 * n_cells / n_all]
cnt_all <- merge(cnt_all, qc[, .(subcluster, qc_eligible)], by = "subcluster",
                 all.x = TRUE)
data.table::fwrite(cnt_all, file.path(TAB_DIR, "S2_proportions_all_clusters.csv"))
cat("\nMAIN denominator: QC-eligible myeloid cells only.\n")
cat("CONTROL written: results/tables/S2_proportions_all_clusters.csv (all 4 clusters)\n")
cat("\nproportion rows (main): ", nrow(cnt),
    " | combos below PB_MIN_CELLS=", PB_MIN_CELLS, ": ", sum(cnt$low_n),
    "\n", sep = "")
cat("n_myeloid per pseudo-sample: min=", min(cnt$n_myeloid),
    " median=", stats::median(cnt$n_myeloid),
    " max=", max(cnt$n_myeloid), "\n", sep = "")

cat("\n--- median subcluster proportion (%) by state ---\n")
piv <- cnt[, .(median_pct = round(stats::median(prop), 2), n_pseudo = .N),
           by = .(subcluster, state_ord)]
print(piv[order(subcluster, state_ord)])

## Main inference: pre-registered one-sided JT permutation trend, BH-FDR.
dat <- cnt[, .(feature = paste0("sub_", subcluster), value = prop,
               state_ord = state_ord)]
tr <- jt_trend_table(dat, B = TREND_NPERM, seed = SEED,
                     alternative = TREND_ALTERNATIVE)
setDT(tr)
tr[, subcluster := sub("^sub_", "", feature)]
data.table::fwrite(tr, file.path(TAB_DIR, T_S2_TREND))
cat("\n--- MAIN JT ordered trend (B=", TREND_NPERM, ", one-sided greater) ---\n", sep = "")
print(tr[, .(subcluster, n_samples, p, FDR, monotone, pass)])

## Robustness: drop noisy low-cell combinations and re-run. This is a
## sensitivity check, NOT the main test and NOT a threshold change.
dat_s <- cnt[low_n == FALSE,
             .(feature = paste0("sub_", subcluster), value = prop,
               state_ord = state_ord)]
tr_s <- jt_trend_table(dat_s, B = TREND_NPERM, seed = SEED,
                       alternative = TREND_ALTERNATIVE)
setDT(tr_s)
tr_s[, subcluster := sub("^sub_", "", feature)]
tr_s[, check := "drop_low_n"]
data.table::fwrite(tr_s, file.path(TAB_DIR, "S2_trend_sensitivity_lowN.csv"))
cat("\n--- sensitivity: drop combos with n_cells < ", PB_MIN_CELLS, " ---\n", sep = "")
print(tr_s[, .(subcluster, n_samples, p, FDR, monotone)])

## Control: full 4-cluster denominator. This is what the previous feasibility
## probe effectively measured, and it shows how much of the apparent drift
## came from the excluded non-myeloid cluster.
dat_all <- cnt_all[, .(feature = paste0("sub_", subcluster),
                       value = prop_all, state_ord = state_ord)]
tr_all <- jt_trend_table(dat_all, B = TREND_NPERM, seed = SEED,
                         alternative = TREND_ALTERNATIVE)
setDT(tr_all)
tr_all[, subcluster := sub("^sub_", "", feature)]
data.table::fwrite(tr_all, file.path(TAB_DIR, "S2_trend_all_clusters_control.csv"))
cat("\n--- CONTROL: all 4 clusters, full denominator ---\n")
print(tr_all[, .(subcluster, n_samples, p, FDR, monotone)])

## Supplementary mixed model. PB_DESIGN is ~ state_ord + batch, but batch is
## a CONSTANT (single platform), so including it makes the design singular.
## Omitting a constant covariate is algebraic simplification, NOT a change of
## the pre-registered design; recorded here for Methods.
cat("\nDESIGN NOTE: batch is constant (", unique(md$batch),
    ") -> ~ state_ord + batch degenerates to ~ state_ord.\n", sep = "")
lmer_rows <- lapply(sort(unique(cnt$subcluster)), function(s) {
  d <- as.data.frame(cnt[subcluster == s])
  out <- data.table(subcluster = s, beta = NA_real_, lrt_chisq = NA_real_,
                    lmer_p = NA_real_, singular = NA, status = "ok")
  full <- try(lme4::lmer(prop ~ state_ord + (1 | donor_base),
                         data = d, REML = FALSE), silent = TRUE)
  null <- try(lme4::lmer(prop ~ 1 + (1 | donor_base),
                         data = d, REML = FALSE), silent = TRUE)
  if (inherits(full, "try-error") || inherits(null, "try-error")) {
    out$status <- "model_failed"; return(out)
  }
  a <- as.data.frame(anova(null, full))
  pc <- grep("^Pr", colnames(a), value = TRUE)[1]
  out$beta <- unname(lme4::fixef(full)["state_ord"])
  out$lrt_chisq <- a[2, "Chisq"]
  out$lmer_p <- if (!is.na(pc)) a[2, pc] else NA_real_
  out$singular <- lme4::isSingular(full)
  out
})
lmer_tab <- rbindlist(lmer_rows)
data.table::fwrite(lmer_tab, file.path(TAB_DIR, "S2_lmer_supplement.csv"))
cat("\n--- supplementary lmer (LRT vs intercept-only) ---\n")
print(lmer_tab)

## ------------------------------------------------------------ step 5: pseudobulk
get_counts <- function(o) {
  tryCatch(Seurat::GetAssayData(o, assay = "RNA", layer = "counts"),
           error = function(e) Seurat::GetAssayData(o, assay = "RNA", slot = "counts"))
}
counts <- get_counts(my)
grp <- paste0(md$pseudo_id, "|", md$subcluster)
lev <- sort(unique(grp))
idx <- match(grp, lev)
M <- Matrix::sparseMatrix(i = seq_along(idx), j = idx, x = 1,
                          dims = c(length(idx), length(lev)))
pb <- counts %*% M
pb <- as.matrix(pb)
colnames(pb) <- lev
rm(counts, M); gc(FALSE)
cat("\npseudobulk matrix: ", nrow(pb), " genes x ", ncol(pb),
    " donor x subcluster columns\n", sep = "")

## Pseudobulk is computed on QC-eligible clusters only; the excluded cluster
## is not myeloid, so treating it as a disease-relevant subcluster would be a
## biological error. Its counts remain in the transparency control table.
info <- cnt[, .(combo = paste0(pseudo_id, "|", subcluster), pseudo_id,
                donor_base, group, state_ord, subcluster, n_cells)]
setorder(info, subcluster, combo)
cat("\npseudobulk computed for QC-eligible subclusters only: ",
    paste(sort(unique(info$subcluster)), collapse = ", "), "\n", sep = "")
stopifnot(all(info$combo %in% colnames(pb)))

deg_list <- list()
pb_rows <- list()
for (s in sort(unique(info$subcluster))) {
  d <- info[subcluster == s & n_cells >= PB_MIN_CELLS]
  n_drop <- info[subcluster == s & n_cells <  PB_MIN_CELLS, .N]
  cat("\n[", s, "] pseudo-samples kept ", nrow(d), " / ",
      info[subcluster == s, .N], " (dropped ", n_drop,
      " below ", PB_MIN_CELLS, " cells)\n", sep = "")
  pb_rows[[s]] <- data.table(subcluster = s, n_pseudo_kept = nrow(d),
                             n_pseudo_dropped = n_drop)
  if (nrow(d) < 6) {
    cat("  SKIPPED: too few pseudo-samples for DESeq2 LRT\n")
    next
  }
  mat <- pb[, d$combo, drop = FALSE]
  keep <- (Matrix::rowSums(mat) >= 10) & (Matrix::rowSums(mat > 0) >= 3)
  mat <- mat[keep, , drop = FALSE]
  cdata <- data.frame(row.names = d$combo,
                      state_ord = as.numeric(d$state_ord),
                      donor_base = as.character(d$donor_base))
  res <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(mat),
                                          colData = cdata,
                                          design = ~ state_ord)
    ## poscounts + fitType=mean avoid the "every gene has a zero" failure
    ## that is common with many small pseudobulk samples.
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq2::estimateDispersions(dds, fitType = "mean")
    dds <- DESeq2::nbinomLRT(dds, reduced = ~ 1)
    as.data.table(as.data.frame(DESeq2::results(dds)), keep.rownames = "gene")
  }, error = function(e) {
    cat("  DESeq2 FAILED: ", conditionMessage(e), "\n", sep = "")
    NULL
  })
  if (!is.null(res) && nrow(res) > 0) {
    res[, subcluster := s]
    deg_list[[s]] <- res
    cat("  genes tested: ", nrow(res),
        " | FDR<", PB_FDR, ": ", sum(res$padj < PB_FDR, na.rm = TRUE),
        "\n", sep = "")
  }
}
pb_tab <- rbindlist(pb_rows)
data.table::fwrite(pb_tab, file.path(TAB_DIR, "S2_pseudobulk_sample_counts.csv"))
if (length(deg_list) > 0) {
  deg <- rbindlist(deg_list, fill = TRUE)
  data.table::fwrite(deg, file.path(TAB_DIR, T_S2_DEG))
  cat("\ntotal pseudobulk DEG rows: ", nrow(deg), "\n", sep = "")
} else {
  cat("\nWARNING: no pseudobulk DEG results produced.\n")
}

## ------------------------------------------------------------ Gate C inputs
g <- tr[, .(subcluster, JT_p = p, JT_FDR = FDR, monotone, trend_pass = pass)]
g <- merge(g, qc[, .(subcluster, qc_eligible, qc_flag, n_epi_top10)],
           by = "subcluster", all.x = TRUE)
g <- merge(g, lmer_tab[, .(subcluster, lmer_p, singular)], by = "subcluster", all.x = TRUE)
g <- merge(g, pb_tab, by = "subcluster", all.x = TRUE)
if (length(deg_list) > 0) {
  nd <- deg[, .(n_DEG_FDR05 = sum(padj < PB_FDR, na.rm = TRUE),
                n_genes_tested = .N), by = subcluster]
  g <- merge(g, nd, by = "subcluster", all.x = TRUE)
} else {
  g[, `:=`(n_DEG_FDR05 = NA_integer_, n_genes_tested = NA_integer_)]
}
g[, gateC_subcluster_pass := (JT_FDR < TREND_FDR) & monotone &
    (qc_eligible %in% TRUE)]
setorder(g, JT_p)
data.table::fwrite(g, file.path(TAB_DIR, "S2_gateC_inputs.csv"))
cat("\n--- Gate C objective inputs (MAIN = QC-eligible clusters only) ---\n")
print(g)
GATE_C_PASS <- any(g$gateC_subcluster_pass, na.rm = TRUE)
cat("\nGATE_C_PASS=", GATE_C_PASS,
    " (criterion: >=1 QC-eligible subcluster with JT BH-FDR<", TREND_FDR,
    " AND monotone direction)\n", sep = "")
cat("If FALSE, apply the frozen downgrade: two-group narrative + retitle.\n")
cat("Do NOT switch compartment, tune parameters, or delete samples.\n")
cat("\nOUTPUT: results/tables/S2_proportions_by_sample.csv\n")
cat("OUTPUT: results/tables/S2_trend_tests.csv\n")
cat("OUTPUT: results/tables/S2_pseudobulk_DEG.csv\n")
cat("OUTPUT: results/tables/S2_gateC_inputs.csv\n")
cat("NO FIGURES WRITTEN\n")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")
cat("=== S2 steps 4-5 DONE ===\n")
