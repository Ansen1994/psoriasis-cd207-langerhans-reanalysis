## S3_function.R - two-group functional analysis after Gate C downgrade.
## Gate C FAILED (three-state progressive trend). The frozen downgrade reframes
## the study as PP (lesional) vs NS (healthy), with PN as a DESCRIPTIVE
## intermediate only. This script therefore builds a fresh two-group
## pseudobulk contrast as the MAIN functional input, and keeps the S2 ordered
## trend LRT genes as an explicitly labelled EXPLORATORY comparison.
##
## Offline branch B (O5): CellChat / Monocle3 are not installed. CellChat and
## trajectory are mutually exclusive in the frozen route; neither is done here.
## Function = enrichment + GSEA + GSVA on sample-level pseudobulk.
##
## All outputs stay inside the project. ASCII only (LC_CTYPE=C).
source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix); library(DESeq2)
  library(clusterProfiler); library(org.Hs.eg.db); library(fgsea)
  library(msigdbr); library(GSVA)
})

logf <- file.path(LOG_DIR, "S3_function.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S3 functional analysis (two-group framing) ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
cat("frozen params: ", PROJ_VERSION, " | PB_MIN_CELLS=", PB_MIN_CELLS,
    " | PB_FDR=", PB_FDR, " | PATHWAY_FDR=", PATHWAY_FDR, "\n", sep = "")
cat("NOTE: CellChat/trajectory NOT run (offline branch B, O5).\n")
set.seed(SEED)
t0 <- Sys.time()

## ---------------------------------------------------------------- 0. inputs
qc <- fread(file.path(TAB_DIR, "S2_subcluster_qc.csv"))
elig <- qc[qc_eligible == TRUE, as.character(subcluster)]
cat("QC-eligible subclusters: ", paste(elig, collapse = ", "), "\n", sep = "")

o <- readRDS(file.path(PROC_DIR, F_S2_OBJECT))
md <- as.data.table(o@meta.data, keep.rownames = "cell")
md[, qc_eligible := subcluster %in% elig]
md[, pseudo_id := paste0(donor_base, "_", group)]
keep_cells <- md[qc_eligible == TRUE, cell]
cnt <- GetAssayData(o, assay = "RNA", layer = "counts")[, keep_cells]
md <- md[cell %in% keep_cells]
stopifnot(identical(colnames(cnt), md$cell))
cat("eligible cells: ", nrow(md), " | pseudo-samples: ", uniqueN(md$pseudo_id), "\n", sep = "")
rm(o); gc(FALSE)

## ---------------------------------------------------------------- 1. pseudobulk
## Identical construction to R/05b_s2_stats.R (O8 donor x state collapse).
grp <- paste0(md$pseudo_id, "|", md$subcluster)
lev <- sort(unique(grp))
M <- sparseMatrix(i = seq_along(grp), j = match(grp, lev), x = 1,
                  dims = c(length(grp), length(lev)))
pb <- as.matrix(cnt %*% M)
colnames(pb) <- lev
rm(cnt, M); gc(FALSE)

info <- md[, .(n_cells = .N),
           by = .(combo = paste0(pseudo_id, "|", subcluster), pseudo_id,
                  donor_base, group, state_ord, subcluster)]
setorder(info, subcluster, combo)
info[, group2 := ifelse(group %in% c("NS", "PP"), group, "PN")]
saveRDS(list(pb = pb, info = info), file.path(PROC_DIR, "S3_pseudobulk.rds"))
cat("pseudobulk matrix: ", nrow(pb), " genes x ", ncol(pb), " combos\n", sep = "")

## ---------------------------------------------------------------- 2. two-group DEG
## MAIN functional input. Wald test on donor-level pseudobulk, NS vs PP only.
## PN is retained in the object for descriptive plotting but NOT tested.
deg2_list <- list()
for (s in sort(unique(info$subcluster))) {
  d <- info[subcluster == s & n_cells >= PB_MIN_CELLS & group2 %in% c("NS", "PP")]
  n_ns <- sum(d$group2 == "NS"); n_pp <- sum(d$group2 == "PP")
  cat("\n[subcluster ", s, "] NS=", n_ns, " PP=", n_pp,
      " (dropped ", sum(info[subcluster == s]$n_cells < PB_MIN_CELLS),
      " low-cell combos, ", sum(info[subcluster == s]$group2 == "PN"),
      " PN kept descriptive only)\n", sep = "")
  if (n_ns < 2 || n_pp < 2) {
    cat("  SKIPPED: fewer than 2 pseudo-samples in one group (no two-group test)\n")
    next
  }
  mat <- pb[, d$combo, drop = FALSE]
  keep <- (rowSums(mat) >= 10) & (rowSums(mat > 0) >= 3)
  mat <- mat[keep, , drop = FALSE]
  cdata <- data.frame(row.names = d$combo,
                      group2 = factor(d$group2, levels = c("NS", "PP")))
  res <- tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = round(mat), colData = cdata,
                                  design = ~ group2)
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- estimateDispersions(dds, fitType = "mean")
    dds <- nbinomWaldTest(dds)
    as.data.table(as.data.frame(results(dds, contrast = c("group2", "PP", "NS"))),
                  keep.rownames = "gene")
  }, error = function(e) {
    cat("  DESeq2 FAILED: ", conditionMessage(e), "\n", sep = ""); NULL
  })
  if (!is.null(res) && nrow(res) > 0) {
    res[, `:=`(subcluster = s, contrast = "PP_vs_NS", n_NS = n_ns, n_PP = n_pp,
               universe_size = nrow(mat))]
    deg2_list[[s]] <- res
    cat("  genes tested: ", nrow(res), " | FDR<", PB_FDR, ": ",
        sum(res$padj < PB_FDR, na.rm = TRUE), "\n", sep = "")
  }
}
if (length(deg2_list)) {
  deg2 <- rbindlist(deg2_list, fill = TRUE)
  fwrite(deg2, file.path(TAB_DIR, "S3_two_group_DEG.csv"))
} else {
  deg2 <- data.table(gene = character(), log2FoldChange = numeric(),
                     padj = numeric(), subcluster = character())
  fwrite(deg2, file.path(TAB_DIR, "S3_two_group_DEG.csv"))
  cat("WARNING: no two-group DEG produced.\n")
}

## ---------------------------------------------------------------- 3. gene sets
## Local MSigDB Hallmark (offline, msigdbr package data). No online KEGG call.
hall <- msigdbr(species = "Homo sapiens", collection = "H")
hsets <- split(hall$gene_symbol, hall$gs_name)
hsets <- hsets[vapply(hsets, length, 1L) >= 10]
cat("\nHallmark gene sets (size >= 10): ", length(hsets), "\n", sep = "")

## Core subcluster signatures for S4 bulk scoring, from frozen S2 markers.
mtop <- fread(file.path(TAB_DIR, "S2_markers_top10.csv"))
sig <- list(
  LC_CD207_signature = mtop[cluster == 2, gene],
  Mac_DC_signature   = mtop[cluster == 0, gene],
  Cycling_signature  = mtop[cluster == 3, gene]
)
## Add two-group up/down genes of the core subcluster when available.
if (nrow(deg2[subcluster == "2" & padj < PB_FDR & log2FoldChange > 1])) {
  u <- deg2[subcluster == "2" & padj < PB_FDR & log2FoldChange > 1]
  setorder(u, padj)
  sig$LC_PP_up_twoGroup <- head(u$gene, 20)
}
sig_tab <- rbindlist(lapply(names(sig), function(n)
  data.table(set_name = n, gene = sig[[n]])), fill = TRUE)
fwrite(sig_tab, file.path(TAB_DIR, "S3_signatures.csv"))
cat("signature sets: ", paste(names(sig), collapse = ", "), "\n", sep = "")

## ---------------------------------------------------------------- 4. enrichment
## enrichGO BP on (a) two-group DEG [MAIN] and (b) S2 ordered-trend LRT genes
## [EXPLORATORY control, explicitly labelled, never mixed in the same claim].
s2deg <- fread(file.path(TAB_DIR, T_S2_DEG))
gene_lists <- list()
if (nrow(deg2)) {
  for (s in sort(unique(deg2$subcluster))) {
    g <- deg2[subcluster == s & padj < PATHWAY_FDR & !is.na(padj), gene]
    if (length(g) >= 5) gene_lists[[paste0("twoGroup_sub", s)]] <- g
  }
}
for (s0 in sort(unique(s2deg$subcluster))) {
  s <- as.character(s0)
  g <- s2deg[as.character(subcluster) == s & padj < PATHWAY_FDR & !is.na(padj), gene]
  if (length(g) >= 5) gene_lists[[paste0("orderedTrend_sub", s)]] <- g
}
cat("\n--- gene lists for enrichment ---\n")
for (n in names(gene_lists)) cat(" ", n, ": ", length(gene_lists[[n]]), " genes\n", sep = "")

enrich_list <- list()
for (n in names(gene_lists)) {
  g <- gene_lists[[n]]
  s <- sub(".*_sub", "", n)
  universe <- unique(c(g, if (nrow(deg2)) deg2[subcluster == s, gene] else character()))
  e <- tryCatch({
    suppressMessages(clusterProfiler::enrichGO(
      gene = g, universe = universe, OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH",
      pvalueCutoff = PATHWAY_FDR, qvalueCutoff = PATHWAY_FDR, readable = FALSE))
  }, error = function(err) { cat("  enrichGO failed for ", n, ": ", conditionMessage(err), "\n", sep = ""); NULL })
  if (!is.null(e)) {
    z <- as.data.table(e@result)
    ## enrichResult@result holds ALL tested terms; filter to the frozen cutoff.
    if (nrow(z)) {
      z[, significant := (pvalue < PATHWAY_FDR) & (p.adjust < PATHWAY_FDR)]
      z[, `:=`(gene_list = n, subcluster = s, method = "enrichGO_BP")]
      enrich_list[[n]] <- z
    }
  }
  cat(" ", n, ": ", if (is.null(enrich_list[[n]])) 0 else
        sum(enrich_list[[n]]$significant, na.rm = TRUE),
      " significant BP terms (of ", if (is.null(enrich_list[[n]])) 0 else
        nrow(enrich_list[[n]]), " tested)\n", sep = "")
}
path_tab <- if (length(enrich_list)) rbindlist(enrich_list, fill = TRUE) else
  data.table(gene_list = character(), ID = character(), Description = character())
fwrite(path_tab, file.path(TAB_DIR, T_S3_PATHWAY))

## ---------------------------------------------------------------- 5. GSEA
## fgsea on full ranked log2FC lists (two-group MAIN, ordered-trend EXPLORATORY).
rank_lists <- list()
if (nrow(deg2)) for (s in sort(unique(deg2$subcluster))) {
  z <- deg2[subcluster == s & !is.na(log2FoldChange) & is.finite(log2FoldChange),
            .(gene, log2FoldChange)]
  z <- z[!duplicated(gene)]
  v <- setNames(z$log2FoldChange, z$gene)
  if (length(v) >= 200) rank_lists[[paste0("twoGroup_sub", s)]] <- v
}
for (s0 in sort(unique(s2deg$subcluster))) {
  s <- as.character(s0)
  z <- s2deg[as.character(subcluster) == s & !is.na(log2FoldChange) & is.finite(log2FoldChange),
             .(gene, log2FoldChange)]
  z <- z[!duplicated(gene)]
  v <- setNames(z$log2FoldChange, z$gene)
  if (length(v) >= 200) rank_lists[[paste0("orderedTrend_sub", s)]] <- v
}
gsea_list <- list()
for (n in names(rank_lists)) {
  set.seed(SEED)
  z <- tryCatch(suppressWarnings(fgsea::fgsea(pathways = hsets, stats = rank_lists[[n]],
                             minSize = 10, maxSize = 500, nperm = 10000,
                             nproc = 1)),
                error = function(err) { cat("  fgsea failed for ", n, ": ", conditionMessage(err), "\n", sep = ""); NULL })
  if (!is.null(z)) {
    z <- as.data.table(z); z[, gene_list := n]
    gsea_list[[n]] <- z
    cat(" ", n, ": ", nrow(z), " sets | padj<", PATHWAY_FDR, " = ",
        sum(z$padj < PATHWAY_FDR, na.rm = TRUE), "\n", sep = "")
  }
}
gsea_tab <- if (length(gsea_list)) rbindlist(gsea_list, fill = TRUE) else data.table(pathway = character())
fwrite(gsea_tab, file.path(TAB_DIR, "S3_gsea.csv"))

## ---------------------------------------------------------------- 6. GSVA
## Sample-level GSVA on log2-CPM pseudobulk, tested PP vs NS per subcluster.
## PN columns are scored but excluded from the two-group test.
gsva_list <- list()
for (s in sort(unique(info$subcluster))) {
  d <- info[subcluster == s & n_cells >= PB_MIN_CELLS]
  if (nrow(d) < 5) { cat("GSVA skipped for subcluster ", s, ": <5 pseudo-samples\n", sep = ""); next }
  mat <- pb[, d$combo, drop = FALSE]
  keep <- rowSums(mat) >= 10
  lib <- colSums(mat); lib[lib == 0] <- 1
  expr <- log2(sweep(mat[keep, , drop = FALSE], 2, lib, "/") * 1e6 + 1)
  res <- tryCatch({
    par <- GSVA::gsvaParam(exprData = expr, geneSets = hsets,
                           kcdf = "Gaussian", minSize = 10, maxSize = 500)
    sc <- as.data.table(as.data.frame(GSVA::gsva(par)), keep.rownames = "pathway")
  }, error = function(err) { cat("  GSVA failed for subcluster ", s, ": ", conditionMessage(err), "\n", sep = ""); NULL })
  if (is.null(res)) next
  long <- melt(res, id.vars = "pathway", variable.name = "combo", value.name = "score")
  long <- merge(long, d[, .(combo, group, group2, pseudo_id, n_cells)], by = "combo")
  long[, subcluster := s]
  out <- long[group2 %in% c("NS", "PP"), .(
    n_NS = sum(group2 == "NS"), n_PP = sum(group2 == "PP"),
    median_NS = median(score[group2 == "NS"]),
    median_PP = median(score[group2 == "PP"]),
    delta_PP_minus_NS = median(score[group2 == "PP"]) - median(score[group2 == "NS"]),
    p_wilcox = if (sum(group2 == "NS") >= 2 && sum(group2 == "PP") >= 2)
      suppressWarnings(wilcox.test(x = score[group2 == "PP"],
                                   y = score[group2 == "NS"],
                                   exact = FALSE)$p.value) else NA_real_),
    by = .(subcluster, pathway)]
  out[, FDR_BH := p.adjust(p_wilcox, method = "BH")]
  gsva_list[[s]] <- out
  cat("GSVA subcluster ", s, ": ", nrow(out), " sets | FDR<", PATHWAY_FDR, " = ",
      sum(out$FDR_BH < PATHWAY_FDR, na.rm = TRUE), "\n", sep = "")
}
gsva_tab <- if (length(gsva_list)) rbindlist(gsva_list, fill = TRUE) else data.table(pathway = character())
fwrite(gsva_tab, file.path(TAB_DIR, "S3_gsva.csv"))

cat("\n=== S3 DONE ===\n")
cat("OUTPUT: results/tables/S3_two_group_DEG.csv\n")
cat("OUTPUT: results/tables/S3_signatures.csv\n")
cat("OUTPUT: results/tables/S3_pathway.csv\n")
cat("OUTPUT: results/tables/S3_gsea.csv\n")
cat("OUTPUT: results/tables/S3_gsva.csv\n")
cat("OUTPUT: data/processed/S3_pseudobulk.rds\n")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")
