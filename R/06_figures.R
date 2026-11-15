## 06_figures.R - S5: ONE-PASS render of Fig 1-4 from frozen S1-S4 outputs.
## Working directory = project root. ASCII ONLY (this machine starts R with
## LC_CTYPE=C). This script does NOT re-run any analysis, change any frozen
## parameter, or compute any new inferential statistic. It only reads the
## already-frozen RDS objects and result tables and draws the four figures.
##
## Frozen figure spec: docs/technical-route doc (v1.1), section 6.
## Mandatory disclosures: logs/gates.md Gate C and Gate D sections.
##   Fig 2: subcluster 1 (648 sweat-gland / secretory epithelial cells) was
##          QC-EXCLUDED; Gate C FAILED; one-sided 'greater' was pre-registered
##          so the strong DECREASE in CD207+ cells is EXPLORATORY ONLY.
##   Fig 3: GSVA 0/50 sets at BH-FDR<0.05 (negative result, must be shown);
##          two-group DEG for subclusters 2/3 were SKIPPED (too few donors);
##          CellChat / trajectory were NOT run (offline branch B, O5).
##   Fig 4: AUC 0.9988 is INTERNAL 5-fold CV, NOT external validation;
##          whole-skin query against a myeloid-only reference is a domain
##          mismatch, transferred proportions are NOT quantitative;
##          mapping-free CD207+ counts are 3 vs 3, descriptive only.
##
## n is ALWAYS the number of donors / samples, NEVER the number of cells.
## Outputs: results/figures/Fig{1..4}_*.pdf (vector) + .png (300 dpi).
## Logs:    logs/S5_figures.log, logs/S5_sessionInfo.txt

source("config/00_paths.R")
source("config/01_params.R")
source("R/figure_cache.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(grid)
})

main <- function() {

set.seed(SEED)
t_start <- Sys.time()
cat("=== S5 one-pass figure rendering ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
cat("frozen params: ", PROJ_VERSION, " | seed=", SEED,
    " | dpi=", FIG_DPI, " | vector=", FIG_VECTOR, "\n", sep = "")
cat("ggplot2 ", as.character(packageVersion("ggplot2")),
    " | patchwork ", as.character(packageVersion("patchwork")),
    " | Seurat ", as.character(packageVersion("Seurat")), "\n", sep = "")

## ------------------------------------------------------------------ 0. setup
## Okabe-Ito palette (colour-blind safe). This overrides FIG_STATE_COLORS at the
## figure layer only; the frozen parameter file is not modified.
state_cols   <- c(NS = "#56B4E9", PN = "#E69F00", PP = "#D55E00")
bulk_cols    <- c(NN = "#56B4E9", PN = "#E69F00", PP = "#D55E00")
lineage_cols <- c(keratinocyte = "#0072B2", fibroblast  = "#E69F00",
                  Tcell        = "#009E73", endothelial = "#56B4E9",
                  myeloid      = "#D55E00", mast        = "#CC79A7")
sub_labs <- c("0" = "0 Mac/DC",
              "1" = "1 QC-EXCLUDED (sweat gland)",
              "2" = "2 CD207+ LC/DC",
              "3" = "3 Cycling myeloid")
sub_cols <- c("0 Mac/DC" = "#0072B2",
              "1 QC-EXCLUDED (sweat gland)" = "#BDBDBD",
              "2 CD207+ LC/DC" = "#D55E00",
              "3 Cycling myeloid" = "#009E73")

## Canvas is now set to the FINAL print width (180 mm = 7.09 in) rather than
## 14-14.5 in, so every size below is the size that reaches the printed page.
## At the old 14.5 in canvas a 10 pt title was rendered at 4.6 pt once the
## journal scaled the figure to 170 mm; that is below the readable floor.
theme_set(theme_bw(base_size = 8) +
          theme(panel.grid = element_blank(),
                strip.background = element_rect(fill = "grey95", colour = NA),
                strip.text = element_text(size = 6.8, face = "bold"),
                plot.tag = element_text(face = "bold", size = 10),
                plot.tag.position = c(0, 1),
                legend.key.size = unit(0.26, "cm"),
                legend.title = element_text(size = 7),
                legend.text = element_text(size = 6.4),
                plot.title = element_text(size = 8.4, face = "bold"),
                plot.subtitle = element_text(size = 6.4),
                plot.caption = element_text(size = 5.8, colour = "grey30")))

read_table <- function(fname) {
  f <- file.path(TAB_DIR, fname)
  if (!file.exists(f)) stop("MISSING frozen input table: ", f)
  x <- data.table::fread(f, showProgress = FALSE)
  cat("  read ", fname, ": ", nrow(x), " rows x ", ncol(x), " cols\n", sep = "")
  x
}
save_fig <- function(p, stem, w, h) {
  if (!fig_should_render(stem)) return(invisible(NULL))
  ## ggsave() 是随机数的消耗点（geom_jitter / geom_text_repel 在渲染阶段抽样），
  ## 不在这里定种就会让 FIG_ONLY 模式产出与完整跑不同的图。固定种子后，
  ## 每张图的输出与其它图是否渲染无关。
  set.seed(SEED)
  fp <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  fn <- file.path(FIG_DIR, paste0(stem, ".png"))
  ggplot2::ggsave(fp, plot = p, width = w, height = h, units = "in",
                  device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(fn, plot = p, width = w, height = h, units = "in",
                  dpi = FIG_DPI, bg = "white")
  cat("  WROTE ", basename(fp), " (", round(file.size(fp) / 1024, 1), " KiB) / ",
      basename(fn), " (", round(file.size(fn) / 1024, 1), " KiB)\n", sep = "")
  invisible(c(fp, fn))
}
text_panel <- function(lines, tag = NULL, size = 2.75) {
  d <- data.frame(x = 0.01, y = 0.98, label = paste(lines, collapse = "\n"))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, label = label)) +
       ggplot2::geom_text(hjust = 0, vjust = 1, size = size,
                          lineheight = 1.28, family = "mono") +
       ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
       ggplot2::theme_void() +
       ggplot2::theme(plot.tag = element_text(face = "bold", size = 12),
                      panel.background = element_rect(fill = "grey98", colour = "grey85"))
  if (!is.null(tag)) p <- p + ggplot2::labs(tag = tag)
  p
}
## Truncate at a word boundary and mark the cut. The old version chopped at a
## fixed character count and appended one dot, so panel B read
## "positive regulation of infla" and looked like a rendering fault.
short_desc <- function(x, n = 40) {
  vapply(x, function(s) {
    if (nchar(s) <= n) return(s)
    cut <- substr(s, 1, n)
    cut <- sub("[[:space:]][^[:space:]]*$", "", cut)
    if (!nzchar(cut)) cut <- substr(s, 1, n)
    paste0(cut, "...")
  }, character(1), USE.NAMES = FALSE)
}
comma <- function(x) format(x, big.mark = ",", scientific = FALSE)

## =============================================================== FIG 1
cat("\n--- Fig 1: global atlas ---\n")
## The 84,440-cell Seurat object takes 33 s to readRDS and only four things are
## used from it, so cache the projection rather than the object (see figure_cache.R).
s1 <- fig_cached("s1_projection",
                 deps  = file.path(PROC_DIR, F_S1_OBJECT),
                 build = function() {
  o <- readRDS(file.path(PROC_DIR, F_S1_OBJECT))
  mk <- c("KRT1", "KRT10", "KRT14", "CD3D", "CD2", "LYZ", "CD14",
          "COL1A1", "DCN", "PECAM1", "LAMP3", "TPSAB1")
  e  <- Seurat::GetAssayData(o, assay = "RNA", layer = "data")
  g  <- intersect(mk, rownames(e))
  list(n_cells   = ncol(o),
       sample_ids = unique(o$sample_id),
       umap      = as.data.frame(Seurat::Embeddings(o, "umap")),
       lineage_marker = as.character(o$lineage_marker),
       expr_markers = as.matrix(e[g, , drop = FALSE]),
       meta = data.frame(sample_id = o$sample_id, group = o$group,
                         lineage_marker = o$lineage_marker,
                         stringsAsFactors = FALSE))
})
cat("  S1 projection: ", s1$n_cells, " cells / ", length(s1$sample_ids),
    " samples\n", sep = "")
s1sum <- read_table(T_S1_SUMMARY)
stopifnot(sum(s1sum$n_cells) == s1$n_cells)
pilot <- read_table("S1_gateB_trend_pilot.csv")

emb1 <- s1$umap
emb1$lineage <- factor(s1$lineage_marker, levels = names(lineage_cols))
## Shuffle so that no lineage is systematically painted over another.
set.seed(SEED)
emb1 <- emb1[sample.int(nrow(emb1)), , drop = FALSE]
p1a <- ggplot(emb1, aes(umap_1, umap_2, colour = lineage)) +
  geom_point(size = 0.14, alpha = 0.85, stroke = 0) +
  scale_colour_manual(values = lineage_cols, name = "Lineage") +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  coord_fixed() + xlab("UMAP 1") + ylab("UMAP 2") +
  labs(tag = "A", title = "Global skin atlas (84,440 cells, 32 samples)",
       subtitle = "Harmony corrected sample_id only; disease group never a covariate") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

mk_genes <- c("KRT1", "KRT10", "KRT14", "CD3D", "CD2", "LYZ", "CD14",
              "COL1A1", "DCN", "PECAM1", "LAMP3", "TPSAB1")
expr <- s1$expr_markers
g <- intersect(mk_genes, rownames(expr))
M <- expr[g, , drop = FALSE]
lin <- s1$lineage_marker
lev <- names(lineage_cols)
avg <- do.call(rbind, lapply(lev, function(l) rowMeans(M[, lin == l, drop = FALSE])))
pct <- do.call(rbind, lapply(lev, function(l) rowMeans(M[, lin == l, drop = FALSE] > 0)))
dot1 <- data.table::data.table(
  gene     = rep(g, times = length(lev)),
  lineage  = rep(lev, each = length(g)),
  avg_expr = as.vector(t(avg)),
  pct_expr = as.vector(t(pct)))
p1b <- ggplot(dot1, aes(gene, lineage, size = pct_expr, colour = avg_expr)) +
  geom_point() +
  scale_size_continuous(range = c(0.6, 5), name = "% expressing",
                        labels = function(x) paste0(round(100 * x), "%")) +
  scale_colour_gradient(low = "grey85", high = "#D55E00", name = "Mean log-expr") +
  labs(tag = "B", title = "Canonical lineage markers", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 7))

md1 <- data.table::as.data.table(s1$meta)
comp <- md1[, .N, by = .(sample_id, group, lineage_marker)]
comp[, prop := N / sum(N), by = sample_id]
comp[, group := factor(group, levels = c("NS", "PN", "PP"))]
comp[, lineage_marker := factor(lineage_marker, levels = names(lineage_cols))]
comp <- comp[order(group, sample_id, lineage_marker)]
nsamp <- md1[, .(n = uniqueN(sample_id)), by = group]
p1c <- ggplot(comp, aes(sample_id, prop, fill = lineage_marker)) +
  geom_col(width = 0.92, position = "stack") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = lineage_cols, name = "Lineage") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(tag = "C", title = "Per-sample lineage composition",
       subtitle = paste0("n = ", nsamp[group == "NS", n], " / ",
                         nsamp[group == "PN", n], " / ",
                         nsamp[group == "PP", n],
                         " samples; no lineage JT trend (min BH-FDR ",
                         sprintf("%.3f", min(pilot$FDR)), ")"),
       x = "individual samples (one bar each)", y = "Fraction of cells") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

qc_flow <- data.table::data.table(
  step = factor(c("1 Raw (33 libraries)", "2 Cell-level QC",
                  "3 Drop low-depth NS-AR008", "4 Remove doublets"),
                levels = c("4 Remove doublets", "3 Drop low-depth NS-AR008",
                           "2 Cell-level QC", "1 Raw (33 libraries)")),
  cells = c(96088L, 89120L, 88961L, 84440L))
p1d <- ggplot(qc_flow, aes(cells, step)) +
  geom_col(fill = "#0072B2", width = 0.66) +
  geom_text(aes(label = comma(cells)), hjust = -0.12, size = 2.4) +
  scale_x_continuous(limits = c(0, 118000), labels = scales::comma_format()) +
  labs(tag = "D", title = "QC funnel (cells)", x = "Cells", y = "",
       subtitle = paste0("nFeature 200-6000; pct.mt <= ", QC_PCT_MT_MAX, "%")) +
  theme(axis.text.y = element_text(size = 6.8))

f1 <- (p1a | p1b) / (p1c | p1d) +
  plot_layout(heights = c(1, 0.92)) +
  plot_annotation(
    title = "Figure 1. Single-cell atlas of lesional, non-lesional and healthy psoriatic skin (GSE173706)",
    caption = paste0("Data: 84,440 QC-passing cells / 32 samples (NS 7, PN 11, PP 14 donors); marker-based lineage assignment, SingleR HPCA cross-check only.\n",
                     "Composition is shown PER SAMPLE (one bar each), not per cell. Sample-level JT trend: no lineage at BH-FDR < 0.05."))
save_fig(f1, "Fig1_global_atlas", 7.09, 5.5)
rm(s1, expr, M, emb1, md1); invisible(gc())

## =============================================================== FIG 2
cat("\n--- Fig 2: myeloid subcluster and reduced CD207+ signature ---\n")
s2 <- readRDS(file.path(PROC_DIR, F_S2_OBJECT))
emb2 <- as.data.frame(Seurat::Embeddings(s2, "umap"))
emb2$sub_lab <- factor(unname(sub_labs[as.character(s2$subcluster)]),
                       levels = names(sub_cols))
n_excl <- sum(s2$subcluster == "1")
cat("  S2 object: ", ncol(s2), " cells; QC-excluded subcluster 1 = ",
    n_excl, " cells\n", sep = "")
p2a <- ggplot(emb2, aes(umap_1, umap_2, colour = sub_lab)) +
  geom_point(size = 0.85, alpha = 0.85, stroke = 0) +
  scale_colour_manual(values = sub_cols, name = "Subcluster") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  coord_fixed() + xlab("UMAP 1") + ylab("UMAP 2") +
  labs(tag = "A", title = "Myeloid/DC compartment (k = 4)",
       subtitle = "Grey = QC-excluded sweat-gland cells",
       caption = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        # ggplot2 right-aligns captions to the PANEL by default, so a caption
        # wider than this narrow coord_fixed() panel overflowed off the left
        # edge. Anchor it to the whole plot and left-align it instead.
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0, size = 7))

top <- read_table("S2_markers_top10.csv")
sel <- top[, head(.SD[order(-avg_log2FC)], 5), by = cluster]
sel[, logFC := avg_log2FC]; sel[, pct_expr := pct.1]
p2b <- ggplot(sel, aes(cluster, gene, size = pct_expr, colour = logFC)) +
  geom_point() +
  scale_size_continuous(range = c(1, 5), name = "pct.1") +
  scale_colour_gradient2(low = "#0072B2", mid = "grey90", high = "#D55E00",
                         midpoint = 2.5, name = "avg log2FC") +
  labs(tag = "B", title = "Top-5 markers per subcluster",
       subtitle = "Adj. p not used for inference",
       x = "Subcluster", y = "Gene") +
  theme(axis.text.y = element_text(size = 7))

## CORRECTED INPUT (revision): the frozen T_S2_PROP long table omitted pseudo-samples
## with zero cells in a subcluster. Panel C and the panel-D statistics now use the
## complete 29 x 3 grid, so structural zeros appear as points at 0% and n = 7/10/12.
prop <- read_table("S9_prop_grid_zero_filled.csv")
prop[, subcluster := as.character(subcluster)]
prop <- prop[subcluster %in% c("0", "2", "3")]
prop[, group := factor(group, levels = c("NS", "PN", "PP"))]
nn <- prop[, .(n = uniqueN(pseudo_id)), by = .(subcluster, group)]
p2c <- ggplot(prop, aes(group, prop, fill = group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.55, linewidth = 0.35) +
  geom_jitter(width = 0.13, size = 0.85, alpha = 0.8) +
  facet_wrap(~ subcluster, ncol = 1,
             labeller = as_labeller(c("0" = "sub 0 Mac/DC",
                                      "2" = "sub 2 CD207+ LC/DC",
                                      "3" = "sub 3 Cycling"))) +
  scale_fill_manual(values = state_cols, guide = "none") +
  scale_x_discrete(labels = c("NS\n(n=7)", "PN\n(n=10)", "PP\n(n=12)")) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
  coord_cartesian(ylim = c(0, 108)) +
  labs(tag = "C", title = "Proportions of QC-eligible subclusters",
       subtitle = "One point per pseudo-sample; zeros kept",
       x = "", y = "% of QC-eligible myeloid cells per donor")

trend <- read_table("S9_prop_trend_zero_filled.csv")
lmer <- read_table("S9_prop_lmer_zero_filled.csv")[dataset == "all"]
gc2 <- read_table("S2_gateC_inputs.csv")
med2 <- prop[subcluster == "2", .(med = round(stats::median(prop), 2)), by = group]
m2 <- setNames(med2$med, as.character(med2$group))
## S9 trend table keys on `feature` ("sub_0"/"sub_2"/"sub_3"), not on a subcluster column
tp <- setNames(trend$p, trend$feature)
tf <- setNames(trend$FDR, trend$feature)
l2 <- lmer[as.character(subcluster) == "2"]
## Panel D: maturation axis (S11), replacing the former statistics text block.
## The Gate-C verdict and the CD207+ medians that used to sit here are in the
## Figure 2 legend and in Table 2. A data panel is worth more than a text panel
## at the same width, and this one carries the "remodelling, not loss" argument.
mat <- read_table("S11_lc_maturation_pseudobulk.csv")
mat <- mat[compartment == "myeloid_qc_eligible" &
           gene %in% c("CD207", "CD1A", "LAMP3", "CD274")]
mat[, group := factor(group, levels = c("NS", "PN", "PP"))]
mat[, gene := factor(gene, levels = c("CD207", "CD1A", "LAMP3", "CD274"))]
medd <- mat[, .(med = 100 * stats::median(detection_fraction)), by = .(gene, group)]
p2d <- ggplot(mat, aes(group, 100 * detection_fraction)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.45,
               linewidth = 0.28, colour = "grey40") +
  geom_point(aes(colour = group),
             position = position_jitter(width = 0.13, seed = 1),
             size = 0.85, alpha = 0.85) +
  geom_line(data = medd, aes(group, med, group = gene), colour = "grey30",
            linewidth = 0.35, inherit.aes = FALSE) +
  geom_point(data = medd, aes(group, med), size = 1.1, colour = "grey10",
             inherit.aes = FALSE) +
  facet_wrap(~ gene, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = state_cols, guide = "none") +
  labs(tag = "D",
       title = "Maturation axis",
       subtitle = "CD207/CD1A fall; LAMP3/CD274 rise",
       x = "", y = "% of cells with a non-zero count")

f2 <- (p2a | p2b) / (p2c | p2d) +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(
    title = "Figure 2. Myeloid/DC subclusters and a reduced, remodelled CD207+ Langerhans/DC compartment in lesional skin",
    caption = paste0(
      "Data: GSE173706, 29 donor x state pseudo-samples (32 samples after O8 collapsing). n = donors, never cells; box = median/IQR, whiskers = 1.5x IQR.\n",
      "Panel A: grey = 648 sweat-gland/secretory epithelial cells (28%), QC-excluded (>= 3 of top-10 markers epithelial; DCD 674.7 vs LYZ 13.6).\n",
      "Gate C (corrected 29 x 3 grid): sub 3 JT p = 0.0105, FDR = 0.031, met in that specification only\n",
      "(FDR 0.131 after dropping three pseudo-samples with <10 myeloid cells), so the three-state claim is not established.\n",
      "CD207+ median NS 9.43 / PN 36.86 / PP 6.81% (non-monotone); count-model rate ratio 0.23 per state step (95% CI 0.17-0.30)."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 5.6)))
save_fig(f2, "Fig2_myeloid_depletion", 7.09, 5.85)
rm(s2, emb2); invisible(gc())

## =============================================================== FIG 3
cat("\n--- Fig 3: two-group function ---\n")
deg <- read_table("S3_two_group_DEG.csv")
deg <- deg[!is.na(pvalue) & !is.na(padj)]
deg[, neglogp := -log10(pvalue)]
deg <- deg[is.finite(neglogp) & is.finite(log2FoldChange)]
deg[, sig := padj < PB_FDR]
nsig <- sum(deg$sig)
lab <- deg[sig == TRUE][order(padj)][seq_len(min(10, .N))]
p3a <- ggplot(deg, aes(log2FoldChange, neglogp)) +
  geom_point(aes(colour = sig), size = 0.55, alpha = 0.55, stroke = 0) +
  scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "#D55E00"),
                      labels = c("FALSE" = "NS", "TRUE" = "FDR < 0.05"),
                      name = "") +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.35,
                           max.overlaps = 12, min.segment.length = 0,
                           segment.size = 0.2, force = 2, box.padding = 0.28) +
  labs(tag = "A", title = "DEG, subcluster 0",
       subtitle = paste0(nsig, " / ", unique(deg$universe_size),
                         " genes at FDR < 0.05; ",
                         unique(deg$n_NS), " vs ", unique(deg$n_PP),
                         " pseudo-samples"),
       x = "log2 fold change (PP vs NS)", y = "-log10 p-value") +
  ## The nine gene labels sit in the upper-right quadrant, where the old
  ## in-panel legend (0.88, 0.88) left S100A8 overlapping "FDR < 0.05".
  ## Top-left of this volcano is empty, so the key moves there.
  theme(legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key = element_blank())

pw <- read_table(T_S3_PATHWAY)
pw <- pw[gene_list == "twoGroup_sub0" & method == "enrichGO_BP" & significant == TRUE]
n_pw_sig <- nrow(pw)
pw <- pw[order(p.adjust)][seq_len(min(15, n_pw_sig))]
pw[, desc := short_desc(Description, 29)]
pw[, neglogp := -log10(p.adjust)]
p3b <- ggplot(pw, aes(FoldEnrichment, reorder(desc, p.adjust),
                      size = Count, colour = neglogp)) +
  geom_point() +
  scale_size_continuous(range = c(1.5, 5), name = "Gene count") +
  scale_colour_gradient(low = "#0072B2", high = "#D55E00", name = "-log10 FDR") +
  labs(tag = "B", title = "GO enrichment",
       subtitle = paste0("Top 15 of ", n_pw_sig, " terms, BH-FDR < 0.05"),
       x = "Fold enrichment", y = "") +
  theme(axis.text.y = element_text(size = 6.0))

gs <- read_table("S3_gsea.csv")
gs <- gs[gene_list == "twoGroup_sub0"]
gs_sig <- sum(gs$padj < PATHWAY_FDR, na.rm = TRUE)
gsel <- gs[order(padj)][seq_len(min(20, nrow(gs)))]
gsel[, set := sub("^HALLMARK_", "", pathway)]
## The 20 set names are drawn as y-axis text. In MSigDB form
## (EPITHELIAL_MESENCHYMAL_TRANSITION) they need ~150 pt of width, which
## patchwork then forces on every panel in that column - it flattened panel A
## from 217 pt to 90 pt once the figure was set to 180 mm. Printed as readable
## labels they need ~110 pt, which also stops "Il6"/"Nfkb" appearing if the
## underscore is simply swapped for a space and the case is lowered blindly.
HALLMARK_LABELS <- c(
  COMPLEMENT                     = "Complement",
  EPITHELIAL_MESENCHYMAL_TRANSITION = "Epithelial-mesenchymal transition",
  IL6_JAK_STAT3_SIGNALING        = "IL-6/JAK/STAT3 signaling",
  INFLAMMATORY_RESPONSE          = "Inflammatory response",
  INTERFERON_ALPHA_RESPONSE      = "IFN-alpha response",
  INTERFERON_GAMMA_RESPONSE      = "IFN-gamma response",
  ADIPOGENESIS                   = "Adipogenesis",
  ALLOGRAFT_REJECTION            = "Allograft rejection",
  HYPOXIA                        = "Hypoxia",
  KRAS_SIGNALING_DN              = "KRAS signaling down",
  IL2_STAT5_SIGNALING            = "IL-2/STAT5 signaling",
  ESTROGEN_RESPONSE_LATE         = "Estrogen response late",
  SPERMATOGENESIS                = "Spermatogenesis",
  COAGULATION                    = "Coagulation",
  OXIDATIVE_PHOSPHORYLATION      = "Oxidative phosphorylation",
  KRAS_SIGNALING_UP              = "KRAS signaling up",
  ESTROGEN_RESPONSE_EARLY        = "Estrogen response early",
  TNFA_SIGNALING_VIA_NFKB        = "TNF-alpha signaling via NF-kB",
  APICAL_JUNCTION                = "Apical junction",
  P53_PATHWAY                    = "p53 pathway")
hallmark_label <- function(x) {
  hit <- unname(HALLMARK_LABELS[x])
  miss <- is.na(hit)
  if (any(miss)) {                       # unknown set: fall back, never drop it
    s <- tolower(gsub("_", " ", x[miss]))
    hit[miss] <- paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
  }
  hit
}
p3c <- ggplot(gsel, aes(NES, reorder(set, NES), fill = padj < PATHWAY_FDR)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#D55E00"),
                    labels = c("FALSE" = "FDR >= 0.05", "TRUE" = "FDR < 0.05"),
                    name = "") +
  scale_y_discrete(labels = hallmark_label) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(tag = "C", title = "GSEA Hallmark",
       subtitle = paste0(gs_sig, " / 50 sets at BH-FDR < 0.05"),
       x = "Normalised enrichment score (NES)", y = "") +
  theme(axis.text.y = element_text(size = 6.8),
        legend.position = c(0.82, 0.15))

gsva <- read_table("S3_gsva.csv")
## Subcluster 2 has only 1 NS vs 3 PP pseudo-samples: Wilcoxon p and BH-FDR
## are NA for all 50 sets. That is NOT a negative result - it is a test that
## could not be computed. Keep the facet and label it explicitly.
gsva[, neglogFDR := ifelse(is.na(FDR_BH), NA_real_, -log10(FDR_BH))]
gsva[, sub_lab := paste0("subcluster ", subcluster)]
n0_sig <- gsva[subcluster == 0 & !is.na(FDR_BH) & FDR_BH < PATHWAY_FDR, .N]
n0_all <- gsva[subcluster == 0 & !is.na(FDR_BH), .N]
n2_na  <- gsva[subcluster == 2 & is.na(FDR_BH), .N]
ann <- data.table::data.table(
  sub_lab = c("subcluster 0", "subcluster 2"),
  txt = c(sprintf("%d/%d FDR<0.05; min attainable 0.51", n0_sig, n0_all),
          sprintf("%d/50 computable (1 vs 3)", n2_na)))
p3d <- ggplot(gsva, aes(delta_PP_minus_NS, neglogFDR)) +
  geom_point(data = gsva[!is.na(neglogFDR)], size = 1.2, alpha = 0.75,
             colour = "#0072B2") +
  geom_hline(yintercept = -log10(PATHWAY_FDR), linetype = "dashed",
             colour = "#D55E00", linewidth = 0.4) +
  facet_wrap(~ sub_lab) +
  labs(tag = "D", title = "GSVA: underpowered",
       subtitle = "sub 0: 0/50 (min attainable FDR 0.51); sub 2: not computable",
       x = "median score difference (PP - NS)", y = "-log10 BH-FDR")

## Column widths must be set INSIDE each row: plot_layout(widths=) applied to
## the stacked (a|b)/(c|d) object only has one column to distribute, so it was
## silently doing nothing and the two rows kept the 1 : 1 default.
f3 <- ((p3a | p3b) + plot_layout(widths = c(1.5, 1))) /
      ((p3c | p3d) + plot_layout(widths = c(1.5, 1))) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 3. Two-group transcriptional and pathway-level evidence (exploratory)",
    caption = paste0(
      "Subclusters 2/3 two-group DEG SKIPPED (1 NS vs 3 PP; 0 NS vs 4 PP after PB_MIN_CELLS = 10); no imputation. All tests use sample-level pseudobulk.\n",
      "CellChat / trajectory NOT run (offline branch B, O5); BayesPrism / CIBERSORTx replaced by ssGSEA (O6). Full GO term names are in the supplementary tables.\n",
      "Panel D: an UNINFORMATIVE GSVA test (sub 0: 0/50, but FDR<0.05 was unattainable at 3 vs 11) beside a test that was NOT COMPUTABLE\n",
      "(sub 2: 1 vs 3). Neither is a negative result."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 5.6)))
save_fig(f3, "Fig3_function", 7.09, 5.7)

## =============================================================== FIG 4
cat("\n--- Fig 4: cross-cohort validation ---\n")
auc <- read_table(T_S4_AUC)
auc_per <- auc[statistic == "per_subcluster"]
mac <- auc[statistic == "gateD_primary", auc_mean]
p4a <- ggplot(auc_per, aes(paste0("sub ", subcluster), auc_mean)) +
  geom_errorbar(aes(ymin = auc_min, ymax = auc_max), width = 0.12,
                colour = "grey40") +
  geom_point(size = 3, colour = "#0072B2") +
  geom_hline(yintercept = mac, linetype = "dashed", colour = "#D55E00",
             linewidth = 0.45) +
  geom_hline(yintercept = GATE_D_MIN_AUC, linetype = "dotted",
             colour = "grey45", linewidth = 0.4) +
  annotate("text", x = 3.45, y = mac - 0.055, label = paste0("macro ", sprintf("%.4f", mac)),
           size = 2.5, colour = "#D55E00", hjust = 1) +
  annotate("text", x = 3.45, y = GATE_D_MIN_AUC + 0.025, label = "0.70 threshold",
           size = 2.5, colour = "grey35", hjust = 1) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
  labs(tag = "A", title = "INTERNAL 5-fold CV AUC (NOT external validation)",
       subtitle = "Label-transfer reproducibility WITHIN the discovery reference only",
       x = "Reference subcluster", y = "AUC")

dm <- read_table("S4_query_direct_markers.csv")
dm[, group := factor(group, levels = c("NS", "PP"))]
p4b <- ggplot(dm, aes(group, pct_CD207_all, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.55, linewidth = 0.35) +
  geom_jitter(width = 0.11, size = 2.2, show.legend = FALSE) +
  scale_fill_manual(values = state_cols, guide = "none") +
  labs(tag = "B", title = "Mapping-free CD207+ cells (GSE162183)",
       subtitle = "3 vs 3; median 0.3376 vs 0.2166%; descriptive, no p",
       x = "", y = "CD207+ cells (% of all cells)")

sc <- read_table(T_S4_SCORES)
lc <- sc[gene_set == "LC_CD207_signature"]
lc[, group := factor(group, levels = c("NN", "PN", "PP"))]
bt <- read_table("S4_bulk_group_tests.csv")
bt <- bt[gene_set == "LC_CD207_signature"]
pt <- read_table("S4_bulk_paired_PN_PP.csv")
pt <- pt[gene_set == "LC_CD207_signature"]
p4c <- ggplot(lc, aes(group, score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.6, linewidth = 0.35) +
  geom_jitter(width = 0.18, size = 0.8, alpha = 0.45, show.legend = FALSE) +
  scale_fill_manual(values = bulk_cols, guide = "none") +
  labs(tag = "C", title = "Bulk LC signature score (GSE13355)",
       subtitle = sprintf("NN %d / PN %d / PP %d samples; PP-NN delta %.4f, p %.2e, BH FDR %.2e",
                          sum(lc$group == "NN"), sum(lc$group == "PN"),
                          sum(lc$group == "PP"),
                          bt$delta_PP_minus_NN, bt$p_wilcox, bt$FDR_BH),
       x = "", y = "ssGSEA score (LC_CD207 signature)")

lc[, group_ch := as.character(group)]
pair <- data.table::dcast(lc[group_ch %in% c("PN", "PP")],
                          donor_id ~ group_ch, value.var = "score")
pair <- pair[!is.na(PN) & !is.na(PP)]
p4d <- ggplot(pair) +
  geom_segment(aes(x = 1, xend = 2, y = PN, yend = PP),
               colour = "grey55", alpha = 0.45, linewidth = 0.3) +
  geom_point(aes(x = 1, y = PN), colour = bulk_cols["PN"], size = 1.3, alpha = 0.75) +
  geom_point(aes(x = 2, y = PP), colour = bulk_cols["PP"], size = 1.3, alpha = 0.75) +
  scale_x_continuous(limits = c(0.85, 2.15), breaks = c(1, 2),
                     labels = c("PN", "PP")) +
  labs(tag = "D", title = "Paired PN -> PP (same donor)",
       subtitle = sprintf("%d pairs, signed-rank p = %.2e (BH FDR %.2e)",
                          nrow(pair), pt$p_signed_rank, pt$FDR_BH),
       x = "", y = "ssGSEA score (LC_CD207 signature)")

lim_a <- c(
  "MANDATORY LIMITATIONS (1/2)",
  "",
  "1. AUC 0.9988 is INTERNAL 5-fold CV within the",
  "   discovery reference: label-transfer reproducibility,",
  "   NOT independent-cohort validation. Do not advertise it.",
  "2. Whole-skin GSE162183 query vs myeloid-only reference",
  "   = DOMAIN MISMATCH: ~98% called macrophage/DC,",
  "   score.max ~1; 0.3/0.5/0.7 filters kept 99.9/99.9/97.7%.",
  "   Transferred proportions are NOT quantitative.",
  "3. Mapping-free CD207+ counts: 3 vs 3 samples; one PP",
  "   sample exceeds two NS samples. Supportive only,",
  "   NOT confirmatory; no significance claim is made.")
lim_b <- c(
  "MANDATORY LIMITATIONS (2/2)",
  "",
  "4. GSE13355 = batch/sex-adjusted log2 RELATIVE expression;",
  "   2^x is a monotone rescaling, not absolute recovery.",
  "   ssGSEA is rank based; the adjustment is irreversible.",
  "5. GSE162183 has none of the 15 MT whitelist symbols:",
  "   no pct.mt filter on that query (22,716 / 24,234 kept).",
  "6. Bulk is the strongest evidence: LC signature",
  "   NN 0.2308 / PN 0.1844 / PP -0.0394. Cycling and",
  "   Mac_DC controls move the OPPOSITE way, so this is not",
  "   a global shift. Gate D PASSES on the frozen criterion.")
## Panels E/F (text limitations) were moved to the Figure 4 legend. At the final
## 180 mm width two text panels render below 5 pt, so the six limitations are now
## carried by the legend and by the caption below; no content was dropped.
f4 <- (p4a | p4b) / (p4c | p4d) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 4. Cross-cohort and bulk evidence for a reduced CD207+ Langerhans/DC signature",
    caption = paste0(
      "Cohorts: GSE162183 (6 whole-skin scRNA samples, 22,716 QC cells) and GSE13355 (180 arrays / 122 subjects; NN 64, PN 58, PP 58).\n",
      "Box = median/IQR, whiskers = 1.5x IQR; panel D connects the SAME donor. Tests: internal 5-fold CV inside the discovery reference (NOT external validation);\n",
      "mapping-free counts with no label transfer; Wilcoxon rank-sum + BH across groups; paired signed-rank + BH for PN-PP.\n",
      "MANDATORY LIMITATIONS (full text in the Figure 4 legend): AUC is internal CV, not external validation; whole-skin query vs myeloid-only reference is a domain mismatch;\n",
      "3-vs-3 counts are descriptive; GSE13355 is batch/sex-adjusted log2 relative expression; GSE162183 lacks the mitochondrial whitelist, so no percent.mt filter;\n",
      "bulk is the strongest evidence and Gate D passes on the frozen criterion."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 5.6)))
save_fig(f4, "Fig4_validation", 7.09, 5.65)

## =============================================================== close
cat("\n=== S5 DONE ===\n")
cat("figures in ", FIG_DIR, ": ",
    paste(list.files(FIG_DIR), collapse = ", "), "\n", sep = "")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
    "\n", sep = "")
invisible(NULL)
}

## ---------------------------------------------------------------- run + log
logf <- file.path(LOG_DIR, "S5_figures.log")
con <- file(logf, open = "wt")
sink(con, split = TRUE)
ok <- tryCatch({ main(); TRUE },
               error = function(e) {
                 cat("\nERROR in S5: ", conditionMessage(e), "\n", sep = "")
                 cat("traceback:\n"); print(traceback())
                 FALSE
               })
sink()
close(con)
writeLines(capture.output(sessionInfo()),
           file.path(LOG_DIR, "S5_sessionInfo.txt"))
if (!ok) quit(save = "no", status = 1)
cat("S5 complete. Log: ", logf, "\n", sep = "")
