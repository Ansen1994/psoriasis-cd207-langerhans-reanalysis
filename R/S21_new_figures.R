## S21_new_figures.R - rendering only for the post hoc supplementary layers
## (S17/S18/S19/S20). Performs NO test, changes NO parameter; every limitation
## a reader needs is printed on the figure itself, as in S13b.
##
## Outputs
##   results/figures/Supplementary_Figure_S7.pdf/.png
##       LC-compartment (subcluster-2) ordered-trend function: fgsea Hallmark
##       NES bars + the GO over-representation null result.
##   results/figures/Supplementary_Figure_S8.pdf/.png
##       NNLS deconvolution (GSE13355) and the GSE226244 treatment layer.
##   results/figures/Graphical_abstract.pdf/.png
##
## ASCII only (LC_CTYPE=C). All outputs stay inside the project.
source("config/00_paths.R")
if (!exists("FIG_DPI") && file.exists("config/01_params.R")) {
  source("config/01_params.R")
}
if (!exists("FIG_DPI")) FIG_DPI <- 600
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

theme_set(theme_bw(base_size = 7) +
            theme(plot.caption = element_text(size = 6.0, colour = "grey35",
                                              hjust = 0)))
th <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.title = element_text(face = "bold", size = 8, hjust = 0),
        plot.subtitle = element_text(size = 6.5, hjust = 0),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7.2),
        legend.text = element_text(size = 7.2),
        plot.caption = element_text(size = 6.2, hjust = 0, colour = "grey35"),
        legend.position = "none")

## ================================================================ Figure S7
g <- fread(file.path(TAB_DIR, "S3_gsea.csv"))
s2 <- g[gene_list == "orderedTrend_sub2"][padj < 0.05][order(NES)]
s2[, pathway := gsub("^HALLMARK_", "", pathway)]
s2[, pathway := gsub("_", " ", pathway)]
s2[, sig_lab := sprintf("BH-FDR %.3g", padj)]
s2[, dir := ifelse(NES > 0, "up", "down")]
pA <- ggplot(s2, aes(x = NES, y = reorder(pathway, NES), fill = dir)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.2) +
  geom_text(data = s2[dir == "up"], aes(x = NES + 0.04, label = sig_lab),
            hjust = 0, size = 2.2, colour = "grey20") +
  geom_text(data = s2[dir == "down"], aes(x = NES - 0.04, label = sig_lab),
            hjust = 0, size = 2.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c("up" = "#D55E00", "down" = "#0072B2"),
                    labels = c("up" = "enriched along NS-PN-PP",
                               "down" = "depleted along NS-PN-PP")) +
  scale_x_continuous(limits = c(-2.6, 3.2)) +
  labs(title = "A | Hallmark sets ranked by the ordered-trend statistic",
       subtitle = paste0("fgsea on all 8,061 subcluster-2 trend-ranked genes",
                         " (log2FC-signed statistic), 50 sets, BH across 50"),
       x = "Normalised enrichment score (NES)",
       y = NULL, fill = NULL, caption = paste0(
         "Post hoc, exploratory. Ranking statistic = signed log2FC of the DESeq2 ordered-trend LRT",
         " (positive = increases NS->PN->PP) on the CD207+ subcluster (sub2);",
         " frozen seed, nperm = 10,000, minSize = 10, maxSize = 500.",
         " Over-representation control: enrichGO BP on the 58 trend genes (FDR < 0.05)",
         " returned 0 significant terms of 134 tested - the ranked test is the",
         " better-powered readout and the GO null is reported deliberately.")) +
  th + theme(legend.position = "top", legend.key.size = unit(2.5, "mm"))
figS7 <- pA + plot_annotation(
  title = "Supplementary Figure S7. Functional shift within the CD207+ Langerhans/DC subcluster (post hoc, exploratory)")
ggsave(file.path(FIG_DIR, "Supplementary_Figure_S7.pdf"), figS7,
       width = 150, height = 90, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Supplementary_Figure_S7.png"), figS7,
       width = 150, height = 90, units = "mm", dpi = FIG_DPI)

## ================================================================ Figure S8
fr <- fread(file.path(TAB_DIR, "S18_deconvolution_fractions.csv"))
t13355 <- fr[cohort == "GSE13355"]
g13355 <- ggplot(t13355[!is.na(share_sub2)],
                 aes(group, share_sub2, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.25,
               colour = "grey25") +
  geom_jitter(width = 0.16, size = 0.5, alpha = 0.55, colour = "grey30") +
  scale_fill_manual(values = c("NN" = "#009E73", "PN" = "#E69F00",
                               "PP" = "#D55E00")) +
  stat_summary(fun = median, geom = "crossbar", width = 0.55, linewidth = 0.2) +
  labs(title = "A | Deconvolved CD207+ share, GSE13355",
       subtitle = sprintf("NN %.3f | PN %.3f | PP %.3f (median share of the four myeloid signatures)",
                          median(t13355$share_sub2[t13355$group == "NN"], na.rm = TRUE),
                          median(t13355$share_sub2[t13355$group == "PN"], na.rm = TRUE),
                          median(t13355$share_sub2[t13355$group == "PP"], na.rm = TRUE)),
       x = NULL, y = "Estimated sub2 share (sum = 1)")
g13355 <- g13355 + th

dif <- t13355[, .(donor_id, group, share_sub2)][group %in% c("PN", "PP")]
w <- dcast(dif, donor_id ~ group, value.var = "share_sub2")
w <- w[!is.na(PN) & !is.na(PP)][, .(donor_id, d = PP - PN)]
gB <- ggplot(w, aes(x = d)) +
  geom_histogram(bins = 24, fill = "#0072B2", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = median(w$d), linetype = 2, linewidth = 0.3, colour = "grey20") +
  labs(title = "B | Within-subject change (58 PN-PP pairs)",
       subtitle = sprintf("PP minus PN; median %.3f | signed-rank p = 3.6e-11", median(w$d)),
       x = "Paired difference (PP - PN)", y = "Donors",
       caption = "NNLS vs the four frozen subcluster mean signatures; 2^x rescaled log2-relative input. Post hoc.") + th

tr <- fread(file.path(TAB_DIR, "S19_treatment_bulk_samples.csv"))
tr[, group := factor(group, levels = c("NN", "PRE", "POST"))]
tr_col <- c("NN" = "#009E73", "PRE" = "#CC79A7", "POST" = "#56B4E9")

pC <- ggplot(tr[!is.na(expr_CD207)], aes(group, expr_CD207, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.25, colour = "grey25") +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.55, colour = "grey30") +
  scale_fill_manual(values = tr_col) +
  labs(title = "C | Single-gene CD207",
       subtitle = "POST vs PRE Wilcoxon FDR 4.8e-05",
       x = NULL, y = "CD207 (2^x of log2 rel.)") + th

pD <- ggplot(tr[!is.na(LC_CD207_signature)], aes(group, LC_CD207_signature, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.25, colour = "grey25") +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.55, colour = "grey30") +
  scale_fill_manual(values = tr_col) +
  labs(title = "D | Frozen LC/CD207 9-gene score",
       subtitle = "PRE 0.589 -> POST 0.709; healthy 0.718",
       x = NULL, y = "ssGSEA score") + th

pE <- ggplot(tr[!is.na(share_sub2)], aes(group, share_sub2, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.25, colour = "grey25") +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.55, colour = "grey30") +
  scale_fill_manual(values = tr_col) +
  labs(title = "E | Deconvolved CD207+ share",
       subtitle = "PRE 0.486 -> POST 0.541; healthy 0.553",
       x = NULL, y = "Estimated sub2 share") + th

cap <- paste0(
  "Panels C-E: GSE226244, GPL570, same platform as GSE13355 (Kim & Krueger, Front Immunol 2023). ",
  "PRE = 33 pretreatment lesional arrays; POST = 28 post-treatment arrays (18 wk12, 9 wk24, 1 wk20 off-drug); ",
  "NN = 8 controls incl. 5 re-analyzed GSE78097 arrays (batch caution). ",
  "PRE vs POST is UNPAIRED (no patient identifier is deposited) and the POST arm mixes time points.\n",
  "Panels A-B: GSE13355 deconvolution. The sub0 (macrophage/DC) signature received zero NNLS weight in every cohort\n",
  "(collinearity; disclosed), so the share is read within the competing-subcluster space, not as an absolute tissue fraction.\n",
  "Direction replicated in GSE13355 only; the two S8 GPL570 cohorts did NOT replicate the share change\n",
  "(paired p = 0.66 / 0.16) and their batch-adjusted log2-relative input is a known limitation (Section 2.12).\n",
  "Single-cell treatment layer (GSE183047): descriptive only, n = 4 paired patients, no test (Supplementary Table).")
figS8 <- (g13355 | gB) / (pC | pD | pE) +
  plot_annotation(
    title = "Supplementary Figure S8. Deconvolution quantification and the treatment layer (post hoc)",
    caption = cap, theme = theme(plot.caption = element_text(size = 6.0, colour = "grey35", hjust = 0)))
ggsave(file.path(FIG_DIR, "Supplementary_Figure_S8.pdf"), figS8,
       width = 190, height = 170, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Supplementary_Figure_S8.png"), figS8,
       width = 190, height = 170, units = "mm", dpi = FIG_DPI)

## ================================================== graphical abstract (P3)
mkbox <- function(lab, x, y, fill = "grey95") {
  annotate("rect", xmin = x - 0.53, xmax = x + 0.53, ymin = y - 0.16, ymax = y + 0.16,
           fill = fill, colour = "grey30", linewidth = 0.3)
}
mktext <- function(lab, x, y, size = 2.6, fontface = "plain", colour = "grey10") {
  annotate("text", x = x, y = y, label = lab, size = size, fontface = fontface,
           colour = colour, lineheight = 0.95)
}
ga <- ggplot(data.frame(x = 0:1, y = 0:1), aes(x, y)) +
  xlim(0, 4) + ylim(-0.06, 1) + theme_void() +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA)) +
  ## header
  mktext("Reduced CD207+ Langerhans-cell signature in lesional psoriasis",
         2, 0.93, size = 3.6, fontface = "bold") +
  mktext("A pre-specified, fully disclosed re-analysis of five public cohorts",
         2, 0.845, size = 2.7, colour = "grey35") +
  ## row 1: the question
  mkbox("", 0.85, 0.62, "#E69F00") +
  mktext("Three-state question\nhealthy -> non-lesional -> lesional\n(sup2, 29 x 3 grid)",
         0.85, 0.62, size = 2.3) +
  mktext("ordered trend NOT established\n(specification-dependent;\nconservative two-group framing retained)",
         0.85, 0.375, size = 2.0, colour = "grey25") +
  annotate("segment", x = 0.85, y = 0.46, xend = 0.85, yend = 0.53,
           arrow = arrow(length = unit(1.6, "mm")), linewidth = 0.3, colour = "grey30") +
  ## row 1 right: main result
  mkbox("", 2.6, 0.62, "#D55E00") +
  mktext("Two-group result\nCD207+ signature REDUCED in lesional skin\n(single-cell proportions + 3 bulk cohorts)",
         2.6, 0.62, size = 2.2) +
  mktext("single-gene CD207 primary readout\nGSE13355 d = -1.55, p = 5.1e-17\n58 within-subject pairs p = 2.8e-10",
         2.6, 0.375, size = 2.0, colour = "grey25") +
  annotate("segment", x = 2.6, y = 0.46, xend = 2.6, yend = 0.53,
           arrow = arrow(length = unit(1.6, "mm")), linewidth = 0.3, colour = "grey30") +
  annotate("segment", x = 1.43, y = 0.62, xend = 2.02, yend = 0.62,
           arrow = arrow(length = unit(1.6, "mm")), linewidth = 0.3, colour = "grey30") +
  mktext("also supported by", 1.72, 0.665, size = 1.9, colour = "grey35") +
  ## row 2: new quantification
  mkbox("", 0.85, 0.14, "#0072B2") +
  mktext("Deconvolution (GSE13355)\nCD207+ share of the myeloid-signature space\n26% -> 4% (PP vs NN FDR 2.3e-20)",
         0.85, 0.14, size = 2.0) +
  mkbox("", 2.6, 0.14, "#009E73") +
  mktext("Treatment (GSE226244, secukinumab)\nCD207 signature recovers toward healthy\nafter effective IL-17A blockade (FDR 4.8e-05)",
         2.6, 0.14, size = 2.0) +
  mktext("negatives, non-computable and contamination-driven results are reported in full; code and tables: Zenodo/GitHub (v1.1.0)",
         2, -0.035, size = 1.85, colour = "grey35")
ggsave(file.path(FIG_DIR, "Graphical_abstract.pdf"), ga,
       width = 180, height = 110, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Graphical_abstract.png"), ga,
       width = 180, height = 110, units = "mm", dpi = FIG_DPI)
cat("figures written.\n")
