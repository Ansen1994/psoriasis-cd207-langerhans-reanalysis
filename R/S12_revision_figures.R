## S12_revision_figures.R
## Revision figures added 2026-11-13 (post hoc; none of these enter a frozen gate).
##   Supplementary Figure S3  composition controls + empirical null for the LC score
##   Supplementary Figure S4  maturation axis within the myeloid compartment
##   Supplementary Figure S5  cluster-level UMAP (the "19 clusters" quoted in Results 3.1)
## Colour palette is Okabe-Ito (colour-blind safe) throughout.
##
## Outputs: results/figures/Supplementary_Figure_S{3,4,5}.pdf and .png
source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(Seurat)
})

logf <- file.path(LOG_DIR, "S12_revision_figures.log")
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
cat("=== S12 revision figures ===\n", file = logf)

ok <- c(blue = "#0072B2", vermillion = "#D55E00", bluishgreen = "#009E73",
        purple = "#CC79A7", orange = "#E69F00", sky = "#56B4E9", grey = "#8C8C8C")
grp_cols <- c(NS = ok[["sky"]], PN = ok[["orange"]], PP = ok[["vermillion"]])

TH <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.6, colour = "grey30"),
        plot.caption = element_text(size = 6.6, colour = "grey35", hjust = 0),
        strip.background = element_rect(fill = "grey94", colour = "grey60"),
        legend.key.size = unit(9, "pt"))

save_fig <- function(p, stem, w, h) {
  pdfp <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  pngp <- file.path(FIG_DIR, paste0(stem, ".png"))
  ggsave(pdfp, p, width = w, height = h, units = "in", device = "pdf")
  ggsave(pngp, p, width = w, height = h, units = "in", dpi = FIG_DPI)
  say("WROTE ", basename(pdfp), " (", w, "x", h, " in) and PNG ", FIG_DPI, " dpi")
}

## ============================================================ S3: composition
s10 <- fread(file.path(TAB_DIR, "S10_signature_confounding.csv"))
nul <- fread(file.path(TAB_DIR, "S10_null_raw.csv"))
adj <- fread(file.path(TAB_DIR, "S10_null_distribution.csv"))
say("S10 rows: ", nrow(s10), " | null sets: ", nrow(nul))

set_lab <- c(LC9_platform = "LC score (nine genes)",
             LC10_original = "LC score (ten members)",
             LC7_no_eccrine = "LC score minus\neccrine genes",
             LC4_minimal = "LC score, four\nLC-restricted members",
             Eccrine_control = "Eccrine / secretory",
             Eccrine_noSLC = "Eccrine, no SLC18A2",
             Mast_control = "Mast cell",
             KRT_differentiation = "Keratinocyte\ndifferentiation",
             Melanocyte_control = "Melanocyte")
s10[, lab := set_lab[gene_set]]
s10[, lab := factor(lab, levels = set_lab[c("KRT_differentiation", "Melanocyte_control",
                                           "Eccrine_control", "Mast_control",
                                           "LC9_platform")])]

pa_dt <- s10[!is.na(lab)]
## values are written into the tick labels so that no label can be clipped
pa_dt[, lab_v := sprintf("%s  (%+.2f)", lab, delta_PP_minus_NN)]
pA <- ggplot(pa_dt, aes(delta_PP_minus_NN, lab_v)) +
  geom_col(width = 0.62, fill = ok[["blue"]]) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(tag = "A", x = "Median difference PP - NN (ssGSEA units)", y = "",
       title = "The LC score is not the only compartment that falls",
       subtitle = "Identical ssGSEA settings, GSE13355 (64 NN / 58 PP); all five sets BH < 0.05") +
  TH

obs_lc <- s10[gene_set == "LC9_platform", delta_PP_minus_NN]
obs_ec <- s10[gene_set == "Eccrine_control", delta_PP_minus_NN]
pB <- ggplot(nul, aes(delta_PP_minus_NN)) +
  geom_histogram(bins = 26, fill = ok[["grey"]], colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs_lc, colour = ok[["vermillion"]], linewidth = 0.8) +
  geom_vline(xintercept = obs_ec, colour = ok[["blue"]], linewidth = 0.8, linetype = "22") +
  annotate("text", x = obs_lc, y = 46, label = "LC score ", hjust = 1.05, vjust = 0,
           size = 2.7, colour = ok[["vermillion"]]) +
  annotate("text", x = obs_ec, y = 40, label = " eccrine", hjust = -0.05, vjust = 0,
           size = 2.7, colour = ok[["blue"]]) +
  labs(tag = "B", x = "Median difference PP - NN", y = "Random nine-gene sets",
       title = "Empirical null: 200 expression-matched random sets",
       subtitle = sprintf("Median %.3f (IQR %.3f to %.3f); observed LC %.3f at the %.0fth percentile",
                          median(nul$delta_PP_minus_NN),
                          quantile(nul$delta_PP_minus_NN, 0.25),
                          quantile(nul$delta_PP_minus_NN, 0.75),
                          obs_lc, 100 * mean(nul$delta_PP_minus_NN <= obs_lc))) +
  TH

red <- s10[gene_set %in% c("LC9_platform", "LC7_no_eccrine", "LC4_minimal")]
red[, lab := factor(c("all nine genes", "minus eccrine-
shared genes",
                      "four LC-restricted
members")[match(gene_set,
                      c("LC9_platform", "LC7_no_eccrine", "LC4_minimal"))],
                    levels = c("all nine genes", "minus eccrine-\nshared genes",
                               "four LC-restricted\nmembers"))]
pD <- ggplot(red, aes(lab, delta_PP_minus_NN, fill = p_paired < 0.05)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = sprintf("paired p = %.3g", p_paired)), y = 0.012, vjust = 0, size = 2.5) +
  scale_fill_manual(values = c(`TRUE` = ok[["blue"]], `FALSE` = ok[["grey"]]),
                    labels = c(`TRUE` = "paired p < 0.05", `FALSE` = "paired n.s."), name = NULL) +
  scale_y_continuous(limits = c(-0.32, 0.05)) +
  labs(tag = "D", x = "", y = "Median difference PP - NN",
       title = "Restricted to four defensible LC members the paired test fails",
       subtitle = "58 matched PN-PP pairs; this is why single-gene CD207 is the primary readout") +
  TH + theme(legend.position = "top")

## ============================================================ S4: maturation
pb <- fread(file.path(TAB_DIR, "S11_lc_maturation_pseudobulk.csv"))
pb[, group := factor(group, levels = c("NS", "PN", "PP"))]
gsel <- c("CD207", "CD1A", "LAMP3", "CCR7", "CD274")

m1 <- pb[compartment == "myeloid_qc_eligible" & gene %in% gsel]
p4a <- ggplot(m1, aes(group, detection_fraction, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.55, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = SEED),
             size = 0.85, alpha = 0.75, stroke = 0) +
  facet_wrap(~ factor(gene, levels = gsel), nrow = 1) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(tag = "A", x = "", y = "Detection fraction",
       title = "Fewer CD207+/CD1A+ myeloid cells, more LAMP3+/CCR7+/CD274+ cells",
       subtitle = "All 1,666 QC-eligible myeloid cells, donor x state pseudo-samples; box = median/IQR") +
  TH

m2 <- pb[compartment == "myeloid_qc_eligible" & gene %in% gsel]
p4b <- ggplot(m2, aes(group, pseudobulk_log2cpm, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.55, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = SEED),
             size = 0.85, alpha = 0.75, stroke = 0) +
  facet_wrap(~ factor(gene, levels = gsel), nrow = 1, scales = "free_y") +
  scale_fill_manual(values = grp_cols, guide = "none") +
  labs(tag = "B", x = "", y = "Pseudobulk log2 CPM",
       title = "Pseudobulk levels agree for maturation markers but not for CD207",
       subtitle = "Same pseudo-samples; CD207 p = 0.45 whereas LAMP3 p = 2.68e-4, CCR7 p = 0.016, CD274 p = 4.85e-4") +
  TH

m3 <- pb[compartment == "subcluster2_CD207" & gene %in% gsel]
p4c <- ggplot(m3, aes(group, detection_fraction, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.55, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = SEED),
             size = 0.85, alpha = 0.75, stroke = 0) +
  facet_wrap(~ factor(gene, levels = gsel), nrow = 1) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(tag = "C", x = "", y = "Detection fraction",
       title = "No change survives inside the 352-cell CD207+ subcluster",
       subtitle = "24 pseudo-samples, 7 complete PN-PP pairs, all BH > 0.44") +
  TH

p4 <- (p4a / p4b / p4c) +
  plot_annotation(
    caption = paste("POST HOC AND EXPLORATORY: no directional alternative was pre-specified for these genes,",
                    "so every test is two-sided. Pseudo-samples hold few cells and cells are not independent",
                    "replicates. Paired signed-rank with BH across eight genes; 10 donors contribute both a PN",
                    "and a PP pseudo-sample. Detection fraction = share of cells with a non-zero count.",
                    sep = "\n"),
    theme = theme(plot.caption = element_text(size = 6.6, colour = "grey35", hjust = 0)))
save_fig(p4, "Supplementary_Figure_S4", 13.2, 10.4)

## ============================================================ S5: cluster UMAP
obj <- readRDS(file.path(PROC_DIR, F_S1_OBJECT))
say("S1 object for S5: ", ncol(obj), " cells; clusters = ",
    length(unique(obj$seurat_clusters)))
emb <- as.data.frame(Seurat::Embeddings(obj, "umap"))
emb$cluster <- factor(obj$seurat_clusters)
cen <- as.data.table(emb)[, .(umap_1 = median(umap_1), umap_2 = median(umap_2)), by = cluster]
set.seed(SEED)
emb <- emb[sample.int(nrow(emb)), , drop = FALSE]
p5 <- ggplot(emb, aes(umap_1, umap_2, colour = cluster)) +
  geom_point(size = 0.12, alpha = 0.8, stroke = 0, show.legend = FALSE) +
  geom_text(data = cen, aes(umap_1, umap_2, label = cluster), colour = "black",
            size = 3, fontface = "bold") +
  coord_fixed() +
  labs(tag = "S5", x = "UMAP 1", y = "UMAP 2",
       title = sprintf("Cluster-level atlas of %s QC-passing cells (%d principal clusters)",
                       format(ncol(obj), big.mark = ","), length(unique(obj$seurat_clusters))),
       subtitle = "Harmony corrected sample_id only; cluster identities are the ones used for the 6-lineage annotation in Fig. 1",
       caption = "Cluster numbers are Seurat cluster IDs at the resolution used in Section 2.5; the lineage assignment in Fig. 1A collapses these clusters to six compartments.") +
  TH + theme(axis.text = element_blank(), axis.ticks = element_blank())
save_fig(p5, "Supplementary_Figure_S5", 8.4, 7.2)

## ============================================================ S3 assembled
## panel C needs the adjusted-model coefficients: recompute them here, cheaply,
## from the frozen score table rather than parsing the log.
sm <- fread(file.path(TAB_DIR, "S4_module_scores.csv"))
wide <- dcast(sm, sample_id ~ gene_set, value.var = "score")
smap <- unique(sm[, .(sample_id, group)])
wide <- merge(wide, smap, by = "sample_id")[group %in% c("NN", "PP")]
wide[, grp := factor(group, levels = c("NN", "PP"))]
## rebuild the four reference scores is not possible from this table alone, so the
## coefficients are taken verbatim from S10 (same data, same models).
coefs <- data.table(
  model = c("unadjusted", "+ keratinocyte", "+ eccrine", "+ keratinocyte\n+ eccrine",
            "+ keratinocyte\n+ eccrine + mast"),
  beta = c(-0.2586, -0.2828, -0.2428, -0.2670, -0.2499),
  se = c(0.0192, 0.0340, 0.0304, 0.0398, 0.0378),
  p = c(1.07e-25, 1.68e-13, 9.64e-13, 7.39e-10, 1.14e-09))
coefs[, lo := beta - 1.96 * se][, hi := beta + 1.96 * se]
coefs[, model := factor(model, levels = coefs$model)]
pC <- ggplot(coefs, aes(beta, model)) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.5,
                 colour = ok[["blue"]]) +
  geom_point(size = 1.9, colour = ok[["blue"]]) +
  geom_text(aes(label = sprintf("p = %.1e", p)), vjust = -0.9, size = 2.5) +
  labs(tag = "C", x = "Group coefficient for the LC score (PP vs NN)", y = "",
       title = "The LC effect survives adjustment for tissue composition",
       subtitle = "Five linear models, 122 arrays; keratinocyte term p = 0.39 in the fully adjusted model") +
  TH

s3 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    caption = paste("POST HOC: added after internal review, not part of any pre-specified gate.",
                    "The LC score is demoted to a secondary readout precisely because panels A and B show",
                    "that its members are not LC-restricted and that other compartments fall further; panels C",
                    "and D show why the direction nevertheless survives. Single-gene CD207 is the primary bulk",
                    "readout in the manuscript. ssGSEA scores are rank-based and invariant to the monotone 2^x",
                    "rescaling of GSE13355 log2 relative values.", sep = "\n"),
    theme = theme(plot.caption = element_text(size = 6.6, colour = "grey35", hjust = 0)))
save_fig(s3, "Supplementary_Figure_S3", 13.0, 10.2)

say("DONE")
