## S13b_tf_figures.R
## Supplementary Figure S6. Rendering only: reads S13 tables, performs NO test,
## changes NO parameter. Every limitation that a reader needs in order to judge
## the panel is printed on the panel itself.
##
## Panels
##   A  TF activity change PP - PN in the pooled myeloid layer, A+B regulon, for
##      the 23 TFs significant on the paired test; filled by whether the change
##      also survives depth adjustment.
##   B  Depth diagnostic: activity change vs correlation of that change with the
##      change in cells per pseudo-sample.
##   C  Specification dependence: how many TFs are robust under each of the four
##      layer x regulon combinations.
##   D  Median activity of the 12 robust TFs across the three states, row-scaled.
##
## Output
##   results/figures/Supplementary_Figure_S6.pdf
##   results/figures/Supplementary_Figure_S6.png   (FIG_DPI dpi)

source("config/00_paths.R")
## FIG_DPI lives in the shared parameter file; this script needs it so the PNG
## rasterises at the same density as the main figures (Dove requires >=600 dpi).
if (!exists("FIG_DPI") && file.exists("config/01_params.R")) {
  source("config/01_params.R")
}
if (!exists("FIG_DPI")) FIG_DPI <- 600
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

t <- fread(file.path(TAB_DIR, "S13_tf_activity_tests.csv"))
t[, class := fifelse(robust, "Robust",
              fifelse(!is.na(paired_FDR) & paired_FDR < 0.05,
                      "Significant but depth-dependent", "Not significant"))]
CLS <- c("Robust" = "#0072B2",
         "Significant but depth-dependent" = "#E69F00",
         "Not significant" = "#BBBBBB")

## ggplot2 4.x stores themes as S7 objects, which patchwork 1.3 cannot pass
## through plot_annotation(theme = ...); setting the global caption theme here
## is what the combined-figure caption actually inherits.
theme_set(theme_bw(base_size = 7) +
            theme(plot.caption = element_text(size = 6.0, colour = "grey35",
                                              hjust = 0)))

th <- theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.title = element_text(face = "bold", size = 8, hjust = 0),
        plot.subtitle = element_text(size = 6, colour = "grey30", hjust = 0),
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 6))

## ---- A ------------------------------------------------------------------
a <- t[layer == "myeloid_all" & network == "AB" & !is.na(paired_FDR) &
         paired_FDR < 0.05][order(paired_delta)]
a[, tf := factor(tf, levels = tf)]
a[, lab := sprintf("%+.2f", paired_delta)]
pA <- ggplot(a, aes(paired_delta, tf, fill = class)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = lab, hjust = ifelse(paired_delta > 0, -0.15, 1.15)),
            size = 2.15, colour = "grey25") +
  scale_fill_manual(values = CLS, drop = TRUE, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.26, 0.26))) +
  theme(axis.text.y = element_text(face = "italic")) +
  labs(title = "A  TF activity change, PP vs PN",
       subtitle = sprintf("Pooled myeloid, A+B. %d/%d significant; %d after depth adjustment",
                          nrow(a), t[layer == "myeloid_all" & network == "AB", .N],
                          a[robust == TRUE, .N]),
       x = "Paired mean change in TF activity score (PP - PN)", y = NULL) +
  th + theme(legend.position = "bottom")

## ---- B ------------------------------------------------------------------
b <- t[layer == "myeloid_all" & network == "AB"]
b[, sig := !is.na(paired_FDR) & paired_FDR < 0.05]
b[, lab := fifelse(robust == TRUE | (sig & abs(rho_delta_depth) >= 0.6), tf, "")]
pB <- ggplot(b, aes(rho_delta_depth, -log10(paired_p))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.25, colour = "grey60") +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = 2, linewidth = 0.25, colour = "grey60") +
  geom_point(aes(colour = class, shape = sig), size = 1.1, stroke = 0.25) +
  geom_text_repel(aes(label = lab), size = 2.12, colour = "grey20", fontface = "italic",
                  min.segment.length = 0, segment.size = 0.15,
                  max.overlaps = 30, box.padding = 0.15) +
  scale_colour_manual(values = CLS, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(title = "B  Depth diagnostic",
       subtitle = "dashed: BH 0.05 / |rho| = 0.6",
       x = "Spearman rho between activity change and change in cells per pseudo-sample",
       y = "-log10 paired p") +
  th

## ---- C ------------------------------------------------------------------
c0 <- t[, .(n_robust = sum(robust),
            n_depth = sum(!robust & !is.na(paired_FDR) & paired_FDR < 0.05),
            n_ns = sum(is.na(paired_FDR) | paired_FDR >= 0.05)),
        by = .(layer, network)]
NLAB <- c(AB = "A+B", ABC = "A+B+C")
c0[, combo := sprintf("%s | %s",
                      fifelse(layer == "myeloid_all", "Pooled myeloid", "Subcluster 0"),
                      NLAB[network])]
c0[, combo := factor(combo, levels = rev(unique(combo)))]
c0[, total := n_robust + n_depth + n_ns]
cm <- melt(c0, id.vars = c("layer", "network", "combo", "total"),
           measure.vars = c("n_robust", "n_depth", "n_ns"),
           variable.name = "cls", value.name = "n")
cm[, cls := factor(cls, levels = c("n_ns", "n_depth", "n_robust"),
                   labels = c("Not significant", "Significant but depth-dependent", "Robust"))]
pC <- ggplot(cm, aes(n, combo, fill = cls)) +
  geom_col(width = 0.62) +
  geom_text(data = c0, aes(x = total, y = combo, label = sprintf("%d robust", n_robust)),
            inherit.aes = FALSE, hjust = -0.12, size = 2.15, colour = "grey20") +
  scale_fill_manual(values = CLS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "C  Specification dependence",
       subtitle = "TF counts per layer x regulon; BH within each combination",
       x = "Number of TFs", y = NULL) +
  th + theme(legend.position = "bottom")

## ---- D ------------------------------------------------------------------
d <- t[layer == "myeloid_all" & network == "AB" & robust == TRUE,
       .(tf, NS = median_NS, PN = median_PN, PP = median_PP)]
dm <- melt(d, id.vars = "tf", variable.name = "state", value.name = "score")
dm[, state := factor(state, levels = c("NS", "PN", "PP"))]
dm[, z := as.numeric(scale(score)), by = tf]
ord <- d[order(NS - PP)]$tf
dm[, tf := factor(tf, levels = ord)]
pD <- ggplot(dm, aes(state, tf, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#56B4E9", mid = "white", high = "#D55E00",
                       midpoint = 0, name = "row z") +
  labs(title = "D  The 12 robust TFs across states",
       subtitle = "median activity score, row-scaled; descriptive only",
       x = NULL, y = NULL) +
  th + theme(panel.grid = element_blank(),
              axis.text.y = element_text(face = "italic"))

fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    caption = paste(
      "TF activity = decoupleR::run_ulm over the dorothea regulon, inferred from mRNA of target genes;",
      "not measured protein activity, not DNA binding, not a perturbation.",
      "\nA TF absent from the A+B regulon is not assessable, not unchanged.",
      "Pseudo-samples hold few cells (median 17 / 34 / 73 for NS / PN / PP), and activity change tracks",
      "\nthat depth for part of the regulon; the third 'robust' criterion therefore also requires |rho| < 0.6.",
      "Exploratory and post hoc; does not enter Gate C/D."))

ggsave(file.path(FIG_DIR, "Supplementary_Figure_S6.pdf"), fig,
       width = 180, height = 176, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Supplementary_Figure_S6.png"), fig,
       width = 180, height = 176, units = "mm", dpi = FIG_DPI)
cat("WROTE Supplementary_Figure_S6.{pdf,png}\n")
