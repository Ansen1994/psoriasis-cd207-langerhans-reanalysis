## S8_figures.R - S8 / P1-4 + P1-5: two SUPPLEMENTARY figures.
##
## NEW deliverables authorised 2026-11-12 (route B: P0 + P1).
## Renders ONLY supplementary figures. Main Fig 1-4 are NOT touched, NOT
## redrawn, NOT overwritten. Reads frozen S4 tables and new S8_ tables only;
## no statistical test is performed here (all tests live in S8_*.R scripts).
##
## Outputs: results/figures/Supplementary_Figure_S1.{pdf,png}
##          results/figures/Supplementary_Figure_S2.{pdf,png}
## ASCII only. 300 dpi PNG, single-page vector PDF.
source("config/00_paths.R")
source("config/01_params.R")
source("config/02_s8_params.R")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

logf <- file.path(LOG_DIR, "S8_figures.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S8 supplementary figures ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")

rd <- function(f) fread(file.path(TAB_DIR, f))
base_theme <- theme_minimal(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "#444444"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7.5),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, face = "bold"))
BLUE <- "#4C78A8"; RED <- "#C44E52"; ORANGE <- "#F58518"; GREY <- "#7F7F7F"

save_fig <- function(p, stem, w, h) {
  pdfp <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  pngp <- file.path(FIG_DIR, paste0(stem, ".png"))
  ggsave(pdfp, p, width = w, height = h, units = "in", device = "pdf")
  ggsave(pngp, p, width = w, height = h, units = "in", dpi = FIG_DPI)
  cat("WROTE ", pdfp, " (", w, "x", h, " in)\n", sep = "")
  cat("WROTE ", pngp, " (300 dpi)\n", sep = "")
}

## ======================================================= S1: inclusion flow ===
cat("\n--- Supplementary Figure S1: data-inclusion flow ---\n")
left <- data.table(
  x = 2.25, w = 3.9, h = 0.95,
  y = c(11, 9.8, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6),
  label = c(
    "GSE173706 DISCOVERY scRNA-seq\n33 libraries | 96,088 cells | 33,538 Ensembl IDs",
    "Cell QC: 200 <= nFeature <= 6,000; pct.mt <= 20%\n89,120 cells retained (92.75%)",
    "Sample-level exclusion: NS-AR008 (low depth)\n88,961 cells / 32 samples",
    "scDblFinder per sample: 4,521 doublets (5.08%)\n84,440 cells | 32 samples | 20 donor base IDs",
    "Gate B: marker lineage, tie-break -> myeloid/DC\n2,314 cells (6 lineages, 19 main clusters)",
    "Myeloid re-clustering k = 4 @ res 0.0178\nsubclusters: 1,186 / 648 / 352 / 128",
    "Contamination QC: subcluster 1 = sweat-gland\n648 cells (28%) EXCLUDED -> 1,666 eligible",
    "CD207+ Langerhans/DC = subcluster 2\n352 cells; inference on 29 donor x state pseudo-samples"),
  fill = c("#DCE6F2", "#DCE6F2", "#DCE6F2", "#C6D9F1",
           "#FCE4D6", "#FCE4D6", "#F8CBAD", "#E2EFDA"))
right <- data.table(
  x = 6.75, w = 3.9, h = 1.55,
  y = c(11, 9.0, 7.0, 5.0),
  label = c(
    "INDEPENDENT sc COHORT\nGSE162183: 3 healthy + 3 lesional\n22,716 QC cells | descriptive 3 vs 3",
    "PRIMARY BULK COHORT\nGSE13355: 180 arrays / 122 subjects\nNN 64 / PN 58 / PP 58; 58 paired PN-PP",
    "S8 EXTERNAL BULK (NEW)\nGSE14905: 82 arrays / 54 subjects\n21 NN, 33 PP, 28 PN; 28 complete pairs",
    "S8 EXTERNAL BULK (NEW)\nGSE30999: 170 arrays / 89 subjects\n81 matched NL/LS pairs + 8 single biopsies; no healthy controls"),
  fill = c("#EDEDED", "#EDEDED", "#FFF2CC", "#FFF2CC"))
boxes <- rbind(left, right, fill = TRUE)
arrows <- data.table(
  x = 2.25, xend = 2.25,
  y = left$y[-nrow(left)] - left$h[1] / 2,
  yend = left$y[-1] + left$h[1] / 2)

p1 <- ggplot() +
  geom_rect(data = boxes,
            aes(xmin = x - w / 2, xmax = x + w / 2,
                ymin = y - h / 2, ymax = y + h / 2, fill = fill),
            colour = "#555555", linewidth = 0.35) +
  geom_text(data = boxes,
            aes(x = x, y = y, label = label), size = 2.55,
            lineheight = 1.15, colour = "#111111") +
  geom_segment(data = arrows,
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
               colour = "#555555", linewidth = 0.4) +
  annotate("text", x = 2.25, y = 12.0, label = "DISCOVERY: cell / sample inclusion",
           size = 3.1, fontface = "bold") +
  annotate("text", x = 6.75, y = 12.0, label = "VALIDATION AND S8 COHORTS",
           size = 3.1, fontface = "bold") +
  annotate("text", x = 4.5, y = 1.15,
           label = paste0("INFERENCE UNIT = donor / sample, NEVER cell.  ",
                          "S8 analyses are POST-HOC supplementary and do NOT\n",
                          "enter the frozen Gate C / Gate D criteria.  ",
                          "Main Fig 1-4 remain unchanged."),
           size = 2.7, colour = "#333333", lineheight = 1.2) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(0, 9), ylim = c(0.4, 12.4), expand = FALSE) +
  labs(tag = "S1", title = "Supplementary Figure S1. Cohort structure and data-inclusion flow",
       subtitle = "Frozen S0-S4 pipeline shown with the two NEW S8 external GPL570 cohorts (yellow)") +
  theme_void() + theme(plot.tag = element_text(size = 11, face = "bold"),
                       plot.title = element_text(size = 10, face = "bold"),
                       plot.subtitle = element_text(size = 8),
                       plot.margin = margin(4, 4, 4, 4))
save_fig(p1, "Supplementary_Figure_S1", 13, 12)

## ============================================== S2: orthogonal S8 evidence ===
cat("\n--- Supplementary Figure S2: orthogonal and external evidence ---\n")

## Panel A: GSE13355 single-gene deltas.
a <- rd("S8_bulk_single_gene_tests.csv")[
  , .(gene, delta = delta_PP_minus_NN, p = FDR_BH_wilcox)]
a[, dir := ifelse(delta < 0, "lower in PP", "higher in PP")]
a[, lab := sprintf("BH %.1e", p)]
pa <- ggplot(a, aes(x = reorder(gene, delta), y = delta, fill = dir)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = lab, y = 3.25), hjust = 0,
            size = 2.3, colour = "#333333") +
  geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.35) +
  coord_flip() +
  scale_y_continuous(limits = c(-1.9, 4.6)) +
  scale_fill_manual(values = c("lower in PP" = BLUE, "higher in PP" = RED)) +
  labs(tag = "A", title = "Single-gene orthogonal check (GSE13355)",
       subtitle = "log2 PP - NN median; 64 NN / 58 PP; Wilcoxon + BH (labels at right); 58 pairs also significant",
       x = NULL, y = "median log2 difference (PP - NN)", fill = "Direction") +
  base_theme + theme(legend.position = "bottom")

## Panel B: LC 10-gene signature across three bulk cohorts.
s4 <- rd("S4_bulk_group_tests.csv")[gene_set == "LC_CD207_signature"]
s8 <- rd("S8_external_bulk_tests.csv")[gene_set == "LC_CD207_signature"]
b <- data.table(
  cohort = c("GSE13355\nPP vs NN", "GSE14905\nPP vs NN", "GSE30999\npaired PN-PP"),
  delta = c(s4$delta_PP_minus_NN,
            s8[cohort == "GSE14905", delta_PP_minus_NN],
            s8[cohort == "GSE30999", delta_paired_PP_minus_PN]),
  p = c(s4$p_wilcox,
        s8[cohort == "GSE14905", p_wilcox_PP_vs_NN],
        s8[cohort == "GSE30999", p_signed_rank_paired]))
b[, lab := sprintf("p = %.1e", p)]
pb <- ggplot(b, aes(x = reorder(cohort, -delta), y = delta)) +
  geom_col(width = 0.6, fill = BLUE) +
  geom_text(aes(label = lab, y = delta - 0.012), size = 2.6, colour = "#333333") +
  geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.35) +
  labs(tag = "B", title = "LC/CD207 signature replicated in two NEW cohorts",
       subtitle = "Same frozen 10-gene signature and ssGSEA settings; all three p < 1e-4",
       x = NULL, y = "median score difference") +
  base_theme

## Panel C: CD207 vs LAMP3 in three cohorts (opposite directions).
sg0 <- rd("S8_bulk_single_gene_tests.csv")[
  gene %in% c("CD207", "LAMP3"),
  .(cohort = "GSE13355\nPP vs NN", gene, delta = delta_PP_minus_NN)]
sge <- rd("S8_external_single_gene_tests.csv")[gene %in% c("CD207", "LAMP3")]
sg1 <- sge[cohort == "GSE14905",
           .(cohort = "GSE14905\nPP vs NN", gene, delta = delta_PP_minus_NN)]
sg2 <- sge[cohort == "GSE30999",
           .(cohort = "GSE30999\npaired PN-PP", gene,
             delta = delta_paired_PP_minus_PN)]
c <- rbind(sg0, sg1, sg2)
c[, cohort := factor(cohort, levels = c("GSE13355\nPP vs NN",
                                        "GSE14905\nPP vs NN",
                                        "GSE30999\npaired PN-PP"))]
pc <- ggplot(c, aes(x = cohort, y = delta, fill = gene)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.58) +
  geom_hline(yintercept = 0, colour = "#333333", linewidth = 0.35) +
  scale_fill_manual(values = c("CD207" = BLUE, "LAMP3" = RED)) +
  labs(tag = "C", title = "Classical LC marker falls while LAMP3 rises",
       subtitle = "Compartment-selective change, NOT a global loss of dendritic cells",
       x = NULL, y = "median log2 difference", fill = "Gene") +
  base_theme + theme(legend.position = "bottom",
                     axis.text.x = element_text(size = 6.8))

## Panel D: descriptive LC subset scores (NO test).
lc <- rd("S8_lc_subset_summary.csv")
ld <- melt(lc[, .(group, LC1 = median_LC1, LC2 = median_LC2,
                  aLC = median_aLC)],
           id.vars = "group", variable.name = "set", value.name = "score")
ld[, group := factor(group, levels = c("NS", "PN", "PP"), ordered = TRUE)]
nl <- lc[, .(group, txt = sprintf("n=%d cells / %d donors", n_cells, n_donors))]
pd <- ggplot(ld, aes(x = group, y = score, colour = set, group = set)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  geom_text(data = nl, aes(x = group, y = -0.72, label = txt),
            inherit.aes = FALSE, size = 2.2, colour = "#333333") +
  scale_colour_manual(values = c("LC1" = BLUE, "LC2" = ORANGE, "aLC" = GREY)) +
  labs(tag = "D", title = "LC-subset scores: DESCRIPTIVE ONLY - NO TEST",
       subtitle = "CD207+ subcluster 2 has 1 NS vs 3 PP donor pseudo-samples; cells are not replicates",
       x = "Disease state", y = "median gene-set z-score", colour = "Signature") +
  base_theme + theme(legend.position = "bottom")

lim_lines <- c(
  "MANDATORY LIMITATIONS FOR THE S8 SUPPLEMENTARY ANALYSES",
  "1. S8 is a POST-HOC layer added after the frozen S1-S4 analysis; it does NOT enter the pre-registered Gate C / Gate D criteria.",
  "2. Single-gene log2 levels are orthogonal to ssGSEA but STILL cannot distinguish fewer CD207+ cells from lower per-cell expression.",
  "3. All three bulk cohorts are GPL570 series matrices of upstream-normalised log2 RELATIVE values; batch/sex adjustment is irreversible.",
  "4. GSE30999 has NO healthy controls (81 paired PN-PP only); GSE14905 has 21 NN / 33 PP / 28 PN and 28 complete pairs.",
  "5. LC-subset scores are DESCRIPTIVE: no p value is reported, no group claim is made, and transcript abundance does not prove protein-level identity.",
  "6. LAMP3 rises in all three cohorts, so the LC signature decrease is compartment-selective rather than a global DC loss.",
  "7. Specificity controls (Cycling, Mac_DC) also change in paired external tests; only the LC signature decreases in every cohort.")
pe <- ggplot() +
  annotate("text", x = 0.01, y = 0.97,
           label = paste(lim_lines, collapse = "\n"),
           hjust = 0, vjust = 1, size = 2.35, lineheight = 1.28,
           colour = "#222222") +
  xlim(0, 1) + ylim(0, 1) + theme_void() +
  theme(panel.background = element_rect(fill = "#F4F4F4", colour = NA))

p2 <- (pa | pb) / (pc | pd) / pe +
  plot_layout(heights = c(1, 1, 0.78)) +
  plot_annotation(
    title = "Supplementary Figure S2. Orthogonal single-gene and two independent bulk-cohort validations",
    subtitle = "Added in S8 (2026-11-12) after authorisation of route B; main Fig 1-4 and all frozen results are unchanged")
save_fig(p2, "Supplementary_Figure_S2", 14, 14)

cat("\n=== S8 figures DONE ===\n")
