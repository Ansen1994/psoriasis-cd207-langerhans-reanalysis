## S23_durability_pasi.R - post hoc descriptive layer III (final pre-submission
## addition, same additive rules as S8/S17-S21: no frozen result, parameter,
## gate or main figure is touched).
##
## A. Recovery durability (GSE226244, Section 3.9 readouts split by treatment
##    timeline: Baseline / on12 / on24 / off20). DESCRIPTIVE medians only:
##    the off20 group holds 1 sample and PRE vs POST is unpaired, so no test
##    is run across timelines.
## B. Residual CD207 compartment vs disease severity (GSE183047): PASI is
##    deposited per sample in the SOFT metadata. Descriptive Spearman between
##    PASI and the per-sample CD207+ fraction among annotated myeloid/DC
##    cells, in LESIONAL psoriasis samples only (controls carry PASI = 0 and
##    would create a composition artefact; reported with this caveat).
##
## Output
##   results/tables/S23_recovery_timeline.csv
##   results/tables/S23_pasi_correlation.csv
##   results/figures/Supplementary_Figure_S9.pdf/.png
##
## ASCII only (LC_CTYPE=C). All outputs stay inside the project.
source("config/00_paths.R")
if (!exists("FIG_DPI") && file.exists("config/01_params.R")) source("config/01_params.R")
if (!exists("FIG_DPI")) FIG_DPI <- 600
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

logf <- file.path(LOG_DIR, "S23_durability_pasi.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S23 recovery durability + PASI correlation (post hoc, descriptive) ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
set.seed(SEED)

## ---------------------------------------------------------------- A. durability
tr <- fread(file.path(TAB_DIR, "S19_treatment_bulk_samples.csv"))
tr[, timeline := sub("^treatment_timeline: *", "", timeline)]  # strip SOFT prefix
tl_levels <- c("Baseline", "on12", "on24", "off20")
tr[, tl := fifelse(group == "NN", "Healthy (NN)",
            fifelse(timeline == "Baseline", "Baseline (PRE)",
            fifelse(timeline == "on12", "Week 12 (on12)",
            fifelse(timeline == "on24", "Week 24 (on24)",
                    "Off drug 20 wk (off20)"))))]
tr[, tl := factor(tl, levels = c("Healthy (NN)", "Baseline (PRE)",
                                 "Week 12 (on12)", "Week 24 (on24)",
                                 "Off drug 20 wk (off20)"))]
tr_col <- c("Healthy (NN)" = "#009E73", "Baseline (PRE)" = "#CC79A7",
            "Week 12 (on12)" = "#56B4E9", "Week 24 (on24)" = "#0072B2",
            "Off drug 20 wk (off20)" = "#E69F00")

tl_tab <- rbindlist(lapply(c("expr_CD207", "LC_CD207_signature", "share_sub2"),
  function(ro) tr[, .(readout = ro, timeline = tl, n = sum(!is.na(get(ro))),
                      median = round(median(get(ro), na.rm = TRUE), 4),
                      q25 = round(quantile(get(ro), 0.25, na.rm = TRUE), 4),
                      q75 = round(quantile(get(ro), 0.75, na.rm = TRUE), 4)),
                  by = tl][order(readout, tl)]))
fwrite(tl_tab, file.path(TAB_DIR, "S23_recovery_timeline.csv"))
cat("\n=== A. medians by timeline ===\n")
print(dcast(tl_tab, readout + timeline ~ ., value.var = "median"))

mkpanel <- function(ro, ylab, ttl, sub) {
  ggplot(tr[!is.na(get(ro))], aes(tl, get(ro), fill = tl)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.25, colour = "grey25") +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.55, colour = "grey30") +
    scale_fill_manual(values = tr_col) +
    labs(title = ttl, subtitle = sub, x = NULL, y = ylab) +
    theme_bw(base_size = 7) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          plot.title = element_text(face = "bold", size = 8, hjust = 0),
          plot.subtitle = element_text(size = 6.5, hjust = 0),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 7.2),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none",
          plot.caption = element_text(size = 6.2, hjust = 0, colour = "grey35"))
}
medlab <- function(ro) paste0("medians: ", paste0(
  c("NN", "PRE", "wk12", "wk24", "off20"), " ",
  sprintf("%.3g", as.numeric(tl_tab[readout == ro][order(timeline), median])),
  collapse = " | "))
pA <- mkpanel("expr_CD207", "CD207 (2^x of log2 rel.)", "A | Single-gene CD207",
              medlab("expr_CD207"))
pB <- mkpanel("LC_CD207_signature", "ssGSEA score", "B | LC/CD207 9-gene score",
              medlab("LC_CD207_signature"))
pC <- mkpanel("share_sub2", "Estimated sub2 share", "C | Deconvolved CD207+ share",
              medlab("share_sub2"))
pC <- pC + labs(caption = paste0("Descriptive only: POST arms are unpaired (GSE226244 deposits no patient ID);",
                                 " the off20 group holds 1 sample, so no across-timeline test is run."))

## ---------------------------------------------------------------- B. PASI
s20 <- fread(file.path(TAB_DIR, "S20_treatment_sc_samples.csv"))
softf <- file.path(RAW_DIR, "PSOR_SC_GSE183047", "GSE183047_family.soft.gz")
sl <- readLines(gzfile(softf, open = "rt"), warn = FALSE)
gsm <- character(0); pasi <- numeric(0)
cur <- NA_character_
for (ln in sl) {
  if (grepl("^\\^SAMPLE", ln)) cur <- trimws(sub("^\\^SAMPLE = ", "", ln))
  if (grepl("area-and-severity", ln)) {
    v <- suppressWarnings(as.numeric(trimws(sub(".*area-and-severity index: *", "", ln))))
    gsm <- c(gsm, cur); pasi <- c(pasi, v)
  }
}
pasi_map <- data.table(sample_id = gsm, pasi = pasi)
cat("\n=== B. PASI parsed for ", nrow(pasi_map), " samples ===\n", sep = "")
s20[, gsm_key := sub("_.*$", "", sample_id)]
m <- merge(s20, pasi_map, by.x = "gsm_key", by.y = "sample_id",
           all.x = TRUE, sort = FALSE)
m[, gsm_key := NULL]
fwrite(m, file.path(TAB_DIR, "S20_treatment_sc_samples_with_pasi.csv"))

les <- m[group %in% c("PreTx", "PostTx_wk12", "PostTx_wk24", "PostTx_wk48") &
         !is.na(pasi) & !is.na(frac_CD207_my)]
ctl <- m[group == "Control" & !is.na(pasi) & !is.na(frac_CD207_my)]
rho_f <- suppressWarnings(cor.test(les$frac_CD207_my, les$pasi,
                                   method = "spearman", exact = FALSE))
res <- data.table(
  comparison = "lesional psoriasis samples (PRE + POST), CD207+ fraction vs PASI",
  n = nrow(les), spearman_rho = round(unname(rho_f$estimate), 4), p = rho_f$p.value)
rho_c <- suppressWarnings(cor.test(c(les$frac_CD207_my, ctl$frac_CD207_my),
                 c(les$pasi, ctl$pasi), method = "spearman", exact = FALSE))
res <- rbind(res, data.table(
  comparison = "incl. 10 controls (composition artefact risk; secondary)",
  n = nrow(les) + nrow(ctl), spearman_rho = round(unname(rho_c$estimate), 4),
  p = rho_c$p.value))
fwrite(res, file.path(TAB_DIR, "S23_pasi_correlation.csv"))
cat("\n=== B. PASI Spearman ===\n"); print(res)

les[, tl2 := factor(timepoint, levels = c("baseline", "week12", "week24", "week48"))]
pD <- ggplot(les, aes(pasi, frac_CD207_my, colour = tl2)) +
  geom_point(size = 1.2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.3, colour = "grey40",
              linetype = 2) +
  scale_colour_manual(values = c("baseline" = "#CC79A7", "week12" = "#56B4E9",
                                 "week24" = "#0072B2", "week48" = "#E69F00"),
                      name = NULL) +
  annotate("text", x = max(les$pasi), y = max(les$frac_CD207_my), hjust = 1,
           vjust = 1.2, size = 2.2,
           label = sprintf("Spearman rho = %.2f | n = %d (descriptive)",
                           unname(rho_f$estimate), nrow(les))) +
  labs(title = "D | Residual CD207+ fraction vs PASI",
       subtitle = "GSE183047 lesional psoriasis samples (controls excluded)",
       x = "PASI (deposited per sample)", y = "CD207+ fraction of myeloid/DC") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 8, hjust = 0),
        plot.subtitle = element_text(size = 6.5, hjust = 0),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7.2),
        legend.text = element_text(size = 7.2),
        legend.position = "bottom",
        plot.caption = element_text(size = 6.2, hjust = 0, colour = "grey35"))

figS9 <- (pA | pB) / (pC | pD) + plot_annotation(
  title = "Supplementary Figure S9. Recovery durability by treatment timeline and the PASI link (post hoc, descriptive)",
  theme = theme(plot.title = element_text(size = 9, face = "bold")))
ggsave(file.path(FIG_DIR, "Supplementary_Figure_S9.pdf"), figS9,
       width = 180, height = 130, units = "mm", device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Supplementary_Figure_S9.png"), figS9,
       width = 180, height = 130, units = "mm", dpi = FIG_DPI)
cat("\nfigures written.\n=== S23 done ===\n")
