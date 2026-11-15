## S8_lc_subsets.R - S8 / P1-3: DESCRIPTIVE LC-subset scoring.
##
## NEW analysis authorised by the user on 2026-11-12 (route B: P0 + P1).
## Additive layer: S8_ outputs only, no frozen table overwritten, no frozen
## parameter changed, no main figure redrawn.
##
## Rationale. Liu et al. (Immunity 2021;54:2305-2320.e11,
## doi 10.1016/j.immuni.2021.08.012) defined human LC1 (CD207-hi / CD1A-hi /
## EPCAM-hi), LC2 (CD1C/CD1B/CD14/ITGAM/CLEC10A/CD63/CD40/CD274) and
## activated/migratory LC (CD83/CCR7/HLA-DRA). LC2 and aLC/migLC were reported
## to be enriched in psoriatic lesions while classical LC1 stayed stable.
## If the frozen 10-gene CD207 signature falls in lesions, an alternative to
## "true depletion" is a RELATIVE loss of the classical LC1 programme with
## expansion of activated/monocyte-derived LC. This script asks whether the
## CD207+ subcluster (subcluster 2) contains such heterogeneity.
##
## HARD RULE: this is DESCRIPTIVE ONLY. Under the frozen plan subcluster 2 has
## 1 NS vs 3 PP donor pseudo-samples (two-group DEG was skipped, all 50 GSVA
## tests NOT COMPUTABLE). NO significance test, NO p value, NO group claim is
## produced here. Cells are not independent replicates.
##
## ASCII only. All outputs stay inside the project.
source("config/00_paths.R")
source("config/01_params.R")
source("config/02_s8_params.R")
suppressPackageStartupMessages({library(data.table)})

logf <- file.path(LOG_DIR, "S8_lc_subsets.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S8 / P1-3 descriptive LC-subset scoring ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
cat("DESCRIPTIVE ONLY. No p values. Cells are NOT independent replicates.\n")
cat("Subcluster 1 (648 sweat-gland contaminants) is excluded as in S2.\n")
set.seed(SEED)

s <- readRDS(file.path(PROC_DIR, F_S2_OBJECT))
meta <- as.data.table(s@meta.data)
cat("\nS2 object: ", ncol(s), " cells | subclusters: ",
    paste(sort(unique(meta$subcluster)), collapse = ","), "\n", sep = "")

a <- Seurat::GetAssayData(s, assay = "RNA", layer = "data")
cat("RNA data layer: ", nrow(a), " genes x ", ncol(a), " cells\n", sep = "")

## QC-eligible cells = all four subclusters EXCEPT the contaminant (1).
keep_cells <- meta$subcluster != 1
cat("QC-eligible myeloid cells: ", sum(keep_cells), " (of ", length(keep_cells),
    ")\n", sep = "")

## z-score each gene across QC-eligible cells, then average per set.
## This is deterministic (no random control-gene sampling) and reproducible.
X <- as.matrix(a[, keep_cells, drop = FALSE])
mu <- rowMeans(X, na.rm = TRUE)
sdv <- apply(X, 1, stats::sd, na.rm = TRUE)
sdv[!is.finite(sdv) | sdv == 0] <- 1
Z <- (X - mu) / sdv
rm(X); gc(FALSE)

sets <- S8_LC_SETS
for (n in names(sets)) {
  hit <- intersect(sets[[n]], rownames(Z))
  cat("set ", n, ": ", length(hit), "/", length(sets[[n]]),
      " genes present -> ", paste(hit, collapse = ","), "\n", sep = "")
  stopifnot(length(hit) >= 2)
}

score_dt <- data.table(cell_id = colnames(Z))
for (n in names(sets)) {
  hit <- intersect(sets[[n]], rownames(Z))
  score_dt[, (n) := colMeans(Z[hit, , drop = FALSE], na.rm = TRUE)]
}
score_dt[, `:=`(
  group = meta$group[keep_cells],
  donor = as.character(meta$donor[keep_cells]),
  donor_base = as.character(meta$donor_base[keep_cells]),
  subcluster = meta$subcluster[keep_cells],
  CD207_expr = as.numeric(Z["CD207", ]),
  CD1A_expr  = as.numeric(Z["CD1A", ])
)]

## Focus: the CD207+ Langerhans/DC subcluster (2).
lc <- score_dt[subcluster == 2]
cat("\nCD207+ subcluster 2: ", nrow(lc), " cells from ",
    uniqueN(lc$donor), " donors / ", uniqueN(lc$donor_base),
    " donor base IDs\n", sep = "")
print(table(lc$group, useNA = "ifany"))

## Dominance label is a per-cell descriptive label, NOT a classifier.
lc[, LC1_minus_LC2 := LC1_classical - LC2_monocytic]
lc[, dominant := ifelse(LC1_minus_LC2 > 0, "LC1-like", "LC2-like")]

fwrite(lc[, .(cell_id, group, donor, donor_base, subcluster,
              LC1_classical, LC2_monocytic, aLC_migLC, LC1_minus_LC2,
              dominant, CD207_expr, CD1A_expr)],
       file.path(TAB_DIR, "S8_lc_subset_cell_scores.csv"))

grp <- lc[, .(
  n_cells = .N,
  n_donors = uniqueN(donor),
  median_LC1 = round(median(LC1_classical), 4),
  median_LC2 = round(median(LC2_monocytic), 4),
  median_aLC = round(median(aLC_migLC), 4),
  median_LC1_minus_LC2 = round(median(LC1_minus_LC2), 4),
  pct_LC1_like = round(100 * mean(dominant == "LC1-like"), 1),
  median_CD207_z = round(median(CD207_expr), 4)
), by = group]
setorder(grp, group)
fwrite(grp, file.path(TAB_DIR, "S8_lc_subset_summary.csv"))
cat("\n--- descriptive group summary (NO test) ---\n"); print(grp)

don <- lc[, .(
  n_cells = .N,
  median_LC1 = round(median(LC1_classical), 4),
  median_LC2 = round(median(LC2_monocytic), 4),
  median_aLC = round(median(aLC_migLC), 4),
  median_LC1_minus_LC2 = round(median(LC1_minus_LC2), 4)
), by = .(group, donor)]
setorder(don, group, donor)
fwrite(don, file.path(TAB_DIR, "S8_lc_subset_donor_summary.csv"))
cat("\n--- per-donor descriptive values (NO test) ---\n"); print(don)

## Within-subcluster correlations, descriptive only.
ct <- suppressWarnings(cor(lc$LC1_classical, lc$CD207_expr,
                           method = "spearman", use = "complete.obs"))
ct2 <- suppressWarnings(cor(lc$LC1_classical, lc$LC2_monocytic,
                            method = "spearman", use = "complete.obs"))
cat("\nSpearman LC1 score vs CD207 z: ", round(ct, 3), "\n", sep = "")
cat("Spearman LC1 score vs LC2 score: ", round(ct2, 3), "\n", sep = "")
fwrite(data.table(
  readout = c("spearman_LC1_vs_CD207_z", "spearman_LC1_vs_LC2"),
  value = round(c(ct, ct2), 4),
  interpretation = c("descriptive within-subcluster co-variation",
                     "descriptive within-subcluster co-variation")),
  file.path(TAB_DIR, "S8_lc_subset_correlations.csv"))

cat("\n=== LIMITATIONS (must be repeated in the manuscript) ===\n")
cat("1. NO group-level test is reported: frozen plan gives 1 NS vs 3 PP donor\n")
cat("   pseudo-samples in subcluster 2, so inference is not computable here.\n")
cat("2. Cells are not independent replicates; group medians are dominated by\n")
cat("   the donors that contributed the most cells.\n")
cat("3. Marker sets come from surface-protein phenotypes (Immunity 2021);\n")
cat("   transcript abundance does not prove protein-level subset identity.\n")
cat("4. A relative LC1/LC2 shift cannot distinguish migration, CD207\n")
cat("   down-regulation, dilution by inflammatory DC, or true depletion.\n")
cat("\n=== S8 LC subsets DONE ===\n")
cat("OUTPUT: S8_lc_subset_cell_scores.csv / S8_lc_subset_summary.csv\n")
cat("OUTPUT: S8_lc_subset_donor_summary.csv / S8_lc_subset_correlations.csv\n")
