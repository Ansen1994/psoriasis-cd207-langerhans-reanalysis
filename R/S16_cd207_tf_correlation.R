## S16_cd207_tf_correlation.R
## Added 2026-11-13 in response to the third-round internal review (R5-1).
##
## Question
##   S13 established that TF activity differs between non-lesional and lesional
##   myeloid pseudobulk, and the manuscript reports this as directional
##   consistency. A reviewer then asks the obvious question: is the TF signal
##   related to CD207 at all, or are they two parallel outputs of the same
##   contrast? This script answers it directly.
##
## Design
##   Unit of observation = the 29 donor x state pseudo-samples (the same ones
##   used everywhere in this project). For each TF we correlate its activity
##   score with two CD207 readouts:
##     (a) pseudobulk log2 CPM  - average expression per pseudo-sample
##     (b) detection fraction   - share of cells with a non-zero CD207 count,
##                                which is closer to "how many CD207+ cells"
##
##   THREE correlations are reported for each TF, in increasing strictness:
##     1. raw Spearman across all 29 pseudo-samples
##     2. partial Spearman controlling for disease state (rank residuals)
##     3. partial Spearman controlling for state AND pseudo-sample depth
##
##   Step 2 is the one that matters. CD207 falls and TF activity rises from
##   healthy to lesional, so the raw correlation is guaranteed to be negative
##   even if no within-state relationship exists. Controlling for state removes
##   that between-group component; whatever survives is a within-state
##   association. If the partial correlation vanishes, the honest reading is
##   that the two layers are parallel descriptions of the same group contrast,
##   not evidence that TF activity drives CD207 loss.
##
## HARD LIMITS (reported in the log and available for the manuscript)
##   1. n = 29, so every correlation is underpowered and the BH family is 111.
##   2. TF activity is inferred from target-gene mRNA, not measured protein
##      activity or DNA binding.
##   3. Correlation is not causation, and nothing here is a perturbation.
##   4. The two readouts share the same cells; they are not independent data.
##
## Outputs
##   results/tables/S16_cd207_tf_correlation.csv
##   logs/S16_cd207_tf_correlation.log
##
## This script does not modify any frozen table, figure or parameter.

source("config/00_paths.R")
source("config/01_params.R")

suppressMessages(library(data.table))

LOG <- file.path(LOG_DIR, "S16_cd207_tf_correlation.log")
if (file.exists(LOG)) file.remove(LOG)
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = LOG, append = TRUE)
}
says <- function(x) say(paste(capture.output(print(x)), collapse = "\n"))
hr <- function(t) say("\n", strrep("-", 72), "\n", t, "\n", strrep("-", 72))

say("S16_cd207_tf_correlation.R  |  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R ", getRversion())

## ------------------------------------------------------------------ 1. load
hr("1. Load S13 sample-level activity and S11 CD207 pseudobulk")

s <- fread(file.path(TAB_DIR, "S13_tf_activity_samples.csv"))
s <- s[layer == "myeloid_all" & network == "AB"]
say("S13 rows (myeloid_all, A+B): ", nrow(s))
say("TFs: ", uniqueN(s$tf), "   pseudo-samples: ", uniqueN(s$pseudo_id))

act <- dcast(s, donor + group + n_cells ~ tf, value.var = "score")
say("activity matrix: ", nrow(act), " pseudo-samples x ", ncol(act) - 3, " TFs")
stopifnot(nrow(act) == 29)

pb <- fread(file.path(TAB_DIR, "S11_lc_maturation_pseudobulk.csv"))
pb <- pb[compartment == "myeloid_qc_eligible" & gene == "CD207"]
say("S11 CD207 rows: ", nrow(pb))
stopifnot(nrow(pb) == 29)

m <- merge(act,
           pb[, .(donor, group, cpm = pseudobulk_log2cpm, det = detection_fraction)],
           by = c("donor", "group"))
stopifnot(nrow(m) == 29)
say("merged: ", nrow(m), " pseudo-samples")
m[, state_ord := fifelse(group == "NS", 0L, fifelse(group == "PN", 1L, 2L))]
says(m[, .(n = .N, cells_min = min(n_cells), cells_med = as.integer(median(n_cells)),
           cells_max = max(n_cells),
           cpm_med = round(median(cpm), 2), det_med = round(median(det), 3)),
       by = group][order(group)])
say("CD207 vs state: rho(cpm,state) = ", round(cor(m$cpm, m$state_ord, method = "spearman"), 3),
    " ; rho(det,state) = ", round(cor(m$det, m$state_ord, method = "spearman"), 3))

## ------------------------------------------------------- 2. correlation fns
## Rank-based partial correlation. Both variables and the covariates are
## replaced by ranks, the covariates are regressed out of each, and the residual
## correlation is computed. With one covariate this is the standard
## Spearman partial correlation; with two it is its natural extension.
pcor_rank <- function(x, y, Z) {
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  X <- cbind(1, Z)
  ex <- resid(lm.fit(X, rx))
  ey <- resid(lm.fit(X, ry))
  r <- cor(ex, ey)
  df <- length(x) - ncol(X) - 1L
  t <- r * sqrt(df / (1 - r^2))
  list(rho = r, p = 2 * pt(-abs(t), df), df = df)
}

tfs <- setdiff(names(act), c("donor", "group", "n_cells"))
ROBUST <- c("IRF1", "STAT1", "STAT2", "STAT3", "RELB", "NFKB1", "NFKB2",
            "JUND", "RARA", "STAT5A", "STAT5B", "NR2F2")

res <- rbindlist(lapply(tfs, function(g) {
  x <- m[[g]]
  c_raw_cpm <- cor.test(x, m$cpm, method = "spearman", exact = FALSE)
  c_raw_det <- cor.test(x, m$det, method = "spearman", exact = FALSE)
  p1 <- pcor_rank(x, m$cpm, cbind(m$state_ord))
  p2 <- pcor_rank(x, m$det, cbind(m$state_ord))
  p3 <- pcor_rank(x, m$cpm, cbind(m$state_ord, m$n_cells))
  p4 <- pcor_rank(x, m$det, cbind(m$state_ord, m$n_cells))
  data.table(
    tf = g,
    robust_12 = g %in% ROBUST,
    rho_cpm = unname(c_raw_cpm$estimate), p_cpm = c_raw_cpm$p.value,
    rho_det = unname(c_raw_det$estimate), p_det = c_raw_det$p.value,
    prho_cpm_state = p1$rho, p_prho_cpm_state = p1$p,
    prho_det_state = p2$rho, p_prho_det_state = p2$p,
    prho_cpm_state_cells = p3$rho, p_prho_cpm_state_cells = p3$p,
    prho_det_state_cells = p4$rho, p_prho_det_state_cells = p4$p)
}))

for (v in c("p_cpm", "p_det", "p_prho_cpm_state", "p_prho_det_state",
            "p_prho_cpm_state_cells", "p_prho_det_state_cells")) {
  set(res, j = sub("^p_", "fdr_", v), value = p.adjust(res[[v]], "BH"))
}
setorder(res, p_prho_det_state)
fwrite(res, file.path(TAB_DIR, "S16_cd207_tf_correlation.csv"))
say("written: ", file.path(TAB_DIR, "S16_cd207_tf_correlation.csv"))

## ------------------------------------------------------------- 3. reporting
hr("2. How many TFs correlate with CD207 (BH < 0.05, family of 111)")
cnt <- res[, .(raw_cpm = sum(fdr_cpm < 0.05), raw_det = sum(fdr_det < 0.05),
               part_state_cpm = sum(fdr_prho_cpm_state < 0.05),
               part_state_det = sum(fdr_prho_det_state < 0.05),
               part_state_cells_cpm = sum(fdr_prho_cpm_state_cells < 0.05),
               part_state_cells_det = sum(fdr_prho_det_state_cells < 0.05))]
say("TFs tested: ", nrow(res), "   (12 of them are the 'robust' set from S13)")
says(cnt)

hr("3. The 12 robust TFs from S13, ordered by within-state association")
says(res[robust_12 == TRUE,
         .(tf,
           rho_cpm = round(rho_cpm, 2), fdr_cpm = signif(fdr_cpm, 2),
           rho_det = round(rho_det, 2), fdr_det = signif(fdr_det, 2),
           prho_det_state = round(prho_det_state, 2),
           fdr_part = signif(fdr_prho_det_state, 2),
           prho_det_state_cells = round(prho_det_state_cells, 2),
           fdr_part2 = signif(fdr_prho_det_state_cells, 2))])

hr("4. Are the 12 robust TFs more strongly (negatively) associated than the other 99?")
## Two-sided Wilcoxon on the within-state partial correlation. A difference here
## would mean the TFs that survive S13 also track CD207 better than average;
## no difference would mean the TF layer is not specifically tied to CD207.
w_det <- wilcox.test(prho_det_state ~ robust_12, data = res, exact = FALSE)
w_cpm <- wilcox.test(prho_cpm_state ~ robust_12, data = res, exact = FALSE)
say("median partial rho (det, control state): robust = ",
    round(median(res[robust_12 == TRUE, prho_det_state]), 3),
    " ; other = ", round(median(res[robust_12 == FALSE, prho_det_state]), 3),
    " ; Wilcoxon p = ", signif(w_det$p.value, 3))
say("median partial rho (cpm, control state): robust = ",
    round(median(res[robust_12 == TRUE, prho_cpm_state]), 3),
    " ; other = ", round(median(res[robust_12 == FALSE, prho_cpm_state]), 3),
    " ; Wilcoxon p = ", signif(w_cpm$p.value, 3))

hr("5. Direction check against the expected biology")
say("If interferon / NF-kB activity tracks CD207 loss, the partial rho should be NEGATIVE.")
say("Negative and BH-significant within state  -> within-state association.")
say("Negative raw but null after controlling for state -> the two layers are")
say("parallel descriptions of the same group contrast, not a within-state link.")
says(res[robust_12 == TRUE,
         .(tf,
           direction_attenuation = fifelse(
             sign(rho_det) != sign(prho_det_state), "SIGN FLIP",
             fifelse(abs(prho_det_state) < 0.5 * abs(rho_det), "ATTENUATED", "PERSISTS")),
           rho_det = round(rho_det, 2),
           prho_det_state = round(prho_det_state, 2))])

hr("6. Limits (must travel with any number quoted from this table)")
say("n = 29 pseudo-samples; 111 TFs in the BH family; no correction across families.")
say("TF activity is inferred from target-gene mRNA, not protein activity or DNA binding.")
say("Both CD207 readouts come from the same cells that produced the TF scores,")
say("so this is an internal association, not independent evidence.")
say("Correlation is not causation; no perturbation was performed.")
say("\nDONE")
