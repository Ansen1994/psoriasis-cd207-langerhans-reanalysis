## S9_prop_zero.R - structural-zero correction for the CD207+ proportion analysis.
##
## WHY THIS SCRIPT EXISTS
##   The frozen pipeline (R/05b_s2_stats.R) built results/tables/S2_proportions_by_sample.csv
##   as a LONG table that only contained pseudo-sample x subcluster rows with >= 1 cell.
##   A pseudo-sample with ZERO cells in a subcluster therefore disappeared instead of
##   contributing a proportion of 0. That is a data-handling bug, not a design choice:
##   the denominator (n_myeloid) is defined for every pseudo-sample, so 0/n_myeloid = 0
##   is a valid and necessary observation.
##
##   This script re-derives the complete 29 x 3 grid with zeros filled, and re-runs the
##   SAME frozen tests (one-sided Jonckheere-Terpstra + lmer) on it. No threshold, no
##   alternative and no model formula is changed. A negative-binomial GLMM on counts is
##   added as an explicitly POST HOC sensitivity analysis, not as a replacement test.
##
## INPUTS  (read only, frozen)
##   results/tables/S2_proportions_by_sample.csv
## OUTPUTS (new, S9_ prefix; nothing frozen is overwritten)
##   results/tables/S9_prop_grid_zero_filled.csv
##   results/tables/S9_prop_group_summary.csv
##   results/tables/S9_prop_trend_zero_filled.csv
##   results/tables/S9_prop_lmer_zero_filled.csv
##   results/tables/S9_prop_zero_counts.csv
##   results/tables/S9_prop_count_glmm.csv
##   logs/S9_prop_zero.log
## ASCII ONLY.

source("config/00_paths.R")
source("config/01_params.R")
source("R/03_trend.R")

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

set.seed(SEED)

TAB <- "results/tables"
LOG <- file.path("logs", "S9_prop_zero.log")
con <- file(LOG, open = "wt")
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  writeLines(txt, con)
}

say("=== S9 structural-zero correction ==========================================")
say("started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("frozen params untouched; seed = ", SEED)

## ---------------------------------------------------------------- 1. full grid
d <- fread(file.path(TAB, "S2_proportions_by_sample.csv"))
say("input rows (as published): ", nrow(d))

base <- unique(d[, .(pseudo_id, donor_base, group, state_ord, batch, n_myeloid)])
say("pseudo-samples in the long table: ", nrow(base))

subs <- sort(unique(d$subcluster))
grid <- CJ(pseudo_id = base$pseudo_id, subcluster = subs, sorted = FALSE)
grid <- merge(grid, base, by = "pseudo_id", all.x = TRUE)
grid <- merge(grid,
              d[, .(pseudo_id, subcluster, n_cells)],
              by = c("pseudo_id", "subcluster"), all.x = TRUE)
setorder(grid, subcluster, pseudo_id)
grid[is.na(n_cells), n_cells := 0L]
grid[, prop := 100 * n_cells / n_myeloid]
grid[, zero_filled := TRUE]
grid[, state_ord := as.integer(state_ord)]
grid[, donor := donor_base]

say("complete grid rows: ", nrow(grid), " (expected ", nrow(base) * length(subs), ")")
say("rows that were MISSING and are now 0: ", nrow(grid) - nrow(d))
print(grid[!(paste(pseudo_id, subcluster) %in% paste(d$pseudo_id, d$subcluster)),
           .(pseudo_id, subcluster, group, n_cells, n_myeloid)],
      nrows = 20)
say("zeros per subcluster x group:")
print(grid[n_cells == 0, .N, by = .(subcluster, group)][order(subcluster, group)])
fwrite(grid, file.path(TAB, "S9_prop_grid_zero_filled.csv"))

## zero-rate by group for the CD207+ subcluster
z <- grid[subcluster == 2, .(n_pseudo = .N, n_zero = sum(n_cells == 0)), by = group]
z[, pct_zero := round(100 * n_zero / n_pseudo, 1)]
say("CD207+ subcluster, pseudo-samples with ZERO cells:")
print(z)
fwrite(z, file.path(TAB, "S9_prop_zero_counts.csv"))
ft <- try(fisher.test(matrix(c(z[group == "NS", n_zero], z[group == "NS", n_pseudo] - z[group == "NS", n_zero],
                               z[group == "PP", n_zero], z[group == "PP", n_pseudo] - z[group == "PP", n_zero]),
                             nrow = 2, byrow = TRUE)), silent = TRUE)
if (!inherits(ft, "try-error")) say("Fisher exact NS vs PP zero-rate p = ", signif(ft$p.value, 4))

## ---------------------------------------------------------------- 2. group summary
s <- grid[, .(n_pseudo = .N,
              n_nonzero = sum(n_cells > 0),
              n_zero = sum(n_cells == 0),
              median_prop = median(prop),
              q1 = quantile(prop, 0.25),
              q3 = quantile(prop, 0.75),
              mean_prop = mean(prop),
              min_prop = min(prop),
              max_prop = max(prop)),
          by = .(subcluster, group)]
setorder(s, subcluster, group)
say("\n=== group summary on the ZERO-FILLED grid (percent of QC-eligible myeloid cells) ===")
print(s)
fwrite(s, file.path(TAB, "S9_prop_group_summary.csv"))

## ---------------------------------------------------------------- 3. frozen JT test
tr <- as.data.table(jt_trend_table(as.data.frame(grid[, .(feature = paste0("sub_", subcluster),
                                                         value = prop,
                                                         state_ord = state_ord)]),
                                   B = TREND_NPERM, seed = SEED,
                                   alternative = TREND_ALTERNATIVE))
say("\n=== pre-registered one-sided JT, zero-filled (same alternative, same B) ===")
print(tr)
tr[, data := "zero_filled_29x3"]
fwrite(tr, file.path(TAB, "S9_prop_trend_zero_filled.csv"))

tr_less <- as.data.table(jt_trend_table(as.data.frame(grid[, .(feature = paste0("sub_", subcluster),
                                                                value = prop,
                                                                state_ord = state_ord)]),
                                        B = TREND_NPERM, seed = SEED, alternative = "less"))
say("POST HOC opposite-direction JT (alternative = 'less'), for description only:")
print(tr_less[, .(feature, p, FDR)])
fwrite(tr_less[, .(feature, p, FDR, test = "post_hoc_less")],
       file.path(TAB, "S9_prop_trend_posthoc_less.csv"))

## ------------------------------------------------- 3b. denominator sensitivity
## Pseudo-samples with very few myeloid cells contribute a 0% that is almost
## guaranteed by chance (e.g. 0/3 with a 10% rate has probability ~0.73).
## Re-run the same JT after removing pseudo-samples with n_myeloid < 10.
thr <- 10L
sub_grid <- grid[n_myeloid >= thr]
say("\n=== SENSITIVITY: JT after removing pseudo-samples with n_myeloid < ", thr, " ===")
say("pseudo-samples retained: ", uniqueN(sub_grid$pseudo_id), " of ", uniqueN(grid$pseudo_id))
print(sub_grid[, .N, by = .(subcluster, group)][order(subcluster, group)])
tr_sub <- as.data.table(jt_trend_table(as.data.frame(sub_grid[, .(feature = paste0("sub_", subcluster),
                                                                  value = prop,
                                                                  state_ord = state_ord)]),
                                       B = TREND_NPERM, seed = SEED,
                                       alternative = TREND_ALTERNATIVE))
print(tr_sub)
tr_sub[, data := paste0("zero_filled_nMyeloid_ge_", thr)]
fwrite(tr_sub, file.path(TAB, "S9_prop_trend_sensitivity_min10.csv"))

## ---------------------------------------------------------------- 4. lmer, identical formula
lm_rows <- list()
for (tag in c("all", "min10")) {
  gd <- if (tag == "all") grid else sub_grid
  for (sc in subs) {
    x <- gd[subcluster == sc]
    f1 <- try(lmer(prop ~ state_ord + (1 | donor), data = x, REML = FALSE), silent = TRUE)
    f0 <- try(lmer(prop ~ 1 + (1 | donor), data = x, REML = FALSE), silent = TRUE)
    if (inherits(f1, "try-error") || inherits(f0, "try-error") || nrow(x) < 8) {
      lm_rows[[paste(tag, sc)]] <- data.table(dataset = tag, subcluster = sc, beta = NA_real_,
                                              lrt_chisq = NA_real_, lmer_p = NA_real_,
                                              singular = NA, n_pseudo = nrow(x))
      next
    }
    an <- anova(f1, f0)
    lm_rows[[paste(tag, sc)]] <- data.table(
      dataset = tag, subcluster = sc,
      beta = unname(lme4::fixef(f1)[2]),
      lrt_chisq = an$Chisq[2],
      lmer_p = an$`Pr(>Chisq)`[2],
      singular = isSingular(f1),
      n_pseudo = nrow(x)
    )
  }
}
lmt <- rbindlist(lm_rows)
say("\n=== exploratory lmer on the zero-filled grid (same formula) ===")
print(lmt)
fwrite(lmt, file.path(TAB, "S9_prop_lmer_zero_filled.csv"))

## ---------------------------------------------------------------- 5. post hoc count model
## Rationale: the per-sample myeloid denominator differs systematically by group
## (median n_myeloid 17 / 33.5 / 72.5 in NS / PN / PP), so a proportion and a
## rate are not interchangeable summaries. A count model with offset(log n_myeloid)
## estimates the rate directly. Poisson is used only where the Pearson dispersion
## is close to 1; over-dispersed fits are flagged and NOT interpreted.
## glmer.nb was tried first and is reported only if it converged; it did not
## (theta hit the boundary at ~651 with an iteration-limit warning), so the
## negative-binomial p-values are NOT used anywhere.
say("\n=== POST HOC sensitivity: Poisson GLMM on counts ===")
say("model: n_cells ~ state_ord + offset(log(n_myeloid)) + (1|donor)")
nb_rows <- list()
for (tag in c("all", "min10")) {
  gd <- if (tag == "all") grid else sub_grid
  for (sc in subs) {
    x <- as.data.frame(gd[subcluster == sc])
    if (nrow(x) < 8) next
    p1 <- try(suppressWarnings(glmer(n_cells ~ state_ord + offset(log(n_myeloid)) + (1 | donor),
                                     family = poisson, data = x)), silent = TRUE)
    p0 <- try(suppressWarnings(glmer(n_cells ~ 1 + offset(log(n_myeloid)) + (1 | donor),
                                     family = poisson, data = x)), silent = TRUE)
    if (inherits(p1, "try-error") || inherits(p0, "try-error")) next
    disp <- sum(residuals(p1, type = "pearson")^2) / df.residual(p1)
    an <- anova(p1, p0)
    b <- unname(lme4::fixef(p1)[2]); se <- sqrt(vcov(p1)[2, 2])
    nb_rows[[paste(tag, sc)]] <- data.table(
      dataset = tag, subcluster = sc, model = "poisson_glmm",
      beta = b, se = se,
      rate_ratio = exp(b),
      rr_lo = exp(b - 1.96 * se), rr_hi = exp(b + 1.96 * se),
      lrt_p = an$`Pr(>Chisq)`[2],
      dispersion = disp,
      dispersion_ok = disp < 1.5,
      singular = isSingular(p1),
      n_pseudo = nrow(x))
  }
}
nbt <- rbindlist(nb_rows, fill = TRUE)
say("rows with dispersion >= 1.5 are over-dispersed: their Poisson p is anti-conservative and is NOT interpreted.")
print(nbt)
fwrite(nbt, file.path(TAB, "S9_prop_count_glmm.csv"))

say("\n=== INTERPRETATION NOTES (for the manuscript) ===")
say("1. The zero-filled medians REPLACE the previously reported 79.59 / 48.72 / 6.95.")
say("2. The frozen one-sided 'greater' alternative is unchanged and still fails.")
say("3. The opposite-direction result remains exploratory in both the lmer and the GLMM.")
say("4. Nothing here touches the bulk analyses.")
say("ended: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
close(con)
cat("\nOUTPUT: results/tables/S9_prop_*.csv\n")
