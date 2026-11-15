#!/usr/bin/env Rscript
## S24_jt_crosscheck.R - independent verification of the self-implemented
## Jonckheere-Terpstra ordered-trend test used in 03_trend.R / S9_prop_zero.R.
##
## clinfun and DescTools are not installed on this machine, so the test could
## not be compared against a package implementation. Four checks that do not
## depend on any package are run instead:
##
##   A. statistic     two implementations written differently (matrix outer()
##                    product vs explicit nested loops over group pairs) must
##                    return the same value on real and synthetic data;
##   B. identities    J(x) + J(-x) must equal the number of between-group pairs,
##                    and under permutation E[J] must equal half that number;
##   C. known values  three hand-computable examples with published answers;
##   D. determinism   the stored p-value must be reproduced exactly from the
##                    frozen seed, and must vary only within Monte-Carlo error
##                    across other seeds / larger B.
##
## ASCII only. Read-only: reads result tables, writes one log. No new analysis.

suppressPackageStartupMessages(library(data.table))
source("R/03_trend.R")

TAB <- file.path("results", "tables")
LOG <- file.path("logs", "S24_jt_crosscheck.log")
fails <- 0L
lines <- character(0)

say <- function(...) {
  msg <- paste0(...)
  lines <<- c(lines, msg)
  cat(msg, "\n", sep = "")
}
ok <- function(pass, what, detail = "") {
  if (!isTRUE(pass)) fails <<- fails + 1L
  say(sprintf("[%s] %-46s %s", if (isTRUE(pass)) "PASS" else "FAIL", what, detail))
}

## ---------------------------------------------------------------- A. stats ---
## Same definition as jt_statistic(), but written from the definition instead of
## with outer(): for every ordered pair of states, loop over the observations.
jt_loop <- function(x, g) {
  states <- sort(unique(g))
  total <- 0
  for (a in seq_along(states)) {
    for (b in seq_along(states)) {
      if (b <= a) next
      lo <- x[g == states[a]]
      hi <- x[g == states[b]]
      for (u in lo) {
        for (v in hi) {
          if (v > u) {
            total <- total + 1
          } else if (v == u) {
            total <- total + 0.5
          }
        }
      }
    }
  }
  total
}

n_between <- function(g) {
  tab <- table(g)
  s <- 0
  for (a in seq_along(tab)) {
    for (b in seq_along(tab)) if (b > a) s <- s + tab[a] * tab[b]
  }
  as.numeric(s)
}

say("=== S24 JT cross-check ===")
say("R ", getRversion(), "; clinfun/DescTools not installed -> independent re-implementation")
say("")

## ------------------------------------------------------------ C. examples ---
say("--- C. hand-computable examples ---")
ex <- list(
  list(x = c(1, 2, 3, 4, 5, 6), g = c(1, 1, 2, 2, 3, 3), want = 12,
       what = "strictly increasing, 3 groups of 2"),
  list(x = c(6, 5, 4, 3, 2, 1), g = c(1, 1, 2, 2, 3, 3), want = 0,
       what = "strictly decreasing, 3 groups of 2"),
  list(x = c(1, 1, 1, 1, 1, 1), g = c(1, 1, 2, 2, 3, 3), want = 6,
       what = "all tied (every pair gets 0.5)"),
  list(x = c(1, 1, 1, 2), g = c(1, 1, 1, 2), want = 3,
       what = "one tie group vs one higher group")
)
for (e in ex) {
  a <- jt_statistic(e$x, e$g)
  b <- jt_loop(e$x, e$g)
  ok(a == e$want && b == e$want, e$what,
     sprintf("J=%.1f loop=%.1f expected=%.0f", a, b, e$want))
}

## ------------------------------------------------------------- B. algebra ---
say("")
say("--- B. algebraic identities ---")
set.seed(11)
syn <- data.table(x = round(rnorm(60), 1), g = rep(1:3, each = 20))
jb <- n_between(syn$g)
ok(jt_statistic(syn$x, syn$g) + jt_statistic(-syn$x, syn$g) == jb,
   "J(x) + J(-x) = number of between-group pairs",
   sprintf("sum=%.1f pairs=%.0f", jt_statistic(syn$x, syn$g) +
             jt_statistic(-syn$x, syn$g), jb))
pm <- jt_perm_p(syn$x, syn$g, B = 4000L, seed = 7L)
## Under random labelling each between-group pair contributes 1/2 in expectation
## (1/2 chance the larger value sits in the higher group; ties always 1/2).
ok(abs(pm$null_mean - jb / 2) < 0.02 * jb,
   "permutation null mean = pairs / 2", sprintf("mean=%.2f expected=%.1f",
                                                pm$null_mean, jb / 2))

## ------------------------------------------------------------- A. real data ---
say("")
say("--- A. real data: 29 x 3 proportion grid ---")
grid <- fread(file.path(TAB, "S9_prop_grid_zero_filled.csv"))
stored <- fread(file.path(TAB, "S9_prop_trend_zero_filled.csv"))
say(sprintf("grid rows: %d (subclusters %s)", nrow(grid),
            paste(sort(unique(grid$subcluster)), collapse = ", ")))

for (k in sort(unique(grid$subcluster))) {
  d <- grid[subcluster == k]
  a <- jt_statistic(d$prop, d$state_ord)
  b <- jt_loop(d$prop, d$state_ord)
  s <- stored[feature == paste0("sub_", k), JT]
  ok(a == b && length(s) == 1L && a == s,
     sprintf("subcluster %d: outer == loop == stored", k),
     sprintf("outer=%.1f loop=%.1f stored=%.1f n=%d", a, b, s, nrow(d)))
}

## --------------------------------------------------------- D. determinism ---
say("")
say("--- D. permutation determinism and Monte-Carlo error ---")
## jt_trend_table() offsets the seed per feature (seed + match(feature, feats)),
## so the seed that produced sub_3 is SEED + 3, not SEED. Reproducing the stored
## number therefore also confirms the seed convention is what the log claims.
feats <- unique(paste0("sub_", grid$subcluster))
SEED <- 20261112L
prim <- grid[subcluster == 3]
seed3 <- SEED + match("sub_3", feats)
p_frozen <- jt_perm_p(prim$prop, prim$state_ord, B = 10000L, seed = seed3)
p_s1 <- jt_perm_p(prim$prop, prim$state_ord, B = 10000L, seed = 1L)
p_s2 <- jt_perm_p(prim$prop, prim$state_ord, B = 10000L, seed = 2L)
p_big <- jt_perm_p(prim$prop, prim$state_ord, B = 200000L, seed = 3L)
stored_p <- stored[feature == "sub_3", p]
ok(isTRUE(all.equal(p_frozen$p_value, stored_p, tolerance = 1e-12)),
   "seed SEED+n reproduces the published p exactly",
   sprintf("seed=%d recomputed=%.8f stored=%.8f", seed3, p_frozen$p_value, stored_p))
spread <- max(abs(c(p_s1$p_value, p_s2$p_value) - p_big$p_value))
ok(spread < 0.01,
   "p varies <0.01 across seeds / B = 200k",
   sprintf("seed1=%.4f seed2=%.4f B200k=%.4f max|diff|=%.4f",
           p_s1$p_value, p_s2$p_value, p_big$p_value, spread))

## ---- donor exchangeability diagnostic (NOT a correction of the frozen test) --
## The pre-registered null permutes state labels freely across all 29 pseudo-
## samples. Ten donors contribute both a non-lesional and a lesional biopsy, so
## those two samples are not independent. Exchanging state labels only WITHIN a
## donor gives the null that respects that pairing - but it also changes the
## question, because the 7 healthy donors and the 9 single-state donors have
## their labels frozen, leaving only the 10 paired donors free to move. The
## number below is therefore a paired-style contrast (PN vs PP within donor),
## not a three-state ordered trend, and it is reported for completeness only.
## The manuscript answers this objection with the donor-clustered mixed model in
## Section 2.10, which is the design-appropriate analysis.
say("")
say("--- D2. donor exchangeability diagnostic (different estimand) ---")
donors <- unique(prim$donor)
p_within <- local({
  set.seed(seed3)
  obs <- jt_statistic(prim$prop, prim$state_ord)
  B <- 10000L
  null <- numeric(B)
  for (b in seq_len(B)) {
    gp <- prim$state_ord
    for (dn in donors) {
      idx <- which(prim$donor == dn)
      if (length(idx) > 1L) gp[idx] <- sample(gp[idx])
    }
    null[b] <- jt_statistic(prim$prop, gp)
  }
  (1 + sum(null >= obs)) / (B + 1)
})
say(sprintf("       free permutation (frozen, 3 states)        p = %.4f", p_frozen$p_value))
say(sprintf("       within-donor exchange (paired PN vs PP)    p = %.4f", p_within))
say("       These test different nulls: the second freezes every single-state")
say("       donor, so it estimates a paired PN-vs-PP contrast with n = 10 pairs")
say("       and is not a re-test of the ordered trend. Not entered into the paper.")

## ------------------------------------------------------------- sensitivity ---
say("")
say("--- E. sensitivity grid reproduced (>= 10 myeloid cells) ---")
min10 <- grid[n_myeloid >= 10]
sens <- fread(file.path(TAB, "S9_prop_trend_sensitivity_min10.csv"))
for (k in sort(unique(min10$subcluster))) {
  d <- min10[subcluster == k]
  a <- jt_statistic(d$prop, d$state_ord)
  s <- sens[feature == paste0("sub_", k), JT]
  ok(a == s, sprintf("min10 subcluster %d matches stored", k),
     sprintf("recomputed=%.1f stored=%.1f n=%d", a, s, nrow(d)))
}

## ------------------------------------------------------------------ close ---
say("")
say(sprintf("RESULT: %s (%d failure%s)", if (fails == 0L) "ALL PASS" else "FAIL",
            fails, if (fails == 1L) "" else "s"))
writeLines(lines, LOG)
