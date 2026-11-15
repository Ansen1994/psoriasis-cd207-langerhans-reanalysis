## 03_trend.R - offline Jonckheere-Terpstra ordered trend test (O1).
## Self-contained because clinfun / DescTools are NOT installed on this machine.
## Pure base R. No data access, no side effects when sourced.
## ASCII ONLY. Used by S2; Methods must state the permutation count and that
## the test is a self-implemented one-sided JT permutation test.

## JT statistic: sum over all ordered group pairs (lower state first) of
## I(x_j > x_i) + 0.5 * I(x_j == x_i). Ties get half credit.
jt_statistic <- function(x, g) {
  stopifnot(length(x) == length(g), length(x) >= 3L)
  d  <- outer(x, x, "-")   # d[i,j]  = x_i - x_j
  gm <- outer(g, g, "-")   # gm[i,j] = g_i - g_j
  sum((gm < 0) * ((d < 0) + 0.5 * (d == 0)))
}

## One-sided permutation p-value (alternative: increasing with state).
## B = number of permutations; +1 correction keeps p strictly > 0.
jt_perm_p <- function(x, g, B = 10000L, seed = 20261112L,
                      alternative = "greater") {
  set.seed(seed)
  obs <- jt_statistic(x, g)
  n <- length(x)
  null <- numeric(B)
  for (b in seq_len(B)) {
    gp <- sample(g, n, replace = FALSE)
    null[b] <- jt_statistic(x, gp)
  }
  if (alternative == "greater") {
    p <- (1 + sum(null >= obs)) / (B + 1)
  } else if (alternative == "less") {
    p <- (1 + sum(null <= obs)) / (B + 1)
  } else {
    stop("only one-sided alternatives are pre-registered")
  }
  list(statistic = obs, p_value = p, B = B, seed = seed,
       null_mean = mean(null), null_sd = stats::sd(null))
}

## Direction check required by TREND_REQUIRE_MONOTONIC: group medians must be
## non-decreasing (or non-increasing for 'less') across the ordered states.
jt_monotone <- function(x, g, alternative = "greater") {
  meds <- stats::aggregate(x, by = list(state = g), FUN = stats::median)$x
  if (alternative == "greater") {
    all(diff(meds) >= 0)
  } else {
    all(diff(meds) <= 0)
  }
}

## Run the test for many subclusters and apply BH-FDR across all of them.
## dat: data.frame with columns (feature, value, state_ord); state_ord is 0/1/2.
jt_trend_table <- function(dat, B = 10000L, seed = 20261112L,
                           alternative = "greater") {
  stopifnot(all(c("feature", "value", "state_ord") %in% colnames(dat)))
  feats <- unique(as.character(dat$feature))
  res <- do.call(rbind, lapply(feats, function(fm) {
    d <- dat[dat$feature == fm, ]
    tt <- jt_perm_p(d$value, d$state_ord, B = B,
                    seed = seed + match(fm, feats),
                    alternative = alternative)
    data.frame(feature = fm,
               n_samples = nrow(d),
               JT = tt$statistic,
               p = tt$p_value,
               monotone = jt_monotone(d$value, d$state_ord, alternative),
               stringsAsFactors = FALSE)
  }))
  res$FDR <- stats::p.adjust(res$p, method = "BH")
  res$pass <- res$FDR < 0.05 & res$monotone
  res[order(res$p), ]
}

## ---------------------------------------------------------------- self test
jt_self_test <- function(verbose = TRUE) {
  set.seed(1)
  ## Perfectly increasing: should give a very small one-sided p.
  g <- rep(0:2, each = 6L)
  x <- rep(c(1, 5, 9), each = 6L) + rnorm(18, 0, 0.01)
  up <- jt_perm_p(x, g, B = 2000L, alternative = "greater")
  ## Random noise: p should not be extreme.
  xr <- rnorm(18)
  rnd <- jt_perm_p(xr, g, B = 2000L, alternative = "greater")
  ## Decreasing data under 'greater' should give a large p.
  dn <- jt_perm_p(rev(x), g, B = 2000L, alternative = "greater")
  ok <- up$p_value < 0.01 && dn$p_value > 0.5 &&
        rnd$p_value > 0.01 && jt_monotone(x, g)
  if (verbose) {
    cat("JT self-test: increasing p=", signif(up$p_value, 4),
        " decreasing p=", signif(dn$p_value, 4),
        " random p=", signif(rnd$p_value, 4),
        " monotone=", jt_monotone(x, g),
        " => ", if (ok) "PASS" else "FAIL", "\n", sep = "")
  }
  invisible(ok)
}
