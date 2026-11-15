## S10_signature_confounding.R
## Reviewer-response analysis (added 2026-11-13). READ ONLY with respect to all
## frozen objects: it re-scores the GSE13355 matrix with additional gene sets and
## never overwrites S1-S4 or S9 tables.
##
## Questions addressed
##   Q1  Is the "10-gene LC/CD207 signature" separable from eccrine / mast /
##       epidermal content? (SLC18A2 and CFAP73, two of its members, are
##       canonical eccrine / ciliated-epithelium genes; CFAP73 is absent from
##       GPL570, so only 9 of 10 members are scored.)
##   Q2  Is a fall of this magnitude unusual, or does any expression-matched
##       random 10-gene set fall in lesional skin? (empirical null)
##   Q3  Does the group effect survive adjustment for keratinocyte, eccrine and
##       mast content?
##
## Outputs
##   results/tables/S10_signature_confounding.csv
##   results/tables/S10_null_distribution.csv
##   logs/S10_signature_confounding.log

source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({library(data.table); library(GSVA)})
set.seed(SEED)
logf <- file.path(LOG_DIR, "S10_signature_confounding.log")
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
cat("=== S10 signature-confounding analysis ===\n", file = logf)

## ---------------------------------------------------------------- expression
f <- file.path(RAW_DIR, "PSOR_BULK_GSE13355", BULK_SOURCE_FILE)
ln <- readLines(gzfile(f, open = "rt"), warn = FALSE)
b <- grep("^!series_matrix_table_begin", ln)[1]
e <- grep("^!series_matrix_table_end", ln)[1]
txt <- paste(ln[(b + 1):(e - 1)], collapse = "\n"); rm(ln); gc(FALSE)
ex <- fread(text = txt, sep = "\t", header = TRUE, check.names = FALSE, showProgress = FALSE)
rm(txt); gc(FALSE)
probe <- gsub('"', "", ex[[1]]); gsms <- gsub('"', "", colnames(ex)[-1])
mat <- as.matrix(ex[, -1, with = FALSE]); rm(ex); gc(FALSE)
storage.mode(mat) <- "double"; colnames(mat) <- gsms
annf <- normalizePath(file.path(RAW_DIR, "..", "..", "..", "data", "ref", "gpl570",
                                "GPL570-55999.txt"), mustWork = TRUE)
ann <- fread(annf, skip = 16, sep = "\t", header = TRUE, check.names = FALSE, quote = "",
             fill = TRUE, select = c("ID", "Gene Symbol"), showProgress = FALSE)
sym <- ann[["Gene Symbol"]][match(probe, ann$ID)]
sym <- sub(" ///.*$", "", sym); sym[is.na(sym) | sym %in% c("", "---", "NA")] <- NA_character_
idx <- !is.na(sym); cnt <- table(sym[idx])
agg <- rowsum(mat[idx, , drop = FALSE], group = sym[idx], reorder = TRUE)
agg <- agg / as.numeric(cnt[rownames(agg)]); rm(mat); gc(FALSE)
expr <- 2^agg; rm(agg); gc(FALSE)
say("gene-level matrix: ", nrow(expr), " symbols x ", ncol(expr), " arrays")

sm <- fread(file.path(TAB_DIR, "S4_module_scores.csv"))
smap <- unique(sm[, .(sample_id, group, donor_id)])

## ---------------------------------------------------------------- gene sets
LC10 <- c("FNDC4", "SLC18A2", "OLFM2", "GTF2IRD1", "CD207",
          "ACOT7", "LIMS2", "CFAP73", "HPGDS", "HLA-DQB2")
S <- list(
  LC10_original   = LC10,
  LC9_platform     = setdiff(LC10, "CFAP73"),
  LC7_no_eccrine   = setdiff(LC9 <- setdiff(LC10, "CFAP73"),
                             c("SLC18A2", "HPGDS")),
  LC4_minimal      = c("CD207", "HLA-DQB2", "ACOT7", "FNDC4"),
  Eccrine_control  = c("DCD", "SCGB2A2", "AQP5", "KRT19", "SLC18A2"),
  Eccrine_noSLC    = c("DCD", "SCGB2A2", "AQP5", "KRT19"),
  Mast_control     = c("CPA3", "TPSAB1", "TPSB2", "KIT", "MS4A2"),
  KRT_differentiation = c("KRT1", "KRT10", "KRT14", "KRT16", "KRT6A", "KRT17"),
  Melanocyte_control  = c("PMEL", "TYR", "MLANA", "DCT", "S100B")
)
mapped <- vapply(S, function(g) length(intersect(g, rownames(expr))), 1L)
say("\n--- mapping of each set onto GPL570 ---")
for (n in names(S)) say(sprintf("  %-20s input %2d  mapped %2d", n, length(S[[n]]), mapped[[n]]))
S <- lapply(S, function(g) intersect(unique(g), rownames(expr)))
S <- S[mapped >= 4]

## Q2: expression-matched random sets
mu <- rowMeans(expr)
lc_mu <- mu[S$LC9_platform]
cand <- setdiff(names(mu)[mu >= min(lc_mu) * 0.7 & mu <= max(lc_mu) * 1.3],
                unlist(S, use.names = FALSE))
say("\nexpression-matched candidate pool: ", length(cand), " genes")
NSET <- 200
rand <- replicate(NSET, sample(cand, 10), simplify = FALSE)
names(rand) <- sprintf("RND%03d", seq_len(NSET))

all_sets <- c(S, rand)
par <- GSVA::ssgseaParam(exprData = expr, geneSets = all_sets, minSize = 4, maxSize = 500,
                         alpha = 0.25, normalize = TRUE)
sc <- as.matrix(GSVA::gsva(par))
long <- as.data.table(as.table(sc)); setnames(long, c("gene_set", "sample_id", "score"))
long <- merge(long, smap, by = "sample_id")

## Q1: group and paired tests
grp <- long[group %in% c("NN", "PP"), {
  a <- score[group == "NN"]; p <- score[group == "PP"]
  .(n_NN = length(a), n_PP = length(p),
    median_NN = median(a), median_PP = median(p),
    delta_PP_minus_NN = median(p) - median(a),
    p_wilcox = suppressWarnings(wilcox.test(p, a)$p.value))
}, by = gene_set]
grp[, FDR_BH := p.adjust(p_wilcox, "BH")]

paired <- long[group %in% c("PN", "PP"), {
  w <- dcast(.SD, donor_id ~ group, value.var = "score")
  w <- w[!is.na(PN) & !is.na(PP)]
  .(n_pairs = nrow(w), median_diff_PP_minus_PN = median(w$PP - w$PN),
    p_paired = suppressWarnings(wilcox.test(w$PP, w$PN, paired = TRUE)$p.value))
}, by = gene_set]
grp <- merge(grp, paired, by = "gene_set", all.x = TRUE)

## Q3: adjustment models on the NN/PP subset
wide <- dcast(long, sample_id ~ gene_set, value.var = "score")
wide <- merge(wide, smap, by = "sample_id")[group %in% c("NN", "PP")]
wide[, grp := factor(group, levels = c("NN", "PP"))]
adj <- data.table()
for (lhs in c("LC9_platform", "LC7_no_eccrine", "LC4_minimal")) {
  forms <- list(
    unadjusted = as.formula(sprintf("%s ~ grp", lhs)),
    plus_KRT   = as.formula(sprintf("%s ~ grp + KRT_differentiation", lhs)),
    plus_ECC   = as.formula(sprintf("%s ~ grp + Eccrine_control", lhs)),
    plus_KRT_ECC = as.formula(sprintf("%s ~ grp + KRT_differentiation + Eccrine_control", lhs)),
    plus_all   = as.formula(sprintf("%s ~ grp + KRT_differentiation + Eccrine_control + Mast_control", lhs))
  )
  for (nm in names(forms)) {
    fit <- lm(forms[[nm]], data = wide); cf <- summary(fit)$coef
    adj <- rbind(adj, data.table(outcome = lhs, model = nm, n = nrow(wide),
                                 beta_group = cf[2, 1], se_group = cf[2, 2],
                                 p_group = cf[2, 4], r_squared = summary(fit)$r.squared))
  }
}
say("\n--- Q1/Q3: group effect before and after covariate adjustment (NN vs PP) ---")
print(adj[order(outcome, model)])

## Q2 summary
rnd <- grp[grepl("^RND", gene_set), delta_PP_minus_NN]
nulltb <- data.table(
  statistic = c("n_random_sets", "median_delta_random", "q025_delta_random",
                "q975_delta_random",
                "delta_LC9", "percentile_LC9", "empirical_p_LC9",
                "delta_KRT", "percentile_KRT",
                "delta_Eccrine", "percentile_Eccrine",
                "delta_Mast", "percentile_Mast"),
  value = c(NSET, median(rnd), quantile(rnd, 0.025), quantile(rnd, 0.975),
            grp[gene_set == "LC9_platform", delta_PP_minus_NN],
            100 * mean(rnd <= grp[gene_set == "LC9_platform", delta_PP_minus_NN]),
            mean(rnd <= grp[gene_set == "LC9_platform", delta_PP_minus_NN]),
            grp[gene_set == "KRT_differentiation", delta_PP_minus_NN],
            100 * mean(rnd <= grp[gene_set == "KRT_differentiation", delta_PP_minus_NN]),
            grp[gene_set == "Eccrine_control", delta_PP_minus_NN],
            100 * mean(rnd <= grp[gene_set == "Eccrine_control", delta_PP_minus_NN]),
            grp[gene_set == "Mast_control", delta_PP_minus_NN],
            100 * mean(rnd <= grp[gene_set == "Mast_control", delta_PP_minus_NN])))
say("\n--- Q2: empirical null of PP-NN median delta (200 expression-matched 10-gene sets) ---")
print(nulltb)

fwrite(grp[!grepl("^RND", gene_set)][order(p_wilcox)],
       file.path(TAB_DIR, "S10_signature_confounding.csv"))
fwrite(grp[grepl("^RND", gene_set), .(gene_set, delta_PP_minus_NN, p_wilcox)],
       file.path(TAB_DIR, "S10_null_raw.csv"))
fwrite(nulltb, file.path(TAB_DIR, "S10_null_distribution.csv"))
say("\nWROTE S10_signature_confounding.csv and S10_null_distribution.csv")
say("DONE")
