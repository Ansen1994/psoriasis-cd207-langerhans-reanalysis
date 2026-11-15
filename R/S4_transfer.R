## S4_transfer.R - independent single-cell cohort (GSE162183) label transfer.
## Frozen route: Seurat label transfer, NO re-integration, NO Harmony (O7).
## Gate D criterion GATE_D_MIN_AUC = 0.70.
##
## AUC definition (pre-specified here, before looking at any result):
##   5-fold stratified cross-validation WITHIN the discovery reference.
##   For each subcluster, one-vs-rest ROC AUC using TransferData's
##   prediction.score.<label> as the continuous score, with pROC direction
##   fixed to "<" (cases are expected to score HIGHER than controls).
##   Macro mean over QC-eligible subclusters is the Gate D statistic.
##   This is an INTERNAL reproducibility check of the transfer procedure,
##   NOT external validation. External evidence is the direction of the
##   transferred labels in the independent cohort and the bulk scores.
##   The one-sided direction is fixed in advance and never flipped.
##
## The query is whole skin (not myeloid-sorted), so every cell receives a
## label. A prediction.score.max filter with three pre-declared thresholds
## (0.3 / 0.5 / 0.7) is reported as a sensitivity analysis; 0.5 is the main
## majority-vote cut. This limitation must be stated in Methods.
##
## ASCII only (LC_CTYPE=C). All outputs stay inside the project.
source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix); library(pROC)
})

logf <- file.path(LOG_DIR, "S4_transfer.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S4 independent scRNA cohort: label transfer ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
cat("GATE_D_MIN_AUC = ", GATE_D_MIN_AUC, " | method = ", VALIDATION_METHOD, "\n", sep = "")
set.seed(SEED)
t0 <- Sys.time()

## ---------------------------------------------------------------- 1. query data
f <- file.path(RAW_DIR, "PSOR_SC_GSE162183",
               "GSE162183_Raw_gene_counts_matrix.tab.gz")
cat("\nreading ", f, "\n", sep = "")
d <- data.table::fread(f, sep = "\t", header = TRUE, showProgress = FALSE,
                       nThread = 4, colClasses = list(character = 1))
gn <- d[[1]]
bc <- colnames(d)[-1]
m <- as.matrix(d[, -1, with = FALSE])
rm(d); gc(FALSE)
sm <- Matrix::Matrix(m, sparse = TRUE)
rm(m); gc(FALSE)
rownames(sm) <- make.unique(gn)
colnames(sm) <- bc
cat("raw validation matrix: ", nrow(sm), " genes x ", ncol(sm), " cells\n", sep = "")
cat("non-zero entries: ", length(sm@x), "\n", sep = "")

q <- CreateSeuratObject(counts = sm, project = "GSE162183",
                        min.cells = 3, min.features = 200)
rm(sm); gc(FALSE)
sid <- sub("_.*$", "", colnames(q))
q$sample_id <- sid
q$group2 <- ifelse(startsWith(sid, "Ctrl"), "NS", "PP")
q$condition_raw <- sid
mt_feats <- intersect(MT_SYMBOL_WHITELIST, rownames(q))
if (length(mt_feats) > 0) {
  q[["percent_mt"]] <- PercentageFeatureSet(q, features = mt_feats)
} else {
  q[["percent_mt"]] <- 0
  cat("WARNING: no MT whitelist genes in query; percent_mt set to 0 (no MT filter applied).\n")
}
cat("samples from barcode prefix:\n"); print(table(q$sample_id, q$group2))
cat("cells before QC: ", ncol(q), "\n", sep = "")
q <- subset(q, subset = nFeature_RNA >= QC_NFEATURE_MIN &
                   nFeature_RNA <= QC_NFEATURE_MAX &
                   percent_mt < QC_PCT_MT_MAX)
cat("cells after QC:  ", ncol(q), "\n", sep = "")
cat("MT whitelist genes found: ", length(mt_feats), "\n", sep = "")

q <- NormalizeData(q, verbose = FALSE)
q <- FindVariableFeatures(q, nfeatures = 2000, verbose = FALSE)
q <- ScaleData(q, verbose = FALSE)
q <- RunPCA(q, npcs = 30, verbose = FALSE)
cat("query PCA computed (30 PCs).\n")

## ---------------------------------------------------------------- 2. reference
qc <- fread(file.path(TAB_DIR, "S2_subcluster_qc.csv"))
elig <- qc[qc_eligible == TRUE, as.character(subcluster)]
o <- readRDS(file.path(PROC_DIR, F_S2_OBJECT))
ref <- subset(o, cells = colnames(o)[o$subcluster %in% elig])
rm(o); gc(FALSE)
cat("\nreference: ", ncol(ref), " QC-eligible myeloid cells; labels: ",
    paste(sort(unique(ref$subcluster)), collapse = ","), "\n", sep = "")
print(table(ref$subcluster))
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
ref <- ScaleData(ref, verbose = FALSE)
ref <- RunPCA(ref, npcs = 30, verbose = FALSE)
cat("reference PCA computed (30 PCs). No Harmony, no integration.\n")

## ---------------------------------------------------------------- 3. CV AUC
## Fixed one-sided direction: cases (cells OF this label) score higher.
labels <- sort(unique(ref$subcluster))
folds <- integer(ncol(ref))
for (cl in labels) {
  ix <- which(ref$subcluster == cl)
  folds[ix] <- sample(rep(1:5, length.out = length(ix)))
}
acc <- list()
for (k in 1:5) {
  tr <- subset(ref, cells = colnames(ref)[folds != k])
  te <- subset(ref, cells = colnames(ref)[folds == k])
  a <- tryCatch(
    FindTransferAnchors(reference = tr, query = te, dims = 1:30,
                        reference.reduction = "pca", verbose = FALSE),
    error = function(e) { cat("fold ", k, " anchors FAILED: ", conditionMessage(e), "\n", sep = ""); NULL })
  if (is.null(a)) next
  p <- tryCatch(TransferData(anchorset = a, refdata = tr$subcluster,
                             dims = 1:30),
                error = function(e) { cat("fold ", k, " transfer FAILED: ", conditionMessage(e), "\n", sep = ""); NULL })
  if (is.null(p)) next
  truth <- te$subcluster[match(rownames(p), colnames(te))]
  for (cl in labels) {
    cn <- paste0("prediction.score.", cl)
    if (!cn %in% colnames(p)) next
    resp <- as.integer(truth == cl)
    if (length(unique(resp)) < 2) next
    z <- tryCatch(as.numeric(pROC::auc(pROC::roc(
      response = resp, predictor = p[[cn]], levels = c(0, 1),
      direction = "<", quiet = TRUE))), error = function(e) NA_real_)
    acc[[length(acc) + 1]] <- data.table(fold = k, subcluster = cl, auc = z,
                                         n_pos = sum(resp), n_neg = sum(1 - resp))
  }
  cat("fold ", k, " done: ", nrow(p), " held-out cells\n", sep = "")
  rm(a, p, tr, te); gc(FALSE)
}
cv <- if (length(acc)) rbindlist(acc) else
  data.table(fold = integer(), subcluster = character(), auc = numeric(),
             n_pos = integer(), n_neg = integer())
cv_sum <- cv[, .(auc_mean = mean(auc, na.rm = TRUE),
                 auc_min = min(auc, na.rm = TRUE),
                 auc_max = max(auc, na.rm = TRUE), n_folds = .N), by = subcluster]
setorder(cv_sum, subcluster)
MACRO <- if (nrow(cv_sum)) mean(cv_sum$auc_mean, na.rm = TRUE) else NA_real_
cat("\n--- per-subcluster CV AUC (one-sided direction fixed) ---\n")
print(cv_sum)
cat("MACRO_AUC = ", round(MACRO, 4), " | threshold = ", GATE_D_MIN_AUC, "\n", sep = "")
auc_tab <- copy(cv_sum)
auc_tab[, statistic := "per_subcluster"]
macro_row <- data.table(subcluster = "MACRO", auc_mean = MACRO,
                        auc_min = NA_real_, auc_max = NA_real_,
                        n_folds = nrow(cv_sum), statistic = "gateD_primary")
auc_tab <- rbind(auc_tab, macro_row, fill = TRUE)
fwrite(auc_tab, file.path(TAB_DIR, T_S4_AUC))
fwrite(cv, file.path(TAB_DIR, "S4_validation_auc_folds.csv"))

## ---------------------------------------------------------------- 4. query transfer
cat("\n--- transfer to independent cohort ---\n")
a <- tryCatch(
  FindTransferAnchors(reference = ref, query = q, dims = 1:30,
                      reference.reduction = "pca", verbose = FALSE),
  error = function(e) { cat("query anchors FAILED: ", conditionMessage(e), "\n", sep = ""); NULL })
if (!is.null(a)) {
  p <- TransferData(anchorset = a, refdata = ref$subcluster,
                    dims = 1:30)
  d <- data.table(cell = rownames(p),
                  sample_id = q$sample_id[match(rownames(p), colnames(q))],
                  group = q$group2[match(rownames(p), colnames(q))],
                  predicted = as.character(p$predicted.id),
                  score_max = p$prediction.score.max)
  cat("transferred cells: ", nrow(d), "\n", sep = "")
  cat("prediction.score.max by group:\n")
  print(d[, .(median = round(median(score_max), 3), q25 = round(quantile(score_max, .25), 3),
              q75 = round(quantile(score_max, .75), 3)), by = group])
  comp <- list()
  for (thr in c(0.3, 0.5, 0.7)) {
    z <- d[score_max >= thr]
    zz <- z[, .(n_cells_kept = .N,
                n_LC = sum(predicted == "2"),
                prop_LC_pct = round(100 * mean(predicted == "2"), 2),
                n_Mac = sum(predicted == "0"),
                prop_Mac_pct = round(100 * mean(predicted == "0"), 2)),
            by = .(sample_id, group)]
    zz[, score_threshold := thr]
    zz[, fraction_kept_of_all := round(n_cells_kept / nrow(d), 4)]
    comp[[length(comp) + 1]] <- zz
  }
  ct <- rbindlist(comp, fill = TRUE)
  setorder(ct, score_threshold, group, sample_id)
  fwrite(ct, file.path(TAB_DIR, "S4_query_composition.csv"))
  cat("\n--- query composition by sample (main threshold 0.5) ---\n")
  print(ct[score_threshold == 0.5][, .(sample_id, group, n_cells_kept,
                                       n_LC, prop_LC_pct)])
  gs <- ct[score_threshold == 0.5,
           .(n_samples = .N, median_LC_pct = round(median(prop_LC_pct), 2)), by = group]
  cat("\n--- group medians at threshold 0.5 (direction check, descriptive) ---\n")
  print(gs)
  fwrite(gs, file.path(TAB_DIR, "S4_query_group_summary.csv"))
  ## Raw score distribution summary for Methods (filtering bias disclosure).
  fs <- d[, .(n = .N, kept_03 = sum(score_max >= .3), kept_05 = sum(score_max >= .5),
              kept_07 = sum(score_max >= .7)), by = group]
  fwrite(fs, file.path(TAB_DIR, "S4_query_score_filter.csv"))
  print(fs)
} else {
  cat("WARNING: query transfer unavailable; Gate D will be judged on AUC + bulk only.\n")
}

cat("\n=== S4 transfer DONE ===\n")
cat("OUTPUT: results/tables/S4_validation_auc.csv\n")
cat("OUTPUT: results/tables/S4_query_composition.csv\n")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")
