## 05_s2_subcluster.R - S2 steps 1-3.
## Subset the Gate B compartment (myeloid/DC), recluster inside it with a
## pre-registered fixed k = 4, and compute display-only marker genes.
## ASCII ONLY. Working directory = project root. NO figures are written.
source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(harmony); library(Matrix)
})

con <- file(file.path(LOG_DIR, "S2_subcluster.log"), open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

set.seed(SEED)
t0 <- Sys.time()
cat("=== S2 steps 1-3: myeloid/DC subcluster (frozen k=4) ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
cat("Seurat: ", as.character(packageVersion("Seurat")), "\n", sep = "")

## ------------------------------------------------------------ step 1: subset
obj <- readRDS(file.path(PROC_DIR, F_S1_OBJECT))
cat("S1 object: ", ncol(obj), " cells x ", nrow(obj), " genes\n", sep = "")
cells <- colnames(obj)[obj$lineage_marker == "myeloid"]
cat("myeloid/DC cells selected: ", length(cells), "\n", sep = "")
stopifnot(length(cells) > 0)
my <- subset(obj, cells = cells)
rm(obj); gc(FALSE)
saveRDS(my, file.path(PROC_DIR, "S2_myeloid.rds"))
cat("saved: data/processed/S2_myeloid.rds\n")

## ------------------------------------------------------------ step 2: recluster
## HVG/PCA/Harmony are recomputed INSIDE the compartment. Reusing the global
## HVGs would bias the split toward genes that separate major lineages.
my <- Seurat::FindVariableFeatures(my, selection.method = "vst",
                                   nfeatures = N_HVG, verbose = FALSE)
cat("HVGs within compartment: ", length(Seurat::VariableFeatures(my)), "\n", sep = "")
my <- Seurat::ScaleData(my, verbose = FALSE)
my <- Seurat::RunPCA(my, npcs = N_PCA, verbose = FALSE)

## O3: Harmony corrects sample_id ONLY. batch is a constant (one platform),
## and disease group is NEVER a covariate.
cat("Harmony covariates: sample_id only (batch constant; disease excluded)\n")
my <- harmony::RunHarmony(my, group.by.vars = "sample_id",
                          theta = HARMONY_THETA, max.iter = HARMONY_MAXITER,
                          verbose = FALSE)
my <- Seurat::FindNeighbors(my, reduction = "harmony", dims = 1:N_DIMS,
                            verbose = FALSE)

## Resolution scan over the frozen range; k = 4 is a Gate B power mitigation,
## NOT a tuning target for significance.
nclust_at <- function(res) {
  tmp <- Seurat::FindClusters(my, resolution = res, verbose = FALSE)
  length(unique(as.character(Seurat::Idents(tmp))))
}
scan <- data.table(resolution = SUB_RES_RANGE)
scan$n_clusters <- vapply(seq_along(scan$resolution),
                          function(i) nclust_at(scan$resolution[i]), integer(1))
cat("\n--- resolution scan (frozen SUB_RES_RANGE) ---\n")
print(scan)

exact <- scan[n_clusters == SUB_CLUST_MIN]
if (nrow(exact) > 0) {
  chosen <- min(exact$resolution)
  how <- "exact k=4 inside frozen SUB_RES_RANGE"
} else {
  ## Bisection targets the pre-registered k=4, never a p-value.
  lo <- 0.01; hi <- 2.0; best_r <- NA_real_; best_k <- NA_integer_
  for (i in seq_len(14)) {
    mid <- (lo + hi) / 2
    k <- nclust_at(mid)
    if (is.na(best_k) || abs(k - SUB_CLUST_MIN) < abs(best_k - SUB_CLUST_MIN)) {
      best_k <- k; best_r <- mid
    }
    if (k == SUB_CLUST_MIN) { best_k <- k; best_r <- mid; break }
    if (k > SUB_CLUST_MIN) hi <- mid else lo <- mid
  }
  chosen <- best_r
  how <- paste0("bisection to reach frozen k=4 (closest k=", best_k, ")")
}
cat("chosen resolution: ", round(chosen, 4), " | ", how, "\n", sep = "")

my <- Seurat::FindClusters(my, resolution = chosen, verbose = FALSE)
my$subcluster <- as.character(Seurat::Idents(my))
k_final <- length(unique(my$subcluster))
cat("final subcluster count: ", k_final, " (target ", SUB_CLUST_MIN, ")\n", sep = "")
if (k_final != SUB_CLUST_MIN) {
  cat("WARNING: k != ", SUB_CLUST_MIN, "; recorded, NOT silently forced.\n", sep = "")
}
cat("subcluster sizes:\n"); print(table(my$subcluster))
cat("subcluster x group:\n"); print(table(my$subcluster, my$group))

my <- Seurat::RunUMAP(my, reduction = "harmony", dims = 1:N_DIMS,
                      verbose = FALSE)
saveRDS(my, file.path(PROC_DIR, F_S2_OBJECT))
data.table::fwrite(scan, file.path(TAB_DIR, "S2_resolution_scan.csv"))

## ------------------------------------------------------------ step 3: markers
## Display only. FindMarkers p-values are NEVER used for inference (O2).
Seurat::Idents(my) <- "subcluster"
mk <- tryCatch(
  Seurat::FindAllMarkers(my, only.pos = TRUE, min.pct = MARKER_PCT1,
                         logfc.threshold = MARKER_LOGFC,
                         test.use = MARKER_TEST, verbose = FALSE),
  error = function(e) { cat("FindAllMarkers failed: ", conditionMessage(e), "\n"); NULL }
)
if (!is.null(mk) && nrow(mk) > 0) {
  setDT(mk)
  data.table::fwrite(mk, file.path(TAB_DIR, T_S2_MARKERS))
  top <- mk[, head(.SD[order(-avg_log2FC)], MARKER_TOPN), by = cluster]
  data.table::fwrite(top, file.path(TAB_DIR, "S2_markers_top10.csv"))
  cat("\nmarker rows: ", nrow(mk), " | top10 rows: ", nrow(top), "\n", sep = "")
  for (s in sort(unique(mk$cluster))) {
    g <- top[cluster == s, head(gene, 8)]
    cat("  subcluster ", s, ": ", paste(g, collapse = ", "), "\n", sep = "")
  }
} else {
  cat("WARNING: no markers produced.\n")
}

cat("\nOUTPUT: data/processed/S2_subcluster.rds\n")
cat("OUTPUT: results/tables/S2_markers.csv, S2_markers_top10.csv, S2_resolution_scan.csv\n")
cat("NO FIGURES WRITTEN\n")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")
cat("=== S2 steps 1-3 DONE ===\n")
