## 03_s1_integrate.R - S1 main pipeline: QC filter -> doublets -> normalization
## -> Harmony (sample/batch ONLY) -> clustering -> annotation -> Gate B inputs.
## ASCII ONLY. Working directory = project root.
## This script computes embeddings but produces NO figure files.
source("config/00_paths.R")
source("config/01_params.R")

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(harmony)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(SingleR)
})

log_file <- file.path(LOG_DIR, "S1_integrate.log")
con <- file(log_file, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

cat("=== S1 integration ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
cat("Seurat: ", as.character(packageVersion("Seurat")), "\n", sep = "")
set.seed(SEED)
t_start <- Sys.time()

dat <- readRDS(file.path(PROC_DIR, "S1_counts.rds"))
counts <- dat$counts
meta <- as.data.frame(dat$meta)
qc <- data.table::fread(file.path(TAB_DIR, "S1_qc_per_cell.csv.gz"))
stopifnot(identical(as.character(qc$cell_id), as.character(meta$cell_id)))
cat("loaded matrix: ", nrow(counts), " genes x ", ncol(counts), " cells\n", sep = "")

## ------------------------------------------------------------ 1. QC filter
keep <- qc$nFeature >= QC_NFEATURE_MIN &
        qc$nFeature <= QC_NFEATURE_MAX &
       !is.na(qc$pct_mt) & qc$pct_mt <= QC_PCT_MT_MAX
cat("QC thresholds used: nFeature in [", QC_NFEATURE_MIN, ",", QC_NFEATURE_MAX,
    "], pct_mt <= ", QC_PCT_MT_MAX, "\n", sep = "")
cat("QC pass: ", sum(keep), " / ", length(keep),
    " (", round(100 * mean(keep), 1), "%)\n", sep = "")

## NS-AR008: low cell count is NOT a quality failure by itself. The decision
## was made from the observed QC distribution and recorded in logs/gates.md.
if (isTRUE(DROP_SMALL_SAMPLE)) {
  sid <- sub("\\.csv\\.gz$", "", FLAG_SMALL_SAMPLE)
  drop <- qc$sample_id == sid
  cat("NS-AR008 explicitly EXCLUDED by DROP_SMALL_SAMPLE=TRUE (", sum(keep & drop),
      " cells)\n", sep = "")
  keep <- keep & !drop
} else {
  cat("NS-AR008 RETAINED: cell count is low but QC metrics are in range; ",
      "decision logged in logs/gates.md\n", sep = "")
}

counts <- counts[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
qc <- qc[keep]
cat("after QC: ", ncol(counts), " cells\n", sep = "")
rm(dat); gc(FALSE)

## ------------------------------------------------------------ 2. doublets
## Doublets can only form WITHIN one 10x lane, so call per sample.
## Cache the calls: scDblFinder is stochastic and takes ~4 min on this cohort.
cat("\n=== scDblFinder (per sample) ===\n")
t0 <- Sys.time()
dbl_file <- file.path(PROC_DIR, "S1_doublet_calls.rds")
if (file.exists(dbl_file)) {
  dbl <- readRDS(dbl_file)
  cat("reusing cached doublet calls: ", length(dbl), "\n", sep = "")
} else {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(sample_id = meta$sample_id)
  )
  sce <- scDblFinder::scDblFinder(sce, samples = sce$sample_id, verbose = FALSE)
  dbl <- as.character(sce$scDblFinder.class)
  saveRDS(dbl, dbl_file)
  rm(sce); gc(FALSE)
}
dbl_tab <- table(dbl)
cat("doublet table: ", paste(names(dbl_tab), as.numeric(dbl_tab), sep = "=", collapse = ", "), "\n", sep = "")
dbl_rate <- mean(dbl == "doublet")
cat("overall doublet rate: ", round(100 * dbl_rate, 2), "%\n", sep = "")
if (dbl_rate > DOUBLET_ABORT) {
  cat("WARNING: doublet rate exceeds ", DOUBLET_ABORT, " -> review required\n", sep = "")
} else {
  cat("doublet rate within expected range [", paste(DOUBLET_EXPECT * 100, collapse = "-"), "%]\n", sep = "")
}
dbl_summary <- data.table::data.table(
  sample_id = meta$sample_id, is_doublet = dbl == "doublet"
)[, .(n_cells = .N, n_doublet = sum(is_doublet),
      pct_doublet = round(100 * mean(is_doublet), 2)), by = sample_id]
data.table::fwrite(dbl_summary, file.path(TAB_DIR, "S1_doublet_summary.csv"))

singlet <- dbl == "singlet"
counts <- counts[, singlet, drop = FALSE]
meta <- meta[singlet, , drop = FALSE]
qc <- qc[singlet]
cat("after doublet removal: ", ncol(counts), " cells\n", sep = "")
cat("doublet_sec: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n", sep = "")
gc(FALSE)

## ------------------------------------------------------------ 3. Seurat object
meta$state_ord <- unname(STATE_ORDER[meta$group])
meta$batch <- "GSE173706"   ## single platform; Harmony corrects sample_id ONLY
meta$percent_mt <- qc$pct_mt
rownames(meta) <- colnames(counts)
obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta,
                                  min.cells = 0, min.features = 0,
                                  project = DISCOVERY_ID)
rm(counts); gc(FALSE)
cat("\nSeurat object: ", ncol(obj), " cells x ", nrow(obj), " genes\n", sep = "")

obj <- Seurat::NormalizeData(obj, normalization.method = NORM_METHOD,
                             scale.factor = 10000, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                    nfeatures = N_HVG, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
obj <- Seurat::RunPCA(obj, npcs = N_PCA, verbose = FALSE)

## O3: NEVER include disease group in Harmony covariates.
cat("Harmony covariates: ", paste(HARMONY_VARS, collapse = ","),
    " (batch is constant here, so sample_id is active)\n", sep = "")
obj <- harmony::RunHarmony(obj, group.by.vars = "sample_id",
                           theta = HARMONY_THETA, max.iter = HARMONY_MAXITER,
                           verbose = FALSE)
obj <- Seurat::FindNeighbors(obj, reduction = "harmony", dims = 1:N_DIMS,
                             verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = CLUST_RES_MAIN, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, reduction = "harmony", dims = 1:N_DIMS,
                       verbose = FALSE)
cat("main clusters at res=", CLUST_RES_MAIN, ": ",
    length(unique(Seurat::Idents(obj))), "\n", sep = "")

## ------------------------------------------------------------ 4. annotation
get_mat <- function(o) {
  tryCatch(Seurat::GetAssayData(o, assay = "RNA", layer = "data"),
           error = function(e) Seurat::GetAssayData(o, assay = "RNA", slot = "data"))
}
expr <- get_mat(obj)
avail <- rownames(expr)
for (nm in names(ANNOT_MARKERS)) {
  g <- intersect(ANNOT_MARKERS[[nm]], avail)
  cat("marker set ", nm, ": ", length(g), "/", length(ANNOT_MARKERS[[nm]]),
      " present\n", sep = "")
}
obj <- Seurat::AddModuleScore(obj, features = ANNOT_MARKERS,
                              name = "Module", ctrl = 100, seed = SEED)
score_cols <- grep("^Module", colnames(obj@meta.data), value = TRUE)
score_cols <- score_cols[order(as.numeric(sub("^Module", "", score_cols)))]
stopifnot(length(score_cols) == length(ANNOT_MARKERS))
cl <- as.character(Seurat::Idents(obj))
sc <- as.data.frame(obj@meta.data[, score_cols, drop = FALSE])
colnames(sc) <- names(ANNOT_MARKERS)
sc$cluster <- cl
cluster_score <- aggregate(. ~ cluster, data = sc, FUN = mean)
best <- names(ANNOT_MARKERS)[apply(cluster_score[, names(ANNOT_MARKERS), drop = FALSE],
                                    1, which.max)]
names(best) <- as.character(cluster_score$cluster)
## A named vector makes Seurat treat the names as CELL IDs, which errors with
## "No cell overlap". AddMetaData with an unnamed vector assigns by position.
obj <- Seurat::AddMetaData(obj, metadata = unname(best[cl]),
                           col.name = "lineage_marker")
cat("marker-based lineage assignment:\n")
print(table(obj$lineage_marker))

## SingleR is an optional cross-check; celldex references need a download.
## Aggregate by cluster (SingleR's clusters= argument) instead of scoring
## 84k cells individually: far faster and more stable on this cohort size.
singler_ok <- FALSE
try({
  options(timeout = 180)
  hpca <- celldex::HumanPrimaryCellAtlasData()
  pred <- SingleR::SingleR(test = expr, ref = hpca,
                           labels = hpca$label.main, clusters = cl)
  cl2label <- stats::setNames(as.character(pred$labels), rownames(pred))
  obj <- Seurat::AddMetaData(obj, metadata = unname(cl2label[cl]),
                             col.name = "singler_label")
  singler_ok <- TRUE
  cat("SingleR HPCA available; cluster-level labels projected to cells.\n")
}, silent = TRUE)
if (!singler_ok) {
  obj <- Seurat::AddMetaData(obj,
                             metadata = rep(NA_character_, ncol(obj)),
                             col.name = "singler_label")
  cat("SingleR SKIPPED: celldex reference not cached / offline. ",
      "marker-based annotation is authoritative; recorded as a limitation.\n", sep = "")
}

## ------------------------------------------------------------ 5. outputs
cl_tab <- as.data.frame(table(cluster = cl,
                              lineage_marker = obj$lineage_marker,
                              group = obj$group))
celltype_summary <- data.table::as.data.table(cl_tab)
data.table::setnames(celltype_summary, "Freq", "n_cells")
data.table::fwrite(celltype_summary, file.path(TAB_DIR, T_S1_SUMMARY))

## Wilcoxon cluster markers are deliberately deferred to S2 (target-compartment
## subclusters), per the frozen 4-figure minimal loop. S1 only needs the
## lineage/cluster composition below for Gate B scoring.

## Gate B inputs: per-lineage cell counts, group composition, marker richness.
lin_tab <- data.table::data.table(
  lineage = obj$lineage_marker, group = obj$group, cluster = cl
)
by_lin <- lin_tab[, .(n_cells = .N, n_clusters = uniqueN(cluster)), by = lineage]
by_state <- lin_tab[, .N, by = .(lineage, group)]
data.table::fwrite(by_lin, file.path(PROC_DIR, "S1_gateB_lineages.csv"))
data.table::fwrite(by_state, file.path(PROC_DIR, "S1_gateB_states.csv"))
cat("\nGate B input by lineage:\n"); print(by_lin)
cat("\nGate B state composition:\n"); print(by_state)

saveRDS(obj, file.path(PROC_DIR, F_S1_OBJECT))
cat("\nOUTPUT: ", file.path(PROC_DIR, F_S1_OBJECT), "\n", sep = "")
cat("OUTPUT: ", file.path(TAB_DIR, T_S1_SUMMARY), "\n", sep = "")
cat("OUTPUT: results/tables/S1_doublet_summary.csv\n")
cat("OUTPUT: data/processed/S1_gateB_lineages.csv + S1_gateB_states.csv\n")
cat("total_sec: ", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), "\n", sep = "")
cat("=== S1 integration DONE; gate B decision in logs/gates.md ===\n")
