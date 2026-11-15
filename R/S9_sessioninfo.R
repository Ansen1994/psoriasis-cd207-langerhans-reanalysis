## S9_sessioninfo.R - record the exact environment used by every analysis stage.
## The original run only captured sessionInfo() for the figure-rendering session
## (logs/S5_sessionInfo.txt), so the versions claimed in Methods for the analysis
## stages had no provenance. This script loads the packages stage by stage and
## writes one file per stage plus a combined file.
## ASCII only. Read-only: it computes nothing and writes only to logs/.

source("config/00_paths.R")

stages <- list(
  S1 = c("Seurat", "harmony", "SingleR", "scDblFinder", "SingleCellExperiment"),
  S2 = c("Seurat", "lme4", "DESeq2", "data.table"),
  S3 = c("clusterProfiler", "fgsea", "GSVA", "msigdbr", "org.Hs.eg.db"),
  S4 = c("Seurat", "pROC", "GSVA", "data.table"),
  S9 = c("data.table", "lme4", "MASS")
)

all_lines <- c("dermatology-sci -- environment snapshot per analysis stage",
               paste0("generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               "")

for (st in names(stages)) {
  cat("\n=== stage", st, "===\n")
  for (p in stages[[st]]) {
    ok <- suppressWarnings(requireNamespace(p, quietly = TRUE))
    cat(sprintf("  %-18s %s\n", p, if (ok) "available" else "NOT INSTALLED"))
    if (ok) suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  f <- file.path(LOG_DIR, paste0("S9_sessionInfo_", st, ".txt"))
  writeLines(capture.output(sessionInfo()), f)
  all_lines <- c(all_lines, paste0("=== stage ", st, " ==="),
                 capture.output(sessionInfo()), "")
  cat("  -> ", f, "\n", sep = "")
}

writeLines(all_lines, file.path(LOG_DIR, "S9_sessionInfo_all.txt"))
cat("\nWROTE ", file.path(LOG_DIR, "S9_sessionInfo_all.txt"), "\n", sep = "")
