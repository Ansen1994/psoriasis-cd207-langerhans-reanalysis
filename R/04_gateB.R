## 04_gateB.R - Gate B objective scoring inputs.
## Sample-level (not cell-level) lineage proportions + ordered trend pilot.
## This is a GATE INPUT, not a result: the main analysis is S2 at subcluster
## level with pseudobulk and donor random effects.
## ASCII ONLY. Working directory = project root. No figures.
source("config/00_paths.R")
source("config/01_params.R")
source("R/03_trend.R")

suppressPackageStartupMessages({ library(data.table); library(Seurat) })

log_file <- file.path(LOG_DIR, "S1_gateB.log")
con <- file(log_file, open = "wt")
sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)

cat("=== Gate B scoring inputs ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")

obj <- readRDS(file.path(PROC_DIR, F_S1_OBJECT))
md <- as.data.table(obj@meta.data)
rm(obj); gc(FALSE)
cat("cells: ", nrow(md), " | samples: ", uniqueN(md$sample_id), "\n", sep = "")

## ---------------------------------------------- annotation cross-check
cat("\n--- SingleR HPCA vs canonical-marker assignment (cells) ---\n")
ct <- as.data.table(table(marker = md$lineage_marker,
                          singleR = md$singler_label))
ct <- ct[N > 0][order(-N)]
print(head(ct, 20))

## ---------------------------------------------- sample-level proportions
props <- md[, .(n_lineage = .N), by = .(sample_id, group, lineage_marker)]
tot <- md[, .(n_total = .N), by = sample_id]
props <- merge(props, tot, by = "sample_id")
props[, prop := 100 * n_lineage / n_total]
props[, state_ord := unname(STATE_ORDER[group])]
setorder(props, lineage_marker, state_ord, sample_id)
data.table::fwrite(props, file.path(PROC_DIR, "S1_gateB_sample_props.csv"))

cat("\n--- sample-level median lineage proportion (%) by state ---\n")
pivot <- props[, .(median_prop = round(stats::median(prop), 3),
                   n_samples = .N,
                   total_cells = sum(n_lineage)),
               by = .(lineage_marker, state_ord)]
print(pivot[order(lineage_marker, state_ord)])

## ---------------------------------------------- ordered trend pilot
cat("\n--- Jonckheere-Terpstra pilot on sample-level proportions ---\n")
dat <- props[, .(feature = lineage_marker, value = prop, state_ord = state_ord)]
tr <- jt_trend_table(dat, B = TREND_NPERM, seed = SEED,
                     alternative = TREND_ALTERNATIVE)
data.table::fwrite(as.data.table(tr), file.path(TAB_DIR, "S1_gateB_trend_pilot.csv"))
print(tr)

## ---------------------------------------------- subclustering feasibility
cat("\n--- S2 feasibility: donor x subcluster cell budget ---\n")
feas <- md[, .(n_cells = .N), by = .(sample_id, lineage_marker)]
feas <- feas[, .(median_cells_per_sample = round(stats::median(n_cells), 1),
                 min_cells_per_sample = min(n_cells),
                 n_samples = .N), by = lineage_marker]
## With k subclusters, expected cells per donor x subcluster.
for (k in c(SUB_CLUST_MIN, 6, SUB_CLUST_MAX)) {
  cn <- paste0("est_per_combo_k", k)
  feas[, (cn) := round(median_cells_per_sample / k, 1)]
}
print(feas[order(lineage_marker)])

cat("\nOUTPUT: data/processed/S1_gateB_sample_props.csv\n")
cat("OUTPUT: results/tables/S1_gateB_trend_pilot.csv\n")
cat("=== Gate B inputs DONE ===\n")
