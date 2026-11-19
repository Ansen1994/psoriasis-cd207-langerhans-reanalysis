## S20_treatment_sc.R - post hoc lean analysis of GSE183047 (immune-cell-
## enriched scRNA-seq of psoriasis lesional skin before and after 12 weeks of
## secukinumab; Kim & Krueger, Front Immunol 2023; 10 repository controls).
##
## QUESTION: does the CD207+ fraction of the myeloid/DC compartment recover
## after effective IL-17A blockade at the single-cell level?
##
## METHOD (lean, mirrors the mapping-free readout of Fig 4B):
##   per sample: 10x MTX -> QC (200 <= nFeature <= 6000, percent.mt <= 20
##   where MT- symbols exist) -> LogNormalize -> AddModuleScore with the
##   SEVEN FROZEN ANNOT_MARKERS sets (config/01_params.R) -> argmax lineage
##   -> myeloid/DC compartment = cells assigned to "myeloid" or "DC".
##   Readout per sample: CD207+ fraction among QC-eligible myeloid/DC cells
##   (detection, count > 0), plus the myeloid/DC share of all cells.
##
## GROUPS: PreTx (n = 10 lesional), postTx week12 (n = 4), week24 (n = 1),
##   week48 (n = 1), one preTx non-lesional sample (PN, descriptive only),
##   healthy controls (n = 10). Paired pre -> week12 exists for patients
##   04/05/07/08 (n = 4 pairs < 6), so NO significance test is run here -
##   descriptive per-patient reporting only, as in the frozen convention.
##
## CAVEATS: immune-cell-ENRICHED protocol - proportions are NOT comparable
##   with the whole-skin discovery cohort; module scores are computed per
##   sample without integration (batch across protocol groups is disclosed);
##   the POST arm mixes weeks 12/24/48 and is reported per time point.
##
## ASCII only (LC_CTYPE=C). All outputs stay inside the project.
source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(Matrix)
})

logf <- file.path(LOG_DIR, "S20_treatment_sc.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== S20 GSE183047 treatment scRNA-seq (post hoc, lean) ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")
set.seed(SEED)
t0 <- Sys.time()

raw <- file.path(RAW_DIR, "PSOR_SC_GSE183047", "RAW")
mtx_files <- list.files(raw, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
ids <- sub("_matrix\\.mtx\\.gz$", "", basename(mtx_files))
cat("samples: ", length(ids), "\n", sep = "")

LINEAGES <- names(ANNOT_MARKERS)
res <- vector("list", length(ids))
for (k in seq_along(ids)) {
  id <- ids[k]
  cat("\n[", k, "/", length(ids), "] ", id, "\n", sep = "")
  m <- ReadMtx(mtx = file.path(raw, paste0(id, "_matrix.mtx.gz")),
               cells = file.path(raw, paste0(id, "_barcodes.tsv.gz")),
               features = file.path(raw, paste0(id, "_features.tsv.gz")))
  obj <- CreateSeuratObject(counts = m, project = id)
  rm(m); gc(FALSE)
  obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  n0 <- ncol(obj)
  obj <- subset(obj, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                       percent_mt <= 20)
  cat("  cells: ", n0, " -> QC: ", ncol(obj), "\n", sep = "")
  if (ncol(obj) < 10) { res[[k]] <- NULL; next }
  obj <- NormalizeData(obj, verbose = FALSE)
  for (ln in LINEAGES) {
    genes <- intersect(ANNOT_MARKERS[[ln]], rownames(obj))
    if (length(genes) < 1) { obj[[paste0("mod_", ln)]] <- NA; next }
    obj <- AddModuleScore(obj, features = list(genes),
                          name = paste0("mod_", ln), verbose = FALSE)
  }
  md <- obj@meta.data
  mod_cols <- paste0("mod_", LINEAGES, "1")   # AddModuleScore appends "1"
  M <- as.matrix(md[, mod_cols, drop = FALSE])
  arg <- apply(M, 1, function(z) LINEAGES[which.max(z)])
  myeloid <- arg %in% c("myeloid", "DC")
  e <- GetAssayData(obj, assay = "RNA", slot = "data")
  cd207_vec <- if ("CD207" %in% rownames(e)) as.numeric(e["CD207", ]) else rep(NA_real_, ncol(e))
  cd1a_vec  <- if ("CD1A"  %in% rownames(e)) as.numeric(e["CD1A", ])  else rep(NA_real_, ncol(e))
  my_idx <- which(myeloid)
  n_my <- length(my_idx)
  cd207_my <- if (n_my) cd207_vec[my_idx] else numeric(0)
  res[[k]] <- data.table(
    sample_id = id,
    title = sub("^GSM[0-9]+_", "", id),
    n_cells_total = n0, n_cells_qc = ncol(obj),
    n_myeloid = n_my,
    frac_myeloid = round(n_my / ncol(obj), 4),
    frac_CD207_my = if (n_my) round(mean(cd207_my > 0), 4) else NA_real_,
    median_CD207_my = if (n_my) round(median(cd207_my), 4) else NA_real_,
    median_CD207_my_pos = if (any(cd207_my > 0)) round(median(cd207_my[cd207_my > 0]), 4) else NA_real_,
    frac_CD207_all = round(mean(cd207_vec > 0), 4),
    median_CD1A_my = if (n_my) round(median(cd1a_vec[my_idx]), 4) else NA_real_)
  rm(obj, e, md, M); gc(FALSE)
}

out <- rbindlist(res, fill = TRUE)
## group assignment from the sample title
out[, group := fifelse(grepl("Control", title), "Control",
        fifelse(grepl("_NL", title), "PN_preTx",
        fifelse(grepl("preTx", title), "PreTx",
        fifelse(grepl("week12", title), "PostTx_wk12",
        fifelse(grepl("week24", title), "PostTx_wk24",
        fifelse(grepl("week48", title), "PostTx_wk48", NA_character_))))))]
out[, patient := sub("^(Psoriasis[0-9]+|Control[0-9]+).*$", "\\1", title)]
out[, timepoint := fifelse(grepl("week12", title), "week12",
       fifelse(grepl("week24", title), "week24",
       fifelse(grepl("week48", title), "week48",
       fifelse(grepl("preTx", title), "baseline", "control"))))]
setorder(out, group, title)
fwrite(out, file.path(TAB_DIR, "S20_treatment_sc_samples.csv"))

cat("\n=== group summary (median [min-max]) ===\n")
print(out[, .(n = .N,
              frac_CD207_my = round(median(frac_CD207_my, na.rm = TRUE), 4),
              frac_myeloid = round(median(frac_myeloid, na.rm = TRUE), 4),
              median_CD1A_my = round(median(median_CD1A_my, na.rm = TRUE), 3)),
          by = group])
cat("\n=== paired patients (pre -> week12), descriptive only ===\n")
pairtab <- dcast(out[patient %in% c("Psoriasis04", "Psoriasis05", "Psoriasis07",
                                    "Psoriasis08") & timepoint %in% c("baseline", "week12")],
                 patient ~ timepoint, value.var = "frac_CD207_my")
print(pairtab)
cat("\nNOTE: n = 4 pairs < 6 -> no test (frozen convention); descriptive only.\n")
cat("\nelapsed: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min\n", sep = "")
cat("=== S20 done ===\n")
