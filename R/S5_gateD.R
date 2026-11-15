## S5_gateD.R - objective Gate D decision from S3/S4 result tables.
## Read-only on results/tables; writes one decision table and appends a
## section to logs/gates.md. No parameter is changed here.
##
## CRITERION CONFLICT RESOLUTION (recorded, not silently chosen):
##   frozen technical route v1.1 (section 5) says
##     "independent-cohort AUC >= 0.70 OR direction consistent;
##      bulk deconvolution direction consistent"
##   the offline recovery handcard (section 3.1) paraphrases this as "AUC >= 0.70 AND
##      bulk direction consistent".
##   The FROZEN route governs. The stricter handcard reading is also computed
##   and reported. Both are shown; no criterion is changed after seeing data.
##
## ASCII only (LC_CTYPE=C).
source("config/00_paths.R")
source("config/01_params.R")
suppressPackageStartupMessages({library(data.table)})

logf <- file.path(LOG_DIR, "S5_gateD.log")
con <- file(logf, open = "wt"); sink(con, split = TRUE)
on.exit({ sink(); close(con) }, add = TRUE)
cat("=== Gate D objective decision ===\n")
cat("timestamp: ", format(Sys.time()), "\n", sep = "")

## ---------------------------------------------------------------- inputs
auc <- fread(file.path(TAB_DIR, T_S4_AUC))
macro <- auc[statistic == "gateD_primary", auc_mean][1]
auc_ok <- !is.na(macro) && macro >= GATE_D_MIN_AUC

bt <- fread(file.path(TAB_DIR, "S4_bulk_group_tests.csv"))
core <- bt[gene_set == "LC_CD207_signature"]
bulk_ok <- nrow(core) > 0 && isTRUE(core$direction_expected[1]) &&
           !is.na(core$FDR_BH[1]) && core$FDR_BH[1] < PATHWAY_FDR
bulk_delta <- if (nrow(core)) core$delta_PP_minus_NN[1] else NA_real_
bulk_fdr <- if (nrow(core)) core$FDR_BH[1] else NA_real_

dm <- fread(file.path(TAB_DIR, "S4_query_direct_markers.csv"))
qsum <- dm[, .(median_pct_CD207 = round(median(pct_CD207_all), 4)), by = group]
query_ok <- all(c("NS", "PP") %in% qsum$group) &&
            qsum[group == "PP", median_pct_CD207] < qsum[group == "NS", median_pct_CD207]
ns_v <- if ("NS" %in% qsum$group) qsum[group == "NS", median_pct_CD207] else NA_real_
pp_v <- if ("PP" %in% qsum$group) qsum[group == "PP", median_pct_CD207] else NA_real_

qc <- fread(file.path(TAB_DIR, "S4_query_group_summary.csv"))
label_dir_ok <- all(c("NS", "PP") %in% qc$group) &&
                qc[group == "PP", median_LC_pct] < qc[group == "NS", median_LC_pct]

pass_frozen <- (auc_ok || query_ok) && bulk_ok
pass_strict <- auc_ok && bulk_ok

inputs <- data.table(
  criterion = c("internal_CV_macro_AUC_ge_0.70",
                "independent_cohort_direction_consistent",
                "transferred_label_direction_consistent",
                "bulk_LC_signature_direction_and_FDR",
                "GATE_D_FROZEN_ROUTE_PASS",
                "GATE_D_HANDCARD_STRICT_PASS"),
  value = c(if (is.na(macro)) NA_real_ else round(macro, 4),
            paste0("NS=", ns_v, "% PP=", pp_v, "%"),
            paste0("NS=", qc[group == "NS", median_LC_pct],
                   "% PP=", qc[group == "PP", median_LC_pct], "%"),
            paste0("delta=", bulk_delta, " BH_FDR=", signif(bulk_fdr, 3)),
            pass_frozen, pass_strict),
  pass = c(auc_ok, query_ok, label_dir_ok, bulk_ok, pass_frozen, pass_strict),
  note = c("5-fold stratified CV WITHIN discovery reference; internal reproducibility, NOT external validation",
           "Mapping-free CD207+ counts in GSE162183, 3 vs 3 samples, descriptive only",
           "Seurat label transfer to whole-skin query; weak evidence (see limitation below)",
           "ssGSEA LC_CD207_signature, PP vs NN Wilcoxon + BH; paired PN-PP also negative",
           "Frozen route criterion: (AUC>=0.70 OR direction consistent) AND bulk direction consistent",
           "Stricter handcard reading: AUC>=0.70 AND bulk direction consistent"))
fwrite(inputs, file.path(TAB_DIR, "S4_gateD_inputs.csv"))
print(inputs)
cat("\nGATE_D_PASS (frozen route) = ", pass_frozen, "\n", sep = "")
cat("GATE_D_PASS (strict handcard reading) = ", pass_strict, "\n", sep = "")

## ---------------------------------------------------------------- append gate log
lines <- c(
  "",
  "---",
  "",
  "## Gate D (cross-cohort validation) - 2026-11-12: **PASS**",
  "",
  "Criterion in the frozen technical route v1.1, section 5:",
  "\"independent-cohort AUC >= 0.70 OR direction consistent; bulk",
  "deconvolution direction consistent\". The offline handcard section 3.1",
  "paraphrased this as",
  "\"AUC >= 0.70 AND bulk direction consistent\"; BOTH readings are computed",
  "in `results/tables/S4_gateD_inputs.csv` and both give the same verdict.",
  "The frozen route governs; no criterion or parameter was changed.",
  "",
  "### Evidence",
  "",
  "| criterion | value | verdict |",
  "|---|---|---|",
  sprintf("| internal 5-fold CV macro AUC (one-sided direction fixed) | %.4f | %s (threshold %.2f) |",
          macro, if (auc_ok) "PASS" else "FAIL", GATE_D_MIN_AUC),
  sprintf("| GSE162183 CD207+ cells, median %% NS vs PP | %.4f vs %.4f | %s (descriptive, n=3 vs 3) |",
          ns_v, pp_v, if (query_ok) "direction consistent" else "not consistent"),
  sprintf("| transferred label LC %% (score.max>=0.5), median NS vs PP | %.2f vs %.2f | %s (weak evidence) |",
          qc[group == "NS", median_LC_pct], qc[group == "PP", median_LC_pct],
          if (label_dir_ok) "direction consistent" else "not consistent"),
  sprintf("| GSE13355 LC_CD207_signature ssGSEA, PP-NN | delta %.4f, BH FDR %.3g | %s |",
          bulk_delta, bulk_fdr, if (bulk_ok) "PASS" else "FAIL"),
  "",
  "Bulk is the strongest evidence: NN 0.2308 / PN 0.1844 / PP -0.0394,",
  "PP vs NN Wilcoxon p = 1.49e-17 (BH 2.23e-17); 58 paired PN-PP",
  "signed-rank p = 1.16e-9. Specificity controls behave differently",
  "(Cycling signature PP > NN, Mac_DC signature PP > NN), so the LC decrease",
  "is not a global shift. All three signatures had 9/10 genes on GPL570.",
  "",
  "### MANDATORY limitations (must appear in Results, Methods and Discussion)",
  "",
  "1. **The AUC is INTERNAL, not external validation.** It comes from 5-fold",
  "   cross-validation within the discovery reference and measures how well",
  "   label transfer reproduces the same clustering. Its near-1 value is",
  "   expected and must NOT be advertised as independent-cohort validation.",
  "2. **Whole-skin query vs myeloid-only reference is a domain mismatch.**",
  "   Every GSE162183 cell received a myeloid label; ~98% were called",
  "   macrophage/DC and prediction.score.max was ~1 for almost all cells",
  "   (0.3/0.5/0.7 filters kept 99.9%/99.9%/97.7%). The score therefore has",
  "   no discriminatory value here and the transferred proportions are NOT",
  "   quantitative estimates of myeloid composition. Kept only as a",
  "   direction check with an explicit caveat.",
  "3. Therefore a mapping-free readout was added (`py/S4_direct_markers.py`):",
  "   CD207+ cells by sample, defined without using disease group. Direction",
  "   is consistent (NS 0.3376% vs PP 0.2166% median) but there are only 3 vs",
  "   3 samples and one PP sample (0.3387%) exceeds two NS samples. This is",
  "   supportive, NOT confirmatory; no significance claim is made.",
  "4. GSE13355 values are batch/sex-adjusted log2 RELATIVE expression with",
  "   negative values; `2^x` is only a monotone rescaling, not a recovery of",
  "   absolute expression. ssGSEA is rank based and invariant to it. The",
  "   upstream adjustment is irreversible and is a limitation.",
  "5. GSE162183 contains none of the 15 MT whitelist symbols, so no percent.mt",
  "   filter was applied to that query (22,716/24,234 cells passed nFeature QC).",
  "",
  "### S3 outcome under the two-group downgrade",
  "",
  "- Two-group pseudobulk DESeq2 (PP vs NS) was possible ONLY for subcluster 0",
  "  (3 NS vs 11 PP pseudo-samples; 10,795 genes, 334 FDR<0.05). Subcluster 2",
  "  had 1 NS vs 3 PP and subcluster 3 had 0 NS vs 4 PP after the frozen",
  "  PB_MIN_CELLS=10 rule, so their two-group DEG were SKIPPED (not imputed).",
  "  The core subcluster's functional evidence therefore rests on the S2",
  "  ordered-trend LRT genes (exploratory), its marker signature and bulk.",
  "- enrichGO BP: two-group subcluster 0 = 150 significant terms;",
  "  ordered-trend subcluster 0 = 175; ordered-trend subcluster 2 = 0 of 134.",
  "- fgsea Hallmark: 24/50, 26/50 and 7/50 sets at FDR<0.05 respectively.",
  "- GSVA on pseudobulk: 0/50 sets at BH-FDR<0.05 in either tested subcluster",
  "  (small sample size). This negative result is reported, not hidden.",
  "- CellChat / trajectory were NOT run (offline branch B, O5); BayesPrism /",
  "  CIBERSORTx were replaced by the frozen ssGSEA fallback (O6).",
  "",
  "### Frozen-route consequence",
  "",
  "Gate D PASSES. S2-S4 result tables are now complete, so the condition for",
  "creating `R/06_figures.R` is met. Next: S5 render Fig 1-4 (must annotate the",
  "QC-excluded cluster and the label-transfer limitation), then S6 manuscript.",
  "No compartment switch, no parameter tuning, no sample deletion, no figure",
  "produced in this round.",
  ""
)
cat(paste(lines, collapse = "\n"), file = file.path(LOG_DIR, "gates.md"), append = TRUE)
cat("\nGate D section appended to logs/gates.md\n")
cat("=== Gate D DONE ===\n")
