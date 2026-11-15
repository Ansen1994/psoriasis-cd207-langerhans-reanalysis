## 01_params.R - frozen analysis parameters for dermatology-sci
## Frozen version: v1.1 (2026-11-12, Gate E passed)
## ASCII only. This machine starts R with LC_CTYPE=C; never put non-ASCII
## characters or non-ASCII paths in any R script.
## All scripts must source("config/00_paths.R") first, then this file.
## Values below are defaults; Gate A recalibrates QC thresholds using the
## observed distribution (median +/- 3*MAD) and records the change in logs/.

## ---------------------------------------------------------------- reproducibility
SEED              <- 20261112
PROJ_VERSION      <- "v1.1"
set.seed(SEED)

## ---------------------------------------------------------------- datasets
## GEO official pages verified 2026-11-12. See the S0 download manifest in docs/.
DISCOVERY_ID      <- "GSE173706"   ## 33 samples: NS=8, PN=11, PP=14
VALIDATION_ID     <- "GSE162183"   ## 6 samples: 3 lesional + 3 healthy
BULK_ID           <- "GSE13355"    ## 180 samples: 58 PP + 58 PN + 64 NN

## Group prefixes are parsed from GSM titles, not from a hard-coded table.
GROUP_PREFIX      <- c(NS = "healthy", PN = "nonlesional", PP = "lesional")
## Ordered state for the trend test (O1). Direction is declared one-sided.
STATE_ORDER       <- c(NS = 0L, PN = 1L, PP = 2L)

## Donor pairing is REAL, not assumed (verified from GSE173706 filelist.txt).
## At least 10 donors have both PN and PP; PP-6215A/6215B are two biopsies of
## one donor; PN/PP-30696 and -30696V3 are two library versions of one donor.
## Consequence: never treat these as independent samples. See O8.
PAIRING_AWARE     <- TRUE
DONOR_RANDOM_EFF  <- TRUE          ## lmer(... + (1|donor))
COLLAPSE_BY_DONOR <- TRUE          ## pseudobulk collapses technical replicates

## Suspect low-cell sample: GSM5277176_NS-AR008.csv.gz is only 0.5 Mb while
## all other GSE173706 CSVs are 2.7-15.9 Mb. Gate A must inspect it.
FLAG_SMALL_SAMPLE <- "GSM5277176_NS-AR008.csv.gz"
## Low cell count alone is not a QC failure. Decision was made from the observed
## QC distribution (see logs/gates.md): median nFeature 305 vs 2152 overall,
## median nCount 513 vs 7785 overall => low-depth library, EXCLUDE the sample.
## NS group still has 7 donors, well above Gate A's minimum of 3.
DROP_SMALL_SAMPLE <- TRUE

## ---------------------------------------------------------------- S1 QC
QC_NFEATURE_MIN   <- 200
QC_NFEATURE_MAX   <- 6000
QC_PCT_MT_MAX     <- 20            ## skin keratinocytes run high; 20% default
QC_PCT_MT_SENS    <- c(10, 15, 20) ## sensitivity check, log the comparison
QC_MAD_K          <- 3             ## recalibration rule at Gate A
## S1 finding (2026-11-12): org.Hs.eg.db 3.20.0 uses NCBI symbols for the
## mitochondrial genome (ND1, COX2, ATP8 ...), NOT HGNC "MT-ND1". A "^MT-"
## pattern matches ZERO symbols in this package, so percent.mt would be all
## zero. Use this explicit whitelist instead and state it in Methods.
## Added after Gate E as a bug fix, not a retrospective threshold change.
MT_SYMBOL_WHITELIST <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
                         "COX1","COX2","COX3","ATP6","ATP8","CYTB",
                         "RNR1","RNR2")
DOUBLET_TOOL      <- "scDblFinder"
DOUBLET_EXPECT    <- c(0.05, 0.08)
DOUBLET_ABORT     <- 0.12          ## >12% triggers a review

NORM_METHOD       <- "LogNormalize" ## standard Seurat log-normalisation (unchanged value)
NORM_SENSITIVITY  <- "SCTransform" ## optional check only, never a second main analysis
N_HVG             <- 3000
N_PCA             <- 50
N_DIMS            <- 30

## Harmony corrects SAMPLE / BATCH only. Correcting by disease group erases
## the real biological effect and destroys the three-state design (O3).
HARMONY_VARS      <- c("sample_id", "batch")
HARMONY_THETA     <- 2
HARMONY_MAXITER   <- 20
HARMONY_FORBIDDEN <- "disease_group"

CLUST_RES_MAIN    <- 0.5           ## aim: 8-12 major lineages
CLUST_RES_SENS    <- c(0.8, 1.0)

ANNOT_TOOL        <- "SingleR"
ANNOT_MARKERS     <- list(
  keratinocyte = c("KRT1", "KRT10", "KRT14"),
  Tcell        = c("CD3D", "CD2"),
  myeloid      = c("LYZ", "CD14"),
  fibroblast   = c("COL1A1", "DCN"),
  endothelial  = c("PECAM1"),
  DC           = c("LAMP3"),
  mast         = c("TPSAB1")
)
## Reference atlas (E-MTAB-8142 / CELLxGENE) is used for LABEL TRANSFER ONLY,
## never integrated into the discovery object.

## ---------------------------------------------------------------- S2 subcluster
SUB_RES_RANGE     <- c(0.4, 0.6, 0.8)
SUB_CLUST_MIN     <- 4
SUB_CLUST_MAX     <- 8
MARKER_TEST       <- "wilcox"
MARKER_LOGFC      <- 0.25
MARKER_PCT1       <- 0.25
MARKER_ADJ        <- "bonferroni"
MARKER_TOPN       <- 10

## O1: ordered trend test replaces pairwise three-group comparison.
TREND_TEST        <- "JonckheereTerpstra"
TREND_NPERM       <- 10000
TREND_ALTERNATIVE <- "greater"     ## declared before seeing the data
TREND_FDR         <- 0.05
TREND_FDR_METHOD  <- "BH"
TREND_REQUIRE_MONOTONIC <- TRUE

## O2: sample-level pseudobulk; cells are NOT independent observations.
PB_MIN_CELLS      <- 10            ## drop sample x subcluster below this
PB_TOOL           <- "DESeq2"
PB_DESIGN         <- ~ state_ord + batch
PB_TEST           <- "LRT"
PB_FDR            <- 0.05
## FindMarkers p-values are for marker display ONLY, never for inference.

## ---------------------------------------------------------------- S3 function
PATHWAY_TOOL      <- "clusterProfiler"
PATHWAY_DB        <- c("GO", "KEGG")
PATHWAY_FDR       <- 0.05
## v2 3.4: trajectory and CellChat are MUTUALLY EXCLUSIVE. Never both.
S3_BRANCH         <- "A"           ## "A" = DEG+pathway+CellChat; "B" = DEG+pathway+Monocle3
CELLCHAT_KEEP     <- 6             ## core subcluster + 4-6 partner types max
CELLCHAT_NPERM    <- 1000
CELLCHAT_MIN_INT  <- 10            ## below this, switch to branch B (O5 fallback)
TRAJECTORY_TOOL   <- "Monocle3"

## ---------------------------------------------------------------- S4 validation
VALIDATION_METHOD <- "SeuratLabelTransfer" ## no re-integration, no Harmony
VALIDATION_AUC_MIN <- 0.70
BULK_SOURCE_FILE  <- "GSE13355_series_matrix.txt.gz"  ## 34 MB, RMA log2 values
BULK_LOG2         <- TRUE          ## MUST apply 2^x before deconvolution
BULK_CAVEAT       <- "values already adjusted for batch and sex effects; record as a limitation"
## O6: tool downgrade chain. Never let an external service block Fig 4.
DECONV_CHAIN      <- c("BayesPrism", "CIBERSORTx", "ssGSEA_module_score")
DECONV_FALLBACK   <- "ssGSEA_module_score"
## GSE13355_RAW.tar (911.2 Mb CEL) is NOT downloaded by default. Only fetch it
## if a reviewer requires a fresh RMA run; then document the extra step.
BULK_CEL_OPTIONAL <- FALSE

## ---------------------------------------------------------------- Gate criteria
GATE_A_MIN_PER_GROUP <- 3
GATE_B_MIN_CELLS     <- 2000
GATE_B_MIN_SUBCLUST  <- 4
GATE_B_TIE_BREAK     <- "myeloid_DC"  ## pre-registered tie-break; avoids fibroblast and T-cell compartments
GATE_C_MIN_SIGNIFICANT <- 1
GATE_D_MIN_AUC       <- 0.70
## Gate E is procedural: this file + the frozen technical-route doc + Figure spec must
## all be frozen before R/06_figures.R is created. NO figures before that.

## ---------------------------------------------------------------- figure spec
FIG_DPI           <- 600
FIG_VECTOR        <- TRUE
FIG_STATE_COLORS  <- c(NS = "#4C78A8", PN = "#F58518", PP = "#C44E52")
FIG_SHOW_N        <- "number of samples (donors), never number of cells"
FIG_STAT_LABEL    <- c("test name", "n", "FDR method", "error bar meaning")

## ---------------------------------------------------------------- output names
F_S1_OBJECT       <- "S1_seurat_integrated.rds"
T_S1_SUMMARY      <- "S1_celltype_summary.csv"
F_S2_OBJECT       <- "S2_subcluster.rds"
T_S2_MARKERS      <- "S2_markers.csv"
T_S2_PROP         <- "S2_proportions_by_sample.csv"
T_S2_TREND        <- "S2_trend_tests.csv"
T_S2_DEG          <- "S2_pseudobulk_DEG.csv"
T_S3_PATHWAY      <- "S3_pathway.csv"
T_S3_COMM         <- "S3_communication.csv"
T_S3_TRAJ         <- "S3_trajectory.csv"
T_S4_AUC          <- "S4_validation_auc.csv"
T_S4_DECONV       <- "S4_deconv_proportions.csv"
T_S4_SCORES       <- "S4_module_scores.csv"
