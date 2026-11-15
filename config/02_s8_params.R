## 02_s8_params.R - ADDITIONAL parameters for the S8 supplementary analyses.
## Created 2026-11-12 after the user authorised route B (P0 + P1).
## IMPORTANT: config/01_params.R is the FROZEN v1.1 parameter file and is
## NOT modified by S8. S8 is a separate, additive layer: all new outputs use
## the S8_ prefix, no frozen S1-S4 table is overwritten, no figure is redrawn,
## and no frozen parameter (resolution, one-sided alternative, gates) changes.
## ASCII only (this machine starts R with LC_CTYPE=C).

S8_VERSION <- "S8.1"

## P1-1: single genes used as an annotation-independent orthogonal readout.
## Chosen because they are canonical LC / DC surface markers present on GPL570,
## not because of any observed group difference.
S8_SINGLE_GENES <- c("CD207", "CD1A", "FCER1A", "LAMP3", "CD1C")

## P1-2: two additional independent GPL570 bulk cohorts downloaded 2026-11-12
## into the shared data catalogue (projects/data/manifest.csv) and linked into
## data/raw/. Both use the SAME 10-gene signatures as the frozen S4 analysis.
S8_EXTERNAL <- list(
  GSE14905 = list(dir = "PSOR_BULK_GSE14905",
                  file = "GSE14905_series_matrix.txt.gz",
                  has_NN = TRUE,
                  note = "82 arrays, 54 subjects: 21 NN, 28 paired PN/PP, 5 PP-only"),
  GSE30999 = list(dir = "PSOR_BULK_GSE30999",
                  file = "GSE30999_series_matrix.txt.gz",
                  has_NN = FALSE,
                  note = "170 arrays, 89 subjects; 81 matched NL/LS pairs (162 arrays) + 8 single biopsies; no healthy controls")
)

## P1-3: LC subset signatures traced to Liu et al., Immunity 2021;54(10):
## 2305-2320.e11 (doi 10.1016/j.immuni.2021.08.012), whose LC1/LC2/aLC/migLC
## surface phenotypes were verified by literature search on 2026-11-12.
## These are DESCRIPTIVE scores only. The CD207+ subcluster has 1 NS vs 3 PP
## donor pseudo-samples, so the frozen plan forbids a group-level test here.
S8_LC_SETS <- list(
  LC1_classical = c("CD207", "CD1A", "EPCAM"),
  LC2_monocytic = c("CD1C", "CD1B", "CD14", "ITGAM", "CLEC10A", "CD63",
                    "CD40", "CD274", "TNFRSF11A"),
  aLC_migLC     = c("CD83", "CCR7", "HLA-DRA", "HLA-DRB1")
)

S8_ALPHA <- 0.05
S8_FDR   <- "BH"
