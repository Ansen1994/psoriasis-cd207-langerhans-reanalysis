# RELEASE NOTES v1.2.0 (2026-11-19)

Post hoc supplementary layers II-III added before journal submission. Additive
only: no frozen table, parameter, decision gate or main figure was modified.

## Added
- `R/S18_deconvolution.R` - NNLS deconvolution of GSE13355/GSE14905/GSE30999
  against the four frozen myeloid-subcluster signatures; estimated CD207+ share
  falls 26.4/26.6 -> 3.7 percent in GSE13355 (PP vs NN BH-FDR 2.3e-20; 58 paired
  PN-PP, signed-rank p = 3.6e-11). Direction replicates in the primary cohort
  only; disclosed.
- `R/S19_treatment_bulk.R` - GSE226244 secukinumab cohort (GPL570): single-gene
  CD207, the frozen LC/CD207 signature, deconvolved shares and LAMP3 all move
  back toward healthy levels after treatment (POST vs PRE Wilcoxon FDR 4.8e-05).
  Unpaired by design (no patient identifier deposited).
- `R/S20_treatment_sc.R` - GSE183047 immune-cell-enriched scRNA-seq, lean
  per-sample pipeline; descriptive only (4 paired donors < 6, no test).
- `R/S23_durability_pasi.R` - recovery durability by treatment timeline (score
  0.589 -> 0.709 (wk12) -> 0.725 (wk24); maintained) and PASI correlation
  (Spearman rho = -0.03, p = 0.91, n = 17 lesional samples; descriptive null).
- `R/S21_new_figures.R` - Supplementary Figures S7-S9 and the graphical abstract.
- `py/S22_package_check.py` - post-rebuild package verification.

## Changed
- `py/S8_build_supplementary.py` - ships S17-S23 tables, S7-S9 figures and logs.
- `.zenodo.json`, `CITATION.cff`, `README.md` - version metadata, seven cohorts.
