# Analysis log (supplementary, English)

Project: dermatology-sci. Public-data re-analysis of psoriasis skin.

## Stage order

| Stage | Script | Output |
|---|---|---|
| S0 | `py/S0_pilot.py` | `S0_*.csv`, Gate A report |
| S1 | `R/01_id_map.R`, `R/02_s1_qc.R`, `R/03_s1_integrate.R` | `S1_*.csv`, integrated object |
| S1 | `R/04_gateB.R` | Gate B selection of the myeloid/DC compartment |
| S2 | `R/05_s2_subcluster.R`, `R/05b_s2_stats.R` | `S2_*.csv`, subclusters, proportions, trend tests |
| S3 | `R/S3_function.R` | `S3_*.csv` |
| S4 | `R/S4_transfer.R`, `R/S4_bulk.R`, `R/S5_gateD.R` | `S4_*.csv`, Gate D |
| S5 | `R/06_figures.R` | Figures 1-4 |
| S8 | `R/S8_bulk_orthogonal.R`, `R/S8_lc_subsets.R`, `R/S8_figures.R` | `S8_*.csv`, Supplementary Figures S1-S2 |
| S9 | `R/S9_prop_zero.R` | `S9_prop_*.csv` (revision) |

Every R script sources `config/00_paths.R` first and uses the single seed
`20261112`. All scripts are ASCII-only.

## Decision gates (frozen before analysis)

- **Gate A** (data usability): PASS. 96,088 cells / 33,538 Ensembl identifiers in
  GSE173706; GSE162183 merged matrix with 24,234 cells; GSE13355 180 arrays.
- **Gate B** (target compartment): PASS. A five-criterion score produced a
  three-way tie at 7 points; the pre-registered tie-break selected myeloid/DC
  (2,314 cells).
- **Gate C** (ordered trend): the frozen criterion is "at least one QC-eligible
  subcluster with Jonckheere-Terpstra BH-FDR < 0.05 and a monotone direction".
  On the originally released proportion table the smallest BH-FDR was 0.222 and
  the gate was recorded as failed. The revision reported here found that the
  table omitted pseudo-samples with zero cells in a subcluster; on the corrected
  29 x 3 grid the smallest BH-FDR is 0.031, for the 128-cell cycling subcluster,
  and 0.131 once pseudo-samples with fewer than 10 myeloid cells are excluded.
  The verdict is therefore specification-dependent and is reported as such; the
  conservative two-group framing is retained.
- **Gate D** (cross-cohort validation): PASS on the frozen criterion
  (internal cross-validation AUC >= 0.70 **or** direction consistent, plus a
  consistent bulk direction). The AUC is internal, not external, validation.
- **Gate E** (route freeze): PASS.

## Corrections recorded during the project

1. A within-myeloid drift at p = 0.0337 seen on coarse S1 clusters was driven by
   a sweat-gland/secretory epithelial contaminant (648 of 2,314 cells, 28%) and is
   **retracted**. It is never cited as supportive evidence.
2. GSVA in subcluster 0 returned 0/50 at BH-FDR < 0.05, but with 3 versus 11
   pseudo-samples the smallest attainable BH-FDR across 50 sets is about 0.51.
   The analysis is **uninformative**, not negative.
3. The proportion analysis omitted structural zeros. `R/S9_prop_zero.R` rebuilds
   the complete 29 x 3 grid and re-runs the same tests without changing any
   threshold, alternative or model formula.

## Stage S13: transcription-factor activity (post hoc, exploratory)

`R/S13_tf_activity.R` scores TF activity with decoupleR 2.12.0 `run_ulm` over the
dorothea human regulon, A+B (118 TFs) as primary and A+B+C (271 TFs) as
sensitivity, at two layers: all QC-eligible myeloid cells pooled, and subcluster
0 alone. Tests are donor-level (Kruskal-Wallis across three states; paired
signed-rank across the 10 donors with both a PN and a PP pseudo-sample; BH within
layer x regulon). Because activity scores correlate with pseudo-sample depth, a
TF is called robust only if paired BH < 0.05, the paired depth-adjusted intercept
BH < 0.05, **and** |rho| < 0.6 for the activity-difference vs depth-difference
correlation.

Result: 23 of 111 TFs are significant on the paired test in the pooled A+B
layer, and 12 of those are robust (IRF1, STAT1, STAT2, STAT3, RELB, NFKB1,
NFKB2, JUND, RARA, STAT5A, STAT5B, NR2F2). At A+B+C confidence 22 are
significant and **none** robust; in subcluster 0 nine are significant and
**none** robust. The layer therefore supports a directional statement
(interferon and NF-kB/AP-1 activity higher in lesional skin) and not a
regulatory mechanism. FOS and RELA pass the pair of statistical criteria but are
removed by the depth screen (both |rho| = 0.64); ARNTL is a control showing why
the screen exists (paired p = 0.11, depth-adjusted p = 0.004, rho = 0.75).

## What is deliberately not claimed

- No external validation from the label-transfer AUC.
- No three-state progressive trajectory.
- No mechanistic or clinical conclusion from bulk signature scores.
- No cell-level significance test; the inference unit is always the donor or
  sample.
- No regulatory mechanism from TF activity: mRNA-inferred, 118-TF coverage,
  specification-dependent.
