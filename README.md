# psoriasis-cd207-langerhans-reanalysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24249682.svg)](https://doi.org/10.5281/zenodo.24249682)

Analysis code for:

> **Reduced CD207+ Langerhans-cell transcriptomic signature in lesional psoriasis:
> a transparency-focused two-group single-cell re-analysis with replication across
> three independent bulk cohorts**
> Huiyu Zhang, Dongmei Shi

The repository contains the R and Python code that reproduces every number, table
and figure of the manuscript from five public GEO series. It does **not** contain
data, intermediate objects or derived result tables; those are downloaded by
`scripts/00_fetch_data.sh` or shipped in the article's supplementary archive.

## What the analysis does

1. **S0** – feasibility pilot on the public matrices (Gate A).
2. **S1** – identifier harmonisation, cell-level QC, doublet removal, Harmony
   integration, lineage annotation, and selection of the myeloid/DC compartment
   (Gate B).
3. **S2** – subclustering of the myeloid/DC compartment, donor-level pseudobulk
   proportions, the pre-specified one-sided Jonckheere–Terpstra permutation test
   with Benjamini–Hochberg correction, and contamination QC.
4. **S3** – two-group pseudobulk differential expression, GO and Hallmark
   enrichment, GSVA.
5. **S4/S5** – label transfer, two independent bulk cohorts, composition
   controls, and the frozen Gate D criterion.
6. **S8–S16** – post hoc supplementary layer (single-gene bulk checks, LC-subset
   scores, structural-zero correction, composition and empirical null,
   maturation axis, transcription-factor activity). These analyses are labelled
   post hoc in the manuscript and do not enter any pre-specified gate.
7. **Figure scripts** render Figures 1–4 and Supplementary Figures S1–S6 from the
   frozen objects; no statistic is computed inside a figure script.

All inference is performed on donor/sample-level pseudo-samples, never on cells.
The analysis plan, thresholds, clustering-resolution range, ordered-trend
alternative and gate criteria were frozen before any biological result was
inspected; every correction applied after that freeze is documented in
`logs/ANALYSIS_SUMMARY.md`.

## Layout

| Path | Contents |
|---|---|
| `R/` | 25 analysis and figure-rendering scripts |
| `py/` | feasibility pilot, mapping-free marker counting, DOI verification, and the Markdown→Word converter used for the submission files |
| `config/` | `00_paths.R` (relative paths only), frozen `01_params.R`, additive `02_s8_params.R` |
| `scripts/00_fetch_data.sh` | downloads the five GEO series, verifies SHA-256 checksums, builds the directory layout the scripts expect |
| `logs/ANALYSIS_SUMMARY.md` | stage order, decision gates, and every correction recorded after the freeze |

## Data

| Series | Role |
|---|---|
| GSE173706 | discovery, three-state single-cell RNA-seq (96,088 cells before QC) |
| GSE162183 | independent whole-skin single-cell RNA-seq (3 healthy vs 3 lesional) |
| GSE13355 | primary bulk microarray cohort (180 arrays, 122 subjects) |
| GSE14905, GSE30999 | two further GPL570 bulk cohorts (post hoc supplementary layer) |

All five series are public and de-identified; no new samples were collected and no
wet-laboratory experiment was performed.

## Reproducing

```bash
bash scripts/00_fetch_data.sh          # about 2 GB; verifies SHA-256
Rscript R/01_id_map.R                  # then run the stages in the order given in
                                       # logs/ANALYSIS_SUMMARY.md
```

Every R script sources `config/00_paths.R` first and uses the single seed
`20261112`. Scripts are ASCII-only and must be run from the repository root.
GPL570 annotation is not redistributed: `scripts/00_fetch_data.sh` prints the
exact GEO annotation file to download.

## Environment

See `ENVIRONMENT.md` (R 4.4.2; Seurat 5.3.0; Bioconductor 3.20).

## Scope and known limits

This is a secondary analysis of public cohorts. The label-transfer AUC is
**internal** cross-validation, not external validation; the 3-versus-3 direct
marker counts are descriptive; GSVA at 3 versus 11 pseudo-samples is
uninformative rather than negative; and the ordered three-state trend is
specification-dependent. These limits are stated in the manuscript and are not
softened here.

## Citation

If you use this code, please cite the manuscript and the archived release of this
repository. See `CITATION.cff` for machine-readable metadata.

- Version DOI (this release, v1.0.0): **10.5281/zenodo.24249682**
- Concept DOI (always resolves to the latest version): 10.5281/zenodo.24249681

```
Zhang H, Shi D. psoriasis-cd207-langerhans-reanalysis (v1.0.0).
Zenodo; 2026. doi:10.5281/zenodo.24249682
```

## License

MIT (see `LICENSE`).
