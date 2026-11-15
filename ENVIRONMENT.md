# Environment

Recorded from the `sessionInfo()` output of the runs that produced the manuscript's
figures (full snapshot: `logs/S9_sessionInfo_all.txt` in the article's supplementary
archive).

## R

```
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64
locale: C, system code page 65001
```

Single seed for the whole project: `20261112`.

## Packages

| Package | Version | Role |
|---|---|---|
| Seurat | 5.3.0 | integration, clustering, UMAP, marker detection |
| SeuratObject | 5.1.0 | object class |
| harmony | 1.2.3 | integration (sample identity only) |
| SingleR | 2.8.0 | label transfer cross-check (HPCA reference) |
| scDblFinder | 1.20.2 | doublet detection |
| SingleCellExperiment | 1.28.1 | container for SingleR / scDblFinder |
| DESeq2 | 1.46.0 | two-group pseudobulk differential expression |
| lme4 | 1.1.36 | donor-level proportion models (LRT) |
| clusterProfiler | 4.14.6 | GO biological-process enrichment |
| fgsea | 1.32.2 | Hallmark GSEA |
| GSVA | 2.0.7 | per-sample pathway scores |
| msigdbr | 24.1.0 | Hallmark gene sets |
| pROC | 1.19.1 | cross-validated AUC |
| decoupleR | 2.12.0 | transcription-factor activity (`run_ulm`) |
| dorothea | 1.18.0 | TF regulons (A+B, A+B+C) |
| org.Hs.eg.db | 3.20.0 | Ensembl to Symbol mapping |
| Matrix | 1.7-1 | sparse matrices |
| data.table | 1.16.4 | table handling |
| ggplot2 | 4.0.3 | figures |
| patchwork | 1.3.0 | multi-panel assembly |

## Python

Python 3.13 with `python-docx`, `Pillow` and `PyMuPDF`. Python is used only for the
feasibility pilot, mapping-free marker counting, DOI verification and the
Markdown-to-Word converter for the submission files; no statistic in the manuscript
is computed in Python.

## System requirements

About 16 GB RAM is needed for the integration step of GSE173706 (96,088 cells) and
about 3 GB of free disk for the five downloaded series plus derived objects.
