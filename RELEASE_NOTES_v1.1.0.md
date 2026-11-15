# v1.1.0 — figure-text corrections

Release date: 2026-11-16

Supersedes v1.0.0. The analysis, the statistics and every number are unchanged;
only the text set inside the figures, and the dpi reported in one console message,
differ. This release exists so that the archived code reproduces the figures that
were actually submitted.

## Changes

1. **Gene symbols inside figures are now italic** — panel axes and labels:
   Fig 1B (lineage-marker dot plot), Fig 2B (top-5 markers), Fig 2D (maturation
   facet strips), Fig 3A (volcano labels), Supplementary Fig S2A (single-gene
   panel), S4 (maturation facet strips), S6 panels A, B and D (transcription
   factors). Gene symbols were already italic in the manuscript body and in the
   figure legends; the raster text now matches.
2. **Figure-level titles removed from inside the images** — "Figure 1." to
   "Figure 4." and "Supplementary Figure S1./S2." were drawn into the artwork and
   duplicated the legends that the journal typesets. Panel titles and panel tags
   (A, B, C, D) are unchanged.
3. **American spelling inside figure text**, matching the manuscript:
   "Grey" to "Gray", "Normalised enrichment score" to "Normalized enrichment score",
   "upstream-normalised" to "upstream-normalized".
4. **`R/S8_figures.R` console message** reported "(300 dpi)" while writing the PNG
   at `FIG_DPI` (600). The message now prints the value actually used. No output
   changed.

## Verification performed

- `py/figure_lint.py`: 10/10 PASS (min font 6.0 pt, widths 179.9-355.6 mm,
  PNG 600-601 dpi).
- All ten PNG renders are 600 dpi; the two patchwork RGBA renders are flattened to
  RGB at submission time (transparent objects are not accepted by the journal).
- No analysis script was re-run: the figure scripts read the frozen RDS objects and
  result tables only, and no number changed.

## Known cosmetic issue

`patchwork` 1.3.0 emits `annotation$theme is not a valid theme` when a caption
theme is passed to `plot_annotation()`; the theme is applied, the warning is
upstream noise, and it is left visible rather than suppressed.

## Files touching line endings

`R/06_figures.R` and `R/S13b_tf_figures.R` previously had mixed CRLF/LF endings;
they are now uniformly CRLF apart from the final line, which is why their diff
appears as a whole-file rewrite. The substantive change is 22 and 6 lines
respectively.
