# Figure contracts

## Figure ID1: simulated identifiability landscape

- Question: how do design-label association, biological signal and technical
  signal change the information available to common validity diagnostics?
- Takeaway: design-label overlap controls estimability; residualisation and
  restriction cease to be interchangeable as overlap collapses, while external
  performance separates biological from design-only signal.
- Forms: line with interval, two heatmaps, coefficient-error line.
- Data: 11,200 simulation replicates over 112 parameter cells.
- Palette: blue for information/external performance, orange for attenuation,
  neutral references; line style and direct labels supplement color.
- Exports: SVG, PDF and 300-dpi PNG.

## Figure ID2: empirical null and composition extension

- Question: are the high design-only AUCs beyond a full-pipeline small-sample
  null, and does a same-cell-set composition representation reproduce the
  residualisation/restriction discrepancy?
- Takeaway: the complete GSE135779 design model is compared with explicit nulls;
  CD4 composition supplies an independent representation family and a direct
  test of the same validity pattern.
- Forms: horizontal violin-plus-observed-dot, paired metric plot, AUC interval
  plot, restriction-minus-matched interval plot.
- Data: GSE135779 44 donors; GSE174188 261 CD4 donors.
- Palette: blue for observed/unadjusted, orange for residualised, grey nulls.
- Exports: SVG, PDF and 300-dpi PNG.
