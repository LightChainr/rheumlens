# PLOS Computational Biology submission snapshot (prepared 2026-09-08, sealed 2026-09-13)

Submission-stage record for:

**Recorded metadata predicts the phenotype label across public single-cell cohorts and changes what internal validation can establish**

This snapshot corresponds to the manuscript prepared for a new PLOS Computational Biology submission after substantial reframing and expansion of PCOMPBIOL-D-26-01856.

## Current scope

- nine phenotype contrasts from five public main-screen datasets (809 donors);
- lupus, COVID-19, influenza, cytomegalovirus infection, and sepsis-versus-COVID-19;
- metadata-only AUC by variable group, the increment of expression over metadata, and cross-fitted residual label variance against cohort-specific permutation nulls;
- free and collection-stratified complete-pipeline permutations, with the null of each stated;
- seven simulation regimes, seven design-diagnosis association levels, two outcome types, and 19,600 simulated cohorts;
- downstream lupus analyses of whether representations contain and classifiers use recorded collection structure;
- explicit limitations for unrecorded structure, sample-quality/demographic association, residualisation, and external validation.

The manuscript no longer claims an identifiability theorem and does not present the checks as a required validity standard; they are complementary diagnostics with stated scope.

## The screen was re-run fold-contained

The permutation pipeline used to fit its top-variance gene filter and its PCA
once on all donors and then cross-fit only the classifier. Neither step reads the
diagnosis, so no label leaked and the observed statistic and its null were always
computed on identical terms - but both saw the held-out donors' expression, so
the reported AUC was transductive rather than out-of-fold, and "the whole
pipeline is refitted on every permutation" overstated what the code did.

Both steps are now inside the training fold, and all 45 (comparison x seed) runs
were regenerated. What this did and did not move:

- **unchanged, bit for bit:** every design-only quantity - `I_D`, `I_D_cv`,
  `design_auc_linear`, `design_auc_rf`, `design_auc_frozen`, `p_design_auc`,
  `n_strata`, `strata_definition` - because they never used that representation.
  The headline counts are the same: recorded metadata significant in 8 of 9,
  collection variables alone in 4 of the 7 that record any, and the CMV negative
  control still at chance. The degenerate and underpowered verdicts are the same
  comparisons as before.
- **changed:** the diagnosis-classifier AUC and the two permutation p-values.
  The free permutation now rejects at the 1/1001 floor in all nine comparisons.
  The collection-preserving permutation still fails to reject on exactly one,
  sepsis-versus-COVID-19 in COMBAT (p = 0.057 to 0.614 across seeds).
- **resolved:** the frozen and the tuned classifier used to disagree by up to
  0.518 AUC on the smallest comparison. They now agree to within 0.066 everywhere.

Before the re-run, the *previous* version of the script was run in the pinned
environment on the build container and reproduced the previously released
`design_screen.tsv` byte for byte, so the differences above are attributable to
the code change and not to the machine.

## Submission PDF

A line-numbered, double-spaced, 60-page PDF with all 15 figures embedded at
315-436 ppi, committed at `01_UPLOAD/Manuscript.pdf`.

SHA256: `e75e5c03173160b911ccc4ee8337a453918227b5a623a792b2f853183b6295d3`

It is built from `01_UPLOAD/Manuscript.md` by
`02_SOURCE/analysis_scripts/build_pdf.py`, which takes its figure map from
`build_html.py` so the PDF, the HTML and the TIFFs cannot disagree about which
file is Figure 3. The previous PDF was assembled from a hand-maintained LaTeX
copy of the body that had already drifted from the markdown, and that copy has
been removed.

Every other file in this directory is covered by `SHA256SUMS`, so there is one place
to check rather than a hash per artefact that can go stale on its own:

```bash
shasum -a 256 -c SHA256SUMS
```

The repository metadata in this branch is prepared for the next versioned archival release (`3.0.0-rc1`).
