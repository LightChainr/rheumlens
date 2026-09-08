# PLOS Computational Biology submission snapshot - 2026-09-08

Submission-stage record for:

**Recorded study metadata predicts the diagnosis in eight of nine public single-cell case-control comparisons**

This snapshot corresponds to the manuscript prepared for a new PLOS Computational Biology submission after substantial reframing and expansion of PCOMPBIOL-D-26-01856.

## Current scope

- nine case-control comparisons from five public main-screen datasets (809 donors);
- lupus, COVID-19, influenza, cytomegalovirus infection, and sepsis-versus-COVID-19;
- design-only AUC and cross-fitted residual label variance against cohort-specific permutation nulls;
- free and collection-preserving complete-pipeline permutations;
- seven simulation regimes, seven design-diagnosis association levels, two outcome types, and 19,600 simulated cohorts;
- downstream lupus analyses of whether representations contain and classifiers use recorded collection structure;
- explicit limitations for unrecorded structure, sample-quality/demographic association, residualisation, and external validation.

The manuscript no longer claims an identifiability theorem and does not present the checks as a required validity standard; they are complementary diagnostics with stated scope.

## Submission PDF

A line-numbered, double-spaced, 72-page PDF with the eight main figures embedded,
generated with XeLaTeX on 2026-09-08 and committed at `01_UPLOAD/Manuscript.pdf`.

SHA256: `3ddddf09e2016386b5a0ccd01c87c2d9ff6d3f3de67e9f71de3435b3d6498eae`

**This PDF is out of date.** It was built before the corrections in the commit that paired every printed AUC with the p-value that tests it, marked the two sample-quality-preserving permutations, and updated the affected numbers in the text. Rebuild it from `01_UPLOAD/Manuscript.md` before submitting, and update this hash and `SHA256SUMS` in the same commit.

Every other file in this directory is covered by `SHA256SUMS`, so there is one place
to check rather than a hash per artefact that can go stale on its own:

```bash
shasum -a 256 -c SHA256SUMS
```

The repository metadata in this branch is prepared for the next versioned archival release (`3.0.0-rc1`).
