# Reference and positioning checklist

This checklist records how the manuscript should be positioned before journal submission.

## Must-cite groups

1. Single-cell foundation models:
   - scGPT.
   - Geneformer.
2. Critical foundation-model evaluation:
   - zero-shot or cross-task evaluations showing limitations of current single-cell foundation models;
   - simple-baseline comparisons in biological prediction tasks.
3. Reporting and validation standards:
   - TRIPOD+AI for prediction-model reporting;
   - DOME for supervised machine-learning validation in biology.
4. Statistical methods:
   - paired bootstrap;
   - permutation tests;
   - Holm correction;
   - Hanley-McNeil AUC interval.
5. Method families and software:
   - kernel mean embeddings;
   - attention multiple-instance learning;
   - scikit-learn;
   - PyTorch;
   - Scanpy/scverse.
6. SLE biology:
   - interferon signature literature;
   - public single-cell lupus cohort papers.

## Positioning language

Use this framing consistently:

- Contribution is the donor-level audit and reproducibility framework, not a claim that a new algorithm dominates.
- Expression baselines are competitive and sometimes statistically stronger in matched-500; across all cohorts, claims should say "competitive" or "numerically strongest" only where supported.
- Covariate sensitivity is a main result. It should not be buried as a secondary caveat.
- GSE285773 has 26 donors and should be described as hypothesis-generating.
- KME and FOCUS are included as structured method families, but the current evidence does not support broad superiority or query-specific mechanism claims.

## Pre-submission manual checks

- Replace placeholder author and affiliation fields.
- Confirm final DOI metadata after Zenodo/figshare deposition.
- Check all DOI/page information in `references.bib`.
- Verify the formal publication status of GSE285773 if a peer-reviewed study becomes available before submission.
