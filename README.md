# Recorded study metadata and patient-level single-cell classification

This repository accompanies the PLOS Computational Biology submission:

**Recorded study metadata predicts the diagnosis in eight of nine public single-cell case-control comparisons**

The current submission snapshot is dated **2026-09-08**. It supersedes the earlier
`v2.0.0-rc3` identifiability-framed manuscript while retaining that historical release
for provenance.

## Current scientific question

Patient-level single-cell classifiers can report very high internal AUC even when the
recorded details of recruitment, processing and sequencing also predict the diagnosis.
The current study asks a narrower empirical question: how strongly does recorded
metadata predict diagnosis in a cohort, which variable groups carry that association,
and what can common validation procedures establish when collection structure and the
label are entangled?

The manuscript does **not** claim a general identifiability theorem and does **not**
present the checks as a required validity standard. They are complementary diagnostics
with explicit scope and failure modes.

## Current study

The main screen contains nine case-control comparisons from five public datasets
(809 donors), spanning lupus, COVID-19, influenza, cytomegalovirus infection, and a
sepsis-versus-COVID-19 contrast.

Key results:

- recorded metadata alone predicts diagnosis significantly in eight of nine comparisons;
- collection variables alone are significant in four of seven comparisons that record a varying collection block;
- a prespecified CMV cohort is the negative control, with metadata-only prediction at chance while the diagnosis classifier remains strong;
- an influenza comparison is not estimable by the collection-preserving permutation because collection strata contain only one label;
- free and collection-preserving complete-pipeline permutations disagree in the sepsis-versus-COVID-19 comparison;
- two COVID-19 cohorts with similar overall metadata predictability are affected through different variable groups, sample quality in one and recruiting hospital/sub-study in the other;
- downstream lupus analyses distinguish representation-level batch information from actual classifier use of that information;
- residualisation and restriction answer different questions and cannot be interpreted as bounds on a biological effect;
- external validation only tests design confounding to the extent that the target acquisition mechanism is genuinely different.

The simulation study covers seven data-generating regimes, seven levels of
design-diagnosis association, binary and continuous outcomes, and **19,600 simulated
cohorts**.

## Submission snapshot

The submission-stage manuscript source and cover letter are under:

`ploscb_submission_2026-09-08/`

The line-numbered, double-spaced submission PDF is generated from that snapshot. The
same snapshot is also packaged for archival release outside GitHub so that figures,
supporting tables, staged scripts, result tables and QA records remain version-locked.

## Reproducibility and scope

All supervised analyses use the donor as the independent unit. Cells from one donor do
not cross training/test partitions. The repository records cohort metadata, design
variable inventories, analysis scripts, simulation outputs, per-seed results and figure
source code.

The checks only see **recorded** structure. A low metadata-only AUC does not rule out
unrecorded collection effects. A significant collection-preserving permutation result
rules out the recorded collection strata used for restriction, not every possible
technical, demographic or sample-quality explanation. None of the internal checks alone
establishes biological attribution or clinical usefulness.

## Historical snapshots

- `release_v2.0.0-rc3/` - previous PLOS Computational Biology review object, retained unchanged for provenance.
- `release_v1.0.2/` - earlier SLE representation benchmark snapshot.

The public package name `rheumlens` is retained for repository continuity; it is not the
name of a method introduced by the current manuscript.
