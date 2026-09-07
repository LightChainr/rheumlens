# Why v3 exists: the PLOS Computational Biology review

Manuscript `PCOMPBIOL-D-26-01856`, "Study design predicts disease and defines the
identifiability boundary for patient-level single-cell classifiers", was rejected
on 2026-09-03 after review by four referees. This note records what the referees
actually said, so that the v3 work does not drift back into the same claims.

## The four convergent objections

**1. The central claim was not supported (editor, R1, R2, R4).**
1-R^2(Y~D) measures how strongly *recorded* design variables correlate with the
label. That is not a formal identifiability limit. R1 asked directly which
estimand is bounded, under what assumptions, and how identifiability differs
from estimability here. R4 added that confounding requires more than D-Y
association: the design must also affect the representation, and the classifier
must exploit that component. 1-R^2 tests one link of that chain.

The project's own pre-submission plan had the correct framing and it was lost in
drafting:

> This is an information/identifiability result, **not a universal theorem about AUC**.

v3 must return to that wording. The quantity is a screening diagnostic —
necessary, not sufficient.

**2. "Five required checks" was not earned (R1, R2, R4).**
No evidence established that the checks are necessary or sufficient for the
failure modes they target. R4 proposed the fix: call them *complementary*
validity checks and state, per check, what it can and cannot establish.

**3. Check 2 was technically wrong (R4).** See `TASK_BRIEF.md`. Standard label
permutation destroys the design-label association, so it tests leakage rather
than design exploitation. The design-preserving null implemented in
`scripts/cohorts/run_design_screen.py` is the repair.

**4. One disease cannot support a general standard (editor, R2, R3).**
All three v2 cohorts were SLE, and coverage was uneven across them (the third
had 26 donors and no processing-batch information). R3 offered the alternative
explicitly: expand validation across diseases and designs, *or* narrow the title
and conclusions to an SLE case study.

## What this repository's v3 work does about it

The multi-disease screen answers objection 4 and repairs objection 3. Objections
1 and 2 are writing problems and are resolved in the manuscript, not in code.

The cohort set is deliberately a *spectrum* rather than "one more disease":
weak (CMV atlas), medium (SLE), strong (two COVID cohorts), plus COMBAT, where
COVID-19, influenza and controls share a single acquisition process so the
disease contrast can change while acquisition is held approximately constant.

If the framework is sound, `CMV_HIHA` should show that design does **not**
predict its label. A screen that flags every cohort is unfalsifiable and would
confirm the referees' suspicion rather than answer it.

## Reviewer-requested items NOT covered by this repository

These are manuscript-side and tracked in the project's response table:

- causal framing (DAG / SEM), and a definition of dataset shift vs confounding
- renaming 1-R^2 so it no longer implies Shannon/Fisher information
- estimand statements for restriction vs residualisation
- withdrawal or proof of the `[residualised, restricted]` "attribution range"
- numbered reference list, figure-number errors in Figures 4-5, panel overlaps,
  line numbers, undefined symbols in section 2.2
- ~5 missing citations named by R4 (ProtoCell4P, scPanel, HiDDEN, and others)
