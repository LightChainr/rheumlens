Dear Editors,

We are submitting **"Recorded study metadata predicts the diagnosis in eight of nine
public single-cell case-control comparisons"** as a new Research Article, following your
decision on PCOMPBIOL-D-26-01856 and your note that a substantially reframed submission
would be welcome.

The manuscript no longer frames the work as an "identifiability boundary", and no longer
presents its checks as a required validity standard. Both framings exceeded the evidence,
and both are gone. In their place is a narrower empirical question: how strongly does
recorded study metadata predict the diagnosis in public patient-level single-cell
comparisons, which recorded variables carry that association, and what can the usual
validation procedures actually establish in that setting?

Four things are substantively new.

**Scope.** The evidence base moves from lupus case studies to nine case-control
comparisons drawn from five public datasets (809 donors), spanning lupus, COVID-19,
influenza, cytomegalovirus infection and a sepsis-versus-COVID-19 contrast. A
cytomegalovirus cohort was chosen in advance as a negative control, and it is the one
comparison the screen does not flag.

**A separation we now keep.** Collection variables are study design in the strict sense.
Demographics are case mix and may be genuine risk factors. Sample-quality summaries sit
downstream of the assay, and can sit downstream of the disease as well. We report the
three groups separately, treat collection as primary for any claim about design, and give
the headline count under both definitions.

**Two permutation questions, not one.** A free complete-pipeline permutation tests
pipeline validity and leakage. A collection-preserving permutation keeps the observed
association between diagnosis and collection strata while removing disease-specific
molecular structure, and asks whether accuracy exceeds what those strata alone produce.
The two disagree on one of our nine comparisons, which is the point of running both.

**The new test's own failure mode, measured rather than caveated.** We simulate a
confounder that survives *inside* the collection strata, with no disease effect anywhere
in the generating model. Both permutation tests are calibrated when that confounder is
absent, and both reject in up to 89% of runs when it is present. A significant
collection-preserving result is therefore evidence against the collection strata, not
evidence of biological signal, and we say so in the Results rather than in a caveat.

The work suits PLOS Computational Biology because the problem sits between computational
genomics, machine learning and study design: a patient-level classifier can reach very
high internal AUC while learning acquisition structure specific to one cohort, and
nothing in a cross-validated number distinguishes the two. The analyses give developers
and readers a cheap way to diagnose that ambiguity before reading leaderboard performance
as disease prediction.

We think it can move the field by changing what is reported alongside a patient-level
classifier: metadata-only predictability by variable group, permutation tests matched to
the question being asked, an explicit statement when the restricted null is undefined
rather than a p-value that reads as "no signal", and external validation on a cohort
whose acquisition process is genuinely distinct from the development cohort.

All data are public. We release the analysis code, the cohort registry with API-verified
donor counts, a machine-generated manifest of every design column entering every model,
all per-seed outputs, the simulation and calibration output, and the checks that verify
the manuscript against those outputs.

We have no competing interests to declare.

Yours sincerely,

Dan Liu, on behalf of the authors
Department of Rheumatology and Immunology, Shanghai Pudong Hospital,
Fudan University Pudong Medical Center
