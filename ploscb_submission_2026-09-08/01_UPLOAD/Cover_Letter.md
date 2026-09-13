Dear Editors,

We are submitting **"Recorded metadata predicts the phenotype label in eight of nine
public single-cell cohort comparisons"** as a new Research Article. It follows your
decision on PCOMPBIOL-D-26-01856 (*Study design predicts disease and defines the
identifiability boundary for patient-level single-cell classifiers*) and the editor's
note that a submission which refined the framework, recalibrated the claims and expanded
the validation would be welcome. A point-by-point account of how each editorial and
reviewer comment was addressed is attached as a separate file.

The manuscript no longer claims an identifiability boundary and no longer proposes a
required validity standard. It asks an empirical question instead: across public
patient-level single-cell comparisons, how strongly do the recorded metadata predict the
phenotype label, which groups of variables carry that association, and what do the usual
validation procedures establish once it is present?

**Scope.** The evidence moves from one disease to nine phenotype contrasts from five
public blood datasets (809 donors): lupus, COVID-19, influenza, cytomegalovirus
serostatus and sepsis versus COVID-19, with two further lupus cohorts for mechanism.

**Findings.** Metadata-only classifiers predicted the phenotype label in eight of the nine
comparisons, in all five split seeds and after Benjamini-Hochberg correction within each
seed. Collection variables alone did so in four of the seven comparisons that record
them, and two COVID-19 cohorts with similar overall metadata AUC were driven by different
variable groups. The cytomegalovirus comparison, specified in advance as a negative
control, was the one exception. Adding expression to the metadata raised AUC by 0.087 to
0.354 in seven comparisons and by nothing measurable in the other two.

**Validation procedures answer different questions.** We run a free label permutation
and a collection-stratified permutation, the latter following Neto et al. (2019), and
state the null each one tests. In sepsis versus COVID-19 the free test rejects in every
seed while the stratified test rejects in none; in influenza every collection stratum is
label-pure, so the stratified test is uninformative rather than non-significant. A
simulation shows that when a recorded variable outside the conditioning set drives both
phenotype and expression, the stratified test correctly rejects its null although the
generating model contains no disease effect, so passing it does not certify disease
biology. At zero strength both tests reject at close to the nominal 5% (6.3% and 7.0% of
300 runs).

**Adjustment is not a contamination measure.** Residualising on a wide, ragged metadata
matrix with no phenotype information costs 0.363 AUC on average in a simulation that
reuses the real cytomegalovirus layout, which contains the 0.326 loss observed in that
negative control.

Every fitted pipeline in the paper now runs through one shared implementation in which
imputation, encoding, feature filtering, PCA and residualisation are fitted on training
donors only, and the hyperparameter table is generated from those pipeline objects rather
than written by hand.

We think the work suits PLOS Computational Biology because it sits between computational
genomics, machine learning and study design, and because its practical output is a
reporting change for patient-level classifiers: metadata-only discrimination by variable
group, the increment of expression over that metadata, permutation tests matched to the
question, and external evaluation on a cohort collected differently. The permutation
methods themselves are established; the contribution is the cross-disease evidence of
how often and in what form the problem appears in public single-cell data.

All data are public. The analysis code, cohort registry, machine-generated metadata
manifest, per-seed outputs, simulation output and the checks that bind manuscript numbers
to those outputs are released with the manuscript.

We have no competing interests to declare.

Yours sincerely,

Dan Liu, on behalf of the authors
Department of Rheumatology and Immunology, Shanghai Pudong Hospital,
Fudan University Pudong Medical Center
