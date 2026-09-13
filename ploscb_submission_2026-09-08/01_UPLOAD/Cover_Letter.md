Dear Editors,

We are submitting **"Recorded metadata predicts the phenotype label across public single-cell
cohorts and changes what internal validation can establish"** as a new Research Article. It follows your
decision on PCOMPBIOL-D-26-01856 (*Study design predicts disease and defines the
identifiability boundary for patient-level single-cell classifiers*) and the editor's
note that a submission which refined the framework, recalibrated the claims and expanded
the validation would be welcome. A point-by-point account of how each editorial and
reviewer comment was addressed is attached as a separate file.

The manuscript no longer claims an identifiability boundary or proposes a required
validity standard. Its central point is that patient-level single-cell validation is not
solely a property of the classifier: it depends jointly on the phenotype-acquisition
structure of the cohort and on which parts of that structure the validation procedure
preserves.

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

**Validation procedures answer different questions.** A free label permutation and a
collection-stratified permutation (Chaibub Neto et al., 2019) test different nulls. The
comparisons fall into three regimes: no detectable metadata association (the negative
control), association with label-mixed strata in which the conditional question can be
asked (seven comparisons; in sepsis versus COVID-19 only the free test rejects), and
structural non-overlap (influenza), where no collection stratum holds both labels and the
conditional question cannot be answered with the available donors. A
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

The contribution is not another patient-level classifier or a new permutation statistic.
It is a systematic characterisation of the validation structure of public single-cell
cohorts: phenotype labels can be encoded in recorded acquisition, demographic and
sample-quality metadata; cohorts with similar overall metadata predictability can be
affected through different mechanisms; and commonly used validation procedures answer
different questions once this structure is present. These results make cohort
acquisition structure an explicit object of classifier evaluation rather than an
unexamined property of the benchmark, which we believe suits the computational biology
readership of PLOS Computational Biology.

All data are public. The analysis code, cohort registry, machine-generated metadata
manifest, per-seed outputs, simulation output and the checks that bind manuscript numbers
to those outputs are released with the manuscript.

We have no competing interests to declare.

Yours sincerely,

Dan Liu, on behalf of the authors
Department of Rheumatology and Immunology, Shanghai Pudong Hospital,
Fudan University Pudong Medical Center
