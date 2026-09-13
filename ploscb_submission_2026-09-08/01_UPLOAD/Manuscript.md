# Recorded metadata predicts the phenotype label in eight of nine public single-cell cohort comparisons

**Short title:** Study metadata and single-cell patient classification

Hongyu Ying<sup>1</sup>, Dandan Yun<sup>1</sup> and Dan Liu<sup>1,\*</sup>

<sup>1</sup> Department of Rheumatology and Immunology, Shanghai Pudong Hospital,
Fudan University Pudong Medical Center, Shanghai 201399, China
<sup>\*</sup> Correspondence: Dan Liu; danliu600@126.com

---

## Abstract

Patient-level classifiers reduce thousands of single cells to one score per donor and
routinely report cross-validated AUC above 0.95. Such a score is a statement about disease
only if the variables recorded about how a donor was recruited, processed and sequenced do
not themselves predict which group that donor is in. We measured that association in nine
phenotype contrasts from five public blood datasets (809 donors), covering lupus,
COVID-19, influenza, cytomegalovirus serostatus and sepsis versus COVID-19. Recorded
metadata were split into collection, demographic and sample-quality variables, and each
group's out-of-fold AUC judged against a permutation null built from that comparison's own
metadata matrix. Two permutation tests were then applied to the expression classifier: one
permuting the phenotype label freely, one permuting it only within strata sharing a
collection configuration.

Metadata alone predicted the phenotype at out-of-fold AUC 0.404 to 1.000, significantly in
eight of nine comparisons in all five split seeds and after Benjamini-Hochberg correction
within each seed; collection variables alone were significant in four of the seven
comparisons recording any. The exception was a cytomegalovirus comparison specified in
advance as a negative control: its metadata model stayed near chance while its expression
classifier reached AUC 0.894 to 0.933. At the other extreme, a 21-donor influenza
comparison had no collection stratum holding both labels, which makes the stratified
permutation uninformative rather than non-significant. Between these ends the variable
groups separated cohorts a single number would have equated: two COVID-19 cohorts with
similar overall metadata AUC were driven by sample quality in one and recruiting site in
the other. Adding expression to the metadata raised AUC by 0.087 to 0.354 in seven
comparisons and by nothing measurable in the other two.

These analyses measure association with recorded variables. They do not establish that a
classifier's discrimination is biological.

**Keywords:** patient-level classification; single-cell RNA sequencing; study metadata;
confounding; permutation tests; external validation; systemic lupus erythematosus;
COVID-19

---

## Author summary

Single-cell sequencing of a blood sample can be turned into a diagnostic score, and
published scores often separate patients from healthy donors almost perfectly. Before
trusting one, a clinician needs to know whether it reads the disease or reads the study.
In most public datasets patients and healthy donors were not recruited the same way: they
may come from different hospitals, years or laboratory runs, and a model can score well
simply by recognising where a sample came from.

We measured how large this problem is across nine comparisons from five public datasets -
lupus, COVID-19, influenza and cytomegalovirus infection - by asking how well the recorded
study details alone predict which group a donor is in. In eight of the nine they
predict it well above chance. In the ninth, chosen in advance as a negative control, they
predict nothing while the expression classifier still separates the groups. Which details
carry the prediction differs between cohorts that look equally affected overall.

The usual test of these models shuffles the patient labels, which destroys the pattern
that needs examining. Shuffling labels only among donors collected the same way asks the
complementary question, and in one comparison the two disagree.

---

## 1. Introduction

Single-cell RNA sequencing supports patient-level classification by capturing
differences in immune-cell composition and molecular state, and published classifiers
routinely separate patients from controls at cross-validated AUC above 0.95. For such a
result to inform clinical research, the discrimination has to extend beyond the
recruitment and processing conditions of the development cohort. That requirement is hard
to assess when the phenotype label is associated with collection site, sequencing batch or
sample quality.

A growing family of methods competes to produce that one number per patient: pseudobulk
aggregation, phenotype-oriented summaries, prototype networks, sparse gene panels,
multicellular factor models, graph representations, and methods that first identify the
disease-relevant cell population [6-11,32-35]. Frozen single-cell foundation models offer
a convenient route - encode each cell once, average, train a small classifier on donors -
though benchmarks report that such embeddings can retain batch structure and can lose to
simple expression baselines [13,14]. Our concern is not which of these wins. It is that
the comparison between them is made on cohorts where the winner may be decided by
something other than biology.

That question is hard to answer in public data because cases and controls are usually
not collected the same way. They may be recruited from different clinics, sequenced in
different years, run on different chemistries, or processed in different batches. When
that happens, a model can reach a high accuracy by learning the collection process.
Batch-label confounding has been known to bias cross-validation for a decade [16], and
a fully confounded single-cell design is not identifiable [17]. Adjusting for technical
covariates does not fix it: regression removes the part of the data that design can
predict, whether that part is technical, biological, or both [18,19]. Work in
predictive modelling has shown the same thing from three directions - adjustment must
be done inside training folds, confound regression can leak information to nonlinear
learners, and residual confounding survives apparently successful correction [27-29].
Benchmark-side work raises the matching concern about how the datasets themselves are
built [12,15].

Methods to quantify confounding and to test conditional association already exist, and
we use them: the permutation approach of Chaibub Neto and colleagues [30], the
model-level confounding statistics of Spisak [29] and the Confounding Index of Ferrari
and colleagues [25] all address the same estimand from different directions. What has not
been established is what those methods return when they are applied across the public
single-cell cohorts on which patient-level methods are benchmarked, where the acquisition
structure differs from cohort to cohort and much of the relevant metadata is missing or
partial. We therefore ask four questions of each cohort. How strongly does the recorded
metadata predict the phenotype label? Which groups of variables carry that prediction?
How does the answer compare with what the same metadata matrix would produce under a
permuted label? And when is the question answerable at all?

A word on "recorded metadata", because the distinction runs through the whole paper. The
released tables mix three kinds of variable. Collection facts are study design in the
strict sense. Demographics are case mix, and may also be genuine risk factors.
Sample-quality summaries sit downstream of the assay, but can sit downstream of the
disease too. Lumping them together and calling the result "study design" would
overstate what a positive result means. We keep them apart, treat collection as primary,
and report the headline count under both definitions.

We analyse nine phenotype contrasts from five public single-cell datasets and examine
what metadata association does to donor-level validation. First, we separate collection
variables, demographics and sample-quality measures, and identify which groups predict
the phenotype. Second, we compare a free label permutation with one confined to
collection strata, to ask whether discrimination exceeds the association those strata
retain. Third, we use simulations and two detailed lupus analyses to examine why
residualisation, within-stratum restriction and external validation can support different
readings of the same cohort.

Three results follow. The association between recorded metadata and the phenotype label
is graded and measurable, and varies across routinely reused cohorts by more than half
the usable AUC scale. It can be traced to an interpretable group of variables, and
cohorts with near-identical overall metadata AUC are driven by different groups. And the
permutation test normally used to validate these pipelines states a different null from
the one at issue; we give a comparison where the two tests point in different directions,
and we measure what the collection-stratified test does and does not cover.

Section 8 states the scope of each quantity and what the approach cannot see.

## 2. Study framework and evaluation measures

### 2.1. The setting

Figure 1 sets out the situation in the notation of graphical causal models [36,37]. A
recruitment process **R** assigns each donor a place in the recorded design **D** -
site, batch, calendar period, chemistry, sample-quality profile - and is also
associated with the diagnosis **Y**, because cases and controls are usually found
through different routes. The measured cells **X** are affected by **D** through
technical effects and by **Y** through biological ones. A classifier maps X to a
prediction Ŷ.

Two further nodes matter in practice.

**C** stands for covariates that are consequences of the disease: severity, WHO score,
comorbidity, medication, time since symptom onset. These sit on the path Y → C → X.
They are not collection variables. Treating them as metadata variables converts "the
disease is predictable" into "the design is predictable", and Section 4.1 quantifies
how badly.

**U** stands for collection structure that was never recorded. Nothing in this paper
can see U. That is a limit of the data rather than of the method, and Section 8 says so.

The graph was specified before the analysis and is released in machine-checkable
`dagitty` syntax [38].

One thing the arrow between D and Y is not. We draw the association as induced by R
rather than as a directed edge. Both directions occur in practice. A patient may be
recruited at a referral centre because they are a patient, and a centre's case mix may
determine which diagnoses appear in a batch. These checks measure the association.
They do not tell you which way it runs.

### 2.2. The primary measure, and the variance summary reported beside it

**The primary measure is the metadata-only AUC.** For each comparison we fit a
classifier on the columns of the recorded metadata matrix **D** and nothing else, and
report its out-of-fold AUC. Significance comes from permuting the phenotype label 200
times and pushing each permuted label through the identical fixed-hyperparameter
pipeline, so the observed statistic and its null are the same quantity computed the same
way. Which columns are in D is stated every time it is used, and three names are kept
distinct throughout: the **metadata-only AUC** uses every recorded variable, the
**collection-only AUC** uses the collection block alone, and the **demographic-only** and
**sample-quality-only** AUCs use those blocks. The same permutation test is run for each
group, and a nonlinear arm is run alongside the linear one so that a low value is not an
artefact of assuming linearity.

**A linear variance summary is reported beside it.** Let *y* be the centred phenotype
vector. Write

  V_D = 1 − R²(y ~ D).

This is the share of phenotype variation that the recorded metadata cannot explain in a
linear model. It is a residual variance, not a Shannon or Fisher information: it is
linear, and it depends on how D is coded. Two things make it readable.

*Cross-fitting.* An in-sample R² inflates with the ratio of metadata columns to donors,
so V_D is also computed out-of-fold with five-fold cross-fitting. The correction is not
cosmetic. The CMV comparison declares 89 metadata columns for 108 donors, so its
in-sample R² is inflated to 0.532 and its in-sample V_D falls to 0.468, which would rank
it the third-most affected of the nine. Its matched in-sample null sits at 0.496: the
observed in-sample V_D is not significantly lower than its own permutation reference
distribution.

*A null built from the same matrix.* The phenotype is permuted 200 times and each
permuted label is pushed through the same five folds, giving a one-sided p-value, p(V_D),
for that specific matrix at that specific donor count. A raw V_D cannot be compared
across cohorts; its position inside its own null can. The in-sample pair is released
alongside with its own in-sample null; p(V_D) refers throughout to the cross-fitted pair,
and p(V_D, in-sample) to the other.

**Why the AUC and not V_D carries the test.** V_D is an *unpenalised* linear R², and
cross-fitting an unpenalised fit loses power as the matrix gets wider. The Ren COVID-19
comparison shows this: its 26 one-hot columns on 182 donors leave the cross-fitted V_D at
its ceiling, because the linear model explains nothing out of fold. A penalised logistic
model on the very same matrix reaches AUC 0.81 to 0.85, and its own permutation null puts
that at p = 0.005 in every seed. Reporting only p(V_D) there would say "no association"
about a comparison in which the recorded metadata plainly predicts the phenotype. So the
primary test is the permutation null for the metadata-only AUC; V_D and p(V_D) are the
linear-variance summary, and where the two disagree we say so.

**What this does and does not cover.** Both quantities look at one arrow in Figure 1, the
one between D and Y. Confounding of a classifier also requires that D affects X and that
the classifier uses that part of X. A high metadata-only AUC is therefore necessary but
not sufficient evidence that a reported expression AUC is metadata-driven. A low one does
not rule out an association through unrecorded structure. Every use of these quantities
below is subject to this.

### 2.3. Two permutation tests that ask different questions

The **free permutation** shuffles the phenotype label without restriction and refits the
whole pipeline. Its null is that the label is exchangeable across all donors, so
rejecting it is evidence of an association between expression and the label - of any
origin. It is the test these pipelines are normally validated with, and we keep it,
because a pipeline that cannot reject it is not measuring anything. What it does not do
is separate the two origins: a free permutation destroys the metadata-label association
along with the biology, so a small p-value is equally consistent with either. Nor does
rejecting it certify that the pipeline is free of leakage. A label-dependent step taken
*before* the permutation - a feature list carried in from a prior analysis, say - is
present in the observed and the permuted runs alike, and this test cannot see it. Every
representation step in the pipeline used here is refitted inside the training fold on
every permutation (Section 9.6), which is the property that has to be arranged in the
code; the permutation test does not establish it.

The **collection-stratified permutation** shuffles the phenotype label only *within*
strata that share a collection configuration - batch, pool, site and, where it varies,
sequencing chemistry - so that each stratum keeps its observed class counts and only the
correspondence between labels and molecular profiles is randomised inside it. This is the
restricted permutation scheme introduced by Chaibub Neto and colleagues [30] and used
since in the confounding literature [25]; it is not new here, and what is new is the
comparison of what the two tests return across cohorts with different acquisition
structures. Its null is conditional exchangeability given those strata, so rejecting it
is evidence of discrimination beyond what the recorded collection strata retain.

**What it conditions on, and what it does not.** The strata are built from collection
variables only. The recorded metadata also contains sample-quality and demographic
columns, and those are *not* conditioned on. So a significant result means the classifier
scored above what these collection strata retain. It does not mean the classifier beat
everything recorded, and it is not evidence that some part of the discrimination is
free of association with the recorded metadata. If a
sample-quality variable is still associated with the phenotype inside a stratum, and
drives expression, then the conditional-independence null is false for a reason that has
nothing to do with disease, and a classifier can reject it with no disease effect present
at all. Section 3.4 measures how often that happens rather than assuming it does not.
Conversely, a non-significant result is a failure to reject under this conditioning set,
not a demonstration that the discrimination is entirely collection-driven.

Conditioning on the full metadata instead would require a conditional randomisation
scheme for a matrix containing continuous columns, which is a different and harder
problem; we state the limit rather than overreach past it.

The two tests are not defined on the same data. The collection-stratified test needs at
least one stratum containing both classes. When every stratum is label-pure the
within-stratum shuffle leaves all labels unchanged, the permutation group contains only
the observed labelling, and the reference distribution is degenerate. We report the test
as **uninformative** and its p-value as NA, because it cannot distinguish the observed
result from a non-trivial reference; the free permutation is still defined in that
situation and is still reported. Section 4.4 gives a comparison where exactly this
happens.

### 2.4. Three verdicts: estimable, underpowered, and not answerable with these data

We separate two failure modes that are easy to run together and that call for
different responses.

A comparison **cannot be answered with these data** when the quantity has no meaning
for it: no collection stratum holds both a case and a control, so the
collection-stratified permutation cannot be formed. These are dropped from the main
ordering and reported separately.

A discrepancy between a tuned and a fixed-hyperparameter classifier is deliberately
**not** used as a feasibility rule. Comparing two algorithms across a cut-point would not
establish that a comparison is unanswerable, and a tuned model outperforming a fixed one
is not the same kind of problem as a comparison with no donors collected under
overlapping conditions. The
fixed-hyperparameter pipeline is the inferential statistic everywhere, because it is what
both permutation nulls are built from. The tuned AUC sits beside it as a description.
Their difference is released as a column, so a reader can judge pipeline sensitivity
directly instead of through a threshold.

A comparison is **underpowered** when the quantity is well defined but unstable: fewer
than 40 donors, or a minority class below 15. These stay in the ordering, flagged.
Deleting them would hide the instability a reader needs to see. The gates themselves are
deterministic - they depend on donor counts and stratum structure, not on the seed - but
what an underpowered comparison then reports is not: Section 4.4 gives one whose
collection-stratified p ranges from 0.057 to 0.614 over five seeds, an order of magnitude
of movement in a quantity that never crosses 0.05.

# Results

Sections 3 to 6 report results. Section 3 establishes how the checks behave on data with
a known generating mechanism. Section 4 applies them to nine phenotype contrasts.
Section 5 asks whether a classifier actually uses the structure the checks detect.
Section 6 states what each check can and cannot support.

## 3. How the checks behave in simulation

A simulation with one binary design variable, an additive design effect and constant
noise would not show whether these checks behave when those choices are relaxed. We
therefore crossed seven data-generating regimes with seven levels of metadata-phenotype
association and two outcome types, binary and continuous, each predicted and scored on
its own scale: 19,600 simulated cohorts of 200 donors and 100 features (Figure 2).
Nothing in the simulation is fitted to the real cohorts.

### 3.1. Calibration of the metadata-association statistic

When design carries no information about the diagnosis, the checks must not fire. They
do not. The fraction of runs with p(V_D) ≤ 0.05 is 0.025 to 0.065 across the seven regimes,
which brackets the nominal 0.05 (Figure 2A). It rises at ρ = 0.25 and reaches one by
ρ = 0.50 in every regime. This is the property the cross-fitted estimator and the
per-cohort null were introduced to obtain; an in-sample estimator judged against an
in-sample null does not have it.

### 3.2. Residualisation and restriction give different answers, in every regime

Two ways of handling design are in common use. Residualisation regresses the data on
the metadata variables and analyses what is left. Restriction compares only donors
collected under matching conditions. In the additive linear regime the two agree when
design carries no diagnosis information and separate as it does (Figure 2B). At ρ = 0
the unadjusted, residualised, restricted and external AUCs are 0.809, 0.809, 0.787 and
0.817. At ρ = 0.98 they are 0.879, 0.538, 0.724 and 0.891. Residualisation has removed
0.34 AUC and restriction 0.16, on data where the biological effect never changed. At
ρ = 1 restriction is not merely worse. It is undefined, because no stratum holds both
labels.

The gap between the two opens in **every** regime (Figure 2C), so it is not an artefact
of the linear setting.

The continuous-outcome arm is a genuine test of outcome type here, not a relabelling.
Its cohorts are generated with a continuous latent outcome, predicted with ridge
regression, and scored by concordance over pairs, which equals AUC when the outcome is
binary and so puts both arms on one scale. The ordering survives: restriction beats residualisation at rho = 0.98 in all seven regimes for both outcome types. The magnitude does not. With a binary outcome the gap is 0.094 to 0.200; with a continuous one it is 0.013 to 0.038. So the qualitative conclusion does not depend on outcome type, and the quantitative one does - a simulation run only on continuous outcomes would have understated the problem by roughly a factor of five.

### 3.3. Effects of class prevalence and of target acquisition structure

**What a rare diagnosis does, and what it does not.** Two things travel together in an
imbalanced cohort and are easy to confuse. One is prevalence itself. The other is that a
collection boundary drawn at the middle of a population no longer lines up with a
phenotype drawn from its top fifth, so the metadata cannot separate the labels even when
the two are perfectly linked in the latent model. We separate them with two arms. In
`imbalanced`, both thresholds sit at the 80th percentile, so prevalence is the only
thing that changed. In `imbalanced_offset`, the diagnosis sits at the 80th percentile
while the collection boundary stays at the median. The result is that prevalence alone changes nothing. With thresholds matched, the imbalanced arm tracks the balanced one almost exactly: metadata-only AUC still reaches 1.000 at rho = 1 and V_D still falls to 0.000. It is the offset arm that blunts the checks, with metadata-only AUC saturating at 0.880 and V_D flooring at 0.728. The reassurance a rare diagnosis appears to give is therefore not about how rare it is. It is about whether the collection boundaries line up with the diagnosis, and a cohort can have that problem at any case rate.

Two of our own comparisons have minority classes of 10 and 11 donors, so this matters
for reading them.

**An external cohort collected the same way will confirm a metadata-driven model.** In
the base regime the external AUC does not fall as confounding rises. It *rises*, from
0.817 at ρ = 0 to 0.904 at ρ = 1, tracking the internal AUC almost exactly (Figure 2D).
Only when the target is collected under a different design mechanism does the external
estimate come apart, flattening at 0.816-0.829 across the whole range while the
internal estimate climbs to 0.899.

The practical reading is uncomfortable. External validation checks design confounding
only to the extent that the second cohort's collection process is genuinely
independent of the first. Being a different dataset does not establish that. Two
cohorts assembled by similar consortia, on similar platforms, through similar
recruitment routes, can agree with each other and both be wrong. This sharpens Check 5
in Section 6 and is, in our reading, the most consequential thing the simulation adds.

### 3.4. A non-disease pathway outside the conditioning set makes the stratified test reject

The collection-stratified permutation conditions on collection strata. The recorded
metadata also holds sample-quality and demographic columns, which it does not condition
on. If one of those is associated with the phenotype *inside* a stratum, and also drives
expression, then the classifier can reject the stratified null with no disease effect
present anywhere. This deserves a number rather than a caveat.

We simulated it directly. Cohorts of 200 donors are assigned to eight collection sites. A
sample-quality variable QC predicts the phenotype within each site, at a strength γ we
vary. Expression depends on the site and on QC, and on nothing else - the generating
model contains no arrow from the phenotype to expression at all. Both permutation tests
then run exactly as they run on the real cohorts (Figure S1, Table S1).

Be precise about what the rejection rate is for each γ. The generating model is QC → Y
and QC → X, so conditioning on the collection strata leaves X and Y dependent: for γ > 0
the conditional-exchangeability null the stratified test states is **false**, and
rejecting it is a correct rejection of a false null, not a type-I error. Only the γ = 0
row measures calibration. At γ = 0 both tests are calibrated: the collection-stratified
permutation rejects in 19 of 300 runs, 6.3% (exact 95% CI 3.9 to 9.7; binomial p = 0.29
against a nominal 5%), and the free permutation in 21 of 300, 7.0% (4.4 to 10.5;
p = 0.11). The classifier's mean AUC over those runs is 0.500.

Above γ = 0 the rate rises with the strength of the within-stratum association and the
two tests move together, neither more protective than the other. At a mean within-stratum
absolute correlation of 0.29 the stratified test rejects in 22.0% of runs and the free
test in 22.7%; at 0.61 they reject in 95.3% and 95.7%. The classifier's own AUC rises
from 0.500 to 0.668 over the same range, with no disease effect anywhere in the
generating model.

That is the useful result, and it is a statement about interpretation rather than about
error rates. A significant collection-stratified result is evidence that the phenotype is
associated with expression beyond what the collection strata retain. It is **not**
evidence that the association is disease biology, because a recorded variable outside the
conditioning set can carry it. Where a sample-quality or demographic block is itself
associated with the phenotype - which Table 1 reports for every comparison - that
possibility has to stay open.

### 3.5. Residualisation loss without any metadata-phenotype association

Residualisation loss is often read as a measure of how contaminated a cohort was. It is
not. The size of the loss that has nothing to do with contamination can be measured, and
it turns out to be the same size as the loss observed in the cohort where nothing is
detected.

We generated cohorts in which the metadata matrix is independent of the phenotype by
construction, with a real biological effect calibrated so the unadjusted AUC lands near
0.90, and residualised on matrices of increasing width. The unadjusted and residualised
arms are the same call to the same function with one argument added, so folds, feature
filter, standardiser, PCA, classifier and scoring rule are shared and the difference
between them is the adjustment alone (Section 9.9). A second condition reuses the exact
level sizes of the real CMV batch-by-pool crossing, which is far more ragged than equal
blocks: 46 levels for 108 donors, 17 of them holding one person.

Loss appears without any metadata-phenotype association, and it is large. Its size tracks
how much of the sample the matrix can address rather than the column count alone. At one
to eight columns residualisation costs almost nothing - mean loss 0.004 to 0.023, because
removing a few directions leaves the biological one intact. Loss then rises steeply,
peaking at 0.333 for 48 columns and 0.302 for 64, and falls back to 0.222 at 80 and 0.159
at 92, the width closest to the CMV comparison's 89: past a point, extra one-hot columns
are increasingly redundant with the ones already there and ridge shrinks them together.
Run-to-run spread is wide throughout; the central 95% of runs at 48 columns spans 0.224
to 0.460.

The ragged CMV-shaped layout is the one to read. It realises 45 columns from the observed
pool sizes and costs **0.363** on average, central 95% 0.272 to 0.464 - more than balanced
blocks of any width, and enough to contain the 0.326 loss the real CMV comparison shows
(Section 5.3, Figure S2). The metadata-only AUC in these simulated cohorts stays at 0.48
to 0.51 throughout, confirming that the matrix carries no phenotype information to
remove. Residualising expression on a wide, ragged matrix over 108 donors therefore costs
about a third of an AUC point's separation whether or not there is anything to remove,
and a loss of that size is not evidence that there was.

Two limits on how far this goes. Column count is not effective dimension: the real CMV
collection block has 80 columns at centred rank 45, and equal blocks of 80 columns span a
different subspace, which is why the ragged layout exists. And the simulation fixes one
adjustment procedure, ridge at a fixed penalty on a training-fold standardised matrix; a
differently penalised residualiser would trade loss against completeness differently.


## 4. Nine phenotype contrasts from five public datasets

### 4.1. Datasets, and what the recorded metadata contains

We assembled five public datasets [2,39-42] through the CELLxGENE Discover curation API
and reduced each to one row per donor: log1p CPM pseudobulk expression, plus variables
built from schema-guaranteed observation fields (`donor_id`, `disease`, `sex`,
`development_stage`, `assay`, `tissue`) and cohort-specific collection variables
recovered from the distributed object.

**These three groups are not the same kind of thing, and we do not treat them as one.**

- **Collection** (site, sub-study, sequencing chemistry, assay, batch, pool) is study
  design in the strict sense: it describes how the material was gathered and processed,
  and nothing about it is caused by the disease.
- **Demographic** (age, sex, ethnicity) is case mix. An imbalance here is a property of
  who was recruited, so it is a design fact; but age and sex are also genuine risk
  factors for most of these diagnoses, so an association with the label can be
  aetiological rather than procedural. The two readings cannot be separated from the
  data.
- **Sample quality** (cells per donor, mean UMI, genes per cell, mitochondrial
  fraction) is measured downstream of everything. It reflects the assay, but it can also
  reflect the disease, through cell composition, treatment or the state of the sample at
  collection. Section 4.1's mediator result shows exactly how much damage that can do.

Because of this, **collection is the primary definition** wherever a result is read as a
statement about study design, and the other two groups are reported beside it as
prespecified sensitivity analyses. The union of all three is what we call the *recorded
metadata*, and every claim about the union is worded as such. Section 4.2 gives the
headline count under both definitions, because they differ.

Covariates that are consequences of the disease are excluded from all three groups by an
explicit list: severity, outcome, comorbidity, medication, symptom timing, diagnosis
fields, stage, grade, WHO score. These are the C node of Figure 1.

The exact raw columns and expanded metadata columns entering every group, for every
comparison, are generated directly from the model-matrix construction code and released
as Table S2. A machine check confirms that the released inventory reproduces the width
of every model matrix actually fitted (34 of 34 blocks) and that no diagnosis label or
label proxy appears in any of them.

This exclusion is not a formality, and we quantified what it is worth. In a simulation
where metadata and phenotype are completely unrelated, admitting a single disease-caused
covariate to the metadata variables raises the metadata-only AUC from 0.488 to 0.793. It also
flags the cohort as significantly affected in **100% of runs**, against 5.5% without it
(Supplementary Figure S3, Supplementary Table S3). Residualising on the enlarged set
then removes 0.12 AUC from a cohort in which nothing was confounded at all. A single
mediator admitted by accident is enough to manufacture the entire finding.

Nine case-control comparisons were formed (Table 1): lupus versus healthy in GSE174188
(CD4 subset and all cells), COVID-19 versus healthy in three independent cohorts,
influenza versus healthy and sepsis versus COVID-19 in COMBAT, and CMV-positive versus
CMV-negative in the HIHA cohort. The five cohorts hold 809 distinct
donors (261, 196, 120, 124 and 108). Comparisons drawn from the same cohort share
donors, so the per-comparison counts in Table 1 do not sum to a donor total.

| Comparison | Donors | Minority | Metadata AUC | p | Collection AUC | p | Expression AUC | + over metadata | p free | p strat. | Verdict |
|---|---:|---:|---:|---|---:|---|---:|---:|---|---|---|
| Influenza · COMBAT | 21 | 10 | 1.000 | 0.0050 | 1.000 | 0.0050 | 0.982-1.000 | +0.000 | 0.0010 | — | **not estimable** |
| Sepsis vs COVID-19 · COMBAT | 111 | 11 | 0.996-0.998 | 0.0050 | 0.982-0.985 | 0.0050 | 0.894-0.962 | -0.001 to +0.001 | 0.0010 | 0.0569-0.6144 | **underpowered** |
| SLE · GSE174188 (CD4) | 261 | 99 | 0.872-0.889 | 0.0050 | 0.762-0.774 | 0.0050 | 0.978-0.985 | +0.087 to +0.101 | 0.0010 | 0.0010 | estimable |
| COVID-19 · COMBAT | 110 | 10 | 0.862-0.884 | 0.0050-0.0100 | 0.500-0.504 | 0.7413-0.8657 | 1.000 | +0.116 to +0.138 | 0.0010 | 0.0010 | **underpowered** |
| COVID-19 · Stephenson | 115 | 29 | 0.845-0.882 | 0.0050 | 0.501-0.594 | 0.0149-0.2338 | 0.972-0.996 | +0.104 to +0.141 | 0.0010 | 0.0010 | estimable |
| COVID-19 · Ren | 182 | 25 | 0.812-0.854 | 0.0050 | 0.729-0.768 | 0.0050 | 0.964-0.982 | +0.100 to +0.164 | 0.0010 | 0.0010 | estimable |
| SLE · GSE174188 (all cells) | 261 | 99 | 0.819-0.832 | 0.0050 | n/a | n/a | 0.992-0.995 | +0.158 to +0.172 | 0.0010 | 0.0010† | estimable |
| COVID-19 · Ren (10x 5' v2) | 161 | 25 | 0.712-0.729 | 0.0050-0.0149 | n/a | n/a | 0.954-0.977 | +0.227 to +0.262 | 0.0010 | 0.0010† | estimable |
| CMV · HIHA | 108 | 45 | 0.404-0.518 | 0.3831-0.8806 | 0.365-0.489 | 0.5274-0.9453 | 0.894-0.933 | +0.306 to +0.354 | 0.0010 | 0.0010 | estimable |

**Table 1. The nine comparisons, ordered by how well the recorded metadata predicts the
phenotype label.** Generated from the screen output by `tools/build_table1.py`. Every
comparison was run at five split seeds; each value is the range over those five, and a
single value means all five agreed to the printed precision. "Minority" is the smaller of
the two label groups. "Metadata AUC" uses all three variable groups together;
"Collection AUC" uses the collection block alone, which is study design in the strict
sense, and n/a means the comparison records no collection variable that varies.
Every AUC is cross-fitted under the prespecified fixed-hyperparameter pipeline, which is
the statistic each p-value beside it was computed against and the one both permutation
tests are built from; the tuned nested-CV AUCs, the two column counts, V_D and p(V_D) are
in Table S4, per seed. "+ over metadata" is the change in AUC when expression is added to
the recorded metadata in one joint model, out of fold (Section 4.5, Table S5). Each p is
one-sided against a permutation null built from that comparison's own metadata matrix. The
two floors, 0.0050 and 0.0010, are 1/201 and 1/1001: the metadata-only AUC and V_D nulls
use 200 permutations, the complete-pipeline tests 1,000. A dash under "p strat." means the
collection-stratified permutation is uninformative because every stratum is label-pure. A
dagger marks the two comparisons that record no collection variable at all, where the
stratified permutation falls back to tertiles of log cells per donor: for those two it is
a sample-quality-stratified permutation, and they are the same two rows showing n/a in the
collection columns. Verdicts are read from the screen's own gates: **not estimable** when
no collection stratum holds both labels, which is the single criterion, and
**underpowered** when the minority group has fewer than 15 donors.

### 4.2. Metadata predictability in the CMV comparison, a negative control specified in advance

We chose the CMV cohort before computing anything, on its recorded structure alone.
Its 108 donors were assayed on one chemistry and one suspension type, and the CMV
status of a donor is a serological result rather than a recruitment criterion, so we
expected the collection labels to be close to independent of who is positive. If these
checks measure what they claim to measure, this cohort has to come out at the bottom.

The structure needs stating precisely, because the phrase "single-centre cohort" invites
the wrong picture. This cohort is **not** one batch: the released metadata records 36 batch
identifiers and 46 pool identifiers, and their crossing gives 46 strata, 17 of them
holding a single donor. Neither is it donor-paired in any sense this analysis uses -
each donor contributes one row with one label, and no within-donor contrast is formed
anywhere in the paper. What makes it a negative control is not a simple design but a
weak association between the collection labels and the phenotype label, which is exactly the
quantity being screened.

It comes out at the bottom. Across five split seeds the metadata-only AUC for the full
recorded set is 0.404 to 0.518, straddling chance, with permutation p 0.383 to 0.881 -
significant in no seed. Collection variables alone give 0.365 to 0.489, p 0.527 to 0.945,
and sample quality 0.453 to 0.507, p 0.478 to 0.716. The same 108 donors support an
expression classifier at AUC 0.894 to 0.933, with both permutation tests at their 1/1001
resolution floor in every seed. That is the prediction, and it holds.

One exception within the comparison. Demographic variables alone reach a metadata-only
AUC of 0.591 to 0.639, with permutation p 0.015 to 0.085; applying Benjamini-Hochberg
[24] within each seed to that seed's 34 variable-group tests leaves q = 0.020 in two of
the five seeds and 0.074 to 0.115 in the other three. So the demographic block is
significant after correction in a minority of seeds, and we report it that way rather
than as a uniform null. Sex and ethnicity are established correlates of cytomegalovirus
seroprevalence in their own right, so an association there is as easily aetiology as
recruitment. This is what the three-way split in Section 4.1 is for: it is a reason to
read that block as case mix, and a reason not to let it stand for study design. It is
also the reason the negative control is stated at the level of the collection block,
where nothing is detected in any seed.

Two further points of interpretation. First, this is the concrete case for building each
statistic's null from the same matrix rather than against a fixed reference. Read without
its own null, the in-sample V_D of 0.468 looks like a strong association - 89 columns fit
108 donors closely whatever the labels are. Against the matched in-sample null it is
0.468 against a null mean of 0.496, with p(V_D, in-sample) 0.284 to 0.353. The
cross-fitted pair, which is what p(V_D) means everywhere else in this paper, sits at its
ceiling of 1.000 with p = 1.000.

Second, no detected association is not the same as no association. The metadata-only test
has whatever power 108 donors and this matrix give it, and the matrix is wide: the
collection block declares 80 one-hot columns but its centred rank is 45, and the full
recorded block declares 89 columns at rank 54. Column count is not effective dimension,
and neither number should be read as the number of independent directions the phenotype
was tested against. What the comparison supports is that the recorded collection labels
carry no detectable information about CMV status at this sample size, which is what a
negative control is for.

### 4.3. Variation in metadata predictability across cohorts

Ordered by the metadata-only AUC, the nine comparisons run from 0.404 to 1.000 (Figure
3A, Table 1), with no gap that would justify a threshold. Positions are stable across
seeds: excluding the CMV control, whose near-chance value moves 0.114 between seeds
precisely because it is near chance, the widest five-seed spread is 0.042. Eight of the
nine have a significant association with the recorded metadata in every seed, and the
count is still eight after Benjamini-Hochberg correction within each seed, whether the
family is that seed's nine all-metadata tests or all 34 of its variable-group tests. The
CMV negative control is the ninth.

Two comparisons show why the metadata-only AUC rather than V_D carries the test
(Section 2.2). In COVID-19 Ren (26 metadata columns, 182 donors) and COVID-19 COMBAT
(7 columns, 110 donors) the cross-fitted V_D is unstable across seeds, running from
0.847 to its ceiling of 1.000 in Ren and 0.906 to 1.000 in COMBAT. When it reaches the
ceiling the matched null reaches it too, so p(V_D) is 1.000; in the other seeds the same
comparison gives p(V_D) = 0.005. The metadata-only AUC on the very same matrices is
stable by contrast, 0.812-0.854 and 0.862-0.884, significant in every seed (p = 0.005 and
0.005-0.010). Counting by p(V_D) alone would give six of nine rather than eight; we
report both counts and use the better-powered one.

Under the primary, collection-only definition the count is smaller, and we give it rather
than leave the union to stand for study design. Seven of the nine comparisons record any
collection variable at all; in four of those seven the collection-only AUC is significant
in all five seeds, and two of those four are the comparisons flagged as not estimable or
underpowered. Among comparisons that are both estimable and record collection variables,
two of four are significant: lupus in GSE174188 CD4 (0.762-0.774) and COVID-19 in Ren
(0.729-0.768). Two comparisons record no collection variable that varies - lupus in
GSE174188 all cells, and the single-assay Ren subset - so their stratified permutation
conditions on tertiles of log cells per donor instead. For those two the test is
sample-quality-stratified, and Table 1 marks them rather than reporting them as
collection-stratified. This difference is the point of separating the groups: much of
what the union detects is carried by sample quality and case mix, which are not
experimental design and which carry their own interpretations (Section 4.1).

Grouping the metadata variables shows that similar overall strength can come from
different parts of the collection process (Figure 3B). The clearest case is a pair of
COVID-19 cohorts of similar size and closely comparable metadata-only AUC (Figure 3C):

- **Stephenson** (n = 115, recorded-metadata AUC 0.845-0.882). Sample-quality variables
  alone reach 0.865-0.884 (p = 0.005 in every seed), while collection variables reach
  only 0.501-0.594 and are not significant in every seed (p = 0.015-0.234). The
  metadata-phenotype association is carried mainly by sample quality.
- **Ren** (n = 182, recorded-metadata AUC 0.812-0.854). The mirror image: collection
  variables reach 0.729-0.768 (p = 0.005 in every seed), while sample quality reaches
  only 0.590-0.632 and is not significant (p = 0.050-0.154). The association is carried
  mainly by recruiting hospital and sub-study.

Similar overall metadata-only AUCs therefore concealed different patterns of association.
The Stephenson result motivates closer examination of library-quality measures, and of
whether harmonised preprocessing changes the estimate; the Ren result motivates
examination of recruiting site and sub-study composition, where the relevant question is
whether donors exist under overlapping collection conditions. A single number would call
the two cohorts equally affected and point both at the same follow-up. This is also the
direct answer to the point that V_D depends on how the metadata matrix is coded: it does,
and the response is to report the grouping rather than to hunt for a coding-free number.

The Ren cohort also shows how an incomplete set of metadata variables misleads. With only
sequencing chemistry recorded as collection metadata, the collection-only AUC is 0.553 -
apparently untroubled by site. Recovering the recorded hospital (12 centres, five of them
case-only) and sub-study identifier (six source studies) raised that group
to 0.750 and the full set from 0.758 to 0.818. The problem was there and invisible.
These checks are only as good as the metadata they are given, and a reassuring result
computed on thin metadata is weak evidence.

### 4.4. Free and collection-stratified permutation results

**Where the stratified test is uninformative.** The COMBAT influenza comparison has 21
donors and metadata that separates the labels exactly (metadata-only AUC 1.000,
V_D 0.000). Every collection stratum is label-pure (Figure S4), so within-stratum
permutation leaves all labels unchanged and the reference distribution is degenerate; we
report the stratified test as uninformative, with p = NA. The free permutation is still
defined and rejects at the 1/1001 floor in all five seeds, behind an expression AUC of
0.982 to 1.000. Reported alone, that reads as a comfortably validated classifier. What
the two results together say is that recorded metadata and the phenotype cannot be
separated in this comparison with the recorded data, and that the test which would have
detected it is the one that has no reference distribution here. We report the comparison
outside the main ordering.

**Where the two tests point in different directions.** In the COMBAT
sepsis-versus-COVID-19 comparison (n = 111, minority class 11), metadata predicts the
phenotype almost perfectly, at 0.996 to 0.998 across seeds. The free permutation yields
p = 0.001 in every seed; the collection-stratified permutation yields p = 0.057 to 0.614,
significant in none (Figure 3D). The classifier therefore shows evidence of an overall
association between expression and the phenotype, but no evidence of discrimination
beyond what the specified collection strata retain. The two results are not contradictory
- the tests state different nulls - but a report carrying only the first would read as a
validated classifier, and that reading is not supported here.

Two cautions on this row. The smaller class holds 11 donors, so the stratified test has
little power and its failure to reject is not a demonstration that the discrimination is
entirely collection-driven. And the fixed-hyperparameter pipeline scores 0.011 to 0.066
AUC above the tuned one here, the largest gap among the eight comparisons that return a
number. That is not a feasibility failure - the fixed pipeline is the inferential
statistic, and it is what both permutation tests are built from - but with 11 donors in
the minority class a reader should be able to see how much of the headline AUC depends on
model selection.

### 4.5. What expression adds over the recorded metadata

Reporting the metadata-only AUC beside the expression AUC says whether each predicts the
phenotype. It does not say whether expression carries information the recorded metadata
does not already have, which is the question a reader of a patient-level classifier has.
We fitted the third model: inside each training fold, expression is reduced to principal
components on the training donors, the metadata matrix is built from the same donors, and
one classifier is fitted on the concatenation. All three arms run the identical
fixed-hyperparameter specification through identical folds, so a difference between them
is a difference of information (Table S5).

In seven of the nine comparisons the increment over metadata alone is large: adding
expression raises the AUC by 0.087 to 0.354 across five seeds, and in each of those the
joint model is within 0.023 of the expression-only model, so the metadata contributes
little once expression is present. Lupus in GSE174188 CD4 gains 0.087 to 0.101 over
metadata alone; the single-assay Ren subset, which records the least collection metadata,
gains the most at 0.227 to 0.262.

The two exceptions are the two comparisons already flagged. In sepsis-versus-COVID-19,
metadata alone reaches 0.996 to 0.998 and adding expression changes the AUC by −0.001 to
+0.001; expression adds nothing measurable to what the recorded metadata already
provides, while metadata adds 0.036 to 0.104 over expression alone. In influenza both
models are at 1.000 and the increment is exactly zero in every seed. The increment
therefore agrees with the permutation results, and agrees with them for a reason that
requires no permutation: these are the two comparisons where the recorded metadata
already accounts for the separation.

CMV runs the other way and is worth stating. Its metadata-only AUC is at chance, its
expression AUC is 0.894 to 0.933, and the joint model is *lower* than expression alone, by
0.069 to 0.136. Appending 89 uninformative one-hot columns to 50 principal components
costs discrimination rather than adding it - the same property that makes residualisation
on that matrix expensive (Sections 3.5, 5.3).

The increment is reported as a descriptive quantity over five seeds, with no p-value
attached. A significance test for it would have to hold the metadata association fixed
while permuting, which is the collection-stratified permutation already reported in
Table 1 on the statistic it is already reported for.

### 4.6. Scope of the ordering

The ordering says how well recorded metadata predicts the phenotype label. It does not
show that any particular published classifier depends on that metadata, because it does
not test whether the classifier uses the metadata-related part of X - the D → X → Ŷ path
in Figure 1. Sections 5.2 to 5.4 test that path directly in the two lupus cohorts, and
reach a more equivocal answer than the ordering alone would suggest.

## 5. Whether the classifier uses the metadata-related structure: two lupus cohorts

The ordering measures the D-Y arrow. This section follows the D → X → Ŷ path in the
two lupus cohorts for which complete donor metadata is held, using frozen
Geneformer [4] embeddings and two expression pseudobulk representations under one
fold-contained pipeline. GSE285773 [3] serves as the independent transfer source.
Figure 4A shows how cases and controls are distributed across the recorded design in
both cohorts; Figure 4B shows how well the recorded metadata alone predicts disease status.

### 5.1. All three representations encode batch

In GSE174188 the three representations classified disease at AUC 0.982, 0.977 and
0.975. The same representations identified processing wave, one against the rest, at
0.9997, 0.999 and 0.997. In GSE135779 disease AUCs were 0.870, 0.933 and 0.945, while
the best single-batch AUCs were 0.993, 0.961 and 0.961 (Figure 5A).

This shows that the representations contain batch information. It shows nothing more.
A representation that encodes batch can still carry disease information, and a
batch-aware representation need not produce a batch-driven decision rule. Batch
predictability is therefore not a validity check on its own, and we report it below only
with that scope attached.

### 5.2. Batch predictability did not imply sensitivity to batch adjustment

GSE135779 provides the counterexample directly. In that cohort, Geneformer embeddings
predicted the most distinguishable batch at AUC 0.993 and both pseudobulk representations
at 0.961, close to perfect batch recognition, whereas batch alone predicted disease
status at AUC 0.499. Batch residualisation then changed disease AUC by +0.002, −0.002 and
−0.009 (Figures 5B, 6C) - three shifts well inside the split-to-split spread, so the
classifier's dependence on batch is small rather than demonstrably zero. Strong batch
information in a representation therefore did not imply a substantial change in disease
discrimination under the batch adjustment evaluated here, and "the representation encodes
batch" does not license "the classifier uses batch" on these data.

### 5.3. Residualisation and restriction answer different questions

Restriction asks whether the separation survives among donors collected under
overlapping conditions. Residualisation deletes the part of the data that a chosen set
of metadata variables can predict. When metadata and phenotype are collinear these are not
the same question, and we compare them to show they are not interchangeable, not to
rank them.

In GSE174188, ridge residualisation on processing-wave fractions reduced the three
representations from 0.982, 0.977 and 0.975 to 0.714, 0.747 and 0.735. Restriction to the
one large near-balanced wave, compared against size- and label-matched random donor
subsets - the comparison that separates the cost of conditioning from the cost of using
fewer donors - cost only 0.014 to 0.049 (Figure 6B). The paired difference between the two
losses was 0.186 to 0.255 across six comparisons, with every descriptive interval above
zero.

In GSE135779, removing the single all-case batch left 36 donors with both labels in
every remaining batch and changed matched performance by at most 0.012. Five ways of
removing the recorded variables give five different answers for the same cohort
(Figure 6A). Projecting out the complete set with ridge collapsed all three
representations to 0.516, 0.522 and 0.520. A random-forest residualiser removed less
(down to 0.748, 0.787 and 0.795), and overlap weighting and batch location adjustment
less still (retaining 0.859 to 0.937). Every full-set removal took away far more than
the batch-only control.

The algebra explains this without needing a technical story. When design predicts the
diagnosis, the part of the data that design can predict necessarily includes
diagnosis-related variation. Removing it removes disease signal regardless of where
that signal came from. A large loss after residualisation is therefore **not** evidence
that the original AUC was spurious. It is equally consistent with removing a shortcut
and with removing real biology, and the size of the loss does not distinguish them.

The CMV comparison settles part of what residualisation loss actually measures. There, no
association is detectable between the full recorded metadata and the phenotype:
metadata-only AUC 0.512 under this pipeline, permutation p 0.383 to 0.881 across five
seeds, and p(V_D) = 1.000. Residualising its representations on those same variables
still reduces classification AUC from 0.925 to 0.599, a loss of 0.326 (Figure 6C). The
discrepancy is worth stating plainly: an adjustment made for a metadata-phenotype
association that no test detects removes a third of the separation anyway.

We tested the obvious explanation rather than asserting it. The matrix declares 89
one-hot columns for 108 donors, at a centred rank of 54, and the collection block alone
declares 80 columns at rank 45; a projection of that reach takes out most of any sample.
The simulation in Section 3.5 generates a metadata matrix independent of the phenotype by
construction and residualises with the identical procedure. A balanced matrix of the same
width (92 columns) removes 0.159, and the loss is not monotone in width, peaking at 0.333
for 48 columns. A matrix built to the ragged shape of the real batch-by-pool crossing -
many levels holding one or two donors - removes **0.363 on average, central 95% 0.272 to
0.464**, which contains the 0.326 observed here (Figure S2).

So the observed CMV loss is the size that this adjustment procedure produces on a matrix
of this shape when there is provably nothing to remove. The conclusion is about
residualisation, not about the cohort: the size of a residualisation loss depends on
properties of the adjustment matrix that have nothing to do with confounding, and cannot
be read as a measure of technical contamination.

It is tempting to read the pair of values from residualisation and restriction as a
range bounding the true biological effect. It is not one. There is no proof that the two
estimates bracket any underlying quantity. They are two estimates made under two
different sets of assumptions, and their disagreement is informative precisely because
neither is identified.

### 5.4. Cell-type composition reproduces the same pattern

Aggregating the same GSE174188 CD4 cells into naive, effector-memory, regulatory and
unassigned fractions gave disease AUC 0.899 to 0.915. The same four-number summaries
predicted processing wave at 0.803 to 0.872. Wave residualisation reduced them to 0.689
to 0.709 (Figure S5A), and within-wave restriction against matched subsets cost 0.027
to 0.062, with every interval crossing zero (Figure S5B). A four-dimensional biological summary therefore reproduces
the whole pattern. That places the problem in the cohort rather than in
foundation-model embeddings, and is consistent with composition being a strong patient
stratification baseline in its own right [26].

### 5.5. Source restriction and cross-cohort discrimination

Restricting the training set to one processing wave is the obvious remedy: train where
the design is uniform, then test elsewhere. It does not work. Training on all 261
GSE174188 donors and testing in the independent 26-donor GSE285773 cohort gave 0.900,
0.944 and 0.969. Restricting the source to the dominant wave lowered every one of them,
and did no better than a random subset of the same size (Figure 7). Within GSE174188,
training on waves 2 and 3 and testing on wave 4 gave 0.842, 0.811 and 0.806, and the
reverse direction 0.919, 0.926 and 0.821. More varied source donors transferred better
than the apparently cleaner restricted source.

Transfer in the other direction, from the 26-donor cohort to the 261-donor target under
strict source-only rules, gave 0.884 [0.843-0.920], 0.919 [0.884-0.946] and
0.926 [0.895-0.953]. Both pseudobulk representations beat Geneformer after paired DeLong
testing [23] and Benjamini-Hochberg adjustment [24] (q = 0.0151, q = 0.000603).

Whether a more expressive donor-level aggregator changes this is a separate question
from the one this paper asks, and we report it in the supplement rather than the main
text. Briefly: three learned pooling methods on the same frozen embeddings never beat
ordinary mean pooling, and every significant paired difference was negative
(Figure S6, Table S6). Adding capacity fitted the source without producing a
representation that travels.

## 6. Integrating metadata assessment and predictive validation

These are not a validity standard. We have not shown that they are sufficient, that they
are minimal, or that passing them licenses a causal claim. They are five checks that look at different arrows in Figure 1
and that fail in different ways. Their value is that they can disagree with each other.

| # | Check | Shows | Does **not** show |
|---|---|---|---|
| 1 | Metadata-only AUC by variable group, cross-fitted V_D, p(V_D) against the cohort's own null, overlap counts | Whether recorded metadata predicts the phenotype label, and which variable groups do | That the classifier uses that metadata; anything about unrecorded structure |
| 2 | Complete-pipeline permutation, both free **and** collection-stratified | Free: evidence of an expression-phenotype association of any origin. Collection-stratified: evidence of discrimination beyond what the collection strata retain | Neither identifies a mechanism, and neither certifies the absence of leakage — a label-dependent step taken before the permutation is present on both sides. The stratified test conditions on collection variables only, so it does not cover sample-quality or demographic association inside a stratum (Section 3.4), and it is uninformative when every stratum is label-pure, where the free test is still defined |
| 3 | How well the representation predicts batch | That the representation encodes collection information | That the phenotype decision rule depends on it — Section 5.2 is a counterexample |
| 4 | Residualisation **paired with** restriction and matched donor controls | How far two non-equivalent adjustments disagree | Which adjustment is right; bounds on a biological effect |
| 5 | Source-only external target | Whether patient ranking survives a change of collection process | Calibration or clinical usefulness; and it carries much less weight when the target was collected the same way, since a shared bias reproduces rather than cancels (Section 3.3) |

Checks 1 and 2 cost minutes of CPU on a donor-level table and are the ones we would run
first on any new cohort. Check 5 is the only one that speaks to transport, and no
combination of the internal checks substitutes for it.

**What a clinician should take from this.** An internal cross-validated AUC measures
discrimination under the conditions the cohort was sampled in, which is a real quantity
and the right starting point. Two further things extend it. A metadata-only model says
whether the recorded collection differences alone could produce that discrimination, and
a conditioned validation - the collection-stratified permutation, or restriction to
overlapping strata - says whether it exceeds them. An independent-cohort evaluation then
addresses transfer to a new setting, and carries much less weight when the second cohort
was assembled the same way as the first, because a shared bias reproduces rather than
cancels. A classifier reported with all three is worth taking seriously as a candidate
diagnostic; one reported with the internal AUC alone has not yet been separated from its
cohort.

Applied to the nine comparisons, the checks give three groups rather than a pass and a
fail. In the first, no association with the recorded metadata is detectable and internal
discrimination is not undermined by anything visible here - which is weaker than saying it
is biological, since unrecorded structure remains possible; only CMV is here. In the
second, the metadata predicts the phenotype strongly and yet the collection-stratified
permutation is still rejected, so the discrimination exceeds what the collection strata
retain; both lupus comparisons and all three COVID-19 cohorts are here. In the third, the
stratified permutation is either not rejected or uninformative, and no internal analysis
settles the question; sepsis-versus-COVID-19 and influenza are here. Only the third group
is a stop condition, and it is a statement about the cohort, not about the method being
evaluated.

### 6.1. Three comparisons traced through the tree

Figure 8 traces four comparisons through the decision tree with the evidence for each
branch stated. Three are given in full here; the fourth, sepsis versus COVID-19, is
described in Section 4.4. The strata, the metadata matrix and the classifier are all
imported from the screen rather than rebuilt for the walkthrough, so Q2 and Q4 ask about
the same partition, the same columns and the same fitted pipeline.

Every branch is written to output what can be estimated and what is still missing, not a
recommended remedy. That distinction is what the previous version of this tree got wrong,
and Section 8 records the correction: a branch may say that no collection association was
detected, and it may say that a conditioned estimate is unavailable, but it cannot on that
basis instruct anyone to residualise. Whether to adjust is decided by the target quantity
and the scientific question, and the evidence here bears on it without settling it.

**COMBAT influenza (n = 21) stops at Q0.** Its stratum variable is the recording
institute, which gives two strata, *neither of which contains both a case and a control*
(Figure S4). Within-stratum permutation therefore leaves every label unchanged and the
stratified test is uninformative. That, and only that, is the stopping condition: the
comparison offers no donors under overlapping collection conditions, so no conditioned
estimate can be formed from these data and no route in the tree is available. Its
metadata-only AUC is also 1.000, which is consistent with the same structure but is not
what triggers the stop - a finite-sample AUC of 1.000 records that the fitted model
separates these donors completely, and is not on its own a proof of non-identifiability.
The fixed and tuned pipelines agree closely here, 0.982 to 1.000 against 0.964 to 1.000
across seeds, so nothing about the stop depends on which is read. Note what the free
permutation would have said alone: p = 0.001 in every seed, which reads as a validated
classifier and is the wrong description of this comparison.

**CMV HIHA (n = 108) reaches Route A.** The stratified permutation is defined - 46 strata
from batch and pool, 21 of them holding both labels - and both permutation tests reject at
the 1/1001 floor. At Q3 no association with the recorded metadata is detected:
metadata-only AUC 0.512 under the walkthrough's pipeline and 0.404 to 0.518 in the main
screen, not significant against its own permutation null in any seed. Route A therefore
outputs *no detected collection association; the unadjusted estimate stands as the
reported quantity*. It does not output a recommendation to residualise, and this
comparison is the reason. Residualising here costs 0.326 AUC (Section 5.3) for an
association no test detects, and the Section 3.5 simulation shows a matrix of this shape
costing 0.363 with the association absent by construction. A block that is not detected is
not thereby shown to be negligible - the test has whatever power 108 donors afford - but
neither is undetectability a reason to adjust. Whether to adjust is a decision about the
target quantity, and Route A reports the evidence rather than making it.

**COVID-19 Ren (n = 182) reaches Route C.** Both permutation tests are rejected and the
association with the recorded metadata is strong: metadata-only AUC 0.832 under the
walkthrough's fixed-penalty pipeline and 0.812 to 0.854 under the main one. Label-mixed
strata exist, so Q4 is reached. Q4 asks whether restriction is *feasible*, and this is the
question the earlier version of the tree left open: it asked only whether a mixed stratum
existed, which let the same comparison be assigned to a restricted estimate at Q4 and then
described as needing redesign a paragraph later. The requirement is now fixed in advance -
a comparison takes the restricted route only if some stratum holds at least five donors of
each class *and* at least 40 donors in total, the same power floor used everywhere else in
this paper - and each comparison therefore has one output.

Ren does not meet it. Exactly one stratum is large enough to cross-validate at all, and it
holds 18 donors; the restricted estimate there is AUC 0.992 against 0.733 in size- and
composition-matched random subsets, but the matched-subset standard deviation is 0.145,
which is the width of the whole effect being estimated. Residualising the 26-column matrix
costs 0.083, from 0.967 down to 0.884. The tree therefore outputs Route C: the association
is present, both adjustments are available, they disagree in direction, and neither is
estimable with useful precision on these donors. Both numbers are reported, no bound is
claimed, and the outstanding requirement is donors collected under overlapping conditions
rather than a further analysis of these ones.

**Restriction is rarely available at these strata, and the tree now says so.** CMV has 21
label-mixed strata and *none* with five donors of each class. Ren has one usable stratum
of 18 donors. Route B is drawn as though restriction were a practical alternative to
residualisation; on real cohorts with realistically fine strata it frequently is not. We
state that rather than coarsening the strata until Route B becomes available, which would
make the tree unfalsifiable.

# Discussion

## 7. Implications for patient-level single-cell evaluation

Across the nine comparisons, phenotype predictability from recorded metadata was
measurable, graded and traceable to specific variable groups. It ran from a metadata-only
AUC of 0.404 to 1.000, more than half the usable scale, and each cohort's position was
stable across split seeds. The *raw* residual label variance, by contrast, is not
comparable across cohorts without its own null: an in-sample R² ranks the least affected
comparison in this set as the third most affected, purely because its metadata matrix is
wide relative to its donor count.

Similar overall metadata-only AUCs did not imply similar underlying associations. In the
Stephenson COVID-19 cohort, sample-quality variables reached AUC 0.865 to 0.884 against
0.501 to 0.594 for collection variables; in Ren, collection variables reached 0.729 to
0.768 against 0.590 to 0.632 for sample quality. The two overall figures differ by about
0.02. The Stephenson result motivates closer examination of library-quality measures and
of whether harmonised preprocessing changes the estimate; the Ren result motivates
examination of recruiting site and sub-study composition, where restriction to balanced
strata is the relevant question and, as Section 6.1 shows, is not available at these
donor counts. Reporting one number hides the part of the answer a reader can act on.

The two permutation tests have the widest reach beyond this dataset. Free label
permutation is standard practice and answers what it is asked: whether expression and the
phenotype are associated at all. It cannot answer whether a classifier exploits the
collection process, because the shuffle destroys that process along with everything else.
Stratified permutation is not new - Chaibub Neto and colleagues introduced it for exactly
this purpose [30] - and what these cohorts add is what the two tests return when the
acquisition structure varies. In sepsis-versus-COVID-19 the free test rejects in every
seed and the stratified test in none; in influenza the stratified test has no reference
distribution at all; in the other seven they agree. A study using label permutation to
argue that a patient-level classifier is valid should report the stratified version
beside it, or state that no label-mixed stratum exists and that the test is uninformative.
That statement is itself the more informative outcome.

The stratified test has a limit of its own, and it is better published as a number than
as a caveat. It conditions on collection strata, not on the whole recorded matrix, so a
sample-quality variable that stays associated with the phenotype inside a stratum can
carry a classifier past it with no disease effect present; Section 3.4 measures how
often. Rejecting it therefore rules out one specific explanation - these strata - and not
the general one. Detecting a conditional association is not a certificate of disease
biology.

The simulation adds a caution about external validation that we did not anticipate.
When the external cohort shares the source's collection mechanism, external AUC rises
with confounding rather than falling. An external result is evidence about transport
only in proportion to how independent the second collection process really is, and
being a different accession number does not establish independence.

For representation research the implication is narrower than a leaderboard result but
harder to work around. In our lupus data pseudobulk beats frozen embeddings in strict
source-only transfer, and added pooling capacity does not change that. The more durable
point is that internal accuracy in these cohorts is not a measurement of the property
that matters, and no amount of model development fixes a cohort in which design predicts
the diagnosis. That is a recruitment and data-release problem. The fix is to record and
distribute collection metadata, and to build cohorts in which cases and controls appear
together within batches - not to build a better encoder.

None of these results supports a diagnostic claim. The external AUCs reported here rank
261 or 26 donors and are not calibrated probabilities. What the analysis supports is a
reporting discipline. A patient-level classifier should be published with three things.
First, the metadata-only AUC of its development cohort **by variable group**, so a reader
can see whether the association runs through collection, case mix or sample quality.
Second, both permutation tests, or a plain statement that the stratified one is
uninformative because every stratum is label-pure. Third, an external result on a cohort
collected through a different route, fitted on the source cohort alone. None of that
requires a new method. It requires releasing the metadata and running two cheap checks.

## 8. Limitations

These analyses address association with recorded metadata rather than causal biological
specificity. Unmeasured collection structure - referral patterns, institutional treatment
protocols, sample handling not captured in released metadata - is invisible to every
quantity reported here, so a low metadata-only AUC is weak evidence of an unaffected
cohort rather than strong evidence of one. What that costs is visible within the Ren
cohort itself: a metadata set that omits the recruiting hospital gives a collection-only
AUC of 0.553 where the completed set gives 0.750, on the same donors.

**The collection-stratified permutation conditions on a subset of the recorded
metadata.** Its strata come from collection variables; sample-quality and demographic
columns are in the metadata matrix but not in the strata. Section 3.4 measures the
consequence directly: with a sample-quality variable associated with the phenotype inside
a stratum and driving expression, the test rejects at rates rising to 95% with no disease
effect anywhere in the generating model. Rejecting it is therefore evidence of
discrimination beyond the collection strata, not evidence of a disease-specific signal,
and not a demonstration that some part of the discrimination is free of association with
the recorded metadata. Building a null
that conditions on the full matrix, including its continuous columns, would require a
conditional randomisation scheme we do not develop here. That is the single largest
methodological gap in this paper.

**The three variable groups are not equivalent, and only one is study design in the
strict sense.** Demographics can be case mix or genuine risk; sample-quality summaries
sit downstream of the assay but can also sit downstream of the disease. We therefore
report the collection block separately everywhere and treat it as primary for any claim
worded as being about study design. The headline count differs between the two
definitions, and both are given.

V_D is linear by construction and depends on how the metadata matrix is coded. Reporting
nonlinear metadata-only classifiers, the variable grouping and the matrix's own null
reduces but does not remove this. The nonlinear arms use two tree ensembles and are not
an exhaustive comparison of learners. Column count is not effective dimension either: the
CMV collection block declares 80 one-hot columns at a centred rank of 45, and no reading
of these results should treat a column count as a number of independent directions.

**Residualisation results are conditional on one adjustment procedure.** Section 3.5
generates a metadata matrix independent of the phenotype and still measures a loss of
0.363 at the CMV matrix's shape, which contains the 0.326 observed. That establishes that
a loss of this size carries no information about contamination; it does not establish
what a different residualiser would do. Ridge at a fixed penalty on a training-fold
standardised matrix is the only adjustment characterised at this depth, and the
random-forest, overlap-weighting and location-adjustment arms in Section 5.3 give
materially different answers on the same cohort.

**The decision tree organises the evidence; it is not a validated instrument.** Two of its
properties are worth stating because they are easy to assume otherwise. A null result at
Q3 reports that no collection association was detected and leaves the unadjusted estimate
standing; it does not recommend an adjustment, and CMV is the case that shows why one
should not be inferred from a null test. And Q4 carries an explicit feasibility
requirement, so a comparison with a label-mixed stratum too small to support an estimate
takes Route C rather than receiving one. Neither the branch structure nor the thresholds
have been validated against an external criterion, and no result in this paper depends on
the tree.

All nine comparisons were run at five seeds, and every value we quote is the range over
those five. Three comparisons have fewer than 15 donors in the minority class. We flag rather than exclude them, but no conclusion should rest on
them.

The nine comparisons come from five datasets and are **not independent**. Comparisons
drawn from the same dataset share donors, so the count of eight in nine is a count of
comparisons, not of independent replications.

Matched random subsets separate the cost of restriction from the cost of using fewer
donors without identifying a causal technical effect, so the restriction results show
non-attributability rather than that retained signal is biological. Coverage is four
diseases from five cohorts, all peripheral blood, all public, all curated under one
schema; tissue-resident and non-blood settings are untested, so the cross-disease claim
here is not a cross-tissue one. Geneformer pretraining overlap with individual benchmark
cells was not tested.

Finally, these are checks. They order cohorts and point at variable groups. They do not
identify a causal effect, and no claim here depends on their doing so.

# Methods

## 9. Materials and methods

### 9.1. Cohorts and the unit of analysis

Seven public datasets were used, in two groups. Table S7 is the registry. It gives, for
each dataset, the identifier, the donor count reported by the source, and the donor, case
and control counts actually modelled - these differ, because donors below the minimum cell
count or carrying a label outside the contrast are dropped. Five entries were verified
against the CELLxGENE Discover API on 2026-09-07; the two used only in Section 5 come from
GEO and are marked as not API-verified.

**Five datasets enter the main screen**, all obtained through the CELLxGENE Discover
curation API: the CD4-positive alpha-beta T-cell subset of GSE174188 (261 donors; 162
cases, 99 controls) [2], the Ren COVID-19 atlas (196 donors) [39], the Stephenson
COVID-19 cohort (120 donors) [40], the COMBAT consortium (124 donors) [41] and the HIHA
cytomegalovirus cohort (108 donors) [42]. Nine case-control comparisons are formed from
these five, so the comparisons are not independent of one another; Table 1 gives the
donor overlap.

**Two further lupus datasets are used only in the downstream analyses of Section 5**,
not in the screen: GSE135779 (44 childhood donors; 33 cases, 11 controls) [1] and
GSE285773 CD4-positive T cells (26 donors; 16 cases, 10 controls) [3], the latter as the
independent transfer target.

Registry entries, dataset identifiers, donor counts and download sizes were verified
against the API on 2026-09-07 and are released in `cohorts/registry.yaml`.

The donor was the independent unit in every supervised analysis. Cells from one donor
never crossed a training/test partition.

### 9.2. Donor-level tables and the metadata variables

Each cohort was streamed from its h5ad in 100,000-cell chunks and reduced to one row per
donor. Expression was aggregated to log1p CPM pseudobulk. The layer used for aggregation
was chosen by sniffing for integer counts rather than assuming a layer name.

Design variables were drawn from schema-guaranteed observation fields (`donor_id`,
`disease`, `sex`, `development_stage`, `assay`, `tissue`) plus cohort-specific collection
variables detected by pattern from the distributed object. Detection patterns cover
batch, pool, run, lane, chip, site, centre, city, region, province, hospital, clinic,
study, sub-study, source, dataset and donor source. Ages given as strings
("25-year-old stage", "fifth decade stage", "90 year-old and over stage") were parsed to
years.

Covariates that are consequences of the disease were excluded by an explicit blocklist
matching severity, outcome, comorbidity, medication, sample timing, symptom, diagnosis,
stage, grade and WHO score. Supplementary Table S2 lists, for every cohort, which columns
entered the metadata matrix, in which group, and which were excluded and why.

Variables were assigned to three groups. **Sample quality**: cells per donor, mean UMI
per cell, mean genes per cell, aggregate mitochondrial percentage, and their logs.
**Demographic**: age in years, sex, ethnicity. **Collection**: every detected
`batch__*` variable, assay and suspension type.

### 9.3. Residual label variance

For centred phenotype-label vector *y* and recorded metadata matrix **D** including an intercept,
V_D = 1 − R²(y ~ D). Two versions are reported. The in-sample version uses ordinary least
squares on all donors and is reported for comparison. The
cross-fitted version replaces the fitted values with out-of-fold predictions from
five-fold `KFold` with shuffling, seeded by the run seed, and is the version used
everywhere a number is interpreted. Negative R² values were clipped at zero before
subtraction.

### 9.4. The metadata matrix's own permutation null

For each design matrix, the phenotype-label vector was permuted 200 times. Each permuted
label was pushed through the **same five folds** used for the observed value, so the null
is a distribution of the cross-fitted statistic, not of a different one. The reported
p-value is one-sided, (1 + #{null ≤ observed}) / (200 + 1), and answers: how often does a
random phenotype label leave at most as much unexplained variation as the observed one? The
null mean and its 2.5th percentile are reported alongside.

The metadata-only AUC is tested the same way, with 200 permutations of the phenotype label
through the identical fixed-hyperparameter pipeline (`p_design_auc` in Table S4); that is the primary
test, for the reason given in Section 2.2. The in-sample V_D is released with its own
in-sample null, as `V_D_insample` and `p_V_D_insample` in Table S4. The two are separate quantities and are never mixed:
p(V_D) always refers to the cross-fitted pair. This null is what makes either version
readable, because its location depends on the number of metadata columns relative to
donors.

### 9.5. The fitted pipelines

**Four analyses fit classifiers, and they are not one pipeline.** The multi-cohort screen,
the decision-tree walkthrough, the calibration simulation and the GSE174188/GSE135779 deep
dive each have their own settings, and Table S8 lists all of them. That table is generated
by `tools/build_hyperparameter_table.py` from the `PipelineSpec` and `ForestSpec` objects
in `scripts/cohorts/pipeline_core.py` that the runners themselves import, so no setting can
be described here as something other than what ran. The paragraphs below state what
differs; Table S8 is the authority on every value.

*Shared by the screen, the walkthrough and the calibration simulation.* One classifier
specification, `SCREEN_FROZEN`, differing only in repeat budget: balanced logistic
regression (liblinear, `max_iter` 5,000) at fixed inverse regularisation 1.0, stratified
five-fold outer splits, the 4,000 highest-variance features selected on the training
donors, and a PCA to 50 components (or fewer when donors are limiting) fitted on those
donors and applied to the held-out ones. The screen runs one repeat, the walkthrough
twenty. A tuned variant of the same specification, with the inverse regularisation
selected from 10⁻⁴ to 10⁴ in decade steps by three-fold resampling inside each outer
training set over five repeats, is reported beside the fixed one as a description; the
fixed one is the inferential statistic, because it is what both permutation nulls refit.

*Metadata matrices.* The metadata matrix is small and dense, so no feature filter and no
PCA are applied to it. Numeric variables are median-imputed from the training donors of
each fold and standardised there; categorical variables are one-hot encoded from the level
set observed in those training donors, dropping the first level, and a level seen for the
first time in a held-out donor is encoded as the reference level. An earlier version of
this screen built the whole matrix once on all donors before cross-fitting. No label is
involved either way, and the change moves no headline count: the metadata-only AUC moves by
at most 0.012 anywhere in the screen, and both the count of eight in nine and the count of
four in seven are unchanged. It is corrected because the description in this section was
otherwise not true of the code. The in-sample V_D is still defined on the whole-cohort
matrix, as an in-sample statistic must be, and `n_design_feature` still reports that
matrix's width.

*The lupus deep dive.* GSE174188 and GSE135779 predate the screen and run a wider budget on
a single dataset: balanced logistic regression (liblinear, `max_iter` 20,000) with the same
inverse-regularisation grid selected by five-fold resampling inside each outer training
set, and 20 repeated stratified five-fold donor splits with integer seeds
20260801-20260820. Out-of-fold predictions are pooled within each repeat; we report the
mean AUC and the 2.5th-97.5th percentile of repeat-level AUCs. That range is a descriptive
split-sensitivity range. It is not a confidence interval, and it does not treat correlated
folds or repeats as independent [31].

**Nonlinear metadata-only arms.** Two tree ensembles rather than a wider search, because
the question is whether a *low-dimensional* metadata-phenotype relationship is being missed
by a linear model, not which learner maximises AUC. The screen's forest uses 500 trees,
unlimited depth, minimum leaf size 2, balanced class weights; the deep dive's uses 256
trees, maximum depth 4, minimum leaf size 3. Histogram gradient boosting, in the deep dive
only, uses 150 iterations, learning rate 0.05, seven terminal leaves, minimum leaf size 5.
All are fitted inside every outer training fold with fixed hyperparameters and are not
selected against held-out donors: a tuned nonlinear arm would not be comparable with the
tuned linear one under the same permutation null.

### 9.6. The two permutation tests

Both permutation tests and the observed statistic they are compared against run one
identical fixed-hyperparameter pipeline, refitted from the raw donor-by-gene matrix on every
permutation: within each training fold, the 4,000 highest-variance genes are selected on
the training donors, a standardiser and a PCA to 50 components (or fewer when donors are
limiting) are fitted on those donors and applied to the held-out ones, and a balanced
logistic regression at fixed inverse regularisation 1.0 is fitted, with a single
cross-fitting repeat. Every step that touches expression is therefore inside the fold, so
the reported AUC is genuinely out-of-fold and the null refits the whole pipeline rather
than reusing one representation. Freezing the regularisation is still necessary: comparing
a tuned observed statistic against an untuned null biases p-values downward.

Every representation step is refitted inside the fold on every permutation because the
code arranges it, not because the permutation test would reveal otherwise: a label-dependent
step taken before the permutation appears identically in the observed and permuted runs.
The property is enforced in one place, `pipeline_core.fold_contained_auc`, which every
analysis in this paper calls, and it is checked by an automated gate that fails the build
if a `fit_transform` appears outside a training fold in the permutation path.

The **free** test permutes the phenotype vector without restriction; it is `p_free` in
the released tables. The **collection-stratified** test permutes it only within collection
strata, preserving each stratum's observed class counts; it is `p_collection_preserving`
there. A stratum is the interaction of every detected `batch__*` variable with assay,
where assay varies. Table 1 names the stratum variables for every comparison, and the
screen writes them into its own output (`strata_definition`) so the name and the
implementation cannot drift apart. The stratified scheme is the restricted permutation of
Chaibub Neto and colleagues [30]; it is applied here, not introduced here.

**The conditioning set is collection variables only.** Sample-quality and demographic
columns are in the metadata matrix but not in the strata, so this test conditions on a
subset of what is recorded. State the two nulls explicitly. The free test's null is that
the phenotype label is exchangeable across all donors; the stratified test's is that it is
exchangeable within each collection stratum, that is, conditional independence of
expression and phenotype given those strata. Section 3.4 shows what follows: when a
sample-quality variable is associated with the phenotype inside a stratum and drives
expression, the conditional-independence null is false even though the generating model
contains no disease effect, and the test rejects it at rates rising to 95%. Those
rejections are correct rejections of a false null, not type-I errors, and they are the
reason a significant stratified result does not certify disease biology. Conditioning on
the full matrix would need a conditional randomisation scheme for continuous columns,
which we do not attempt here.

Where a cohort records no collection variable at all, the fallback is tertiles of log
cells per donor. That is a sample-quality variable, so for those comparisons the test is a
sample-quality-stratified permutation and is labelled as such rather than being folded in
silently.

Where every stratum is label-pure, within-stratum permutation leaves all labels unchanged
and the permutation group contains only the observed labelling. The reference distribution
is degenerate, and the test is reported as uninformative with p = NA. That is a statement
about the permutation group, not a claim that no null distribution exists in any
mathematical sense; the free permutation remains defined and is still reported.

A permutation contributes to the collection-stratified null only if at least one stratum
contains more than one donor and both labels; if no permutation qualifies, the p-value is
undefined and reported as such rather than as a number. Main runs used 1,000
permutations, so the smallest value it can return is 1/1001 ≈ 0.0010. The V_D and
metadata-only AUC nulls use 200, so their floor is 1/201 ≈ 0.0050. Every seed used the same
counts; the two floors in Table 1 are these two tests, not two different runs.

### 9.7. Feasibility and power checks

A comparison is reported as **not answerable with these data** on one criterion only: no
collection stratum holds both classes, so the stratified permutation has no non-trivial
reference distribution and no conditioned estimate can be formed. A metadata-only AUC of
1.000 is recorded alongside as corroborating structure, but it is not a second trigger. A
finite-sample AUC of 1.000 states that the fitted model separated these donors
completely; with 21 donors and seven one-hot columns that can happen without the metadata
determining the label in the population, so it is not on its own a proof of
non-identifiability. In the one comparison where both hold, the label-pure strata are the
stated reason.

A comparison is reported as **underpowered**, and kept in the main ordering with a flag,
if it has fewer than 40 donors or fewer than 15 in the minority class. The two gates are
separated because they call for different responses: the first cannot be fixed by more
careful analysis, and the second can be fixed by more donors.

Restriction has its own feasibility requirement, fixed in advance and applied uniformly:
a comparison takes the restricted route only if some collection stratum holds at least
five donors of each class and at least 40 donors in total. Where no stratum qualifies the
restricted estimate is not computed, and the comparison is reported as needing donors
under overlapping conditions rather than being given an estimate from a stratum too small
to support one.

The difference between the tuned and the fixed-hyperparameter pipeline AUC is reported as a column
(`frozen_minus_tuned_auc`, Table S4) and is deliberately **not** used as a feasibility
rule. A gap between two different algorithms does not establish that a comparison is
unanswerable. Every inferential statement in this paper is made about the
fixed-hyperparameter pipeline, which is the statistic both permutation nulls are built
from; the tuned AUC is descriptive. "Fixed-hyperparameter" refers to the statistical
pipeline throughout, and "frozen" only to a foundation model whose weights are not
updated.

### 9.8. Seed stability

All nine comparisons were re-run at five seeds (20260907-20260911), varying the seed for
fold construction, permutation draws, and the cross-fitting split used for V_D and its own
permutation null. All five per-seed outputs are released rather than summarised (Figure S7, Table S9), and
every range quoted in Table 1 and in the text is the minimum and maximum over those five
runs. Five master seeds control fold assignment and permutation draws throughout; the
seed is read at call time in every estimator, so each seed reaches the V_D null as well as
the classifier.

### 9.9. Simulation

Cohorts of 200 donors and 100 features were generated. A latent variable *z* produced a
site indicator and a continuous quality proxy; the diagnosis was generated as
ρ·z + √(1−ρ²)·e and then either used continuously or dichotomised. Biological and
technical loading vectors overlapped partially (technical = 0.7·random + 0.3·biological,
renormalised). Seven regimes were crossed with ρ ∈ {0, 0.25, 0.5, 0.75, 0.9, 0.98, 1}
and both label types, at 200 replicates per cell, giving 19,600 simulated cohorts. The seven
regimes are:

1. **additive linear**, the base case;
2. **nonlinear design effect**: tanh plus a quadratic term in the quality proxy;
3. **heteroscedastic noise**: the standard deviation depends on the site;
4. **design-by-cell-type interaction**: two feature blocks whose mixing weight depends on
   the site, with the biological effect confined to one block;
5. **20% case rate, matched thresholds**: the collection boundary is moved to the same
   quantile as the diagnosis threshold;
6. **20% case rate, offset boundary**: the collection boundary stays at the median, so
   the threshold and the boundary are offset;
7. **different acquisition in the target**: the external cohort is acquired under a
   different design mechanism. The two 20% arms are counted separately
because they behave differently, and separating them is what showed the earlier
"rare diagnosis blunts the checks" result to be a threshold artefact (Section 3.3).

Critically, the external target reuses the **source's** biological loading vector.
Regenerating it makes every transfer chance-level by construction and tests nothing;
only the design mechanism differs between source and target in the shift arm.

A supplementary arm quantifies what admitting a disease-caused covariate does. A severity
variable was generated as 1.2·(y − ȳ) + noise and added to the metadata matrix. Everything
else was held fixed.

The extended simulation's cohorts are 200 donors by 100 features, so it applies no feature
filter and no PCA; its classifier settings are listed in Table S8 under
`extended_simulation` and are deliberately not derived from the screen's, because the
screen's 4,000-feature filter and PCA to 50 would be a different reduction on a
100-column matrix rather than the same step. Its unadjusted, residualised and restricted
arms differ from each other in exactly one operation, which is the property the
calibration study below was missing.

**Calibration arms.** Two further studies use `pipeline_core.SIMULATION`, which *is* the
screen's fixed-hyperparameter specification. Arm A places 200 donors in eight collection
sites, generates a sample-quality variable that predicts the phenotype within site at
strength γ ∈ {0, 0.3, 0.6, 1.0, 1.5, 2.0}, and lets expression depend on the site and on
that variable only - there is no phenotype-to-expression arrow. Both permutation tests
then run at 200 draws each, over 300 replicates per γ. Only γ = 0 measures calibration;
for γ > 0 the stratified test's conditional-exchangeability null is false by construction,
because the sample-quality variable lies outside the conditioning set, so the rejection
rate there is the probability of detecting a non-disease association and is reported under
that name.

Arm B generates 108 donors and 4,000 features with a real biological effect and a metadata
matrix assigned independently of the phenotype, at widths 1 to 92 and in one further layout
that reuses the observed level sizes of the CMV batch-by-pool crossing. Unadjusted and
residualised AUCs come from the same function with one argument added, so the two arms share
folds, feature filter, standardiser, PCA, classifier and scoring rule and differ only in
whether the training-fold ridge fit on the metadata matrix is subtracted first. An earlier
version fitted the unadjusted arm's PCA on all donors and the residualised arm's inside the
fold; the measured gap then contained the change of representation as well as the cost of
adjustment, and came out negative at small widths, which is the signature of that
asymmetry. Sixty replicates per cell. The sweep is shardable across machines - each cell
seeds itself from its own coordinates - and a sharded run reproduces an unsharded one.

### 9.10. Decision-tree walkthrough

Three cohorts were traced through the tree with all evidence recomputed, and every branch
value is released in Table S10. The collection strata, the metadata matrix and the
classifier are imported from the screen rather than reimplemented, so Q2 and Q4 refer to
the same partition, the same columns and the same fitted pipeline; the walkthrough is
`SCREEN_FROZEN` at twenty repeats. An earlier version standardised the whole-cohort
metadata matrix before handing it to a fold-contained ridge, and ridge shrinkage depends
on scale, so a step described as inside the fold was partly outside it; the matrix is now
built and standardised on the training donors of each split.

Restriction was evaluated in the largest stratum meeting the feasibility requirement in
Section 9.7 - at least five donors of each class and at least 40 in total. Where no
stratum qualifies, that is reported and the comparison takes Route C, rather than being
worked around by coarsening the strata or by computing an estimate from a stratum too
small to support one. Matched random subsets of identical size and case/control
composition were drawn 20 times from the full cohort for comparison, which is what
separates the cost of conditioning from the cost of using fewer donors.

### 9.11. Recovery of GSE135779 collection metadata

Batch, collection year, age, sex, race, ethnicity and clinical variables were read from
the original Supplementary Table 1b, and per-library sequencing metrics from
Supplementary Table 1c [1]. Study names were linked to analysis donor identifiers
through the GEO title/accession map; all 44 donors mapped uniquely. The restoration
script writes the donor table, batch-by-label, year-by-label and batch-by-year
cross-tabulations, and SHA256 hashes of every source file. For GSE174188, per-cell
`Processing_Cohort` was read from the distributed CELLxGENE object and reduced to donor
fractions over four waves. The dominant wave is the maximum-fraction wave; a pure-wave
donor has at least 99% of analysed cells in one wave.

### 9.12. Representations

**Frozen Geneformer.** Cell embeddings used Geneformer V2-316M at repository revision
`04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5` [4]; checkpoint, configuration, token
dictionary and gene-median dictionary hashes are in the methods manifest. Gene
identifiers were version-trimmed and mapped to the V2 vocabulary. Within each cell the
ranking value was `raw count / total count × 10,000 / Genecorpus-104M gene median`, and
genes were sorted descending. Sequences used V2 start/end token IDs 2 and 3, padding ID
0 and maximum length 4096, leaving at most 4094 ranked genes. The model output was the
final `last_hidden_state`; token position 0 was kept as the 1,152-dimensional cell
vector. The encoder was frozen. Donor representations are coordinate-wise means of
sampled cell vectors, sampled without replacement with integer seed 1, capped at 500
cells per donor for GSE135779 and 1,000 for GSE174188 and GSE285773. The embedding
environment used Python 3.10, PyTorch 2.5.1 and Transformers 4.46.3.

**Pseudobulk.** We use pseudobulk rather than a learned latent space as the donor-level
representation, so that the comparison between representations is not itself mediated by
a model fitted with batch covariates [5]. Stored donor-by-gene values are `log1p(CPM)`. Inside every outer
training fold, gene variance was computed on training donors, the 4,000
highest-variance genes were selected, and coordinates were standardised by
training-donor mean and standard deviation. HVG pseudobulk uses these values directly.
PCA pseudobulk fits up to 30 components on the training values and applies the fitted
map to held-out donors. For cross-cohort transfer, trimmed gene identifiers were
intersected before source-only HVG selection; feature means, variances, selected genes,
scaler, PCA, classifier and regularisation were all fitted on source donors and applied
unchanged to the target.

### 9.13. Fold-contained residualisation and its variants

Unadjusted and residualised models share outer folds, feature construction, scaling,
PCA and classifier selection. Inside each outer training fold, metadata variables were
encoded as in Section 9.5 and a multivariate ridge regression (`alpha` 1, intercept
fitted and **not** penalised) mapped design to every raw representation coordinate.
Training and held-out coordinates were replaced by observed minus predicted values,
after which the identical pipeline was fitted to training residuals and applied to
held-out residuals. The diagnosis was never given to the residualiser.

Three further arms were run on GSE135779. A multivariate random-forest residualiser
(96 trees, depth 4, minimum leaf 3) fitted only on outer-training donors, with
pseudobulk genes pre-selected by training-donor variance. A batch-location arm
subtracting training-batch mean offsets from training and held-out features; because it
applies no empirical-Bayes scale shrinkage we describe it as ComBat-style location
adjustment rather than ComBat. And an overlap-weighting arm fitting a training-only
logistic diagnosis propensity model from the full design, clipping propensities to
0.025-0.975 and training the representation classifier with overlap weights.

For the central GSE174188 comparison, residualisation loss was paired by split seed. We
bootstrapped the 20 paired repeat-level losses 5,000 times and subtracted the locked
matched-restriction discrepancy. These intervals quantify split sensitivity conditional
on the analysed donors; they are not population confidence intervals.

To probe confound leakage with a nonlinear downstream learner [28], 100 Gaussian
simulations used no biological effect, a technical effect of one, and metadata-phenotype
association 0.75. Random-forest prediction was compared on raw features, globally
residualised features and fold-contained residualised features.

### 9.14. Restriction and matched donor controls

GSE174188 was restricted to dominant wave 4 and separately to pure wave 4. GSE135779
excluded the all-case batch B1 and retained B2 to B6, each containing both labels. The
complete internal pipeline was re-run after restriction. For every observed stratum, 20
donor subsets were drawn without replacement from the full cohort with identical donor
count and case/control count, each using one prespecified split seed. The
observed-minus-matched difference asks whether the stratum is harder than expected from
reduced donor number and label composition alone.

### 9.15. Cell-type composition

The GSE174188 object was restricted to the same CD4-positive alpha-beta T-cell
population and 261 donors as the molecular analysis. Per-donor counts were formed for
`T4_naive`, `T4_em`, `T4_reg` and all other or unassigned labels. We evaluated raw
closed proportions and centred-log-ratio coordinates; for CLR, 0.5 was added to each
donor-category count before closure, the donor mean log abundance was subtracted and the
final coordinate omitted. Each representation was evaluated with and without log total
CD4-cell yield, under the same repeated donor pipeline, with processing wave predicted
one against the rest.

### 9.16. Transfer and learned pooling

Cross-batch analyses trained on waves 2 and 3 and evaluated wave 4, then reversed.
Cross-cohort transfer fitted every statistic on the source cohort; target labels were
accessed only after prediction. Target AUC intervals used 5,000 case/control-stratified
donor bootstrap resamples and the percentile method. Paired AUC comparisons used DeLong
tests on identical target donors [23] with Benjamini-Hochberg correction inside declared
method families [24].

DeepSets, gated-attention multiple-instance learning and pooling by multi-head attention
operated on cap-500 frozen cell embeddings [20-22]. Cell standardisation, a whitened
source PCA32 projection, network weights and hyperparameters were all source-fitted.
Hidden widths were 16 or 32, weight decay 0.01 or 0.1, and Adam used learning rate 0.02
for 150 epochs; five-fold source-donor cross-validation selected the configuration,
three initialisations were fitted on all source donors, and target probabilities were
averaged. An ordinary donor mean in the identical PCA32 space isolates the contribution
of learned pooling. The PCA32 projection used for hyperparameter selection was fitted
once on all source cells before source-donor folds and never used target cells or
labels. Source-internal pooling AUC is therefore not an unbiased generalisation
estimate, and only independent-target and paired target-donor comparisons support the
learned-pooling conclusion.

### 9.17. Software and reproducibility

The donor-level analyses used Python 3.11-3.13 with NumPy 2.3-2.5, pandas 2.3-3.0,
scikit-learn 1.9.0, PyArrow and joblib; exact per-run versions are recorded in each
output manifest. Figures were produced in R 4.6.0 with ggplot2 4.0.3 and patchwork
1.3.2. The screen, the extended simulation, the permutation calibration and the
decision-tree walkthrough were run on a 16-vCPU ARM container; all other analyses were
run locally. Input hashes, seeds, parameter grids, donor-level predictions, complete null
distributions and output manifests accompany the analysis.

Random seeds are derived with `hashlib.blake2b` from a recorded master seed and the cell
identity. Python's built-in `hash` is not used: it is salted per process for strings, so
it would not reproduce across interpreter sessions. The derived seed for every cell is
written into the released result tables. No script contains a machine-specific path: the
figure scripts locate the project root from their own file position or from
`RHEUMLENS_ROOT`, and the analysis scripts take input and output paths as arguments.

Everything downstream of the compute-heavy runs is rebuilt by one command,
`bash tools/build_all.sh`: the supplementary tables, Table 1, all fifteen figures, the
typeset manuscript, and two checks that fail with a non-zero exit status. The first
re-derives numbers quoted in the text from the result tables and requires each to appear
verbatim. The second checks the length of the abstract and author summary, that the abstract is
unstructured, that figures and tables are cited in the order they are declared, and that
each declared item is cited at least once. It also refuses any placeholder, retired term
or withdrawn phrase, and it looks inside the rendered figures and the released table
headers, where a check on the manuscript text alone cannot see them. The submission package is assembled by a further script from
the same outputs, so a figure number cannot differ between the manuscript and the files.

### 9.18. Use of generative AI

Generative AI coding assistants were used during this work: Anthropic Claude (models
Claude Opus 4.1 and Claude Opus 5) through the Claude Code command-line interface, and
OpenAI Codex through its command-line interface. They were used in the following places.

- **Code review and refactoring** of the analysis scripts in `scripts/cohorts/`, `sim/`
  and `walkthrough/`. Every statistical definition, gate and threshold in those files was
  specified by the authors; the model's contribution was implementation, review of
  implementation, and the identification of two defects that the authors then confirmed
  by inspection and by rerunning the affected analyses.
- **Figure code.** All figures were drawn by R scripts written with model assistance
  (`figures/src/*.R`). Every number plotted was read from a released result table, and
  the correspondence between figure and table was checked by the authors panel by panel.
- **Manuscript language and structure**: rewriting for length and readability,
  reorganisation of sections, and consistency checking of cross-references and reference
  numbering.

No text, number, figure or citation was accepted without author verification against the
underlying result tables. Cohort selection, the statistical definitions, every gate and
threshold, the decision of which analyses to report, and the interpretation of results
were made by the authors. No data of any kind were generated by a model: every value in
this paper traces to a released result table produced by the analysis scripts. The
authors take full responsibility for the content of this publication.


## 10. Data and code availability

All cohorts are public. The lupus datasets are at GEO accessions GSE135779, GSE174188
and GSE285773; the COVID-19, influenza and CMV cohorts were obtained through the
CELLxGENE Discover curation API and their dataset identifiers are in
`cohorts/registry.yaml`. The MIT-licensed repository is at
https://github.com/LightChainr/rheumlens with a permanent Zenodo record at
https://doi.org/10.5281/zenodo.20813922. The versioned research object accompanying
this manuscript contains the cohort registry with API-verified donor counts, the
donor-level interface tables, every result table underlying Figures 2 to 8, all five
per-seed outputs rather than a summary, the 19,600-cohort simulation output, analysis
scripts, locked environments, `REPRODUCE.md` and `SHA256SUMS`. Public raw single-cell
matrices remain at their original accessions.

## 11. Declarations

**Author contributions.** Conceptualization, H.Y. and D.L.; methodology, H.Y.; software,
H.Y.; validation, H.Y. and D.Y.; formal analysis, H.Y.; data curation, H.Y.;
writing - original draft preparation, H.Y.; writing - review and editing, D.Y. and D.L.;
visualization, H.Y.; supervision, D.L.; project administration, D.L. All authors have
read and agreed to the submitted version of the manuscript.

**Funding.** This research received no external funding.

**Competing interests.** The authors declare no competing interests.

**Ethics approval.** Not applicable. This study analysed only publicly available,
de-identified single-cell datasets released by their original investigators under the
consent and approvals obtained at each contributing site. No new human or animal data
were collected.

**Informed consent.** Not applicable, for the reason given above.

**Use of AI tools.** See Section 9.18, which states what was used, where, and how the
output was checked.

**Acknowledgments.** The authors thank the investigators of GSE135779, GSE174188 and
GSE285773, and the contributors to the COMBAT, Stephenson, Ren and CMV cohorts
distributed through CELLxGENE Discover, for releasing donor-level data that made this
re-analysis possible.

## Figure legends

**Figure 1. What is being measured.**
A recruitment process R places each donor in the recorded design D and is also
associated with the diagnosis Y; the measured cells X are affected by D (technical)
and by Y (biological); a classifier maps X to Ŷ. C denotes covariates caused by the
disease - severity, medication, symptom timing - which lie on Y → C → X and are
excluded from D by blocklist. U denotes collection structure that was never recorded
and which nothing in this paper can see. The dotted line marks the association these
checks quantify. They measure it; they do not tell you which way it runs.

**Figure 2. Simulation across seven data-generating regimes.**
19,600 simulated cohorts (200 donors, 100 features): seven regimes × seven levels of
metadata-phenotype association ρ × two label types × 200 replicates.
**(A)** Fraction of runs with p(V_D) ≤ 0.05, binary labels. At ρ = 0 this is 0.025 to
0.065 across regimes, bracketing the nominal 0.05.
**(B)** Unadjusted, external, restricted and residualised AUC against ρ in the additive
linear regime. Restriction is undefined at ρ = 1 because no stratum holds both labels.
**(C)** Restricted minus residualised AUC for all seven regimes. The gap opens in every
one. It is blunted in the `imbalanced_offset` regime, where the collection boundary and
the phenotype threshold sit at different quantiles, and not in `imbalanced`, where the
case rate is the same 20% but the two are matched (Section 3.3).
**(D)** Internal (solid) and external (dashed) AUC when the target shares the source's
collection mechanism versus when it does not. External AUC rises with confounding in
the shared case.

**Figure 3. Nine case-control comparisons.**
**(A)** Comparisons ordered by cross-validated metadata-only AUC (filled, coloured by
disease) with the corresponding diagnosis-classifier AUC (open) on the same axis. Bars
give the range over five random seeds. Dashed line is chance.
Italic labels mark the two comparisons that the feasibility check removes from or flags
within the ordering.
**(B)** Metadata-only AUC by variable group, rows aligned to (A). A dark outline marks
p ≤ 0.05 for the metadata-only AUC against that group's own permutation null; "n/a" marks a group absent for
that comparison.
**(C)** The two COVID-19 cohorts with almost equal overall metadata-only AUC are affected
through different variable groups: Stephenson through sample quality, Ren through
recruiting hospital. Solid points are significant.
**(D)** Free and collection-stratified permutation p-values, one point per seed, for all
nine comparisons. Dashed line is 0.05. Only the sepsis-versus-COVID-19
comparison separates the two tests.

**Figure 4. Design and diagnosis in the two lupus cohorts.**
**(A)** Donors per collection stratum, split by case and control, for processing wave in
GSE174188 and for sequencing batch and collection year in GSE135779. A triangle marks a
stratum holding only one kind of donor; such strata cannot be permuted within, and
cannot be restricted to.
**(B)** Cross-validated AUC for predicting the diagnosis from recorded metadata alone, by
variable group, with 2.5th-97.5th percentiles over 20 repeated donor splits. GSE174188
has no recorded collection year, so its line is broken there. In GSE135779 batch alone
sits at chance (0.499) while the other groups do not, so the association runs through
collection year, sample quality and demographics rather than through batch.

**Figure 5. Containing batch information is not the same as using it.**
**(A)** For each representation, the cross-validated AUC for the diagnosis (filled) and
for the most predictable single recorded batch, one against the rest (open), in both
cohorts. Every representation recognises a batch about as well as it recognises the
diagnosis.
**(B)** GSE135779 negative control. The representation recovers the batch at AUC 0.993,
but the batch predicts the diagnosis at 0.499.
**(C)** Change in disease AUC after fold-contained batch residualisation in the same
cohort: +0.002, −0.002 and −0.009, all inside the split-to-split spread. Batch is
predictable from the representation, unpredictive of disease, and its removal changes
little, so batch predictability on its own supports no conclusion about what the
classifier uses.

**Figure 6. Residualisation and restriction give different answers, and neither is
identified.**
**(A)** GSE135779 diagnosis AUC under five ways of removing the recorded metadata, for
three representations, with 2.5th-97.5th percentiles over 20 repeated splits. The same
cohort scores 0.94 or 0.52 depending only on which adjustment is chosen.
**(B)** GSE174188 CD4, restricted to the one large near-balanced processing wave
(n = 66), against size-matched random donor subsets. The separation largely survives
inside a single wave.
**(C)** Negative control. In the CMV comparison, where the metadata-only AUC is 0.512 and
no metadata-phenotype association is detectable in any seed, the same ridge residualisation
still costs 0.326 AUC. A simulated matrix of the same ragged shape, with the phenotype
independent of the metadata by construction, costs 0.363 on average (central 95%
0.272-0.464; Figure S2). Residualisation loss is a function of the shape of the adjustment
matrix and cannot be read as a measure of contamination.

**Figure 7. Restricting the training set to one processing wave does not improve
transfer.**
AUC in GSE285773 (26 donors) for classifiers trained on GSE174188, by source set: all
261 donors, a random subset of 89, and the dominant processing wave (n = 89). The
wave-restricted source is no better than a random subset of the same size and worse than
using every donor. Bars are DeLong intervals, except for the random subset, where they
are 2.5th-97.5th percentiles over 20 draws.


**Figure 8. A decision tree for applying the checks.**
The feasibility question comes first (Q0): a comparison is set aside when no collection
stratum holds both labels, so that no conditioned estimate can be formed. The permutation
step is then split into two questions that state different nulls - the free test (Q1) and
the collection-stratified test (Q2). Q3 reads the metadata-only AUC against its own
permutation null rather than against a fixed cut-point; a null result there reports that
no collection association was detected and leaves the unadjusted estimate standing, and
does not recommend an adjustment. Q4 requires restriction to be feasible - some stratum
with at least five donors of each class and at least 40 in total - so each comparison has
one output; where it is not, the comparison takes Route C. Every terminal box states what
can be estimated and what is still missing, not a remedy to apply. Coloured paths trace
four real comparisons. The strata are the ones the collection-stratified permutation uses,
not a separate definition.

## Supporting information

**Figure S1. Rejection rate of both permutation tests under a non-disease within-stratum
association.**
**(A)** Fraction of runs rejecting at p ≤ 0.05 for both permutation tests, as a function
of how strongly a sample-quality variable predicts the phenotype *within* a collection
stratum. The generating model contains no disease effect at all. Only the leftmost point
(γ = 0) measures calibration; to its right the conditional-exchangeability null is false
by construction, because the sample-quality variable lies outside the conditioning set,
so the rate is the probability of detecting a non-disease association rather than a
type-I error rate. **(B)** The mean absolute within-stratum correlation that conditioning
on the collection strata leaves in place.

**Figure S2. Residualisation loss with no metadata-phenotype association present.**
AUC lost to fold-contained ridge residualisation in simulated cohorts where the metadata
matrix is independent of the phenotype by construction, against the width of that matrix,
at 108 donors. The unadjusted and residualised arms are the same call with one argument
added, so they share folds, feature filter, standardiser, PCA and classifier and differ
only in the adjustment. The biological effect is calibrated so the unadjusted AUC lands
near 0.90. The ribbon is the 2.5th-97.5th percentile over 60 draws. The triangle reuses
the level sizes of the real CMV batch-by-pool crossing, which is far more ragged than
equal blocks; it loses 0.363 on average. The dashed line marks the 0.326 observed in that
comparison, which falls inside the triangle's interval.

**Figure S3. Admitting one disease-caused covariate manufactures the finding.**
Metadata-only AUC, V_D, the fraction of runs flagged significant, and residualisation
loss, with and without a severity variable added to the metadata matrix, across five
levels of true metadata-phenotype association. At ρ = 0, where metadata and phenotype are
unrelated by construction, admitting the covariate raises metadata-only AUC from 0.488 to
0.793 and
flags the cohort in 100% of runs against 5.5% without it.

**Figure S4. Design composition of each cohort.**
Case and control counts by collection stratum for all nine comparisons, showing which
strata are label-pure and how many donors sit in mixed strata.

**Figure S5. Cell-type composition reproduces the whole pattern.**
**(A)** Diagnosis AUC in GSE174188 from four cell-type fractions alone, under four
encodings, before and after removing the processing wave. Removing the wave costs about
0.20 AUC in every encoding.
**(B)** The same encodings restricted to the dominant wave and to the pure wave, each
against a size- and label-matched random donor subset, 20 repeated splits per box.
Restriction costs 0.03 to 0.04 AUC on average - far less than removal costs in (A). A
four-number biological summary shows the same behaviour as a 1152-dimensional
foundation-model embedding, which places the problem in the cohort rather than in the
embedding.

**Figure S6. Learned pooling never beats averaging the cells.**
Paired differences in AUC between three learned pooling methods (DeepSets,
gated-attention multiple-instance learning, pooling by multi-head attention) and four
baselines, in both transfer directions, with paired DeLong intervals. Filled points are
significant after Benjamini-Hochberg adjustment; every significant difference is
negative.

**Figure S7. Seed stability.**
Per-seed values of metadata-only AUC and its permutation p, V_D and p(V_D), diagnosis AUC
and both permutation p-values, for all nine comparisons at five seeds. The p-value panels
use a log axis; 0.005 is the smallest value 200 permutations can return.

**Table S1.** Calibration study: summary output of both arms, the permutation rejection-rate
error of Section 3.4 and the residualisation-width study of Section 3.5.
**Table S2.** Design manifest, generated directly from the model-matrix construction
code: for every comparison and every variable group, the raw columns consumed and the
number of expanded columns produced. **Table S2b** lists every column present in a
covariate file that enters no block, with the reason. **Table S2c** is a machine check that the
manifest reproduces the width of every fitted model matrix (34 of 34 blocks agree).
**Table S3.** Simulation summary: all seven regimes × seven ρ × two label types, plus the
disease-caused-covariate arm.
**Table S4.** Full screen output: every comparison × every variable group, with
metadata-only AUC (linear and random forest), V_D in-sample and cross-fitted, null mean,
p(V_D), diagnosis AUC, both permutation p-values and stratum counts.
**Table S5.** Metadata-only, expression-only and joint AUC for every comparison at each
of five split seeds, with the two increments. All three arms use the same
fixed-hyperparameter specification and the same folds; expression is reduced to principal
components fitted on the training donors before the metadata columns are appended.
**Table S6.** Learned pooling: per-method target metrics and all paired DeLong tests,
both transfer directions.
**Table S7.** Cohort registry, generated from the registry file and the donor tables:
for all seven datasets, the identifier, citation, the donor count reported by the source,
the download size, the case and control labels, the donor, case and control counts
actually modelled in each comparison, and how and when each entry was verified.
**Table S8.** Every fitted pipeline in the paper and its settings, generated from the
specification objects in `pipeline_core.py` that the analysis scripts import. The four
analyses are separate pipelines and the table keeps them separate.
**Table S9.** Seed-stability table, all five seeds, all nine comparisons.
**Table S10.** Decision-tree walkthrough: all evidence for every branch, three cohorts.

## References


1. Nehar-Belaid D, Hong S, Marches R, et al. Mapping systemic lupus
   erythematosus heterogeneity at the single-cell level. *Nature Immunology*.
   2020;21:1094-1106. doi:10.1038/s41590-020-0743-0.
2. Perez RK, Gordon MG, Subramaniam M, et al. Single-cell RNA-seq reveals
   cell type-specific molecular and genetic associations to lupus. *Science*.
   2022;376:eabf1970. doi:10.1126/science.abf1970.
3. NCBI Gene Expression Omnibus. GSE285773: Single cell RNA profiling of
   blood CD4+ T cells identifies distinct helper and dysfunctional regulatory
   clusters in children with SLE. Public 2025-09-09.
4. Theodoris CV, Xiao L, Chopra A, et al. Transfer learning enables predictions
   in network biology. *Nature*. 2023;618:616-624.
   doi:10.1038/s41586-023-06139-9.
5. Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for
   single-cell transcriptomics. *Nature Methods*. 2018;15:1053-1058.
   doi:10.1038/s41592-018-0229-2.
6. He B, Thomson M, Subramaniam M, Perez R, Ye CJ, Zou J. CloudPred:
   predicting patient phenotypes from single-cell RNA-seq. *Pacific Symposium
   on Biocomputing*. 2022;27:337-348.
7. Martorell-Marugan J, Lopez-Dominguez R, Villatoro-Garcia JA, et al.
   Explainable deep neural networks for predicting sample phenotypes from
   single-cell transcriptomics. *Briefings in Bioinformatics*.
   2025;26:bbae673. doi:10.1093/bib/bbae673.
8. Liu T, De Brouwer E, Verma A, et al. Learning multi-cellular
   representations of single-cell transcriptomics data enables characterization
   of patient-level disease states. *Cell Systems*. 2026;17:101570.
   doi:10.1016/j.cels.2026.101570.
9. Wu Q, Ding J, He R, et al. Exploring phenotype-related single-cells through
   attention-enhanced representation learning. *Genome Medicine*. 2026;18:21.
   doi:10.1186/s13073-026-01598-x.
10. Perez M, Hong J, Zweig A, Azizi E. Domain-invariant feature learning for
    patient-level phenotype prediction from single-cell data. *bioRxiv*.
    2025. doi:10.1101/2025.09.22.677881.
11. Clemente-Larramendi I, Hillion S, Cornec D, Jamin C, Foulquier N.
    FloREN: Decoding immune regulatory networks through interpretable graph
    transformer patient representations. *bioRxiv*. Posted 2026-07-20.
    doi:10.64898/2026.07.12.738088.
12. Muhire J. Donor-aware scRNA-seq benchmarks for IBD classification.
    *arXiv*. 2026;2605.03281.
13. Kedzierska KZ, Crawford L, Amini AP, Lu AX. Zero-shot evaluation reveals
    limitations of single-cell foundation models. *Genome Biology*.
    2025;26:101. doi:10.1186/s13059-025-03574-x.
14. Wu J, Ye Q, Wang Y, et al. Biology-driven insights into the power of
    single-cell foundation models. *Genome Biology*. 2025;26:334.
    doi:10.1186/s13059-025-03781-6.
15. Ali S. Auditing pretraining contamination in single-cell foundation model
    benchmarks. *arXiv*. 2026;2607.20572.
16. Soneson C, Gerster S, Delorenzi M. Batch effect confounding leads to strong
    bias in performance estimates obtained by cross-validation. *PLOS ONE*.
    2014;9:e100335. doi:10.1371/journal.pone.0100335.
17. Song F, Chan GMA, Wei Y. Flexible experimental designs for valid
    single-cell RNA-sequencing experiments allowing batch effects correction.
    *Nature Communications*. 2020;11. doi:10.1038/s41467-020-16905-2.
18. Risso D, Perraudeau F, Gribkova S, Dudoit S, Vert JP. A general and
    flexible method for signal extraction from single-cell RNA-seq data.
    *Nature Communications*. 2018;9:284.
    doi:10.1038/s41467-017-02554-5.
19. Chaibub Neto E. Causality-aware counterfactual confounding adjustment as an
    alternative to linear residualization in anticausal prediction tasks based on
    linear learners. *Proceedings of Machine Learning Research*.
    2021;139:8034-8044.
20. Zaheer M, Kottur S, Ravanbakhsh S, et al. Deep Sets. *Advances in Neural
    Information Processing Systems*. 2017;30.
21. Ilse M, Tomczak JM, Welling M. Attention-based deep multiple instance
    learning. *Proceedings of Machine Learning Research*. 2018;80:2127-2136.
22. Lee J, Lee Y, Kim J, Kosiorek A, Choi S, Teh YW. Set Transformer: a
    framework for attention-based permutation-invariant neural networks.
    *Proceedings of Machine Learning Research*. 2019;97:3744-3753.
23. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or
    more correlated receiver operating characteristic curves: a nonparametric
    approach. *Biometrics*. 1988;44:837-845. doi:10.2307/2531595.
24. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical
    and powerful approach to multiple testing. *Journal of the Royal
    Statistical Society Series B*. 1995;57:289-300.
    doi:10.1111/j.2517-6161.1995.tb02031.x.
25. Ferrari E, Retico A, Bacciu D. Measuring the effects of confounders in
    medical supervised classification problems: the Confounding Index (CI).
    *Artificial Intelligence in Medicine*. 2020;103:101804.
    doi:10.1016/j.artmed.2020.101804.
26. Halter C, Andreatta M, Carmona SJ. Cell type composition drives patient
    stratification in single-cell RNA-seq cohorts. *bioRxiv*. 2026.
    doi:10.64898/2026.03.27.714811.
27. Chyzhyk D, Varoquaux G, Milham M, Thirion B. How to remove or control
    confounds in predictive models, with applications to brain biomarkers.
    *GigaScience*. 2022;11:giac014. doi:10.1093/gigascience/giac014.
28. Hamdan S, Love BC, von Polier GG, et al. Confound-leakage: confound
    removal in machine learning leads to leakage. *GigaScience*.
    2023;12:giad071. doi:10.1093/gigascience/giad071.
29. Spisak T. Statistical quantification of confounding bias in machine
    learning models. *GigaScience*. 2022;11:giac082.
    doi:10.1093/gigascience/giac082.
30. Chaibub Neto E, Pratap A, Perumal TM, Tummalacherla M, et al. A permutation
    approach to assess confounding in machine learning applications for digital
    health. In: *Proceedings of the 25th ACM SIGKDD International Conference on
    Knowledge Discovery and Data Mining*. 2019:54-64.
    doi:10.1145/3292500.3330903.
31. Zeng T, Li H, Zhang S, Tan J. Spurious model comparisons are widespread in
    biomedical artificial intelligence. *bioRxiv*. Posted 2026-05-20.
    doi:10.64898/2026.05.17.724301.
32. Xiong G, Bekiranov S, Zhang A. ProtoCell4P: an explainable prototype-based
    neural network for patient classification using single-cell RNA-seq.
    *Bioinformatics*. 2023;39:btad493. doi:10.1093/bioinformatics/btad493.
33. Xie Y, Yang J, Ouyang JF, Petretto E. scPanel: a tool for automatic
    identification of sparse gene panels for generalizable patient
    classification using scRNA-seq datasets. *Briefings in Bioinformatics*.
    2024;25:bbae482. doi:10.1093/bib/bbae482.
34. Wagle MM, Wang Y, Samanta S, et al. Deep interpretable learning of sample
    representations for characterizing disease states in single-cell
    transcriptomics. *bioRxiv*. Posted 2026-07-22.
    doi:10.64898/2026.07.21.738207.
35. Goeva A, Dolan MJ, Luu J, et al. HiDDEN: a machine learning method for
    detection of disease-relevant populations in case-control single-cell
    transcriptomics data. *Nature Communications*. 2024;15:9468.
    doi:10.1038/s41467-024-53666-8.
36. Pearl J. *Causality: Models, Reasoning and Inference*. 2nd ed. Cambridge
    University Press; 2009.
37. Peters J, Janzing D, Scholkopf B. *Elements of Causal Inference:
    Foundations and Learning Algorithms*. MIT Press; 2017.
38. Textor J, van der Zander B, Gilthorpe MS, Liskiewicz M, Ellison GTH.
    Robust causal inference using directed acyclic graphs: the R package
    'dagitty'. *International Journal of Epidemiology*. 2016;45:1887-1894.
    doi:10.1093/ije/dyw341.
39. Ren X, Wen W, Fan X, et al. COVID-19 immune features revealed by a large-scale
    single-cell transcriptome atlas. *Cell*. 2021;184:1895-1913.e19.
    doi:10.1016/j.cell.2021.01.053.
40. Stephenson E, Reynolds G, Botting RA, et al. Single-cell multi-omics analysis of
    the immune response in COVID-19. *Nature Medicine*. 2021;27:904-916.
    doi:10.1038/s41591-021-01329-2.
41. COMBAT Consortium. A blood atlas of COVID-19 defines hallmarks of disease severity
    and specificity. *Cell*. 2022;185:916-938.e58. doi:10.1016/j.cell.2022.01.012.
42. Yazar S, Alquicira-Hernandez J, Wing K, et al. Single-cell eQTL mapping identifies
    cell type-specific genetic control of autoimmune disease. *Science*.
    2022;376:eabf3041. doi:10.1126/science.abf3041.
