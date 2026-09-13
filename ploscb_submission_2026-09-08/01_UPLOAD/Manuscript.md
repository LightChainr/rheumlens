# Recorded metadata predicts the phenotype label across public single-cell cohorts and changes what internal validation can establish

**Short title:** Recorded metadata and single-cell classifier validation

Hongyu Ying<sup>1</sup>, Dandan Yun<sup>1</sup> and Dan Liu<sup>1,\*</sup>

<sup>1</sup> Department of Rheumatology and Immunology, Shanghai Pudong Hospital,
Fudan University Pudong Medical Center, Shanghai 201399, China
<sup>\*</sup> Correspondence: Dan Liu; danliu600@126.com

---

## Abstract

Patient-level single-cell classifiers routinely report cross-validated AUC above 0.95. A
high AUC supports disease-related discrimination only if the metadata recorded about
recruitment, processing and sequencing do not themselves predict the phenotype label. We
measured this association in nine phenotype contrasts from five public blood datasets
(809 donors), covering lupus, COVID-19, influenza, cytomegalovirus serostatus and sepsis
versus COVID-19. Recorded metadata were split into collection variables (study design), demographics (case mix or risk
factors) and sample-quality summaries (assay or disease state). Each group's out-of-fold
AUC was tested against a permutation null built from that comparison's own metadata.
Expression classifiers were then tested with a free label permutation, asking whether
expression is associated with the phenotype at all, and a collection-stratified
permutation, asking whether discrimination exceeds the recorded collection strata.

Metadata alone predicted the phenotype label in eight of nine comparisons, in all five
seeds and after Benjamini-Hochberg correction (out-of-fold AUC 0.404 to 1.000);
collection variables alone did so in four of the seven comparisons that record them.
Similar overall metadata AUCs masked different sources of predictability: sample quality
in one COVID-19 cohort, recruiting site in another. The comparisons fell into three
validation regimes. A prespecified cytomegalovirus negative control showed no
detectable metadata association. Seven comparisons were associated with metadata but had
label-mixed strata, so the conditional question could be asked; in sepsis versus COVID-19
only the free test rejected. In influenza no collection stratum held both labels, so it could
not be answered with these donors. Adding expression to metadata raised AUC by 0.087 to 0.354 in seven
comparisons and by nothing measurable in the other two.

Metadata-phenotype association is therefore a measurable property of a benchmark cohort
that changes what internal validation can support. These analyses quantify association
with recorded variables; they do not establish that discrimination is biological.

**Keywords:** patient-level classification; single-cell RNA sequencing; study metadata;
confounding; permutation tests; external validation; systemic lupus erythematosus;
COVID-19

---

## Author summary

Single-cell sequencing of a blood sample can be turned into a diagnostic score, and
published scores often separate patients from healthy donors almost perfectly. Before
trusting such a score, one needs to know whether it reflects the disease or the way the
study was run. In many public datasets, patients and healthy donors were recruited at
different hospitals, sequenced at different times or processed in different batches, and
a model can score well by recognising where a sample came from.

We measured how often this happens in nine comparisons from five public datasets covering
lupus, COVID-19, influenza and cytomegalovirus infection. Using only the recorded study
details, without any gene expression, we could predict which group a donor belonged to in
eight of the nine. The ninth, chosen in advance as a negative control, showed no such
pattern. Cohorts that looked equally affected were affected through different details.

The usual validation test shuffles patient labels and asks whether expression predicts
them at all. Shuffling labels only among donors collected the same way asks whether
expression predicts them beyond the collection process. The two tests can disagree, and
in one cohort the second cannot be run.

---

## 1. Introduction

Single-cell RNA sequencing supports patient-level classification by capturing differences
in immune-cell composition and molecular state, and published classifiers routinely
separate patients from controls at cross-validated AUC above 0.95. For such a result to
inform clinical research, the discrimination has to extend beyond the recruitment and
processing conditions of the development cohort. This is hard to assess when the
phenotype label is associated with collection site, sequencing batch or sample quality.

Many methods now produce one representation per patient: pseudobulk aggregation,
phenotype-oriented summaries, prototype networks, sparse gene panels, multicellular factor
models, graph representations, and methods that first identify disease-relevant cell
populations [6-11,32-35]. Frozen single-cell foundation models offer another route -
encode each cell, average per donor, train a small classifier - although benchmarks report
that such embeddings can retain batch structure and can underperform simple expression
baselines [13,14]. These methods are compared on public cohorts, and the comparison is
informative only if those cohorts do not let a model succeed by recognising how samples
were collected.

In public data, cases and controls are often collected differently: they come from
different clinics, years, chemistries or processing batches. Batch-label confounding has
been known to bias cross-validation for a decade [16], and a fully confounded single-cell
design is not identifiable [17]. Regression on technical covariates does not resolve the
problem, because it removes whatever part of the data the covariates predict, technical or
biological [18,19]. Work in predictive modelling has shown that adjustment must be fitted
inside training folds, that confound regression can leak information to nonlinear
learners, and that residual confounding can survive apparently successful correction
[27-29]. Benchmark studies raise the related question of how the datasets themselves are
assembled [12,15].

Existing work provides statistical tools for detecting confounding and testing conditional
association [25,29,30]. What remains poorly characterised in patient-level single-cell prediction is the benchmark cohort itself. It is not known how often recorded acquisition structure predicts the phenotype, which parts of the metadata carry that association, or whether standard validation procedures preserve or destroy it. Nor is it known when the available donors contain enough design overlap for a conditional question to be asked at all. We therefore
treat metadata-phenotype association as a measurable property of the cohort rather than as
a post hoc explanation for classifier performance, and examine its consequences across
nine phenotype contrasts from five public single-cell datasets.

We organise the analysis around four questions that any patient-level evaluation answers,
explicitly or not. **Presence:** does the recorded metadata predict the phenotype?
**Source:** which group of variables carries that prediction? **Validation:** does
expression discriminate beyond the acquisition structure preserved by the relevant null?
**Transport:** does the discrimination survive a genuinely different collection process?

Recorded metadata are not one kind of variable. Collection variables describe study
design in the strict sense. Demographics describe case mix and may also be risk factors.
Sample-quality summaries are measured downstream of the assay and can also be affected by
disease. We analyse the three groups separately, treat collection as primary for any
statement about study design, and report headline counts under both definitions.

The study makes three contributions. The first is empirical: across diseases, we quantify
phenotype predictability from recorded metadata as a cohort-level property of public
single-cell benchmarks and decompose it into collection, demographic and sample-quality
components. The second concerns validation: free and collection-stratified permutations
test different nulls, they give different answers in real cohorts, and in a cohort
without design overlap the stratified test has no informative reference distribution. The third concerns interpretation: two quantities often read as evidence of robustness have no single interpretation. The AUC lost after metadata adjustment depends on the shape of the adjustment matrix, and agreement in a second cohort can reflect a shared acquisition mechanism. Our contribution
is therefore not a new classifier or a new permutation test, but an empirical
characterisation of benchmark validity in patient-level single-cell prediction.

## 2. Study framework and evaluation measures

### 2.1. The setting

Figure 1 sets out the setting in the notation of graphical causal models [36,37]. A
recruitment process **R** places each donor in the recorded design **D** - site, batch,
calendar period, chemistry, sample-quality profile - and is also associated with the
phenotype label **Y**, because cases and controls are usually found through different
routes. The measured cells **X** are affected by **D** through technical effects and by
**Y** through biological ones. A classifier maps X to a prediction Ŷ.

Two further nodes matter in practice. **C** denotes covariates caused by the disease, such
as severity, WHO score, comorbidity, medication and time since symptom onset. They lie on
the path Y → C → X and are not collection variables; including them in D would turn
disease predictability into apparent design predictability, and Section 3.1 quantifies the
effect. **U** denotes collection structure that was never recorded, which no analysis here
can detect (Section 8). The graph was specified before the analysis and is released in dagitty syntax [38].

We draw the D-Y association as induced by R rather than as a directed edge, because both
directions occur: a patient may be recruited at a referral centre because of their
disease, and a centre's case mix may determine which diagnoses appear in a batch. The
checks below measure the association, not its direction.

### 2.2. The primary measure and a secondary variance summary

**Metadata-only AUC.** For each comparison we fit a classifier on the columns of the
recorded metadata matrix **D** alone and report its out-of-fold AUC. Its null distribution
comes from 200 permutations of the phenotype label, each passed through the same
fixed-hyperparameter pipeline. Four names are used consistently: the **metadata-only AUC**
uses every recorded variable, and the **collection-only**, **demographic-only** and
**sample-quality-only** AUCs use one group each. Each group has its own permutation test,
and a nonlinear arm is run beside the linear one.

**Residual label variance.** With *y* the centred phenotype vector, we also report

  V_D = 1 − R²(y ~ D),

the share of phenotype variation that the metadata cannot explain in a linear model. V_D is
a residual variance rather than an information measure, and it depends on how D is coded.
It is computed out of fold with five-fold cross-fitting, because in-sample R² inflates with
the ratio of columns to donors. The CMV comparison has 89 metadata columns for 108 donors;
its in-sample V_D of 0.468 would rank it the third most affected of the nine, although its
matched in-sample null mean is 0.496. Each V_D is compared with a null from 200 label
permutations through the same five folds, giving p(V_D). Unless marked in-sample, p(V_D)
refers to the cross-fitted statistic.

**Why the AUC carries the test.** V_D is an unpenalised linear R², so its cross-fitted
version loses power as D widens. In the Ren COVID-19 comparison, 26 one-hot columns on 182
donors leave the cross-fitted V_D at its ceiling, whereas a penalised logistic model on the
same matrix reaches AUC 0.81 to 0.85 (p = 0.005 in every seed). The metadata-only AUC is
therefore the primary test and V_D a secondary summary, and disagreements between them are
reported.

**Scope.** Both quantities concern the D-Y association in Figure 1. For a classifier to be
driven by metadata, D must also affect X and the classifier must use that part of X. A high
metadata-only AUC is therefore necessary but not sufficient evidence of metadata-driven
classification, and a low value does not exclude association through unrecorded structure.

### 2.3. Two permutation tests ask two validation questions

The **free permutation** shuffles the phenotype label across all donors and refits the
whole pipeline. Its null is exchangeability of the label across donors, so rejection
indicates an association between expression and the label, whatever its origin. Because
the shuffle also destroys the metadata-label association, a small p-value cannot
distinguish disease signal from collection structure. Rejection also does not certify the
absence of leakage: a label-dependent step performed before permutation, such as a feature
list carried over from an earlier analysis, affects observed and permuted runs alike. In
our pipeline every representation step is refitted inside the training fold on every
permutation (Section 9.6); that property is established by the code, not by the test.

The **collection-stratified permutation** shuffles labels only within strata that share a
collection configuration - batch, pool, site and, where it varies, sequencing chemistry.
Each stratum keeps its class counts, and only the pairing of labels with molecular
profiles is randomised. This restricted permutation scheme was introduced by Chaibub Neto
and colleagues [30] and has since been used in the confounding literature [25]. Its null
is exchangeability within the collection strata, so rejection indicates discrimination
beyond what those strata retain.

The strata use collection variables only; sample-quality and demographic columns are not
conditioned on. If a sample-quality variable remains associated with the phenotype within
a stratum and also affects expression, the stratified null is false for a reason unrelated
to disease; Section 4.4 quantifies this case. Conversely, a non-significant stratified test
is a failure to reject under this conditioning set, not evidence that discrimination is
entirely collection-driven. Conditioning on the full metadata, including continuous
columns, would require a conditional randomisation scheme that we do not develop here.

The stratified test also requires at least one stratum containing both classes. When every
stratum is label-pure, within-stratum shuffling leaves the labels unchanged and the
reference distribution is degenerate. We then report the test as **uninformative**, with
p = NA. The free permutation remains defined and is still reported.

### 2.4. Three validation regimes

Together, the metadata-only AUC and the stratum structure place each comparison in one of
three regimes. The regimes are defined by the cohort, not by the classifier.

1. **Low recorded association.** The metadata-only AUC is not distinguishable from its
   null. No recorded variable undermines the unadjusted estimate, although unrecorded
   structure remains possible.
2. **Associated but conditionally evaluable.** The metadata predicts the phenotype, and at
   least one stratum holds both labels, so the stratified test can be run. Its result
   indicates whether discrimination exceeds what the collection strata retain.
3. **Structurally non-overlapping.** No collection stratum holds both labels. The
   conditional question cannot be answered with these donors by any analysis; answering it
   requires donors collected under overlapping conditions. Such comparisons are labelled
   **not estimable** in Table 1 and reported outside the main ordering.

Separately, a comparison is flagged **underpowered** when it has fewer than 40 donors or a
minority class below 15. Flagged comparisons stay in the ordering because their
instability is informative: Section 3.4 reports one whose collection-stratified p ranges
from 0.057 to 0.614 across five seeds. These gates depend on donor counts and stratum
structure, not on the seed. The fixed-hyperparameter pipeline is the inferential statistic
throughout, because both permutation nulls are built from it; the tuned AUC is released
beside it as a description, and their difference is not used as a gate.

# Results

Section 3 describes the nine comparisons, Section 4 uses simulation to interpret the
patterns they show, Section 5 asks whether classifiers in two lupus cohorts use the
detected structure, and Section 6 brings the checks together.

## 3. Metadata-phenotype association across nine public comparisons

### 3.1. Datasets and recorded metadata

We assembled five public datasets [2,39-42] through the CELLxGENE Discover curation API and
reduced each to one row per donor: log1p CPM pseudobulk expression, plus variables built
from standard CELLxGENE schema fields (donor identifier, disease, sex, developmental stage, assay and tissue) and cohort-specific collection variables recovered from the distributed object. The datasets were chosen to span different
acquisition structures rather than sampled systematically; Section 9.1 gives the selection
rationale and a retrospective count of eligible datasets.

The metadata variables fall into three groups.

- **Collection** (site, sub-study, sequencing chemistry, assay, batch, pool) describes how
  material was gathered and processed, and is not caused by the disease.
- **Demographic** (age, sex, ethnicity) describes case mix. An imbalance is a fact about
  recruitment, but age and sex are also risk factors for most of these diagnoses, and the
  two readings cannot be separated with these data.
- **Sample quality** (cells per donor, mean UMI, genes per cell, mitochondrial fraction)
  reflects the assay, but can also reflect disease through cell composition, treatment or
  the state of the sample at collection.

Collection is the primary definition for any statement about study design, and the other
two groups are prespecified sensitivity analyses. The union of all three is called the
*recorded metadata*. Covariates caused by the disease - severity, outcome, comorbidity,
medication, symptom timing, diagnosis fields, stage, grade, WHO score - are excluded from
all groups (node C in Figure 1). Table S1 lists the raw and expanded columns entering every
group; we verified that it reproduces the width of all 34 fitted metadata matrices and that no label or label proxy enters any of them.

The exclusion matters. In a simulation where metadata and phenotype are unrelated,
admitting one disease-caused covariate raises the metadata-only AUC from 0.488 to 0.793 and
flags the cohort as significant in 100% of runs, against 5.5% without it (Figure S1,
Table S2). Residualising on the enlarged set then removes 0.12 AUC from a cohort with no
confounding.

Nine phenotype contrasts were formed (Table 1): lupus versus healthy in GSE174188 (CD4
subset and all cells), COVID-19 versus healthy in three cohorts, influenza versus healthy
and sepsis versus COVID-19 in COMBAT, and CMV-positive versus CMV-negative in HIHA. The five
cohorts contain 809 distinct donors (261, 196, 120, 124 and 108). Comparisons from the same
cohort share donors, so the per-comparison counts in Table 1 do not sum to 809.

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

**Table 1. The nine comparisons, ordered by metadata-only AUC.** Each comparison was run at five split seeds; values are
ranges over the five seeds, and a single value means that all five agreed at the printed
precision. "Minority" is the smaller label group. "Metadata AUC" uses all three variable
groups; "Collection AUC" uses the collection block alone, and n/a means that no varying
collection variable is recorded. All AUCs are cross-fitted with the prespecified
fixed-hyperparameter pipeline, which is the statistic tested by each adjacent p-value and
the one both permutation tests use. Tuned nested-CV AUCs, column counts, V_D and p(V_D)
are given per seed in Table S3. "+ over metadata" is the change in out-of-fold AUC when
expression is added to the metadata in one joint model (Section 3.5, Table S4). Each p is
one-sided against a permutation null built from that comparison's own metadata matrix. The
floors 0.0050 and 0.0010 are 1/201 and 1/1001: the metadata-only AUC and V_D nulls use 200
permutations and the complete-pipeline tests 1,000. A dash under "p strat." marks an
uninformative stratified test, because every stratum is label-pure. A dagger marks the two
comparisons that record no collection variable; their stratification uses tertiles of log
cells per donor, so the test is sample-quality-stratified. Verdicts come from the screen's
gates: **not estimable** when no collection stratum holds both labels (regime 3,
Section 2.4), and **underpowered** when the minority group has fewer than 15 donors.

### 3.2. Recorded metadata predicts the phenotype label in eight of nine comparisons

Ordered by metadata-only AUC, the nine comparisons range from 0.404 to 1.000 (Figure 2A,
Table 1), with no gap that would justify a threshold. Positions are stable across seeds:
excluding the CMV control, whose near-chance value varies by 0.114, the widest five-seed
spread is 0.042. Eight of the nine comparisons show significant metadata association in
every seed. The count remains eight after Benjamini-Hochberg correction [24] within each
seed, whether the family is that seed's nine all-metadata tests or all 34 of its
variable-group tests. The ninth is the CMV negative control (Section 3.3).

Counting by p(V_D) alone would give six rather than eight. In COVID-19 Ren (26 columns, 182
donors) and COVID-19 COMBAT (7 columns, 110 donors), the cross-fitted V_D ranges across
seeds from 0.847 to its ceiling of 1.000 and from 0.906 to 1.000. At the ceiling the matched
null reaches the ceiling too, so p(V_D) is 1.000, whereas other seeds give 0.005. The
metadata-only AUCs on the same matrices are stable, 0.812-0.854 and 0.862-0.884, and
significant in every seed (p = 0.005 and 0.005-0.010).

The count is smaller under the collection-only definition. Seven comparisons record a
varying collection variable, and in four of them the collection-only AUC is significant in
all five seeds. Two of those four are the comparisons flagged as not estimable or
underpowered. Among comparisons that are estimable and record collection variables, two of
four are significant: lupus in GSE174188 CD4 (0.762-0.774) and COVID-19 in Ren
(0.729-0.768). The two comparisons without collection variables, GSE174188 all cells and
the single-assay Ren subset, are stratified on tertiles of log cells per donor and marked
in Table 1. Much of the association detected by the full metadata is therefore carried by
sample quality and case mix rather than by study design.

**Similar overall association, different sources.** Two COVID-19 cohorts of similar size
have closely comparable metadata-only AUCs but different drivers (Figure 2B, C).

- **Stephenson** (n = 115; metadata-only AUC 0.845-0.882). Sample-quality variables reach
  0.865-0.884 (p = 0.005 in every seed), whereas collection variables reach 0.501-0.594 and
  are not consistently significant (p = 0.015-0.234).
- **Ren** (n = 182; metadata-only AUC 0.812-0.854). Collection variables reach 0.729-0.768
  (p = 0.005 in every seed), whereas sample quality reaches 0.590-0.632 and is not
  significant (p = 0.050-0.154).

The appropriate follow-up differs: library quality and harmonised preprocessing in
Stephenson, recruiting site and sub-study composition in Ren. A single metadata AUC would
treat the two cohorts as equally affected. Reporting by variable group is also our
response to the dependence of V_D on how D is coded.

**Incomplete metadata can hide association.** With only sequencing chemistry recorded as
collection metadata, the Ren collection-only AUC is 0.553. Recovering the recruiting
hospital (12 centres, five of them case-only) and the sub-study (six source studies) raises
it to 0.750, and raises the metadata-only AUC from 0.758 to 0.818. A reassuring result
computed on sparse metadata is weak evidence.

### 3.3. A negative control specified in advance

We selected the CMV cohort as a negative control before computing any statistic. Its 108
donors were assayed with one chemistry and one suspension type, and CMV status is a
serological result rather than a recruitment criterion, so the collection labels were
expected to be nearly independent of serostatus. The cohort is not a single batch: its
metadata record 36 batch and 46 pool identifiers, which cross into 46 strata, 17 of them
holding one donor. Each donor contributes one row and one label.

As predicted, CMV shows the weakest metadata association. Across five seeds the
metadata-only AUC is 0.404 to 0.518 (p 0.383 to 0.881). Collection variables give
0.365 to 0.489 (p 0.527 to 0.945), and sample quality 0.453 to 0.507 (p 0.478 to 0.716). The
expression classifier on the same donors reaches AUC 0.894 to 0.933, with both permutation
tests at the 1/1001 floor in every seed.

Demographic variables are the exception, with AUC 0.591 to 0.639 and p 0.015 to 0.085.
After Benjamini-Hochberg correction across each seed's 34 variable-group tests, q = 0.020 in
two seeds and 0.074 to 0.115 in the other three. Sex and ethnicity are established
correlates of CMV seroprevalence, so this association may reflect aetiology as much as
recruitment. We therefore state the negative control at the level of the collection block,
where no seed shows association.

Two further points concern interpretation. Read without a null, the in-sample V_D of 0.468
would suggest strong association; its matched in-sample null mean is 0.496
(p(V_D, in-sample) 0.284 to 0.353), and the cross-fitted V_D is at its ceiling of 1.000
with p = 1.000. And no detected association does not mean no association. The collection
block has 80 one-hot columns but a centred rank of 45, and the full block 89 columns at rank
54, so column counts overstate the number of independent directions tested. What the
comparison shows is that the recorded collection labels carry no detectable information
about CMV status at this sample size.

### 3.4. Free and stratified permutations answer different validation questions

The two permutation tests agree in seven comparisons. The two exceptions illustrate
regimes 2 and 3 of Section 2.4.

**Sepsis versus COVID-19: the tests disagree.** In COMBAT (n = 111, minority class 11) the
metadata-only AUC is 0.996 to 0.998. The free permutation gives p = 0.001 in every seed,
whereas the collection-stratified permutation gives p = 0.057 to 0.614 and is significant
in no seed (Figure 2D). Expression is associated with the phenotype, but there is no
evidence of discrimination beyond the collection strata. With 11 donors in the minority
class the stratified test has little power, so this is not evidence that discrimination is
entirely collection-driven. The fixed-hyperparameter pipeline also scores 0.011 to 0.066
AUC above the tuned one here, the largest gap among the eight comparisons with a
stratified p-value.

**Influenza: the conditional question cannot be asked.** The COMBAT influenza comparison has
21 donors, and metadata separates the labels exactly (metadata-only AUC 1.000, V_D 0.000).
Every collection stratum is label-pure (Figure S2), so the stratified test is
uninformative. The free permutation rejects at the 1/1001 floor in all five seeds, with an
expression AUC of 0.982 to 1.000. Reported alone, the free test would suggest a
well-validated classifier. Together, the two results show that recorded collection and
phenotype cannot be separated with these donors, and that the test designed for this
question has no reference distribution.

All three regimes therefore occur in public data: low recorded association (CMV),
association with conditionally evaluable strata (seven comparisons), and structural
non-overlap (influenza). In the third regime no choice of statistic recovers the missing
comparison; only donors collected under overlapping conditions would.

### 3.5. What expression adds over the recorded metadata

The metadata-only and expression AUCs show whether each predicts the phenotype, but not
whether expression adds information beyond the metadata. We therefore fitted a joint
model. Within each training fold, expression is reduced to principal components fitted on
the training donors, the metadata matrix is built from the same donors, and one classifier
is fitted on both. The three models share the fixed-hyperparameter specification and the
folds (Table S4).

In seven comparisons, adding expression raises the AUC by 0.087 to 0.354 across seeds, and
the joint model is within 0.023 of the expression-only model. Lupus in GSE174188 CD4 gains
0.087 to 0.101; the single-assay Ren subset, which records the least collection metadata,
gains most (0.227 to 0.262).

The two exceptions are the comparisons discussed in Section 3.4. In sepsis versus COVID-19,
metadata alone reaches 0.996 to 0.998 and adding expression changes the AUC by −0.001 to
+0.001, whereas metadata adds 0.036 to 0.104 to expression alone. In influenza both models
are at 1.000 and the increment is zero in every seed. The increment thus agrees with the
permutation results without requiring a permutation.

In CMV the joint model is 0.069 to 0.136 below expression alone. Appending 89
uninformative one-hot columns to 50 principal components reduces discrimination, the same
property that makes residualisation on this matrix costly (Sections 4.5 and 5.3).

The increment is descriptive and carries no p-value; a test would have to hold the metadata
association fixed during permutation, which the stratified test in Table 1 already does.
Neither the ordering nor the increment shows that a particular classifier uses
metadata-related structure. Section 5 examines that path directly.

## 4. Simulation: what these patterns can and cannot mean

The simulations address the patterns seen in the real data: whether the metadata statistic
is calibrated, why residualisation and restriction disagree, when an external cohort tests
confounding, and what the stratified test and residualisation loss can establish. We
crossed seven data-generating regimes with seven levels of metadata-phenotype association
and two outcome types, each predicted and scored on its own scale: 19,600 simulated cohorts
of 200 donors and 100 features (Figure 3). No simulation parameter was fitted to the real
cohorts.

### 4.1. Calibration of the metadata-association statistic

When design carries no phenotype information, p(V_D) ≤ 0.05 occurs in 0.025 to 0.065 of
runs across the seven regimes, close to the nominal 0.05 (Figure 3A). The rate rises at
ρ = 0.25 and reaches one by ρ = 0.50 in every regime. Cross-fitting and a per-cohort null
provide this calibration; an in-sample estimator judged against an in-sample null does not.

### 4.2. Residualisation and restriction diverge in every regime

Residualisation regresses the data on metadata variables and analyses the residuals;
restriction compares donors collected under matching conditions. In the additive linear
regime the two agree when design carries no phenotype information and diverge as the
association increases (Figure 3B). At ρ = 0 the unadjusted, residualised, restricted and
external AUCs are 0.809, 0.809, 0.787 and 0.817; at ρ = 0.98 they are 0.879, 0.538, 0.724
and 0.891. Residualisation removes 0.34 AUC and restriction 0.16, although the biological
effect is unchanged. At ρ = 1 restriction is undefined, because no stratum contains both
labels.

The gap appears in every regime (Figure 3C). The continuous-outcome arm uses a continuous
latent outcome, ridge regression and pairwise concordance, which equals AUC for a binary
outcome. Restriction outperforms residualisation at ρ = 0.98 in all seven regimes for both
outcome types, but the magnitude differs: 0.094 to 0.200 with binary outcomes and 0.013 to
0.038 with continuous ones. A simulation run only on continuous outcomes would understate
the gap about five-fold.

### 4.3. Threshold alignment and external validation

Class imbalance changes two things at once: prevalence, and the alignment between
collection boundaries and the phenotype threshold. In the imbalanced regime both thresholds are at the 80th percentile; in the imbalanced-offset regime the phenotype threshold is at
the 80th percentile and the collection boundary at the median. With matched thresholds,
imbalance has little effect: the metadata-only AUC reaches 1.000 at ρ = 1 and V_D reaches
0.000. With offset thresholds, the metadata-only AUC saturates at 0.880 and V_D at 0.728.
The apparent reassurance of a rare phenotype therefore comes from threshold alignment, not
from prevalence. This matters for the two real comparisons with minority classes of 10 and
11 donors.

In the base regime, external AUC rises with confounding, from 0.817 at ρ = 0 to 0.904 at
ρ = 1, tracking internal AUC (Figure 3D). When the target cohort has a different design
mechanism, external AUC stays at 0.816-0.829 while internal AUC rises to 0.899. External
validation therefore tests metadata-driven discrimination only to the extent that the
second collection process differs from the first. A different accession does not ensure
this: cohorts assembled by similar consortia, on similar platforms and through similar
recruitment routes can reproduce the same bias.

### 4.4. A non-disease pathway outside the conditioning set makes the stratified test reject

The stratified permutation does not condition on sample-quality or demographic variables.
We simulated cohorts of 200 donors across eight collection sites in which a sample-quality
variable QC predicts the phenotype within sites at strength γ, and expression depends only
on site and QC. The generating model contains no effect of the phenotype on expression.
Both permutation tests were run as for the real cohorts (Figure S3, Table S5).

Because QC → Y and QC → X, conditioning on the collection strata leaves X and Y dependent
whenever γ > 0. The stratified null is then false, and its rejections are correct
rejections rather than type-I errors; only γ = 0 measures calibration. At γ = 0 the
collection-stratified permutation rejects in 19 of 300 runs, 6.3% (exact 95% CI 3.9 to 9.7;
binomial p = 0.29 against a nominal 5%), and the free permutation in 21 of 300, 7.0% (4.4 to
10.5; p = 0.11), with a mean AUC of 0.500.

As γ increases the two tests move together. At a mean within-stratum absolute correlation of
0.29 the stratified and free tests reject in 22.0% and 22.7% of runs; at 0.61, in 95.3% and
95.7%. Classifier AUC rises from 0.500 to 0.668 without any disease effect. A significant
stratified result therefore shows association beyond the collection strata, not disease
biology, because a recorded variable outside the conditioning set can produce it. Where
sample-quality or demographic variables are themselves associated with the phenotype
(Table S3), this possibility remains open.

### 4.5. Residualisation loss without metadata-phenotype association

Residualisation loss is often read as a measure of contamination. We generated cohorts in
which the metadata matrix is independent of the phenotype, calibrated the biological effect
to an unadjusted AUC near 0.90, and residualised on matrices of increasing width. The
unadjusted and residualised models are the same function call with one added argument, so
they share folds, feature filter, standardiser, PCA, classifier and scoring (Section 9.9).
A second layout reproduces the level sizes of the real CMV batch-by-pool crossing: 46
levels for 108 donors, 17 of them holding one donor.

Loss occurs without any association. With one to eight columns it is small (mean 0.004 to
0.023). It peaks at 0.333 for 48 columns and 0.302 for 64, then falls to 0.222 at 80 and
0.159 at 92, the width closest to CMV's 89, as additional one-hot columns become redundant
and ridge shrinks them together. Variation between runs is wide; at 48 columns the central
95% of runs spans 0.224 to 0.460.

The CMV-shaped layout realises 45 columns and costs **0.363** on average (central 95% 0.272
to 0.464), more than balanced blocks of any width and enough to contain the 0.326 observed
in CMV (Section 5.3, Figure S4). The metadata-only AUC in these cohorts stays at 0.48 to
0.51. A loss of this size is therefore not evidence that anything was removed. The result
applies to this adjustment procedure - ridge at a fixed penalty on a training-fold
standardised matrix - and to this layout. Column count is not effective dimension, and
other residualisers would trade loss against completeness differently.

## 5. Whether classifiers use metadata-related structure: two lupus cohorts

Section 3 measures the D-Y association. This section examines the D → X → Ŷ path in two
lupus cohorts with complete donor metadata, using frozen Geneformer [4] embeddings and two
pseudobulk expression representations in one fold-contained pipeline. GSE285773 [3] is the
independent transfer cohort. Figure 4A shows how cases and controls are distributed across
the recorded design, and Figure 4B how well metadata alone predicts disease status.

### 5.1. All three representations encode batch

In GSE174188 the three representations classified disease at AUC 0.982, 0.977 and 0.975,
and identified processing wave, one against the rest, at 0.9997, 0.999 and 0.997. In
GSE135779 the disease AUCs were 0.870, 0.933 and 0.945, and the best single-batch AUCs
0.993, 0.961 and 0.961 (Figure 5A). The representations contain batch information; this
alone does not show that a disease classifier uses it.

### 5.2. Batch predictability did not imply sensitivity to batch adjustment

In GSE135779, Geneformer embeddings predicted the most distinguishable batch at AUC 0.993
and both pseudobulk representations at 0.961, whereas batch alone predicted disease at AUC
0.499. Batch residualisation changed disease AUC by +0.002, −0.002 and −0.009 (Figures 5B,
6C), within split-to-split variation. Dependence on batch was therefore small, though not
shown to be zero, and batch information in a representation did not imply batch use by the
classifier.

### 5.3. Residualisation and restriction answer different questions

Restriction asks whether separation persists among donors collected under overlapping
conditions; residualisation removes the part of the data predictable from chosen metadata.
When metadata and phenotype are collinear these are different questions, and we compare
them to show that they are not interchangeable.

In GSE174188, ridge residualisation on processing-wave fractions reduced the three
representations from 0.982, 0.977 and 0.975 to 0.714, 0.747 and 0.735. Restriction to the
one large, near-balanced wave, compared with size- and label-matched random subsets, cost
0.014 to 0.049 (Figure 6B). The paired difference between the two losses was 0.186 to
0.255 across six comparisons, with every descriptive interval above zero.

In GSE135779, removing the single all-case batch left 36 donors with both labels in every
remaining batch and changed matched performance by at most 0.012. Five ways of removing the
recorded variables gave five answers (Figure 6A). Ridge on the complete set reduced the
representations to 0.516, 0.522 and 0.520, a random-forest residualiser to 0.748, 0.787
and 0.795, and overlap weighting and batch location adjustment retained 0.859 to 0.937.
Every full-set removal cost more than the batch-only control.

When design predicts the phenotype, the part of the data predictable from design
necessarily includes phenotype-related variation, so removing it removes disease signal
whatever its origin. A large residualisation loss is therefore consistent both with
removing a shortcut and with removing biology.

CMV shows this without any detected association: metadata-only AUC 0.512 under this
pipeline, permutation p 0.383 to 0.881 across seeds, and p(V_D) = 1.000. Residualising on
the same variables reduces AUC from 0.925 to 0.599, a loss of 0.326 (Figure 6C). The matrix
has 89 one-hot columns for 108 donors at a centred rank of 54. In simulation (Section 4.5),
a phenotype-independent matrix with the same ragged layout loses 0.363 on average, which
contains the observed value (Figure S4). The size of a residualisation loss therefore
depends on the adjustment matrix and cannot be read as a measure of contamination.

Nor do residualised and restricted estimates bound a biological effect. They are estimates
under different assumptions, neither is identified, and their disagreement is the
informative result.

### 5.4. Cell-type composition reproduces the pattern

Summarising the same GSE174188 CD4 cells as naive, effector-memory, regulatory and
unassigned fractions gave disease AUC 0.899 to 0.915, and the same four fractions predicted
processing wave at 0.803 to 0.872. Wave residualisation reduced disease AUC to 0.689 to
0.709 (Figure S5A), whereas within-wave restriction against matched subsets cost 0.027 to
0.062, with every interval crossing zero (Figure S5B). A four-dimensional biological
summary thus reproduces the pattern, which places it in the cohort rather than in
foundation-model embeddings and is consistent with composition being a strong
patient-stratification baseline [26].

### 5.5. Source restriction and cross-cohort discrimination

Restricting training to one processing wave did not improve transfer. Training on all 261
GSE174188 donors and testing in the independent 26-donor GSE285773 cohort gave 0.900, 0.944
and 0.969. Training on the dominant wave lowered every value and performed no better than a
random subset of the same size (Figure 7). Within GSE174188, training on waves 2 and 3 and
testing on wave 4 gave 0.842, 0.811 and 0.806, and the reverse direction 0.919, 0.926 and
0.821. More varied source donors transferred better than a restricted source.

In the other direction, from the 26-donor cohort to the 261-donor target under strict
source-only rules, AUCs were 0.884 [0.843-0.920], 0.919 [0.884-0.946] and
0.926 [0.895-0.953]. Both pseudobulk representations outperformed Geneformer after paired
DeLong testing [23] and Benjamini-Hochberg adjustment [24] (q = 0.0151, q = 0.000603). Three
learned pooling methods on the same frozen embeddings never outperformed mean pooling, and
every significant paired difference was negative (Figure S6, Table S6). Increasing model
capacity improved fitting to the source cohort but did not improve cross-cohort transfer.

## 6. Integrating metadata assessment and predictive validation

The five checks below examine different arrows in Figure 1 and fail in different ways. They
are not a validity standard: we have not shown that they are sufficient or minimal, or that
passing them supports a causal claim. Their value is that they can disagree.

| # | Check | Shows | Does **not** show |
|---|---|---|---|
| 1 | Metadata-only AUC by variable group, cross-fitted V_D, p(V_D) against the cohort's own null, overlap counts | Whether recorded metadata predicts the phenotype label, and which variable groups do | That the classifier uses that metadata; anything about unrecorded structure |
| 2 | Complete-pipeline permutation, both free **and** collection-stratified | Free: an expression-phenotype association of any origin. Collection-stratified: discrimination beyond what the collection strata retain | A mechanism, or the absence of leakage. The stratified test does not cover sample-quality or demographic association within strata (Section 4.4), and is uninformative when every stratum is label-pure |
| 3 | How well the representation predicts batch | That the representation encodes collection information | That the decision rule depends on it (Section 5.2) |
| 4 | Residualisation **paired with** restriction and matched donor controls | How far two non-equivalent adjustments disagree | Which adjustment is right; bounds on a biological effect |
| 5 | Source-only external target | Whether patient ranking survives a change of collection process | Calibration or clinical usefulness; little about confounding when the target was collected the same way (Section 4.3) |

Checks 1 and 2 take minutes on a donor-level table. Check 5 alone addresses transport, and
no combination of internal checks substitutes for it.

The four questions of Section 1 map onto these checks: presence and source onto check 1,
validation onto check 2 with support from checks 3 and 4, and transport onto check 5. For a
clinical reader, an internal cross-validated AUC measures discrimination under the
development cohort's sampling conditions. A metadata-only model shows whether recorded
collection differences could produce that discrimination. A conditioned validation - the
stratified permutation, or restriction to overlapping strata - shows whether the
discrimination exceeds those differences. An independent cohort collected differently
tests transfer. Together, these evaluations provide stronger evidence for a candidate
diagnostic classifier than an internal AUC alone.

Applied to the nine comparisons, the checks reproduce the three regimes of Section 2.4. CMV
shows no detected metadata association. In both lupus comparisons and all three COVID-19
cohorts, the metadata predicts the phenotype but the stratified test still rejects, so
discrimination exceeds what the strata retain. In sepsis versus COVID-19 the stratified test
does not reject, and in influenza it is uninformative; for these two, no internal analysis
settles the question, and the limitation lies in the cohort rather than in the classifier.

### 6.1. Three comparisons traced through the tree

Figure 8 traces four comparisons through the decision tree. Three are described here, and
sepsis versus COVID-19 in Section 3.4. The walkthrough imports the strata, metadata matrix
and classifier from the screen, so Q2 and Q4 use the same partition, columns and pipeline.
Each branch reports what can be estimated and what is missing (Table S7). It does not
prescribe an adjustment, which depends on the target quantity and the scientific question.

**COMBAT influenza (n = 21) stops at Q0.** Its stratum variable, the recording institute,
gives two strata, neither containing both a case and a control (Figure S2). The stratified
test is uninformative, and no conditioned estimate can be formed. This is the only stopping
criterion. The metadata-only AUC of 1.000 is consistent with the same structure but does not
trigger the stop, because a finite-sample AUC of 1.000 records complete separation of these
donors rather than proving non-identifiability. The fixed and tuned pipelines agree closely
(0.982 to 1.000 and 0.964 to 1.000), and the free permutation alone gives p = 0.001 in every
seed.

**CMV HIHA (n = 108) reaches Route A.** The stratified permutation is defined - 46
batch-by-pool strata, 21 of them holding both labels - and both tests reject at the 1/1001
floor. At Q3 no metadata association is detected: metadata-only AUC 0.512 in the
walkthrough pipeline and 0.404 to 0.518 in the screen, not significant in any seed. Route A
reports that no collection association was detected and that the unadjusted estimate
stands. It does not recommend residualisation, which here costs 0.326 AUC (Section 5.3),
close to the 0.363 lost by a phenotype-independent matrix of the same shape (Section 4.5).
An undetected association is not thereby shown to be negligible, but non-detection is not a
reason to adjust.

**COVID-19 Ren (n = 182) reaches Route C.** Both tests reject, and the metadata association
is strong: metadata-only AUC 0.832 in the walkthrough and 0.812 to 0.854 in the screen.
Label-mixed strata exist, so Q4 asks whether restriction is feasible, which requires at
least one stratum with five donors of each class and 40 donors in total, the power floor
used throughout. Ren does not meet this requirement. Its only stratum large enough for
cross-validation holds 18 donors; the restricted AUC there is 0.992 against 0.733 in size-
and composition-matched random subsets, with a matched-subset standard deviation of 0.145.
Residualising the 26-column matrix reduces AUC from 0.967 to 0.884, a loss of 0.083. Route C
therefore reports that association is present, that both adjustments are available but
disagree, and that neither is estimable with useful precision; donors collected under
overlapping conditions are needed.

**Restriction is rarely feasible at these strata.** CMV has 21 label-mixed strata but none
with five donors of each class, and Ren has one usable stratum of 18 donors. Route B assumes
that restriction is a practical alternative to residualisation, but in real cohorts with
fine strata it often is not. We do not coarsen the strata to make Route B available,
because that would change the conditioning set used by the stratified test.

# Discussion

## 7. Implications for patient-level single-cell evaluation

Patient-level single-cell validation is not solely a property of the classifier. It
depends jointly on the phenotype-acquisition structure of the cohort and on which parts of
that structure the validation procedure preserves. The nine comparisons bear on each of the
four questions set out in Section 1.

**Presence.** Phenotype predictability from recorded metadata was measurable, graded and
stable across seeds, with metadata-only AUCs from 0.404 to 1.000. It has to be judged
against a null built from the same matrix: an in-sample R² ranks CMV, the least affected
comparison, as the third most affected, because its metadata matrix is wide relative to its
donor count. The count of eight in nine describes these comparisons, which are neither
independent nor a random sample of public cohorts (Sections 8 and 9.1); it is not an
estimate of prevalence.

**Source.** Similar overall metadata AUCs concealed different associations. In Stephenson,
sample-quality variables reached 0.865 to 0.884 and collection variables 0.501 to 0.594; in
Ren, collection variables reached 0.729 to 0.768 and sample quality 0.590 to 0.632. The
overall AUCs differ by about 0.02, yet the follow-up differs: preprocessing and library
quality in one cohort, recruiting site and sub-study in the other. A combined metadata AUC
alone would hide this.

**Validation.** Free label permutation tests for any expression-phenotype association. It
cannot test whether a classifier exploits collection structure, because it destroys that
structure. The collection-stratified permutation [30] asks the conditional question, and
across these cohorts the two tests gave three different pictures: agreement in seven
comparisons, disagreement in sepsis versus COVID-19, and no informative reference
distribution in influenza. The last case is a substantive finding rather than a technical
failure. Some validation questions require design overlap, and without it no statistical
method can recover the missing comparison. The stratified test has its own limit: it
conditions on collection strata only, and a sample-quality variable associated with the
phenotype within strata can make it reject without any disease effect (Section 4.4). A
detected conditional association is therefore not evidence of disease biology. Studies that
use label permutation to support a patient-level classifier should report the stratified
version beside it, or state that no label-mixed stratum exists. The increment of expression
over metadata points the same way: 0.087 to 0.354 in seven comparisons, and essentially zero
in the two where the stratified test did not support discrimination beyond collection.

**Transport.** Agreement in an external cohort is weaker evidence than it appears when that
cohort shares the source's acquisition mechanism; in simulation, external AUC then rises with
confounding. Its value depends on how different the second collection process is, not on
whether the data carry a different accession. Residualisation loss has a similar limit: it
depends on the shape of the adjustment matrix, and in CMV it reached 0.326 with no detected
association.

For representation research, pseudobulk outperformed frozen embeddings in strict source-only
transfer in our lupus data, and learned pooling did not change this. More generally, internal
accuracy in these cohorts does not measure transportable discrimination, and no model can
correct a cohort in which design predicts the phenotype. That requires recording and
releasing collection metadata, and building cohorts in which cases and controls are
collected together within batches.

None of these results supports a diagnostic claim; the external AUCs rank 261 or 26 donors
and are not calibrated probabilities. The contribution is an evaluation layer rather than
another representation learner, and it addresses questions a model leaderboard cannot. We
suggest that patient-level classifiers be reported with three additions: the metadata-only
AUC of the development cohort by variable group; both permutation tests, or a statement that
the stratified test is uninformative; and an external result on a cohort collected through a
different route, fitted on the source cohort alone. These additions make cohort acquisition
structure an explicit part of classifier evaluation rather than an unexamined property of
the benchmark.

## 8. Limitations

These analyses address association with recorded metadata, not causal biological
specificity. Unmeasured collection structure - referral patterns, institutional treatment
protocols, sample handling absent from the released metadata - is invisible to every
quantity reported here, so a low metadata-only AUC is weak evidence that a cohort is
unaffected. Within Ren, omitting the recruiting hospital gives a collection-only AUC of 0.553,
against 0.750 with the completed metadata.

The collection-stratified permutation conditions on a subset of the recorded metadata
(Section 4.4). A null that conditions on the full matrix, including continuous columns,
would require a conditional randomisation scheme that we do not develop; this is the main
methodological gap of the study.

Only the collection block represents study design in the strict sense, and demographic and
sample-quality associations have other interpretations, so both definitions are reported.
V_D is linear and depends on coding; the nonlinear arms use two tree ensembles rather than
an exhaustive set of learners; and column counts overstate effective dimension.

Residualisation results are specific to ridge at a fixed penalty on a training-fold
standardised matrix. The simulation shows that, under this procedure, a loss of the size seen
in CMV carries no information about contamination, but not what other residualisers would do;
the random-forest, overlap-weighting and location-adjustment arms in Section 5.3 already
differ materially.

The decision tree organises evidence but has not been validated against an external
criterion, and no result depends on it.

The nine comparisons come from five datasets, share donors within datasets, and are not
independent replications. The datasets were selected to span acquisition structures, not
sampled from all eligible public cohorts (Section 9.1), so eight of nine is not a prevalence
estimate. Three comparisons have fewer than 15 donors in the minority class; they are
flagged, and no conclusion rests on them alone. All datasets are peripheral blood from public
repositories curated under one schema, so the results do not establish cross-tissue
behaviour. Matched random subsets separate the cost of restriction from the cost of using
fewer donors but do not identify a technical effect, so signal retained after restriction is
not shown to be biological. Geneformer pretraining overlap with the benchmark cells was not
tested.

# Methods

## 9. Materials and methods

### 9.1. Cohorts and the unit of analysis

Seven public datasets were used, in two groups. Table S8 is the registry. It gives, for
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
against the API on 2026-09-07 and are released in the repository's cohort registry.

**Dataset selection.** The five screen datasets were chosen before the screen was run, to span different acquisition structures rather than to sample every eligible cohort. They are a standardised single-assay atlas expected to show weak association (HIHA), a lupus cohort with processing waves (GSE174188), a three-site COVID-19 cohort (Stephenson), a multi-centre atlas with two sequencing chemistries (Ren), and a consortium cohort in which several phenotypes share one acquisition structure (COMBAT). Each had to provide human
peripheral blood, donor identifiers, a phenotype contrast within the dataset, and at least
200 cells per retained donor.

To show how large the pool was, we counted, after the analysis, the CELLxGENE Discover datasets published by 2026-09-07 that meet the screen's dataset-level requirements (API response retrieved 2026-09-13). A dataset qualifies if it is human, carries a blood tissue label and a single-cell suspension, has a "normal" label alongside at least one disease label, and has at least 40 donors. Of 2,226
datasets, 35 in 15 collections qualify, and the five analysed collections are among them
(Table S9). Dataset-level metadata do not give donor counts per class, so some of the
other ten collections may fail the per-class floors. The count of eight in nine therefore
describes a purposive selection covering a third of the eligible collections; it is not a
prevalence estimate.

The donor was the independent unit in every supervised analysis. Cells from one donor
never crossed a training/test partition.

### 9.2. Donor-level tables and the metadata variables

Each cohort was streamed from its h5ad in 100,000-cell chunks and reduced to one row per
donor. Expression was aggregated to log1p CPM pseudobulk. Raw counts were identified by checking which stored layer held integer values.

Metadata variables were drawn from standard schema fields (donor identifier, disease, sex, developmental stage, assay, tissue) plus cohort-specific collection
variables detected by pattern from the distributed object. Detection patterns cover
batch, pool, run, lane, chip, site, centre, city, region, province, hospital, clinic,
study, sub-study, source, dataset and donor source. Ages given as strings
("25-year-old stage", "fifth decade stage", "90 year-old and over stage") were parsed to
years.

Covariates that are consequences of the disease were excluded by an explicit blocklist
matching severity, outcome, comorbidity, medication, sample timing, symptom, diagnosis,
stage, grade and WHO score. Supplementary Table S1 lists, for every cohort, which columns
entered the metadata matrix, in which group, and which were excluded and why.

Variables were assigned to three groups. **Sample quality**: cells per donor, mean UMI
per cell, mean genes per cell, aggregate mitochondrial percentage, and their logs.
**Demographic**: age in years, sex, ethnicity. **Collection**: every detected batch-type variable, assay and suspension type.

### 9.3. Residual label variance

For centred phenotype-label vector *y* and recorded metadata matrix **D** including an intercept,
V_D = 1 − R²(y ~ D). Two versions are reported. The in-sample version uses ordinary least
squares on all donors and is reported for comparison. The
cross-fitted version replaces the fitted values with out-of-fold predictions from
shuffled five-fold cross-validation, seeded by the run seed, and is the version used
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
through the identical fixed-hyperparameter pipeline (Table S3); that is the primary
test, for the reason given in Section 2.2. The in-sample V_D is released with its own
in-sample null in Table S3. The two are separate quantities and are never mixed:
p(V_D) always refers to the cross-fitted pair. This null is what makes either version
readable, because its location depends on the number of metadata columns relative to
donors.

### 9.5. The fitted pipelines

**Four analyses fit classifiers, and they are not one pipeline.** The multi-cohort screen,
the decision-tree walkthrough, the calibration simulation and the GSE174188/GSE135779 deep
dive each have their own settings, and Table S10 lists all of them. Table S10 was generated directly from the analysis code, and the paragraphs below state what differs between the pipelines.

*Shared by the screen, the walkthrough and the calibration simulation.* One classifier specification is used, differing only in repeat budget. It is a balanced logistic regression (liblinear solver, at most 5,000 iterations) at fixed inverse regularisation 1.0 with stratified five-fold outer splits. The 4,000 highest-variance features are selected on the training donors, and a PCA to 50 components (or fewer when donors are limiting) is fitted on those donors and applied to the held-out ones. The screen runs one repeat, the walkthrough
twenty. A tuned variant selects the inverse regularisation from 10⁻⁴ to 10⁴ in decade steps by three-fold resampling inside each outer training set, over five repeats. It is reported beside the fixed pipeline as a description; the fixed pipeline is the inferential statistic, because both permutation nulls refit it.

*Metadata matrices.* The metadata matrix is small and dense, so no feature filter and no
PCA are applied to it. Numeric variables are median-imputed from the training donors of each fold and standardised there. Categorical variables are one-hot encoded from the levels observed in those training donors, dropping the first level; a level seen for the first time in a held-out donor is encoded as the reference level. The in-sample V_D is defined on the whole-cohort matrix, as an in-sample statistic must be.

*The lupus deep dive.* GSE174188 and GSE135779 predate the screen and run a wider budget on a single dataset. The classifier is a balanced logistic regression (liblinear solver, at most 20,000 iterations) with the same inverse-regularisation grid, selected by five-fold resampling inside each outer training set, over 20 repeated stratified five-fold donor splits with integer seeds 20260801-20260820. Out-of-fold predictions are pooled within each repeat; we report the
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

Both permutation tests and the observed statistic run one fixed-hyperparameter pipeline, refitted from the raw donor-by-gene matrix on every permutation. Within each training fold, the 4,000 highest-variance genes are selected on the training donors. A standardiser and a PCA to 50 components (or fewer when donors are limiting) are fitted on those donors and applied to the held-out ones. A balanced logistic regression at fixed inverse regularisation 1.0 is then fitted, with a single cross-fitting repeat. Every step that touches expression is therefore inside the fold, so
the reported AUC is genuinely out-of-fold and the null refits the whole pipeline rather
than reusing one representation. Freezing the regularisation is still necessary: comparing
a tuned observed statistic against an untuned null biases p-values downward.

Every representation step is refitted inside the fold on every permutation because the
code arranges it, not because the permutation test would reveal otherwise: a label-dependent
step taken before the permutation appears identically in the observed and permuted runs.
All analyses in this paper use one shared implementation of these fold-contained steps.

The **free** test permutes the phenotype vector without restriction. The **collection-stratified** test permutes it only within collection strata, preserving each stratum's observed class counts. A stratum is the interaction of every detected batch-type variable with assay, where assay varies, and Table S3 lists the stratum variables for every comparison. The stratified scheme is the restricted permutation of
Chaibub Neto and colleagues [30]; it is applied here, not introduced here.

**The conditioning set is collection variables only.** Sample-quality and demographic
columns are in the metadata matrix but not in the strata, so this test conditions on a
subset of what is recorded. The free test's null is that
the phenotype label is exchangeable across all donors; the stratified test's is that it is
exchangeable within each collection stratum, that is, conditional independence of
expression and phenotype given those strata. Section 4.4 shows what follows: when a
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
(Table S3) and is deliberately **not** used as a feasibility
rule. A gap between two different algorithms does not establish that a comparison is
unanswerable. Every inferential statement in this paper is made about the
fixed-hyperparameter pipeline, which is the statistic both permutation nulls are built
from; the tuned AUC is descriptive. "Fixed-hyperparameter" refers to the statistical
pipeline throughout, and "frozen" only to a foundation model whose weights are not
updated.

### 9.8. Seed stability

All nine comparisons were re-run at five seeds (20260907-20260911), varying the seed for
fold construction, permutation draws, and the cross-fitting split used for V_D and its own
permutation null. All five per-seed outputs are released rather than summarised (Figure S7, Table S11), and
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
because together they separate the effect of prevalence from that of threshold
alignment (Section 4.3).

Critically, the external target reuses the **source's** biological loading vector.
Regenerating it makes every transfer chance-level by construction and tests nothing;
only the design mechanism differs between source and target in the shift arm.

A supplementary arm quantifies what admitting a disease-caused covariate does. A severity
variable was generated as 1.2·(y − ȳ) + noise and added to the metadata matrix. Everything
else was held fixed.

The extended simulation's cohorts are 200 donors by 100 features, so it applies no feature filter and no PCA. Its classifier settings are listed in Table S10. They are not derived from the screen's, because a 4,000-feature filter and a PCA to 50 would be a different reduction on a 100-column matrix. Its unadjusted, residualised and restricted arms differ from each other in exactly one operation.

**Calibration arms.** Two further studies use the screen's fixed-hyperparameter specification. Arm A places 200 donors in eight collection
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
whether the training-fold ridge fit on the metadata matrix is subtracted first. Sixty replicates were run per cell, each seeded from its own coordinates.

### 9.10. Decision-tree walkthrough

Three cohorts were traced through the tree with all evidence recomputed, and every branch
value is released in Table S7. The collection strata, the metadata matrix and the
classifier are imported from the screen rather than reimplemented, so Q2 and Q4 refer to
the same partition, the same columns and the same fitted pipeline; the walkthrough uses the screen's fixed-hyperparameter specification with twenty repeats. The metadata matrix is built and
standardised on the training donors of each split, because ridge shrinkage depends on scale.

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
through the GEO title/accession map; all 44 donors mapped uniquely. The restored donor table and its batch-by-label, year-by-label and batch-by-year cross-tabulations are released with checksums of every source file. For GSE174188, the per-cell processing-cohort field was read from the distributed CELLxGENE object and reduced to donor
fractions over four waves. The dominant wave is the maximum-fraction wave; a pure-wave
donor has at least 99% of analysed cells in one wave.

### 9.12. Representations

**Frozen Geneformer.** Cell embeddings used Geneformer V2-316M at repository revision 04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5 [4]; checkpoint, configuration, token
dictionary and gene-median dictionary checksums are released with the code. Gene
identifiers were version-trimmed and mapped to the V2 vocabulary. Within each cell the
ranking value was raw count / total count × 10,000 / Genecorpus-104M gene median, and
genes were sorted descending. Sequences used V2 start/end token IDs 2 and 3, padding ID
0 and maximum length 4096, leaving at most 4094 ranked genes. The model output was the
final hidden layer; token position 0 was kept as the 1,152-dimensional cell
vector. The encoder was frozen. Donor representations are coordinate-wise means of
sampled cell vectors, sampled without replacement with integer seed 1, capped at 500
cells per donor for GSE135779 and 1,000 for GSE174188 and GSE285773. The embedding
environment used Python 3.10, PyTorch 2.5.1 and Transformers 4.46.3.

**Pseudobulk.** We use pseudobulk rather than a learned latent space as the donor-level
representation, so that the comparison between representations is not itself mediated by
a model fitted with batch covariates [5]. Stored donor-by-gene values are log1p(CPM). Inside every outer
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
encoded as in Section 9.5 and a multivariate ridge regression (penalty 1, intercept
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
population and 261 donors as the molecular analysis. Per-donor counts were formed for the naive (T4_naive), effector-memory (T4_em) and regulatory (T4_reg) labels and for all other or unassigned labels. We evaluated raw
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
scikit-learn 1.9.0, PyArrow and joblib; exact per-run versions are recorded with each
output. Figures were produced in R 4.6.0 with ggplot2 4.0.3 and patchwork 1.3.2. The
screen, the extended simulation, the permutation calibration and the decision-tree
walkthrough were run on a 16-vCPU ARM container; all other analyses were run locally.

Random seeds are derived deterministically from a recorded master seed, and every derived
seed is written into the released result tables, together with input checksums, parameter
grids, donor-level predictions and complete null distributions. All tables and figures are
regenerated from the released result tables by the released code, and every number quoted
in the text was checked against those tables.

### 9.18. Use of generative AI

Generative AI tools were used during this work: Claude Code (Anthropic) with the Claude
Opus 5 model, and Codex (OpenAI) with the GPT 6 Astra model. They were used for:

- **Code**: implementation and review of analysis and figure scripts written to
  statistical definitions, gates and thresholds specified by the authors, and
  identification of implementation defects that the authors then confirmed by inspection
  and by rerunning the affected analyses.
- **Manuscript preparation**: language editing, reorganisation of sections, and
  consistency checking of cross-references, figure and table numbering, and reference
  numbering.

Cohort selection, the statistical definitions, the choice of analyses to report and the
interpretation of results were made by the authors. No data were generated by these tools;
every reported value comes from the analysis outputs. The authors reviewed all AI-assisted
code and text, verified it against the underlying results, and take full responsibility
for the content of this article.

## 10. Data and code availability

All cohorts are public. The lupus datasets are at GEO accessions GSE135779, GSE174188
and GSE285773; the COVID-19, influenza and CMV cohorts were obtained through the
CELLxGENE Discover curation API and their dataset identifiers are listed in Table S8 and in the repository's cohort registry. The MIT-licensed repository is at
https://github.com/LightChainr/rheumlens. Its releases are archived on Zenodo under the
concept DOI https://doi.org/10.5281/zenodo.20813922, which resolves to the latest version;
this manuscript corresponds to release v3.0.0. The versioned research object accompanying
this manuscript contains the cohort registry with API-verified donor counts, the
donor-level interface tables, every result table underlying Figures 2 to 8, all five
per-seed outputs rather than a summary, the 19,600-cohort simulation output, analysis scripts, locked software environments, reproduction instructions and file checksums. Public raw single-cell
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
A recruitment process R places each donor under the recorded metadata D and is also
associated with the phenotype label Y; the measured cells X are affected by D (technical)
and by Y (biological); a classifier maps X to Ŷ. C denotes covariates caused by the
disease - severity, medication, symptom timing - which lie on Y → C → X and are
excluded from D. U denotes collection structure that was never recorded
and which these checks cannot detect. The dotted line marks the association the
checks quantify; they measure it but do not orient it.

**Figure 2. Nine phenotype contrasts.**
**(A)** Comparisons ordered by cross-validated metadata-only AUC (filled, coloured by
disease) with the corresponding expression-classifier AUC (open) on the same axis. Bars
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

**Figure 3. Simulation across seven data-generating regimes.**
19,600 simulated cohorts (200 donors, 100 features): seven regimes × seven levels of
metadata-phenotype association ρ × two label types × 200 replicates.
**(A)** Fraction of runs with p(V_D) ≤ 0.05, binary labels. At ρ = 0 this is 0.025 to
0.065 across regimes, bracketing the nominal 0.05.
**(B)** Unadjusted, external, restricted and residualised AUC against ρ in the additive
linear regime. Restriction is undefined at ρ = 1 because no stratum holds both labels.
**(C)** Restricted minus residualised AUC for all seven regimes. The gap opens in every
one. It is blunted in the imbalanced-offset regime, where the collection boundary and
the phenotype threshold sit at different quantiles, and not in the imbalanced regime, where the
case rate is the same 20% but the two are matched (Section 4.3).
**(D)** Internal (solid) and external (dashed) AUC when the target shares the source's
collection mechanism versus when it does not. External AUC rises with confounding in
the shared case.

**Figure 4. Collection strata and disease status in the two lupus cohorts.**
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
0.272-0.464; Figure S4). Residualisation loss is a function of the shape of the adjustment
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

**Figure S1. Admitting one disease-caused covariate produces a spurious metadata association.**
Metadata-only AUC, V_D, the fraction of runs flagged significant, and residualisation
loss, with and without a severity variable added to the metadata matrix, across five
levels of true metadata-phenotype association. At ρ = 0, where metadata and phenotype are
unrelated by construction, admitting the covariate raises metadata-only AUC from 0.488 to
0.793 and
flags the cohort in 100% of runs against 5.5% without it.

**Figure S2. Collection-stratum composition of each cohort.**
Case and control counts by collection stratum for all nine comparisons, showing which
strata are label-pure and how many donors sit in mixed strata.

**Figure S3. Rejection rate of both permutation tests under a non-disease within-stratum
association.**
**(A)** Fraction of runs rejecting at p ≤ 0.05 for both permutation tests, as a function
of how strongly a sample-quality variable predicts the phenotype *within* a collection
stratum. The generating model contains no disease effect at all. Only the leftmost point
(γ = 0) measures calibration; to its right the conditional-exchangeability null is false
by construction, because the sample-quality variable lies outside the conditioning set,
so the rate is the probability of detecting a non-disease association rather than a
type-I error rate. **(B)** The mean absolute within-stratum correlation that conditioning
on the collection strata leaves in place.

**Figure S4. Residualisation loss with no metadata-phenotype association present.**
AUC lost to fold-contained ridge residualisation in simulated cohorts where the metadata
matrix is independent of the phenotype by construction, against the width of that matrix,
at 108 donors. The unadjusted and residualised arms share folds, feature filter, standardiser, PCA and classifier and differ
only in the adjustment. The biological effect is calibrated so the unadjusted AUC lands
near 0.90. The ribbon is the 2.5th-97.5th percentile over 60 draws. The triangle reuses
the level sizes of the real CMV batch-by-pool crossing, which is far more ragged than
equal blocks; it loses 0.363 on average. The dashed line marks the 0.326 observed in that
comparison, which falls inside the triangle's interval.

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
Per-seed values of metadata-only AUC and its permutation p, V_D and p(V_D), expression-classifier AUC
and both permutation p-values, for all nine comparisons at five seeds. The p-value panels
use a log axis; 0.005 is the smallest value 200 permutations can return.

**Table S1.** Metadata manifest: for every comparison and every variable group, the raw columns consumed and the
number of expanded columns produced. **Table S1b** lists every column present in a
covariate file that enters no variable group, with the reason. **Table S1c** confirms that the
manifest reproduces the width of every fitted model matrix (34 of 34 variable groups agree).
**Table S2.** Simulation summary: all seven regimes × seven ρ × two label types, plus the
disease-caused-covariate arm.
**Table S3.** Full screen output: every comparison × every variable group, with
metadata-only AUC (linear and random forest), V_D in-sample and cross-fitted, null mean,
p(V_D), expression-classifier AUC, both permutation p-values and stratum counts.
**Table S4.** Metadata-only, expression-only and joint AUC for every comparison at each
of five split seeds, with the two increments. All three arms use the same
fixed-hyperparameter specification and the same folds; expression is reduced to principal
components fitted on the training donors before the metadata columns are appended.
**Table S5.** Calibration study: summary output of both arms, the permutation rejection rates of Section 4.4 and the residualisation-width study of Section 4.5.
**Table S6.** Learned pooling: per-method target metrics and all paired DeLong tests,
both transfer directions.
**Table S7.** Decision-tree walkthrough: all evidence for every branch, three cohorts.
**Table S8.** Cohort registry, compiled from the dataset records and the donor tables:
for all seven datasets, the identifier, citation, the donor count reported by the source,
the download size, the case and control labels, the donor, case and control counts
actually modelled in each comparison, and how and when each entry was verified.
**Table S9.** Retrospective eligibility audit of CELLxGENE Discover: every collection whose datasets meet the screen's dataset-level requirements, with donor counts, disease labels, assays and whether it was analysed. Performed after the analysis and not used for selection.
**Table S10.** Every fitted pipeline in the paper and its settings. The four analyses are separate
pipelines and are listed separately.
**Table S11.** Seed-stability table, all five seeds, all nine comparisons.

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
