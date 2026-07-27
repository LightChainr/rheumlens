# When study design predicts disease: an identifiability limit for patient-level single-cell classifiers

Hongyu Ying<sup>1</sup>, Dandan Yun<sup>1</sup> and Dan Liu<sup>1,*</sup>

<sup>1</sup> Department of Rheumatology and Immunology, Shanghai Pudong
Hospital, Fudan University Pudong Medical Center, Shanghai 201399, China  
<sup>*</sup> Correspondence: Dan Liu; danliu600@126.com

## Abstract

Single-cell RNA sequencing resolves disease-associated immune states at cellular
resolution, but clinical classification concerns patients. Patient-level models
compress thousands of cells into one donor representation and can report
cross-validated discrimination above 0.95. Such performance is ambiguous when
recruitment, processing or sequencing design also predicts the disease label.

We define a design-label estimability scale and test it in three systemic lupus
erythematosus cohorts. In an additive model, the fraction of label information
remaining after design adjustment is \(1-R^2(Y\sim D)\). Across 11,200 simulations,
increasing design-label association depleted this information, destabilised
biological-effect recovery and created internal-to-external gaps. A
low-entanglement positive control passed all five proposed validity checks.

The empirical cohorts occupied the failure regime. In GSE135779, recorded design,
sequencing and demographic variables predicted disease at AUC 0.952 and exceeded
1,000 complete-pipeline label permutations (\(p=0.000999\)). Nonlinear design-only
models also discriminated disease (AUC 0.815-0.824). Yet batch provided an internal
negative control: representations identified batch at AUC up to 0.993, while batch
alone predicted disease at 0.499 and batch residualisation changed disease AUC by at
most 0.009. In GSE174188, design restriction preserved most discrimination whereas
linear residualisation removed 0.230-0.268 AUC. Random-forest residualisation in
GSE135779 also attenuated all representations, while batch location adjustment and
overlap weighting largely preserved them. Strict source-only transfer was weaker
than internal validation; pseudobulk remained strongest in the 261-donor target and
learned pooling did not repair the gap.

We convert these results into five executable checks spanning estimability,
full-pipeline nulls, design exposure, overlap restriction and external transfer.
High internal AUC establishes separability; attribution requires design overlap, and
transport requires an independent acquisition process.

**Keywords:** patient-level classification; single-cell RNA sequencing;
identifiability; study design; batch confounding; residualisation; external
validation; Geneformer; pseudobulk; systemic lupus erythematosus

## 1. Introduction

Systemic lupus erythematosus is a heterogeneous autoimmune disease in which blood
interferon activity and immune-cell states vary substantially between patients.
Single-cell RNA sequencing can resolve those states, but diagnosis concerns a
patient while the measurements are made in individual cells. Patient-level modelling
must therefore compress a variable-sized cell population into one representation
without confusing biological signal with the way the cohort was assembled.

Several solutions now compete for this role. Pseudobulk aggregates expression,
CloudPred and singleDeep learn phenotype-oriented summaries, and recent methods such
as PaSCient, scPhase and FloREN learn multicellular or graph-based patient
representations [6-11]. Frozen single-cell foundation models offer a particularly
attractive shortcut: encode each cell once, average the embeddings, and train a small
donor classifier. This strategy is computationally convenient, but high within-cohort
accuracy does not establish that the representation will survive a change of
hospital, population, acquisition period or sequencing workflow.

The hidden difficulty is that many public disease cohorts are not random samples from
one acquisition process. Cases and controls may come from different recruitment
sources, years, multiplexing pools or processing waves. Batch-label confounding has
long been known to bias cross-validation [16], and a completely confounded
single-cell design is statistically non-identifiable [17]. This is not merely a
preprocessing problem. If disease and design are collinear, the observed data cannot
distinguish a biological disease effect from a technical design effect without
additional overlap or external information. Regressing a representation on technical
covariates does not restore that missing contrast: it removes the disease-aligned
component whether that component is technical, biological or both [18,19].

This problem has become more urgent, not less. Broad benchmarks show that frozen
single-cell foundation-model embeddings can retain batch structure and can lose to
simple expression baselines [13,14]. Patient-level methods already report high
internal performance in SLE, including external validation between GSE174188 and
GSE135779 [7]. FloREN now evaluates both cohorts with supervised graph
representations [11], while pretraining-contamination work argues that benchmark
validity must accompany model rank [15]. What remains missing is a patient-level
framework that connects the estimability limit to executable diagnostics, then tests
whether design predictability, residualisation, design restriction and external
transfer tell the same story on the same donors. They do not.

Predictive-modelling research has already shown that confound adjustment must be
contained within training folds, that featurewise confound regression can leak or
reveal information to nonlinear learners, and that residual confounding can persist
after apparently successful mitigation [27-29]. Permutation-based confounding tests
and classifier-dependent indices provide complementary diagnostics [25,30]. These
results motivate caution, but they do not quantify how much disease-label contrast
remains estimable from a patient-level single-cell cohort, nor compare adjustment
with observed-overlap restriction and source-only transfer. A recent unsupervised
single-cell benchmark further established cell-type composition as a strong
patient-stratification baseline [26], making composition an essential rather than
optional control.

Here we define a design-adjusted label-information measure, construct observational
equivalence at perfect design-label collinearity and map the resulting
failure regime in simulation. We then recover donor design metadata for GSE135779,
combine it with cell-distributed processing metadata in GSE174188, and test frozen
Geneformer plus two expression pseudobulk representations. Design-only prediction,
complete-pipeline permutation, representation-to-design prediction, fold-contained
residualisation, overlap restriction, matched donor controls, cross-batch transfer,
cross-cohort transfer and learned pooling are evaluated under one donor-matched
framework. This moves beyond another representation leaderboard: it defines when a
patient-level disease signal can be separated, when it can be attributed and when it
can be expected to transport.

## 2. Results

### 2.1. Study design predicts disease in two widely reused SLE cohorts

GSE174188 contains 261 analysed CD4-positive donors: 162 SLE cases and 99 controls
[2]. The distributed single-cell object records four processing waves. In the
CD4 subset, their control/case counts were 32/0, 17/97, 7/19 and 43/46
(Figure 1A). Wave 1 contains no cases, and wave 2 is dominated by cases. A logistic
model using only donor processing-wave fractions achieved AUC 0.922
(2.5th-97.5th percentile over 20 repeated five-fold donor analyses:
0.916-0.928). Sequencing/QC variables reached 0.781, and the complete measured
design block reached 0.953 (Table 1; Figure 1D).

For GSE135779, we restored batch, collection year, demographics, clinical fields and
sequencing QC for all 44 childhood-cohort donors from the original Supplementary
Tables 1b and 1c and mapped them to the analysis identifiers through GEO [1].
No donor or selected field was missing. The six sequencing batches contained
case/control counts of 8/0, 3/1, 6/2, 6/2, 6/2 and 4/4 (Figure 1B). Collection
years 2014 and 2016 contained only cases, whereas 2018 was balanced at 8/8
(Figure 1C). Collection year alone predicted disease at AUC 0.721, demographics at
0.715 and sequencing/QC at 0.707. Their combination with batch reached AUC 0.952
(0.893-0.992).

These two cohorts encode design-label entanglement differently. GSE174188 has one
dominant processing coordinate that is itself a strong disease classifier.
GSE135779 distributes the label channel across year, sequencing and demographic
variables. A batch table alone would expose the first design but miss the second.

**Table 1. Cohort design structure and design-only discrimination.**

| Cohort | Donors (case/control) | Recorded design | Complete-design AUC | Most informative single block |
|---|---:|---|---:|---|
| GSE135779 | 44 (33/11) | 6 batches; 4 collection years; sequencing QC; demographics | 0.952 [0.893-0.992] | collection year: 0.721 |
| GSE174188 CD4 | 261 (162/99) | 4 processing-wave fractions; sequencing QC; age; sex | 0.953 [0.945-0.961] | processing wave: 0.922 |
| GSE285773 CD4 | 26 (16/10) | donor labels and cell profiles; no recoverable donor processing batch | not estimated | not estimated |

### 2.2. A quantitative estimability scale for design-label association

Let \(Y\) denote the centred donor disease label, \(D\) the recorded design matrix
including an intercept, and \(M_D\) the projection onto the orthogonal complement of
the column space of \(D\). The fraction of label information remaining after linear
design adjustment is

\[
\mathcal{I}_D=\frac{Y^\top M_DY}{Y^\top M_1Y}=1-R^2(Y\sim D).
\]

This scale is one when design is unrelated to the label and zero when design
perfectly spans it. In the additive model, relative variance for estimating a
disease coefficient from the remaining contrast increases as
\(1/\mathcal{I}_D\). At perfect collinearity, a biological-only mechanism
\(X=\beta Y+\epsilon\) and a technical-only mechanism
\(X=\delta D+\epsilon\) can produce exactly the same observed \((X,D,Y)\)
distribution. The construction is a direct consequence of collinearity; its value
here is to turn that familiar fact into a cohort-level scale and a testable
patient-representation protocol.

We mapped this boundary in 11,200 donor-level simulations crossing four biological
effects, four technical effects and seven levels of design-label association
(Figure 2). Mean \(\mathcal{I}_D\) fell monotonically from 0.996 at zero
association to 0.189 at \(\rho=0.90\), 0.039 at \(\rho=0.98\) and zero at perfect
collinearity. Median variance inflation reached 5.27 at \(\rho=0.90\) and 30.18 at
\(\rho=0.98\). With a purely technical signal
(\(\beta=0,\delta=2,\rho=0.90\)), internal unadjusted AUC was 0.936, while
residualised, overlap-restricted and independent-target AUCs were 0.499, 0.501 and
0.499. With mixed biological and technical signal
(\(\beta=1,\delta=1,\rho=0.90\)), internal AUC was 0.947, residualisation reduced it
to 0.573, overlap restriction retained 0.794 and independent transfer retained
0.814. Residualisation therefore measured removal of design-aligned variation, not
the technical share of the original classifier.

The same simulation provided a positive control. At lower design-label association
(\(\beta=1,\delta=0.5,\rho=0.25\)), mean information fraction was 0.933,
design-only AUC was 0.626, and the largest loss under residualisation or restriction
was 0.028. External AUC was not lower than internal AUC, and a complete-pipeline
permutation test gave \(p=0.003984\). All five checks therefore passed in a setting
where biological signal was recoverable and portable.

### 2.3. The GSE135779 design signal survives complete-pipeline null tests

The high complete-design AUC in GSE135779 could, in principle, reflect a flexible
pipeline applied to only 44 donors. We therefore repeated the full nested
imputation, encoding, scaling, regularisation selection and outer prediction
pipeline under 1,000 label permutations per design block. The observed complete
design AUC of 0.952 exceeded every permutation except the finite-sample correction
required for an empirical upper-tail \(p=0.000999\) (Figure 3A). An independent
null that shuffled each design column separately while preserving its marginal
distribution gave \(p=0.003984\) over 250 repetitions. Collection year, sequencing
QC and demographics individually produced weaker evidence
(\(p=0.017\), 0.040 and 0.032 under label permutation); batch alone was null
(AUC 0.499, \(p=0.470\)). The combined label channel was distributed across
variables rather than carried by one categorical batch.

Model flexibility did not explain the result. Random-forest and gradient-boosting
design-only classifiers achieved mean AUC 0.824 and 0.815, respectively, under the
same repeated outer donor splits. Both exceeded chance but remained below the
regularised logistic model, indicating that the strongest recorded design-label
structure was predominantly low-dimensional rather than an artefact of a flexible
nonlinear learner.

The information calculation reached the same conclusion from a different direction.
Batch left 87.9% of centred label variation available after adjustment. The complete
21-column design matrix left only 20.3%, corresponding to a relative variance
inflation of 4.93 before any molecular feature was considered (Figure 3B).
Thirty-six donors remained in mixed-label sequencing batches and 32 in mixed-label
collection years, so overlap-based analyses were still possible; adjustment of the
complete design nevertheless approached the low-information regime identified by
simulation.

### 2.4. The representations remember the design more accurately than the disease

We next asked whether the donor representations retained the recorded design
identity. Under the same fold-contained pipeline, each representation predicted
case/control status and every batch one-versus-rest (Figure 4).

In GSE174188, Geneformer, HVG pseudobulk and PCA pseudobulk classified disease at
AUC 0.982, 0.977 and 0.975. Their best processing-wave AUCs were 0.9997, 0.999
and 0.997. In GSE135779, disease AUCs were 0.870, 0.933 and 0.945, while the
best individual-batch AUCs were 0.993, 0.961 and 0.961. Near-perfect batch
recognition was therefore not a peculiarity of one representation family or one
cohort.

This result establishes exposure, not attribution. A representation that predicts
batch can still contain disease information, and a batch-aware classifier need not
use batch as its disease decision rule. The distinction is resolved by restriction,
not by correlation alone.

### 2.5. An internal negative control separates design exposure from shortcut use

GSE135779 supplied a direct counterexample to the common inference that “batch is
predictable from the representation” means “the disease classifier uses batch.”
The best one-versus-rest batch AUC was 0.993 for Geneformer and 0.961 for both
pseudobulk representations. Yet batch alone predicted disease at AUC 0.499, and
fold-contained batch residualisation changed disease AUC by only +0.002, -0.002 and
-0.009. Thus representation-to-design predictability documents exposure, but
dependence of the disease decision must be demonstrated separately.

### 2.6. Cell-type composition reproduces the design-adjustment contradiction

We next asked whether the identifiability problem extended beyond high-dimensional
molecular representations. For every GSE174188 donor, we aggregated the same
CD4-positive cells into naive, effector-memory, regulatory and
other-or-unassigned fractions. Raw proportions, centred-log-ratio composition and
each representation augmented with log cell yield were evaluated with the locked
20-repeat donor pipeline.

Composition alone was a strong disease classifier. Mean AUC was 0.899 for raw
proportions, 0.903 for centred-log-ratio composition and 0.915 after adding cell
yield (Figure 3C). Yet each composition vector also predicted processing-wave
identity: one-versus-rest wave AUCs reached 0.803-0.872 for wave 1 and
0.817-0.847 for wave 4. Fold-contained processing-wave residualisation reduced all
four disease models to AUC 0.689-0.709.

Restriction again gave a different answer. Within dominant wave 4, observed
composition AUC was only 0.027-0.062 below size- and label-matched donor subsets;
within pure wave 4, the difference was 0.034-0.037 (Figure 3D). Every
2.5th-97.5th percentile interval crossed zero. Thus a simple biological aggregate
reproduced the same empirical pattern as Geneformer and expression pseudobulk:
design was readily recoverable, residualisation removed a large disease-aligned
component, and observed overlap preserved substantially more discrimination than
residualisation implied.

### 2.7. Disease discrimination survives design restriction

We first restricted GSE174188 to dominant wave 4, the only large near-balanced wave
(n=89; 46 cases and 43 controls). Because reducing 261 donors to 89 inevitably lowers
performance, we compared the observed stratum with 20 random donor subsets having
the same donor count and case/control composition (Figure 5A). Relative to those
matched subsets, AUC changed by -0.014 for Geneformer, -0.028 for HVG pseudobulk and
-0.049 for PCA pseudobulk.

The stricter analysis retained only donors with at least 99% of analysed cells in
wave 4 (n=66; 25 cases and 41 controls). Relative to size- and label-matched random
subsets, AUC changed by -0.026, -0.044 and -0.044 (Figure 5B). These are detectable
losses, but they are far smaller than the raw full-cohort-to-stratum difference. Most
of that raw difference is explained by using fewer donors.

GSE135779 provided an independent design. Removing the all-case B1 batch left 36
donors from batches B2-B6, with both labels represented in every retained batch.
Geneformer AUC was 0.846 versus 0.858 in matched random subsets (difference -0.012).
Both pseudobulk representations were 0.919 versus 0.916 (difference +0.003;
Figure 5C). Thus the disease signal persisted after the most obvious label-exclusive
batch was removed.

Across both cohorts, design restriction weakens the claim that internal AUC is a pure
technical shortcut. The representations carry design information, but the measured
disease separation does not disappear when the clearest design channels are held
fixed.

### 2.8. Residualisation can erase the signal that restriction preserves

Residualisation produced a sharply different conclusion. All unadjusted and adjusted
models were recomputed on identical outer folds with the same training-only HVG
selection, scaling, PCA construction, regularisation grid and classifier.

In GSE174188, ridge residualisation on processing-wave fractions reduced Geneformer
AUC from 0.982 to 0.714, HVG pseudobulk from 0.977 to 0.747 and PCA pseudobulk from
0.975 to 0.735 (losses 0.268, 0.230 and 0.241; Figure 6A). The same design, tested by
dominant- or pure-wave restriction with matched donor controls, cost only
0.014-0.049. Residualisation therefore removed five to nineteen times more AUC than
the matched restriction discrepancy, depending on representation and stratum. The
paired split-bootstrap difference between residualisation loss and matched
restriction loss was 0.186-0.255 across the six comparisons; every descriptive
2.5th-97.5th percentile interval remained above zero (Figure 6F).

GSE135779 revealed why. Batch alone had a cross-validated disease AUC of 0.499, and
batch residualisation changed representation AUC by +0.002, -0.002 and -0.009
(Figure 6B). Collection year had AUC 0.721 and its removal cost 0.005-0.039.
The complete measured design block had AUC 0.952; projecting it out collapsed
Geneformer, HVG pseudobulk and PCA pseudobulk to 0.516, 0.526 and 0.515, removing
0.354, 0.407 and 0.430 AUC (Figure 6C-D). Yet exclusion of the all-case batch had
changed matched performance by at most 0.012.

The conclusion was not specific to a linear residualiser. A training-fold random
forest fitted from the complete design to each representation retained AUC 0.748,
0.787 and 0.795, corresponding to losses of 0.122-0.150 (Figure 6E). In contrast,
training-only batch location adjustment retained AUC 0.864-0.911, and full-design
overlap weighting retained 0.859-0.937. The method matters quantitatively, but every
full-design feature-removal arm attenuated discrimination much more than the
batch-only negative control.

The projected covariate block predicts how much discrimination residualisation will
erase. This is expected algebraically: when design predicts the label, the component
of a representation predictable from design includes disease-aligned variation.
Residualisation cannot determine whether that variation is technical, biological or
both. It answers how much discrimination survives deleting design-aligned signal,
not how much of the original discrimination was caused by technical batch.

This empirical discrepancy is the paper's central finding. A large
post-residualisation loss is not evidence that the original AUC was spurious. Without
overlap-stratum evaluation, it is equally compatible with successful removal of a
shortcut and destructive removal of outcome-aligned biology.

In a separate synthetic sensitivity analysis, linear confound regression followed by
a random-forest classifier did not inflate null performance in our Gaussian
construction: fold-contained residualisation gave mean AUC 0.500 and global
residualisation 0.438. This does not rule out confound leakage under skewed or
limited-precision features [28]; it shows that the empirical attenuation reported
here was not accompanied by that failure mode in the tested construction.

### 2.9. Transfer degrades across processing waves more than across cohorts

Design restriction preserved strong within-cohort discrimination, but transfer
weakened it at every level tested (Figure 7). Training on processing waves 2+3 and
testing on wave 4 in GSE174188 produced AUC 0.842, 0.811 and 0.806 for Geneformer,
HVG pseudobulk and PCA pseudobulk. Reversing the direction produced 0.919, 0.926 and
0.821. Wave 1 could not be evaluated as a target because it contained no cases.

Strict source-only transfer from the 26-donor GSE285773 cohort to the independent
261-donor GSE174188 target yielded AUC 0.884 [0.843-0.920] for Geneformer,
0.919 [0.884-0.946] for HVG pseudobulk and 0.926 [0.895-0.953] for PCA
pseudobulk. On the same 261 target donors, both pseudobulk methods exceeded
Geneformer after paired DeLong testing and Benjamini-Hochberg adjustment
(q=0.0151 and q=0.000603). The reverse direction retained the same ranking but was
underpowered at 26 target donors.

Restricting the GSE174188 source to dominant wave 4 did not improve transfer to
GSE285773. Target AUC fell to 0.794, 0.894 and 0.888 and was no better than
size- and label-matched random source subsets. More varied source donors transferred
better than the apparently cleaner restricted source.

The cross-cohort result has a direct translational interpretation without implying
clinical readiness: a molecular classifier is useful beyond its development cohort
only if patient ranking survives changes in age distribution, recruitment and
technical context. Internal AUC did not predict that property here.

### 2.10. Learned pooling does not repair the transfer gap

Mean pooling might suppress informative cell heterogeneity, so we tested DeepSets,
gated-attention multiple-instance learning and pooling by multi-head attention
[20-22]. All operators used the same frozen cell embeddings, cap-500 cells per donor,
a source-fitted whitened PCA32 projection and source-only hyperparameter selection.
Independent target data were never used to select a model (Figure 8).

In the 261-donor target, none of the learned operators beat ordinary mean pooling.
Relative to the 1,152-dimensional Geneformer donor mean, their paired AUC differences
were -0.133, -0.055 and -0.042. Relative to a donor mean in the identical PCA32
space, the differences were -0.077, +0.001 and +0.014; none of the positive
differences was significant. All learned operators remained below both pseudobulk
representations. In the 26-donor target, no learned-pooling comparison was
significant.

This is a useful negative method result. The external gap was not repaired by adding
permutation-invariant capacity to the frozen embedding. A more expressive pooling
function can fit a source cohort without creating a transportable patient
representation.

### 2.11. An executable validity standard

The results support five minimum requirements for patient-level single-cell
classification:

1. **Quantify design-label estimability.** Report design-only disease AUC,
   \(1-R^2(Y\sim D)\), design rank and overlap counts before interpreting molecular
   discrimination.
2. **Test the complete pipeline under null data.** Label permutation must repeat
   imputation, encoding, feature handling, hyperparameter selection and outer
   prediction rather than test only a final fixed score.
3. **Report representation-to-design predictability.** A representation that retains
   design identity at near-perfect AUC has documented exposure to cohort structure,
   even if attribution remains unresolved.
4. **Pair residualisation with design restriction and matched donor controls.**
   Residualisation alone is uninterpretable under label-design collinearity.
5. **Require a source-only external target.** Target labels must not choose feature
   transformations, and internal performance must be reported beside external
   performance.

These checks are inexpensive relative to embedding generation or deep pooling. They
convert “the model achieved high AUC” into five separable statements: the design-label
contrast is estimable, the pipeline exceeds its null, the design is visible, the
signal survives observed overlap, and the ranking transfers.

## 3. Discussion

This study identifies an estimability boundary that representation leaderboards do
not measure. When study design predicts the disease label, the biological contrast
available after design adjustment shrinks as \(1-R^2(Y\sim D)\); at perfect
collinearity, biological and technical explanations are observationally equivalent.
No increase in model capacity can recover information absent from the acquisition
design. The simulations establish this result across controlled mechanisms, while
the two SLE cohorts show that commonly reused patient-level benchmarks can occupy the
low-information regime.

The decisive empirical comparison is between two adjustment strategies. Restriction
asks whether discrimination survives among donors observed under the same or
overlapping design. Residualisation deletes coordinates predictable from a chosen
design block. When design and disease are collinear, these operations target
different quantities. The two-cohort result makes the distinction visible:
mixed-batch restriction preserved GSE135779 performance, while residualising the
complete design block reduced every representation to chance. In GSE174188, wave
restriction caused modest matched losses, while wave residualisation removed roughly
one quarter of the AUC scale. Random-forest residualisation attenuated less than
ridge but preserved the same ordering, whereas overlap weighting and batch location
adjustment retained substantially more discrimination. CD4-subtype composition
independently reproduced the discrepancy, showing that the problem is not peculiar
to foundation-model embeddings.

The result is more consequential than a contest between Geneformer and pseudobulk.
Recent benchmarks already show that simple expression baselines can outperform
frozen foundation-model embeddings [13,14], and patient-level representation learning
is an active field [8-11]. Our contribution is to formalise and measure why an
internal comparison can remain scientifically ambiguous even when donor leakage is
prevented, folds are strict and adjustment is performed inside each fold. The
design-only permutation test confirms that the empirical signal is not an artefact
of a small-sample modelling pipeline; the information fraction measures how little
outcome contrast remains after adjustment. Correct computation does not repair a
non-identifiable design.

The external results identify what can be concluded. Strong discrimination persisted
inside design-restricted strata, supporting a reproducible disease-aligned signal.
The same representations lost performance across processing waves and cohorts,
showing that the signal was not sufficiently invariant. Pseudobulk converted
additional source donors into transferable information more effectively than frozen
Geneformer in the evaluated CD4 setting, and learned pooling did not reverse the
ordering. These are not contradictory findings: internal separability can be real
while external generalisation remains inadequate.

For clinical translation, this reframes representation selection as a cross-centre
validity problem. A classifier that ranks patients inside one retrospective cohort
but changes under a new acquisition design is not ready for probability
interpretation or individual decisions. The present external AUCs quantify ranking,
not calibrated risk. The immediate deliverable is therefore not a diagnostic model
but a validation protocol that future diagnostic models should pass.

## 4. Limitations

This study is confined to three public SLE cohorts, frozen Geneformer final-layer CLS
embeddings, fixed expression baselines and a small set of learned pooling operators.
The estimability scale is derived for an additive model and does not claim that AUC
must change monotonically under every data-generating process. The nonlinear
sensitivity arm used one tree-ensemble residualiser and two design-only classifiers;
it is not an exhaustive comparison of adjustment algorithms. The simulation
isolates one binary design variable, whereas real cohorts can contain interacting,
continuous and unrecorded design channels. The CD4 composition extension uses the
cell labels distributed with GSE174188 and does not re-annotate subtypes
independently.
The overlap strata are smaller than the full cohorts, and matched random subsets
separate sample-size loss from design restriction without identifying a causal
technical effect. Recorded batch can combine recruitment source, age, ancestry,
disease severity and laboratory processing; our result is non-attributability, not
proof that the retained signal is purely biological. GSE285773 lacks recoverable
donor processing-batch metadata. External bootstrap intervals condition on one fitted
source model, the 26-donor target has limited power, and external probabilities are
not clinically calibrated. Geneformer pretraining overlap with individual benchmark
cells was not tested. These boundaries do not weaken the central empirical finding:
residualisation, restriction and transfer produced different answers on the same
donors and cannot substitute for one another.

## 5. Methods

### 5.1. Cohorts, outcomes and statistical unit

We analysed GSE135779 (44 childhood donors; 33 SLE, 11 controls), the
CD4-positive alpha-beta T-cell subset of GSE174188 (261 donors; 162 SLE,
99 controls) and GSE285773 CD4-positive T cells (26 donors; 16 SLE,
10 controls) [1-3]. The donor was the independent unit in every supervised analysis.
Cells from one donor never crossed training and test partitions.

### 5.2. Recovery of GSE135779 design metadata

GSE135779 batch, collection year, age, sex, race, ethnicity and clinical variables
were read from the original Supplementary Table 1b; per-library sequencing metrics
were read from Supplementary Table 1c [1]. Study names were linked to analysis donor
identifiers through the GEO title/accession map. All 44 donors mapped uniquely.
The restoration script writes the donor table, batch-by-label, year-by-label and
batch-by-year cross-tabs, and SHA256 hashes of every source file.

For GSE174188, per-cell `Processing_Cohort` was read from the distributed CELLxGENE
object and reduced to donor fractions over four waves. Dominant wave was the
maximum-fraction wave; a pure-wave donor had at least 99% of analysed cells in one
wave.

### 5.3. Design-label information and observational equivalence

For centred binary label vector \(Y\), intercept-only projection \(M_1\), recorded
design matrix \(D\) including an intercept and residual projection \(M_D\), we
defined the design-adjusted label-information fraction as
\(\mathcal{I}_D=(Y^\top M_DY)/(Y^\top M_1Y)=1-R^2(Y\sim D)\). Values below
floating-point tolerance were set to zero. Relative variance inflation was
\(1/\mathcal{I}_D\), reported as infinite at zero information. Design rank and
numbers of mixed-label levels and donors were reported beside this measure.

Observational equivalence was demonstrated constructively. At \(D=Y\), the models
\(X=\beta Y+\epsilon\) and \(X=\delta D+\epsilon\) generate identical observations
when \(\beta=\delta\), despite assigning the entire effect to different mechanisms.
The numerical test used identical Gaussian noise and required maximum absolute
difference zero.

### 5.4. Identifiability simulation

The additive donor simulation generated 240 source and 2,000 independent target
donors with 20 features. The first coordinate received biological effect
\(\beta\in\{0,0.5,1,2\}\) and technical effect
\(\delta\in\{0,0.5,1,2\}\). A binary design variable was generated with
association \(\rho\in\{0,0.25,0.5,0.75,0.9,0.98,1\}\) to the balanced disease
label; target design was independent of disease. Each of 112 parameter cells used
100 deterministic replicates.

Balanced logistic regression compared unadjusted source cross-validation,
training-fold residualisation, maximum-overlap design restriction and source-only
target prediction. Restriction was evaluated in one design stratum selected for
maximum label overlap; probabilities from separately fitted strata were not pooled.
Biological-coefficient recovery used ordinary least squares on disease and design
and was undefined at rank deficiency. All 11,200 replicate rows and the exact
observational-equivalence check were retained.

The prespecified low-entanglement positive control used
\(\beta=1,\delta=0.5,\rho=0.25\). Five illustrative pass criteria were evaluated:
information fraction at least 0.80; complete-pipeline permutation
\(p\leq0.01\); design-only AUC at most 0.70; maximum residualisation or
restriction loss at most 0.05; and internal-minus-external AUC at most 0.05. The
thresholds define an interpretable demonstration of the protocol rather than a
universal acceptance rule.

### 5.5. Complete-pipeline design nulls

For every GSE135779 design block, 1,000 label permutations repeated the entire
nested donor pipeline, including fold construction, imputation, one-hot encoding,
scaling, inner regularisation selection and outer prediction. A second null used 250
independent shuffles of every design column while keeping disease labels fixed and
preserving each column's marginal distribution. Empirical upper-tail probabilities
used \((1+\sum I[T_b\ge T_{obs}])/(B+1)\). Null draws used eight parallel workers;
no molecular feature entered these tests.

### 5.6. CD4-subtype composition

The GSE174188 CELLxGENE object was restricted to the same CD4-positive alpha-beta
T-cell population and 261 donors used in the molecular analysis. Per-donor counts
were formed for `T4_naive`, `T4_em`, `T4_reg` and all other or unassigned labels.
We evaluated raw closed proportions and centred-log-ratio coordinates. For CLR,
0.5 was added to each donor-category count before closure; donor mean log abundance
was subtracted and the final coordinate omitted. Each representation was evaluated
with and without log total CD4-cell yield.

The four composition variants used the locked repeated donor classifier. Processing
wave was predicted one-versus-rest. Fold-contained ridge residualisation used
processing-wave fractions. Dominant-wave-4 and at-least-99%-pure-wave-4 restrictions
were compared with donor subsets matched exactly for size and case/control count
using the same 20 seeds as the main analysis.

### 5.7. Frozen Geneformer embeddings

Cell embeddings used Geneformer V2-316M at repository revision
`04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5` [4]. The archived checkpoint,
configuration, token dictionary and gene-median dictionary SHA256 hashes are provided
in the methods manifest. Gene identifiers were version-trimmed and mapped to the V2
vocabulary. Within each cell, the ranking value was
`raw count / total count * 10,000 / Genecorpus-104M gene median`; genes were sorted
descending. Sequences used V2 start/end token IDs 2/3, padding ID 0 and maximum
length 4096, leaving at most 4094 ranked genes.

The model output was the final `last_hidden_state`; token position 0 was retained as
the 1,152-dimensional cell CLS vector. The encoder was frozen. Donor representations
were coordinate-wise means of sampled cell vectors. Sampling without replacement used
integer seed 1 and caps of 500 cells per donor for GSE135779 and 1,000 for GSE174188
and GSE285773. The embedding environment used Python 3.10, PyTorch 2.5.1 and
Transformers 4.46.3.

### 5.8. Pseudobulk representations

Stored donor-by-gene values were `log1p(CPM)`. Within every outer training fold, gene
variance was computed on training donors, the 4,000 highest-variance genes were
selected, and coordinates were standardised by training-donor mean and standard
deviation. HVG pseudobulk used these values directly. PCA pseudobulk fitted up to
30 components on the training values and applied the fitted map to held-out donors.

For cross-cohort transfer, trimmed gene identifiers were intersected before
source-only HVG selection. Feature means, variances, selected genes, scaler, PCA,
classifier and regularisation were fitted on source donors and applied unchanged to
the target.

### 5.9. Repeated donor-level evaluation

Internal analyses used 20 repeated stratified five-fold donor splits with integer
seeds 20260801-20260820. The classifier was balanced logistic regression with the
liblinear solver and `max_iter=20,000`. Inverse regularisation strength was selected
from 10^-4 through 10^4 in decade steps using five-fold resampling within each outer
training set. The full outer-training representation was used for this inner
regularisation selection; no outer test donor entered HVG selection, scaling, PCA or
classifier fitting.

Out-of-fold donor predictions were pooled within each repeat. We report mean AUC and
the 2.5th-97.5th percentile of repeat-level AUCs. This is a descriptive
split-sensitivity range, not a population confidence interval or a test based on
treating correlated folds or repeats as independent [31].

### 5.10. Design-only and representation-to-design prediction

Numeric design variables were median-imputed and standardised inside the outer
training fold. Categorical variables were most-frequent-imputed and one-hot encoded,
with held-out unknown levels ignored. The same logistic model and regularisation grid
were used. GSE135779 blocks were batch, collection year, sequencing QC, demographics
and their union. GSE174188 blocks were processing-wave fractions, sequencing QC,
demographics and their union.

Each molecular representation also predicted every observed batch one-versus-rest
under five repeated donor-level cross-validations. The number of folds was capped by
the minority batch count.

As a model-flexibility sensitivity analysis, the complete GSE135779 design block was
also evaluated with random forests (256 trees, maximum depth 4, minimum leaf size 3)
and histogram gradient boosting (150 iterations, learning rate 0.05, seven terminal
leaves, minimum leaf size 5). Encoding was fitted inside every outer training fold;
these arms used fixed hyperparameters and were not selected against held-out donors.

### 5.11. Fold-contained residualisation

Unadjusted and residualised models shared outer folds, feature construction, scaling,
PCA and classifier selection. Within each outer training fold, design variables were
encoded as above and a multivariate ridge regression (`alpha=1`, intercept included)
mapped design to every raw representation coordinate. Training and held-out
coordinates were replaced by observed minus predicted values, after which the
identical representation-specific pipeline was fitted to training residuals and
applied to held-out residuals. The disease label was not included in the residualiser.

### 5.12. Nonlinear adjustment, weighting and attenuation uncertainty

The GSE135779 complete-design sensitivity was repeated with a multivariate
random-forest residualiser fitted only on outer-training donors (96 trees, maximum
depth 4, minimum leaf size 3). Pseudobulk genes were selected by training-donor
variance before this high-dimensional fit. A separate batch-location arm subtracted
training-batch mean offsets from training and held-out features; because it did not
apply empirical-Bayes scale shrinkage, we interpret it as ComBat-style location
adjustment rather than ComBat. A third arm fitted a training-only logistic disease
propensity model from the complete design, clipped propensities to 0.025-0.975 and
trained the representation classifier with overlap weights.

For the central GSE174188 comparison, residualisation loss was paired by split seed.
We bootstrapped the 20 paired repeat-level losses 5,000 times and subtracted the
locked matched-restriction discrepancy for dominant or pure wave 4. These intervals
quantify split sensitivity conditional on the analysed donors and matched-control
estimate; they are not population confidence intervals.

To probe confound leakage with a nonlinear downstream learner [28], 100 Gaussian
simulations used no biological effect, a technical effect of one and
design-label association 0.75. Random-forest disease prediction was compared on raw
features, features residualised globally and features residualised inside each
training fold.

### 5.13. Design restriction and matched controls

GSE174188 was restricted to dominant wave 4 and separately to pure wave 4.
GSE135779 excluded all-case batch B1 and retained B2-B6, each of which contained both
labels. The complete internal pipeline was rerun after restriction.

For every observed stratum, 20 donor subsets were drawn without replacement from the
full cohort with the identical donor count and case/control count. Each draw used one
prespecified split seed. The observed-minus-matched difference quantifies whether the
stratum is harder than expected from reduced donor number and label composition.

### 5.14. Cross-batch and cross-cohort transfer

Cross-batch analyses trained on waves 2+3 and evaluated wave 4, then reversed the
direction. Cross-cohort transfer fitted every statistic on the source cohort. Target
labels were accessed only after prediction.

Target AUC intervals used 5,000 case/control-stratified donor bootstrap resamples and
the percentile method. Paired AUC comparisons used DeLong tests on the identical
target donors [23] and Benjamini-Hochberg correction within declared method families
[24].

### 5.15. Learned pooling

DeepSets, gated-attention multiple-instance learning and pooling by multi-head
attention operated on cap-500 frozen cell embeddings [20-22]. Cell standardisation, a
whitened source PCA32 projection, network weights and hyperparameters were
source-fitted. Hidden widths were 16 or 32, weight decay 0.01 or 0.1, and Adam used a
learning rate of 0.02 for 150 epochs. Five-fold source-donor cross-validation selected
the configuration; three initialisations were fitted on all source donors and target
probabilities were averaged. An ordinary donor mean in the identical PCA32 space
isolated the contribution of learned pooling.

The PCA32 projection used for hyperparameter selection was fitted once on all source
cells before source-donor folds. It never used target cells or labels. Consequently,
source-internal pooling AUC was not used as an unbiased generalisation estimate; only
independent-target and paired target-donor comparisons support the learned-pooling
conclusion.

### 5.16. Software and provenance

The locked donor-level reanalysis used Python 3.13.7, NumPy 2.4.6, pandas 3.0.3,
SciPy 1.18.0, scikit-learn 1.9.0, PyArrow 21.0.0 and joblib 1.5.3. Input hashes,
seeds, parameter grids, donor-level predictions, complete null distributions and
output manifests accompany the analysis. The identifiability extension contains 18
automated validation checks covering simulation cardinality, exact observational
equivalence, information monotonicity, empirical donor counts, composition
attenuation and publication-file completeness.

## 6. Data and code availability

The source datasets are publicly available through GEO accessions GSE135779,
GSE174188 and GSE285773. The MIT-licensed public repository is available at
https://github.com/LightChainr/rheumlens, and its permanent Zenodo record is
https://doi.org/10.5281/zenodo.20813922. The complete versioned research object
submitted for review accompanies this manuscript. It contains restored GSE135779
metadata, donor-level inputs and predictions, analysis scripts, two locked software
environments, 11,200 simulation replicates, complete null distributions, robustness
outputs, figure source tables, `REPRODUCE.md`, `SHA256SUMS`,
`RELEASE_MANIFEST.json` and an executable package validator. Public raw
single-cell matrices remain at their GEO accessions.

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
19. Chaibub Neto E. Using linear residualization to remove a potential
    confounder from a machine learning model's predictions. *Proceedings of
    Machine Learning Research*. 2021;139:1361-1371.
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
30. Chaibub Neto E, Pratap A. A permutation approach to assess confounding in
    machine learning applications for digital health. In: *Proceedings of the
    25th ACM SIGKDD International Conference on Knowledge Discovery and Data
    Mining*. 2019:54-64. doi:10.1145/3292500.3330790.
31. Zeng T, Li H, Zhang S, et al. Widespread use of invalid statistical tests
    in biomedical machine learning. *bioRxiv*. 2026.
    doi:10.64898/2026.05.17.724301.

## Figure files

1. `figures/main/Figure_1_two_cohort_design_channels.svg`
2. `figures/main/Figure_2_simulated_identifiability_landscape.svg`
3. `figures/main/Figure_3_empirical_null_and_composition.svg`
4. `figures/main/Figure_4_disease_and_design_information.svg`
5. `figures/main/Figure_5_design_restriction_matched_controls.svg`
6. `figures/main/Figure_6_residualisation_failure_modes.svg`
7. `figures/main/Figure_7_generalisation_ladder.svg`
8. `figures/main/Figure_8_learned_pooling.svg`
