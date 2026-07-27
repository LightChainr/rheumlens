# When design predicts disease: residualisation failure and generalisation gaps in donor-level single-cell classifiers

Hongyu Ying<sup>1</sup>, Dandan Yun<sup>1</sup> and Dan Liu<sup>1,*</sup>

<sup>1</sup> Department of Rheumatology and Immunology, Shanghai Pudong
Hospital, Fudan University Pudong Medical Center, Shanghai 201399, China  
<sup>*</sup> Correspondence: Dan Liu; danliu600@126.com

## Abstract

Single-cell RNA sequencing can resolve disease-associated immune states at cellular
resolution, but clinical classification is performed at the patient level. A rapidly
growing literature therefore converts thousands of cells into one donor
representation and reports cross-validated discrimination, often above ROC-AUC 0.95.
Whether these values reflect disease biology, study design, or a mixture of both is
rarely tested.

We reconstructed donor-level design metadata for two widely reused systemic lupus
erythematosus (SLE) cohorts and evaluated frozen Geneformer donor means, highly
variable gene (HVG) pseudobulk and PCA pseudobulk under one locked donor-level
pipeline. In GSE174188, processing-wave composition alone predicted disease at AUC
0.922, while in GSE135779 the complete measured design, sequencing and demographic
block reached 0.952. Every representation recovered recorded batch identity with
AUC up to 0.9997. Yet the inference that internal discrimination was merely a batch
shortcut was contradicted by design restriction: limiting GSE174188 to dominant or
pure wave 4 reduced AUC by only 0.014-0.049 relative to size- and label-matched donor
subsets, and removing the all-case B1 batch from GSE135779 changed AUC by
-0.012 to +0.003.

Conventional residualisation gave the opposite answer. Removing processing-wave
coordinates in GSE174188 reduced representation AUC by 0.230-0.268. In GSE135779,
residualising batch alone changed almost nothing, whereas projecting out the complete
design block collapsed all three representations to AUC 0.515-0.526. The magnitude of
this collapse tracked how accurately the projected covariates predicted the outcome,
not how much performance disappeared under design restriction. Residualisation
therefore cannot identify technical confounding when design and disease are
collinear; it removes disease-aligned signal by construction.

Independent transfer exposed the failure that internal adjustment could not resolve.
Moving across processing waves reduced AUC to 0.806-0.926, and strict source-only
cross-cohort transfer favoured expression pseudobulk over frozen Geneformer in the
261-donor target. Three learned set-pooling architectures did not recover the gap.
We introduce an executable validity standard for patient-level single-cell
classification: report design-only AUC, representation-to-design predictability,
design-restricted performance with matched donor controls, and source-only external
performance. High internal AUC is evidence of separability; only this chain tests
whether the separation is attributable and transportable.

**Keywords:** patient-level classification; single-cell RNA sequencing; study design;
batch confounding; residualisation; external validation; Geneformer; pseudobulk;
systemic lupus erythematosus

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
single-cell design is statistically non-identifiable [17]. Nevertheless, most
patient-level single-cell studies treat this problem as a preprocessing question:
regress a representation on technical covariates and re-evaluate the residuals.
When those covariates are also proxies for disease, residualisation removes the
disease-aligned component whether it is technical or biological [18,19]. A lower
post-residualisation AUC can therefore be mistaken for proof that the original model
used a shortcut.

This problem has become more urgent, not less. Broad benchmarks show that frozen
single-cell foundation-model embeddings can retain batch structure and can lose to
simple expression baselines [13,14]. Patient-level methods already report high
internal performance in SLE, including external validation between GSE174188 and
GSE135779 [7]. FloREN now evaluates both cohorts with supervised graph
representations [11], while pretraining-contamination work argues that benchmark
validity must accompany model rank [15]. What remains missing is a direct empirical
test of whether design predictability, residualisation, design restriction and
external transfer tell the same story on the same donors. They do not.

Here we recover donor design metadata for GSE135779 and combine it with the
cell-distributed processing metadata in GSE174188. Across frozen Geneformer and two
expression pseudobulk representations, we place design-only prediction,
representation-to-design prediction, fold-contained residualisation,
design-restricted evaluation, matched donor controls, cross-batch transfer,
cross-cohort transfer and learned pooling under one donor-matched framework. To our
knowledge, this is the first patient-level single-cell study to quantify this complete
validity chain across two disease cohorts. The result is not another leaderboard. It
is a demonstration that separability, adjustment and transportability are different
claims and require different experiments.

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

### 2.2. The representations remember the design more accurately than the disease

We next asked whether the donor representations retained the recorded design
identity. Under the same fold-contained pipeline, each representation predicted
case/control status and every batch one-versus-rest (Figure 2).

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

### 2.3. Disease discrimination survives design restriction

We first restricted GSE174188 to dominant wave 4, the only large near-balanced wave
(n=89; 46 cases and 43 controls). Because reducing 261 donors to 89 inevitably lowers
performance, we compared the observed stratum with 20 random donor subsets having
the same donor count and case/control composition (Figure 3A). Relative to those
matched subsets, AUC changed by -0.014 for Geneformer, -0.028 for HVG pseudobulk and
-0.049 for PCA pseudobulk.

The stricter analysis retained only donors with at least 99% of analysed cells in
wave 4 (n=66; 25 cases and 41 controls). Relative to size- and label-matched random
subsets, AUC changed by -0.026, -0.044 and -0.044 (Figure 3B). These are detectable
losses, but they are far smaller than the raw full-cohort-to-stratum difference. Most
of that raw difference is explained by using fewer donors.

GSE135779 provided an independent design. Removing the all-case B1 batch left 36
donors from batches B2-B6, with both labels represented in every retained batch.
Geneformer AUC was 0.846 versus 0.858 in matched random subsets (difference -0.012).
Both pseudobulk representations were 0.919 versus 0.916 (difference +0.003;
Figure 3C). Thus the disease signal persisted after the most obvious label-exclusive
batch was removed.

Across both cohorts, design restriction weakens the claim that internal AUC is a pure
technical shortcut. The representations carry design information, but the measured
disease separation does not disappear when the clearest design channels are held
fixed.

### 2.4. Residualisation can erase the signal that restriction preserves

Residualisation produced a sharply different conclusion. All unadjusted and adjusted
models were recomputed on identical outer folds with the same training-only HVG
selection, scaling, PCA construction, regularisation grid and classifier.

In GSE174188, ridge residualisation on processing-wave fractions reduced Geneformer
AUC from 0.982 to 0.714, HVG pseudobulk from 0.977 to 0.747 and PCA pseudobulk from
0.975 to 0.735 (losses 0.268, 0.230 and 0.241; Figure 4A). The same design, tested by
dominant- or pure-wave restriction with matched donor controls, cost only
0.014-0.049. Residualisation therefore removed five to nineteen times more AUC than
the matched restriction discrepancy, depending on representation and stratum.

GSE135779 revealed why. Batch alone had a cross-validated disease AUC of 0.499, and
batch residualisation changed representation AUC by +0.002, -0.002 and -0.009
(Figure 4B). Collection year had AUC 0.721 and its removal cost 0.005-0.039.
The complete measured design block had AUC 0.952; projecting it out collapsed
Geneformer, HVG pseudobulk and PCA pseudobulk to 0.516, 0.526 and 0.515, removing
0.354, 0.407 and 0.430 AUC (Figure 4C-D). Yet exclusion of the all-case batch had
changed matched performance by at most 0.012.

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

### 2.5. Generalisation, rather than internal separability, is the reproducible failure

Design restriction preserved strong within-cohort discrimination, but transfer
weakened it at every level tested (Figure 5). Training on processing waves 2+3 and
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

### 2.6. Learned pooling does not repair the transfer gap

Mean pooling might suppress informative cell heterogeneity, so we tested DeepSets,
gated-attention multiple-instance learning and pooling by multi-head attention
[20-22]. All operators used the same frozen cell embeddings, cap-500 cells per donor,
a source-fitted whitened PCA32 projection and source-only hyperparameter selection.
Independent target data were never used to select a model (Figure 6).

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

### 2.7. An executable validity standard

The results support four minimum requirements for patient-level single-cell
classification:

1. **Report design-only disease AUC.** Batch, year, centre, cell yield, sequencing
   quality and demographic variables should be evaluated as a label-only baseline.
2. **Report representation-to-design predictability.** A representation that retains
   design identity at near-perfect AUC has documented exposure to cohort structure,
   even if attribution remains unresolved.
3. **Pair residualisation with design restriction and matched donor controls.**
   Residualisation alone is uninterpretable under label-design collinearity.
4. **Require a source-only external target.** Target labels must not choose feature
   transformations, and internal performance must be reported beside external
   performance.

These checks are inexpensive relative to embedding generation or deep pooling. They
convert “the model achieved high AUC” into four separable statements: the label is
predictable, the design is visible, the signal survives observed overlap, and the
ranking transfers.

## 3. Discussion

This study identifies a validity failure that representation leaderboards do not
measure. In two widely reused SLE cohorts, measured study design can classify disease
almost as accurately as the molecular representations, and every representation can
recover batch identity with near-perfect accuracy. That combination should trigger a
confounding investigation. It does not, however, tell us whether the disease
classifier is a batch classifier.

The decisive comparison is between two adjustment strategies. Restriction asks
whether discrimination survives among donors observed under the same or overlapping
design. Residualisation deletes every coordinate predictable from a chosen design
block. When design and disease are collinear, these operations target different
quantities. The two-cohort result makes the distinction visible: mixed-batch
restriction preserved GSE135779 performance, while residualising the complete
design block reduced every representation to chance. In GSE174188, wave restriction
caused modest matched losses, while wave residualisation removed roughly one quarter
of the AUC scale.

The result is more consequential than a contest between Geneformer and pseudobulk.
Recent benchmarks already show that simple expression baselines can outperform
frozen foundation-model embeddings [13,14], and patient-level representation learning
is an active field [8-11]. Our contribution is to show why an internal comparison can
remain scientifically ambiguous even when donor leakage is prevented, folds are
strict and adjustment is performed inside each fold. Correct computation does not
repair a non-identifiable design.

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

### 5.3. Frozen Geneformer embeddings

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

### 5.4. Pseudobulk representations

Stored donor-by-gene values were `log1p(CPM)`. Within every outer training fold, gene
variance was computed on training donors, the 4,000 highest-variance genes were
selected, and coordinates were standardised by training-donor mean and standard
deviation. HVG pseudobulk used these values directly. PCA pseudobulk fitted up to
30 components on the training values and applied the fitted map to held-out donors.

For cross-cohort transfer, trimmed gene identifiers were intersected before
source-only HVG selection. Feature means, variances, selected genes, scaler, PCA,
classifier and regularisation were fitted on source donors and applied unchanged to
the target.

### 5.5. Repeated donor-level evaluation

Internal analyses used 20 repeated stratified five-fold donor splits with integer
seeds 20260801-20260820. The classifier was balanced logistic regression with the
liblinear solver and `max_iter=20,000`. Inverse regularisation strength was selected
from 10^-4 through 10^4 in decade steps using five-fold resampling within each outer
training set. The full outer-training representation was used for this inner
regularisation selection; no outer test donor entered HVG selection, scaling, PCA or
classifier fitting.

Out-of-fold donor predictions were pooled within each repeat. We report mean AUC and
the 2.5th-97.5th percentile of repeat-level AUCs. This is a split-sensitivity range,
not a population confidence interval.

### 5.6. Design-only and representation-to-design prediction

Numeric design variables were median-imputed and standardised inside the outer
training fold. Categorical variables were most-frequent-imputed and one-hot encoded,
with held-out unknown levels ignored. The same logistic model and regularisation grid
were used. GSE135779 blocks were batch, collection year, sequencing QC, demographics
and their union. GSE174188 blocks were processing-wave fractions, sequencing QC,
demographics and their union.

Each molecular representation also predicted every observed batch one-versus-rest
under five repeated donor-level cross-validations. The number of folds was capped by
the minority batch count.

### 5.7. Fold-contained residualisation

Unadjusted and residualised models shared outer folds, feature construction, scaling,
PCA and classifier selection. Within each outer training fold, design variables were
encoded as above and a multivariate ridge regression (`alpha=1`, intercept included)
mapped design to every raw representation coordinate. Training and held-out
coordinates were replaced by observed minus predicted values, after which the
identical representation-specific pipeline was fitted to training residuals and
applied to held-out residuals. The disease label was not included in the residualiser.

### 5.8. Design restriction and matched controls

GSE174188 was restricted to dominant wave 4 and separately to pure wave 4.
GSE135779 excluded all-case batch B1 and retained B2-B6, each of which contained both
labels. The complete internal pipeline was rerun after restriction.

For every observed stratum, 20 donor subsets were drawn without replacement from the
full cohort with the identical donor count and case/control count. Each draw used one
prespecified split seed. The observed-minus-matched difference quantifies whether the
stratum is harder than expected from reduced donor number and label composition.

### 5.9. Cross-batch and cross-cohort transfer

Cross-batch analyses trained on waves 2+3 and evaluated wave 4, then reversed the
direction. Cross-cohort transfer fitted every statistic on the source cohort. Target
labels were accessed only after prediction.

Target AUC intervals used 5,000 case/control-stratified donor bootstrap resamples and
the percentile method. Paired AUC comparisons used DeLong tests on the identical
target donors [23] and Benjamini-Hochberg correction within declared method families
[24].

### 5.10. Learned pooling

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

### 5.11. Software and provenance

The locked donor-level reanalysis used Python 3.13.7, NumPy 2.4.6, pandas 3.0.3,
SciPy 1.18.0, scikit-learn 1.9.0, PyArrow 21.0.0 and joblib 1.5.3. Input hashes,
seeds, parameter grids, donor-level predictions and output manifests accompany the
analysis.

## 6. Data and code availability

The source datasets are publicly available through GEO accessions GSE135779,
GSE174188 and GSE285773. The current public benchmark archive is available at
https://github.com/LightChainr/rheumlens and concept DOI
https://doi.org/10.5281/zenodo.20813922. The design-validity reconstruction,
restored GSE135779 metadata, locked donor predictions and new figure source tables
will be deposited as a new major scientific release after final manuscript and
author-metadata validation.

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
14. Wu Y, et al. A comprehensive benchmark of single-cell foundation models.
    *Genome Biology*. 2025;26:334. doi:10.1186/s13059-025-03781-6.
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

## Figure files

1. `figures_locked/Figure_1_two_cohort_design_channels.svg`
2. `figures_locked/Figure_2_disease_and_design_information.svg`
3. `figures_locked/Figure_3_design_restriction_matched_controls.svg`
4. `figures_locked/Figure_4_residualisation_failure_modes.svg`
5. `figures_v2_refreshed/Figure_5_generalisation_ladder.svg`
6. `figures_v2_refreshed/Figure_6_learned_pooling.svg`
