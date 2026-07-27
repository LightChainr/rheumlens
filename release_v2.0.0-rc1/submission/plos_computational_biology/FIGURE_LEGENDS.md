# Figure legends

## Figure 1. Outcome labels are entangled with study design in two widely reused SLE cohorts

(A) Case and control counts across the four processing waves represented among the
analysed GSE174188 CD4-positive alpha-beta T-cell donors. Wave 1 contains no SLE
cases, whereas wave 2 is dominated by cases. (B) Case and control counts across the
six restored GSE135779 sequencing batches. Batch B1 contains only cases; batches
B2-B6 contain both labels. (C) Case and control counts by GSE135779 collection year.
(D) Repeated donor-level cross-validated ROC-AUC obtained using recorded design
variables alone. Points denote mean AUC across 20 repeated stratified five-fold
analyses and bars denote the 2.5th-97.5th percentile of repeat-level AUCs. The
vertical dashed line marks chance performance. Purple and ochre indicate GSE174188
and GSE135779, respectively.

## Figure 2. High disease discrimination coexists with near-perfect design recognition

(A-B) Mean repeated-cross-validation ROC-AUC for case-control discrimination (filled
points) and for the most predictable recorded batch in a one-versus-rest analysis
(open points) using frozen Geneformer donor means, HVG pseudobulk and PCA
pseudobulk in GSE174188 (A) and GSE135779 (B). Horizontal segments connect the two
tasks for each representation. (C) Matrix of disease AUC, best-batch AUC and their
difference. Values are donor-level estimates; batch predictability documents design
information retained by each representation but does not by itself identify that
information as technical or biological.

## Figure 3. Disease discrimination survives design restriction after controlling for lost sample size

(A) GSE174188 donors assigned to dominant processing wave 4 (n=89). (B) The stricter
GSE174188 pure-wave-4 subset, requiring at least 99% of analysed cells in wave 4
(n=66). (C) GSE135779 donors from mixed-label batches B2-B6 after excluding the
all-case B1 batch (n=36). Coloured points show mean repeated-cross-validation AUC in
the observed stratum. Open points and grey bars show the mean and 2.5th-97.5th
percentile from 20 random donor subsets matched to the observed donor count and
case-control count. Labels report observed minus matched mean AUC. Restriction
produced only small performance changes after accounting for the smaller analysis
set.

## Figure 4. Residualisation is a sensitivity operation, not a causal batch-effect estimate

(A) GSE174188 AUC before and after fold-contained residualisation on processing-wave
fractions; the projected design block itself predicts disease at AUC 0.922.
(B) GSE135779 AUC before and after residualisation on batch; batch alone is
uninformative for disease at AUC 0.499. (C) Change in GSE135779 AUC after projecting
out batch, collection year or the complete measured design block. (D) Across
prespecified design blocks, the mean AUC removed from the molecular representations
increases with the disease predictability of the projected block. All adjusted and
unadjusted estimates use identical outer donor folds and training-only preprocessing.
Residualisation therefore quantifies sensitivity to deleting design-aligned
coordinates; under label-design collinearity it cannot identify a causal technical
batch effect.

## Figure 5. Generalisation, not internal separability, is where donor-level models break

(A) Donor-level ROC-AUC for within-cohort evaluation, bidirectional cross-wave
evaluation in GSE174188 and strict source-only transfer from GSE285773 to GSE174188.
Horizontal lines are 95% donor-bootstrap intervals for transfer estimates; the
within-cohort points summarize 20 repeated stratified five-fold analyses. (B) Target
ROC curves for the 261-donor GSE174188 cohort after training on GSE285773, with
pointwise bands from 200 stratified donor bootstrap resamples. (C) Overlap among
target donors whose scores fell on the wrong side of target prevalence for the three
representations, illustrating shared rather than method-specific failures.
(D) Transfer to GSE285773 after restricting the GSE174188 source to dominant wave 4,
compared with matched random source subsets and all 261 source donors. (E) AUC across
the sequence from within-cohort evaluation to pure-wave, cross-batch and cross-cohort
evaluation. Every transformation and classifier in cross-cohort analyses was fitted
on source donors only.

## Figure 6. Learned pooling does not rescue the frozen representation

(A) Independent-target AUC and 95% stratified donor-bootstrap intervals for DeepSets,
gated-attention multiple-instance learning and pooling by multi-head attention in
both transfer directions. Colours identify the source-target direction. (B) Paired
target-donor AUC differences for the 261-donor target relative to the
1,152-dimensional frozen Geneformer donor mean and to a donor mean in the identical
source-fitted PCA32 space. Red bars denote Benjamini-Hochberg-adjusted q<0.05.
(C) Distribution across source donors of the share of gated-attention mass assigned
to the top 10% of cells, stratified by processing wave; the dashed line is the value
under uniform attention. (D) Spearman correlations among independent-target donor
scores from deterministic and learned pooling operators. (E) Target AUC for learned
operators and fixed mean pooling in the 261-donor target. (F) Strict source-only
cross-cohort AUC for fixed Geneformer mean pooling, the best learned operator and the
two expression pseudobulk baselines. Learned pooling was selected using source donors
only and did not recover the external performance gap.
