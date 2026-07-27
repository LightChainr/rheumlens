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
vertical dashed line marks chance performance.

## Figure 2. Design-label association determines how much disease contrast remains estimable

(A) Mean design-adjusted label-information fraction and 2.5th-97.5th percentile
across 100 replicates at each level of design-label association. (B) Change in
internal AUC after fold-contained residualisation when the biological effect equals
one, across technical-effect magnitudes and design-label associations. (C)
Independent-target minus internal AUC under the same parameter grid; target design
was independent of disease. (D) Median error in recovery of the biological
coefficient. The simulation crossed four biological effects, four technical effects
and seven design-label associations for 11,200 replicates. At perfect collinearity,
biological-only and technical-only mechanisms are observationally equivalent.

## Figure 3. Empirical design diagnostics reproduce the simulated identifiability failure

(A) Observed GSE135779 design-only AUCs over full-pipeline label-permutation
distributions. The complete design reached AUC 0.952 with empirical upper-tail
p=0.000999. (B) Design-only disease AUC and the fraction of label information
remaining after linear design adjustment. (C) GSE174188 disease AUC for raw and
centred-log-ratio CD4-subtype composition, with or without log cell yield, before
and after fold-contained processing-wave residualisation. Points denote means and
bars the 2.5th-97.5th percentile across 20 repeated donor splits. (D) Difference
between observed composition AUC within dominant or pure wave 4 and AUC from donor
subsets matched for size and case-control count. Intervals crossing zero indicate
that the restricted stratum was not consistently harder than its matched controls.

## Figure 4. High disease discrimination coexists with near-perfect design recognition

(A-B) Mean repeated-cross-validation ROC-AUC for case-control discrimination
(filled points) and for the most predictable recorded batch in a one-versus-rest
analysis (open points) using frozen Geneformer donor means, HVG pseudobulk and PCA
pseudobulk in GSE174188 (A) and GSE135779 (B). Horizontal segments connect the two
tasks for each representation. (C) Matrix of disease AUC, best-batch AUC and their
difference. Batch predictability documents design information retained by each
representation but does not by itself identify that information as technical or
biological.

## Figure 5. Disease discrimination survives design restriction after controlling for lost sample size

(A) GSE174188 donors assigned to dominant processing wave 4 (n=89). (B) The stricter
GSE174188 pure-wave-4 subset, requiring at least 99% of analysed cells in wave 4
(n=66). (C) GSE135779 donors from mixed-label batches B2-B6 after excluding the
all-case B1 batch (n=36). Coloured points show mean repeated-cross-validation AUC in
the observed stratum. Open points and grey bars show the mean and
2.5th-97.5th percentile from 20 random donor subsets matched to donor count and
case-control count. Labels report observed minus matched mean AUC.

## Figure 6. Residualisation measures sensitivity, not a causal batch effect

(A) GSE174188 AUC before and after fold-contained residualisation on processing-wave
fractions; the projected design block itself predicts disease at AUC 0.922.
(B) GSE135779 internal negative control: representations predict batch at AUC up to
0.993, but batch alone is uninformative for disease at AUC 0.499 and batch
residualisation changes disease AUC by at most 0.009. (C) Change in GSE135779 AUC
after projecting out batch, collection year or the complete measured design block.
(D) Across prespecified design blocks, the mean AUC removed from molecular
representations increases with disease predictability of the projected block.
(E) GSE135779 AUC after training-only batch location adjustment, complete-design
overlap weighting, random-forest residualisation or ridge residualisation. Points
denote means and bars denote 2.5th-97.5th percentiles across 20 repeated donor
splits. (F) Residualisation loss minus the matched restriction loss for dominant and
pure wave 4. Bars are descriptive 2.5th-97.5th percentile intervals from 5,000
paired resamples of repeat-level residualisation losses, conditional on the locked
matched-restriction estimate. Adjusted and unadjusted estimates use identical outer
donor folds and training-only preprocessing.

## Figure 7. Generalisation, not internal separability, is where donor-level models break

(A) Donor-level ROC-AUC for within-cohort evaluation, bidirectional cross-wave
evaluation in GSE174188 and strict source-only transfer from GSE285773 to GSE174188.
(B) Target ROC curves for the 261-donor GSE174188 cohort after training on
GSE285773, with pointwise bands from 200 stratified donor bootstrap resamples.
(C) Overlap among target donors whose scores fell on the wrong side of target
prevalence for the three representations. (D) Transfer to GSE285773 after
restricting the GSE174188 source to dominant wave 4, compared with matched random
source subsets and all 261 source donors. (E) AUC across the sequence from
within-cohort evaluation to pure-wave, cross-batch and cross-cohort evaluation.
Every transformation and classifier in cross-cohort analyses was fitted on source
donors only.

## Figure 8. Learned pooling does not rescue the frozen representation

(A) Independent-target AUC and 95% stratified donor-bootstrap intervals for DeepSets,
gated-attention multiple-instance learning and pooling by multi-head attention in
both transfer directions. (B) Paired target-donor AUC differences for the 261-donor
target relative to the 1,152-dimensional frozen Geneformer donor mean and a donor
mean in the identical source-fitted PCA32 space. Red bars denote
Benjamini-Hochberg-adjusted q&lt;0.05. (C) Distribution across source donors of the
share of gated-attention mass assigned to the top 10% of cells, stratified by
processing wave. (D) Spearman correlations among independent-target donor scores.
(E) Target AUC for learned operators and fixed mean pooling in the 261-donor target.
(F) Strict source-only cross-cohort AUC for fixed Geneformer mean pooling, the best
learned operator and the two expression pseudobulk baselines. Learned pooling was
selected using source donors only.
