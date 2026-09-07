# Technical audit of the reconstructed manuscript line

**Audit date:** 2026-07-27  
**Status:** closed for the design-validity release candidate.

## Findings requiring correction

### P0 — archived residualisation was not matched to the new classifier pipeline

The draft compared an archived residualisation analysis (fixed `C=1`, no common
representation standardisation) with a new stratified analysis (standardisation and
inner `C` selection), while describing them as the same data, features and folds.

**Resolution:** `scripts/locked_validity_audit.py` recomputes unadjusted and
batch-residualised predictions on identical outer folds with the same HVG/PCA,
standardisation, `C` grid and classifier. The old 5–10-fold headline is withdrawn
unless reproduced by this output.

### P0 — GSE135779 metadata availability was stated incorrectly

The draft said only GSE174188 distributed auditable design metadata. The local public
evidence package contains the original GSE135779 Supplementary Tables 1b and 1c,
including donor batch, collection year, demographics and sequencing QC.

**Resolution:** all 44 childhood-cohort analysis donors map without missing values.
The restoration, cross-tabs and source SHA256 records are under
`results/gse135779_metadata/`. The manuscript must now say two cohorts are auditable;
GSE285773 remains unaudited for donor-level processing batch.

### P1 — learned-pooling source-CV AUC was overstated as an unbiased internal estimate

The learned-pooling target predictions are valid source-only predictions. However,
the source PCA32 projection used in network hyperparameter selection was fitted on all
source cells before source-donor CV. This is unsupervised source-side reuse, not target
leakage, but it makes the numerical source-CV-to-target gap unsuitable as a headline
generalisation estimate.

**Resolution:** remove the 0.10–0.19 source-CV inflation claim. Figure 6A now uses
external target AUC and target-donor bootstrap intervals. Preserve the independent
target and paired pooling comparisons.

### P1 — "pure wave 4, n=89" was a factual mislabel

The n=89 source restriction in `deconfounded_transfer.py` uses dominant wave 4. The
pure-wave-4 set (wave fraction at least 0.99) contains n=66.

**Resolution:** all text and figure captions must call n=89 "dominant wave 4" or
"wave-4 source"; reserve "pure wave 4" for n=66.

### P1 — residualisation and design restriction do not estimate the same causal object

Even under a matched implementation, residualisation removes all linearly
batch-predictable coordinates, including disease-aligned biology when batch and label
are collinear. Restriction estimates performance in observed overlap strata and loses
power. Their difference is informative but not a formal estimate of residualisation
bias without additional causal assumptions.

**Resolution:** use "empirical discrepancy" or "downward distortion under
collinearity", not "true batch effect" or "causal bias estimate". The covariate-only
AUC calibrates non-identifiability.

### P1 — clinical probability interpretation remains unsupported

The 261-to-26 transfer direction has poor Brier and ECE values despite high AUC.

**Resolution:** external AUC is ranking discrimination only. Do not claim clinical
risk prediction, threshold utility or individual decision support. No DCA should be
promoted as clinical evidence without calibrated prospective probabilities.

### P2 — pretraining overlap is not audited

Geneformer V2 may have encountered related public data distributions during
pretraining. The July 2026 scContam work makes this an increasingly visible benchmark
validity issue.

**Resolution:** record the checkpoint revision and data publication dates, state that
instance-level pretraining overlap was not tested, and avoid calling frozen transfer
"zero-shot clinical generalisation". This limitation does not explain pseudobulk's
advantage or the observed batch-label design structure.

## Closed reproducibility details

- Final hidden layer and exact CLS token position: closed.
- Embedding dimension: 1152.
- V2 special tokens, token limit and ranking normalization: closed.
- Checkpoint, config, token dictionary and median dictionary SHA256: closed.
- Cell cap and sampling seed by cohort: closed.
- Donor-level split count, repeat seeds, stratification and class weighting: closed.
- HVG count, PCA dimension, standardisation scope and `C` grid: closed.
- Bootstrap count, stratification and percentile method: closed.
- Paired DeLong target-donor definition and BH families: closed.
- GSE135779 batch provenance and donor mapping: closed.

## Closed in the 2026-07-27 integration pass

1. The locked two-cohort run completed, including GSE135779 batch, collection-year
   and complete-design residualisation.
2. The abstract, Results, Discussion, Limitations and Methods were rewritten in
   `manuscript/design_validity_draft_v2.md`.
3. Figures 1-4 were regenerated from the locked tables as editable SVG/PDF/PNG.
   Figure 6 no longer uses source-internal learned-pooling AUC as an unbiased
   generalisation estimate. The n=89 source is labelled dominant wave 4.
4. The reference and innovation review now covers singleDeep, PaSCient, scPhase,
   FloREN, donor-aware IBD benchmarking, scFM benchmarks, Soneson 2014 and
   residualisation.
5. `scripts/validate_design_validity_release.py` passes and writes a hash-locked
   numerical validation record.

## Repository packaging closure

- Creator metadata now names Hongyu Ying, Dandan Yun and Dan Liu.
- The manuscript, figures, source tables, donor-level inputs and release metadata are
  synchronised under one SHA256 manifest.
- The reconstruction is prepared as `v2.0.0-rc1`, not patch version `v1.0.3`.
- A clean-copy rerun reproduced the donor-level audit and transfer outputs.
- The remaining checks concern author declarations and final DOI/tag substitution,
  not analytical validity.
