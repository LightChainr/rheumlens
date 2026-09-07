# Literature and novelty audit

**Search lock:** 2026-07-26  
**Scope:** donor/patient-level classification from scRNA-seq, patient representations,
single-cell foundation-model benchmarks, design/batch confounding, residualisation,
cross-cohort evaluation, and SLE studies using GSE135779 or GSE174188.

## Search method

Primary publication records and full texts were checked through PubMed/Europe PMC,
OpenAlex, Crossref-resolved DOI pages, bioRxiv and arXiv. Searches combined:
`patient-level`, `donor-level`, `sample representation`, `single-cell`, `classification`,
`foundation model`, `Geneformer`, `pseudobulk`, `multiple-instance learning`,
`cross-cohort`, `batch confounding`, `residualization`, `GSE135779`, `GSE174188`, and
`SLE`. Reference lists and "similar works" records were followed for direct
comparators. The final date filter included records public by 2026-07-26; this matters
because FloREN and scContam appeared during the week of the audit.

## Executive verdict

The paper has a defensible and timely novelty position, but it is not "Geneformer versus
pseudobulk" and it is not the discovery that batch confounding can bias cross-validation.
Both spaces are already occupied.

The strongest claim is narrower and more consequential:

> This is a two-cohort empirical validity audit of donor-level single-cell disease
> classification that places design-only predictability, representation-to-design
> predictability, fold-contained residualisation, design-restricted evaluation,
> size/label-matched controls, cross-batch transfer, cross-cohort transfer, and learned
> pooling under one donor-matched protocol.

The literature search did not identify a prior donor-level scRNA-seq paper that reports
that complete chain. The claim should be made as a conjunction, not by declaring any one
component unprecedented.

### Locked-result addendum (2026-07-27)

The two-cohort rerun strengthened the novelty position. Complete measured design
predicts disease at AUC 0.953 in GSE174188 and 0.952 in GSE135779, yet matched
design restriction causes only small performance changes. In contrast,
residualising the highly label-predictive design block removes 0.230-0.268 AUC in
GSE174188 and 0.354-0.430 in GSE135779. The second cohort also supplies a
negative control: batch alone predicts disease at AUC 0.499 and batch
residualisation changes almost nothing, despite near-perfect representation-to-batch
prediction. This makes the paper a two-cohort demonstration that design exposure,
outcome attribution and transportability are distinct empirical questions.

## What is already occupied

### Patient-level representations and learned pooling

- [CloudPred](https://pubmed.ncbi.nlm.nih.gov/34890161/) already established supervised
  patient phenotype prediction from variable-size single-cell sets and included lupus.
- [ProtoCell4P](https://pmc.ncbi.nlm.nih.gov/articles/PMC10444962/) and
  [CELLECTION](https://pmc.ncbi.nlm.nih.gov/articles/PMC12424641/) expanded supervised
  population-level prediction and interpretable cell selection.
- [singleDeep](https://pmc.ncbi.nlm.nih.gov/articles/PMC11735047/) trained on GSE174188
  and externally evaluated on GSE135779; it is therefore a direct SLE comparator, not
  merely background.
- [PaSCient](https://pubmed.ncbi.nlm.nih.gov/41923638/) is now a final Cell Systems
  paper (2026;17:101570; DOI 10.1016/j.cels.2026.101570), not just a preprint. It
  explicitly learns multicellular patient representations and compares pooling and
  pseudobulk strategies.
- [scPhase](https://doi.org/10.1186/s13073-026-01598-x) combines attention-enhanced
  patient representations with adversarial domain adaptation across five datasets.
- [FloREN](https://www.biorxiv.org/content/10.64898/2026.07.12.738088v1.full) appeared
  on 2026-07-20 and is the closest new collision. It uses both the Nehar-Belaid
  (GSE135779) and Perez (GSE174188) cohorts, supervised patient graph embeddings,
  five-fold cross-validated AUC, pseudobulk/CloudPred/scFeatures/SampleCLR/MrVI
  comparators, and a metric that rewards batch-information removal.

**Consequence:** learned patient representation, SLE patient classification, attention
pooling, and a high within-cohort AUC are not novelty claims.

### Donor-aware and external evaluation

- [Donor-Aware scRNA-seq Benchmarks for IBD Classification](https://arxiv.org/abs/2605.03281)
  uses donor-aware cross-validation, two independent IBD cohorts and cross-dataset
  transfer.
- [singleDeep](https://pmc.ncbi.nlm.nih.gov/articles/PMC11735047/) already uses one of
  our SLE cohorts for training and the other for external validation.

**Consequence:** donor-aware splitting and external transfer are essential strengths,
but neither is individually first.

### Simple baselines can beat foundation models

- [Kedzierska et al.](https://doi.org/10.1186/s13059-025-03574-x) showed that zero-shot
  scGPT and Geneformer embeddings can lose to HVG, scVI or Harmony baselines and that
  Geneformer may retain batch structure.
- [Wu et al.](https://doi.org/10.1186/s13059-025-03781-6) found strong task dependence
  across six single-cell foundation models rather than a universal winner.
- A patient-level cancer benchmark
  ([bioRxiv 2025.10.31.685892](https://www.biorxiv.org/content/10.1101/2025.10.31.685892v1))
  compared nine foundation models with simple baselines and found limited universal
  advantage.
- [Parameter-free representations outperform single-cell foundation models](https://arxiv.org/abs/2602.16696)
  makes the negative-baseline result still less claimable as a discovery.

**Consequence:** "pseudobulk beats a frozen foundation model" is corroborating evidence,
not the paper's main innovation.

### Batch confounding and over-correction

- [Soneson, Gerster and Delorenzi](https://doi.org/10.1371/journal.pone.0100335)
  demonstrated in 2014 that label-batch confounding can strongly bias
  cross-validation and that batch removal need not repair the estimate.
- [Flexible experimental designs for valid scRNA-seq experiments](https://doi.org/10.1038/s41467-020-16905-2)
  formally states that a completely confounded design is non-identifiable.
- [ZINB-WaVE](https://pmc.ncbi.nlm.nih.gov/articles/PMC5773593/) showed that ComBat or
  explicit batch adjustment can remove patient effects when biology and batch are
  confounded.
- [Performance assessment of scRNA-seq normalization](https://pmc.ncbi.nlm.nih.gov/articles/PMC6544759/)
  explicitly warns that regressing on batch without preserving the covariate of
  interest can remove the donor effect.
- [Chaibub Neto 2021](https://proceedings.mlr.press/v139/neto21a.html) treats linear
  residualisation as a potentially unstable confounding adjustment for anticausal
  prediction.
- The [Confounding Index](https://pubmed.ncbi.nlm.nih.gov/32143800/) already provides a
  general supervised-learning diagnostic for misleading covariates.

**Consequence:** the mathematical mechanism is not new. Our novelty is the direct
measurement of how conventional residualisation and design-restricted evaluation
diverge in widely reused donor-level single-cell benchmarks.

## The opening created by the 2026 literature

The newest papers make the proposed audit more, not less, relevant.

1. **FloREN is an immediate case study.** It evaluates GSE135779 and GSE174188 by
   five-fold within-cohort prediction and separately rewards low batch retention. It
   does not report design-only disease AUC, mixed-batch restrictions with matched
   donor controls, or cross-cohort SLE transfer. Our paper tests whether those commonly
   reported quantities certify transportable disease information.
2. **Pretraining-contamination auditing is becoming expected.** The July 2026
   [scContam preprint](https://arxiv.org/abs/2607.20572) argues that benchmark validity,
   not leaderboard rank alone, must accompany scFM evaluation. Our work addresses a
   distinct validity threat: design-label entanglement after the cell embeddings have
   been constructed.
3. **Composition and simple baselines are now accepted as strong comparators.**
   [scECODA](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1) and the
   general scFM benchmarks reduce the editorial value of another representation
   leaderboard. A validity protocol with empirical failure cases has a clearer reason
   to be cited.

## Defensible novelty claims

The following language is aggressive but supportable if the locked two-cohort results
remain positive:

> To our knowledge, this is the first donor-level scRNA-seq study to quantify, on the
> same donors and under one evaluation pipeline, the discrepancy between
> representation residualisation and design-restricted classification under
> label-design collinearity.

> Across two heavily reused SLE benchmark cohorts, we show that a high within-cohort
> AUC can coexist with strong design predictability, that residualisation can
> misattribute disease-aligned signal to the batch, and that neither a frozen
> foundation model nor learned set pooling repairs the resulting generalisation gap.

> We convert these observations into an executable reporting standard: design-only
> AUC, design-restricted performance with matched donor controls, and source-only
> external performance should accompany every donor-level single-cell classifier.

Avoid:

- first patient-level single-cell classifier;
- first donor-aware benchmark;
- first external SLE validation;
- first demonstration that batch confounding biases cross-validation;
- first evidence that pseudobulk can beat a foundation model;
- proof that the retained within-batch signal is causal or purely biological.

## Recommended editorial positioning

**Preferred paper type:** validation protocol / critical benchmark / problem-solving
protocol, not a model-development paper.

**Recommended title:**

> **When design predicts disease: residualisation failure and generalisation gaps
> in donor-level single-cell classifiers**

Alternative, more technical:

> **Residualisation misreads design-entangled disease signal in donor-level
> single-cell classification**

The title should lead with the failure mode and the audit, not with Geneformer,
pseudobulk, SLE, or the name of a software package.

## Release implication

This is a major scientific repositioning. It should not be published as a patch-level
`v1.0.3` unless that tag is explicitly limited to correcting the existing archive.
The reconstructed analysis warrants a new minor or major line (`v1.1.0` prerelease or
`v2.0.0`) after the two-cohort audit, methods lock, manuscript references, author
metadata, and figure-source tables pass the release gate.
