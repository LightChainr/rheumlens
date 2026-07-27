# Identifiability upgrade: evidence, novelty and release decision

Date: 2026-07-27

## Decision

The new analysis supports a scientific reframe, not a patch-level cosmetic update.
Keep `15_design_entanglement_20260726/release/v2.0.0-rc1` frozen. Develop the present
workspace as `v2.0.0-rc2` internally, then release `v2.0.0` only after figure
renumbering, reference verification, manuscript rendering and repository alignment.
Do not publish this work as `v1.0.3`: that version number understates the change in
scientific claim and would confuse the provenance of the locked design analysis.

No additional GPU server is required for the priority analysis. All P0 computations
completed locally.

## What changed scientifically

The predecessor argued that residualisation and restriction can disagree. The new
work explains why and defines the boundary at which the disagreement becomes
unavoidable:

1. The design-adjusted label-information fraction is
   `I_D = 1 - R2(Y ~ D)`.
2. Relative variance inflation is `1 / I_D`.
3. At perfect design-label collinearity, biological-only and technical-only
   mechanisms can induce the same observed `(X, D, Y)` distribution.
4. Therefore, residualisation cannot by itself estimate the technical share of
   disease discrimination. It measures discrimination left after deleting
   design-aligned variation.

This is an estimability result. It is intentionally not presented as a universal
theorem that AUC must decrease monotonically for every nonlinear model.

## New evidence

### Simulation

- 11,200 complete replicates.
- Four biological effects, four technical effects, seven design-label associations,
  and 100 replicates per parameter cell.
- Mean information fraction fell from 0.996 at `rho=0` to 0.189 at `rho=0.90`,
  0.039 at `rho=0.98`, and zero at `rho=1`.
- Median variance inflation reached 5.27 at `rho=0.90` and 30.18 at `rho=0.98`.
- Pure technical example (`beta=0`, `delta=2`, `rho=0.90`):
  internal AUC 0.936; residualised 0.499; restricted 0.501; external 0.499.
- Mixed example (`beta=1`, `delta=1`, `rho=0.90`):
  internal AUC 0.947; residualised 0.573; restricted 0.794; external 0.814.
- Exact observational-equivalence test: maximum difference 0.

### GSE135779 full-pipeline null audit

- Complete recorded design AUC: 0.952.
- 1,000 label permutations repeated the nested pipeline:
  empirical upper-tail `p=0.000999`.
- 250 independent column-shuffle nulls:
  empirical upper-tail `p=0.003984`.
- Complete design left 20.3% of centred label variation after adjustment.
- Relative variance inflation: 4.93.
- Batch alone was null (AUC 0.499); the strong label channel was distributed across
  collection year, sequencing/QC and demographics.
- Removing all-case batch B1 retained 36 donors in mixed-label batches. Observed
  minus matched AUC was -0.012 for Geneformer and +0.003 for both pseudobulk
  representations; exact sign-flip tests were non-significant.

### GSE174188 CD4 composition extension

- Same 261 donors and same CD4 cell set as the molecular benchmark.
- Raw, CLR and cell-yield-augmented composition AUCs: 0.899-0.915.
- Wave-adjusted AUCs: 0.689-0.709.
- Composition predicted processing-wave identity with one-versus-rest AUC up to
  0.872.
- Within dominant or pure wave 4, observed-minus-size/label-matched AUC was
  approximately -0.027 to -0.062 and -0.034 to -0.037, with all 95% split ranges
  crossing zero.
- The residualisation/restriction contradiction is therefore not specific to
  Geneformer or high-dimensional expression.

## Innovation audit

### Prior work that must be acknowledged

- Song, Chan and Wei (Nature Communications, 2020) established that completely
  confounded single-cell batch/cell-type designs are non-identifiable and developed
  valid overlap designs for BUSseq.
- Soneson, Gerster and Delorenzi (PLOS ONE, 2014) showed that batch confounding can
  bias cross-validated genomic classifiers.
- Confounding-index work in biomedical machine learning quantifies how a candidate
  confounder affects classification.
- PaSCient, CloudPred, singleDeep and related methods already perform patient-level
  prediction from single-cell data.
- Halter, Andreatta and Carmona (bioRxiv, 2026) show that CLR cell-type composition
  is a strong baseline for unsupervised patient stratification.
- Recent foundation-model benchmarks already show that simple expression baselines
  can outperform frozen single-cell embeddings.

### Defensible novelty

No located work combines all of the following:

1. a donor-label information fraction tied to study-design estimability;
2. an exact observational-equivalence construction for patient-level disease
   attribution;
3. full-pipeline design-only permutation tests;
4. residualisation versus overlap restriction with size/label-matched controls;
5. source-only external transfer on the same donors and representations;
6. replication across frozen Geneformer, pseudobulk and cell-type composition.

The claim should therefore be:

> We translate experimental-design identifiability into a quantitative and
> executable validity standard for supervised patient-level single-cell
> classification, and show that high-AUC SLE benchmarks can enter a regime in which
> residual adjustment cannot attribute the learned signal.

Do not claim to be the first paper on single-cell non-identifiability, batch
confounding, patient aggregation or composition.

## Innovation assessment

- Conceptual contribution: 8/10
- Empirical contribution: 8/10
- Algorithmic novelty: 4/10
- Reusable methodological value: 9/10
- Desk-editor visibility after reframe: materially stronger than RC1

The principal citation hook is no longer “pseudobulk can beat Geneformer.” It is the
five-part validity protocol and the information fraction that later patient-level
single-cell studies can report.

## Release gates for RC2

1. Renumber the two new figure groups into the main sequence and update all calls.
2. Add the BUSseq, Confounding Index and scECODA citations to the manuscript.
3. Render the manuscript to PDF and inspect equations, tables, references and
   figure legibility.
4. Copy scripts, tests, manifests, source tables and validated SVG/PDF/PNG files into
   a candidate release tree.
5. Re-run the 18 automated checks from the candidate tree.
6. Align README, CITATION, release notes, Zenodo metadata and manuscript numbers.
7. Only then create the GitHub/Zenodo release.

