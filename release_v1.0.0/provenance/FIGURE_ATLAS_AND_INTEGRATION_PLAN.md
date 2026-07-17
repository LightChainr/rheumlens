# Editable SVG figure atlas and manuscript integration plan

## Editorial structure

The atlas contains eleven editable multi-panel SVG drafts. Seven are proposed as main figures and four as supplementary consolidation figures. This preserves the manuscript's clinical-to-computational narrative while moving historical engineering detail out of the main line.

| Proposed role | SVG draft | Manuscript role | Evidence status |
|---|---|---|---|
| Main Figure 1 | `FIG01_study_design_and_evidence_ladder.svg` | Replaces the current workflow figure | Primary, matched design |
| Main Figure 2 | `FIG02_repeated_within_cohort_benchmark.svg` | Replaces current Figure 2 and absorbs Figure S14 | Primary, repeated donor-level CV |
| Main Figure 3 | `FIG03_cross_cohort_transfer_and_calibration.svg` | Replaces current Figure 3 and absorbs Figures S17/S22 | Primary transfer plus protocol sensitivity |
| Main Figure 4 | `FIG06_pooling_stress_test.svg` | Replaces current Figure 5 | Direct cap500 source-only stress test |
| Main Figure 5 | `FIG05_data_scale_and_regularization_map.svg` | Replaces current Figure 6 and absorbs Figure S18 | Primary budget results plus sensitivity |
| Main Figure 6 | `FIG07_interferon_and_matched_module_attribution.svg` | Promotes and consolidates Figures S7/S8/S20 | Primary repeated IFN analysis; matched-module interpretation |
| Main Figure 7 | `FIG11_pseudobulk_volcano_and_gene_overlap.svg` | Replaces current Figure 7 and Figure S19 | Donor-pseudobulk biological interpretation |
| Supplementary Figure S29 | `FIG04_donor_embedding_geometry.svg` | New geometric interpretation | Label-blind 5090 analysis; descriptive |
| Supplementary Figure S30 | `FIG08_source_only_fusion_sensitivity.svg` | Replaces current Figure 4 or remains supplementary if six main figures are required | Source-only sensitivity |
| Supplementary Figure S31 | `FIG09_validity_controls_and_historical_landscape.svg` | Consolidates Figures S25-S27 | Archived, explicitly exploratory |
| Supplementary Figure S32 | `FIG10_archived_resource_and_aggregation_atlas.svg` | Consolidates Figures S1/S23/S24 | Archived engineering sensitivity |

## Captions

**Figure 1. Study design and evidence ladder for donor-level representation benchmarking.** (A) Donor composition of the three public SLE cohorts. (B) The patient-level representation problem: single-cell measurements are encoded, aggregated to donors, and evaluated against donor disease status. (C) Matched representation families. (D) The analysis hierarchy, spanning repeated internal cross-validation, reciprocal source-only transfer, biological-program attribution, and resource-budget analyses. All supervised evaluations use donor-disjoint splits.

**Figure 2. Repeated within-cohort discrimination of frozen Geneformer and expression pseudobulk.** (A) ROC AUC distributions over 20 repeated stratified five-fold donor splits. (B) Within-split AUC differences for HVG and PCA pseudobulk relative to frozen Geneformer. (C) Mean AUC across repeats. (D) AUC standard deviation across splits. All methods share donor assignments within each repeat; source-fitted feature operations are contained within training donors.

**Figure 3. Strict cross-cohort transfer separates discrimination from probability reliability.** (A) Target-cohort ROC AUC and donor-bootstrap 95% confidence intervals after source-only feature harmonization, scaling, dimensionality reduction, and regularization selection. (B) Sensitivity to sharing source-internal cross-validation folds across methods. (C) Paired target-donor probabilities from Geneformer and PCA pseudobulk in the 261-donor target. (D) Brier score and 10-bin expected calibration error. The shared-fold sensitivity preserves the principal method ordering; calibration is descriptive and no target recalibration is applied.

**Figure 4. Distribution-preserving pooling stress test under reciprocal external transfer.** (A) Target ROC AUC for coordinate mean and seven fixed alternatives at 500 cells per donor. (B) Candidate-minus-mean paired AUC differences. (C) Discrimination-calibration trade-off across target cohorts. (D) Performance relative to the pre-specified bidirectional promotion gate. No alternative provides a stable improvement in both directions.

**Figure 5. Donor sample size, cell budget, and regularization define the representation decision map.** (A) Independent-target learning curves as the small source cohort is repeatedly subsampled. (B) External AUC across a fixed logistic-regularization path for 500- and 1000-cell Geneformer donor means. (C) Matched Geneformer within-cohort AUC at 500 and 1000 cells. (D) Independent repeated scGPT low-cell sensitivity. The analyses show a widening pseudobulk advantage as source donors accumulate and diminishing gains as per-donor cell counts increase.

**Figure 6. Geneformer discrimination is associated with a broad interferon-related expression axis.** (A) Repeated donor-level AUC for Geneformer, the fixed ISG score, IFN-residual Geneformer, and their recombination. (B) Distribution of mean AUC decreases after residualizing expression-, variance-, and coexpression-matched random modules, with the observed IFN decrease marked. (C) Standardized expression of the fixed interferon panel across donors ordered by phenotype and IFN score. (D) IFN and matched-module decreases across repeated split sets. The stable IFN association is not uniquely stronger than other structured expression modules.

**Figure 7. Donor pseudobulk recovers recurrent interferon-associated SLE expression.** (A-C) Cohort-specific SLE-minus-control donor-pseudobulk effects and BH-adjusted Welch-test evidence. Prespecified interferon-associated genes are highlighted. (D) Exact intersections among the 200 largest absolute-effect genes in each cohort. The figure connects the transportability of expression summaries to recurrent, interpretable SLE-associated expression programs.

**Supplementary Figure S29. Geometry of frozen donor-mean Geneformer embeddings.** (A-B) Label-blind PCA coordinates of donor means in the reciprocal-transfer cohorts. (C) Entropy effective rank of donor means compared with sampled cell embeddings. (D) Preservation of the donor geometry between 500- and 1000-cell representations. Mean pooling contracts the cell-level spectrum to a low-dimensional donor subspace while the two cell caps yield nearly coincident donor geometry.

**Supplementary Figure S30. Source-only early-fusion sensitivity.** (A) Target ROC AUC for Geneformer, pseudobulk, and concatenated representations. (B) Paired AUC changes for hybrids relative to their component references. (C) Target-donor score correlations. (D) Discrimination-calibration trade-offs. Simple concatenation does not add stable transportable signal.

**Supplementary Figure S31. Validity controls and historical multicohort landscape.** (A) Archived fold-contained residualization for available technical and demographic covariates. (B) Accepted 1000-permutation intervals for expression PCA and multiscale scGPT kernel embeddings. (C) Kernel effective rank and feature-norm association with cell yield. (D) A traceable but heterogeneous historical multicohort method matrix. Panels A and D combine archived branches and are interpreted as sensitivity evidence, not matched primary comparisons.

**Supplementary Figure S32. Archived resource and aggregation sensitivities.** (A) Observed Geneformer model-size and sequence-length configurations. (B) Repeated scGPT aggregation results. (C) Twenty-two fixed-split aggregation strategies in GSE174188. (D) Geneformer 1000-cell and all-available-cell probes. These analyses define an engineering performance plateau but do not replace the reciprocal external-transfer estimand.

## Text integration

Add the following concise findings to the Results after the strict transfer paragraph:

> A shared-fold source-internal sensitivity analysis preserved the external ordering of the three primary representations. In the 261-donor target, Geneformer, HVG pseudobulk, and PCA pseudobulk achieved AUCs of 0.890, 0.916, and 0.916; both paired pseudobulk comparisons remained significant after BH correction (q=0.029). The sensitivity therefore indicates that the larger-target pseudobulk advantage was not created by method-specific source-fold assignments.

Add the following paragraph after the resource-scaling results:

> Donor-mean Geneformer geometry was markedly lower-dimensional than the sampled cell space. Entropy effective rank was 6.8 and 6.7 for the two donor-mean matrices, compared with 51.4 and 31.8 in the sampled cell spectra. Nevertheless, cap500 and cap1000 donor geometries were nearly identical: centered cosine medians were 0.998 and 0.993, and donor-distance-matrix correlations were 0.999 and 0.997. These findings distinguish a stable cell-budget plateau from evidence that the mean representation retains the full cell-level geometry.

Add one sentence to the Limitations and scope paragraph:

> The low-rank analysis is descriptive and does not establish that spectral contraction itself causes the transfer gap.

## Manual SVG refinement notes

- Keep canvas proportions and panel lettering so the manuscript references remain stable.
- Preserve text as text; do not outline fonts until final export.
- Keep the method colors consistent across all figures: Geneformer blue, HVG gold, PCA coral, hybrid green/purple.
- Retain the “archived” qualifier in Figures S31-S32.
- Export final print files as SVG and PDF; produce a 600-dpi RGB TIFF only if the submission portal rejects vector files.
