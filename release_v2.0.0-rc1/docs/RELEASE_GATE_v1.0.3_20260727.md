# Release gate: design-validity reconstruction

**Decision date:** 2026-07-27  
**Decision:** **Do not publish this reconstruction as `v1.0.3`.**  
**Recommended line:** prepare `v2.0.0-rc1`, then publish `v2.0.0` after repository
packaging and author-metadata synchronisation.

## Why this is not a patch release

The reconstruction changes the scientific object, not just its presentation:

- the primary question moves from Geneformer-versus-pseudobulk benchmarking to
  validity under design-label entanglement;
- GSE135779 becomes a second design-auditable cohort after restoration of batch,
  collection year, sequencing QC and demographic metadata;
- the central result becomes the disagreement among residualisation, design
  restriction and external transfer;
- the old 5-10-fold adjustment headline is replaced by an identical-pipeline
  two-cohort analysis;
- Figures 1-4, the title, abstract, Results, Discussion and Methods are new;
- the learned-pooling claim is restricted to independent-target evidence.

Calling this `v1.0.3` would hide a major scientific change behind a patch number and
make the Zenodo version history harder to interpret.

## Scientific gate

| Gate | Status | Evidence |
|---|---|---|
| GSE135779 donor mapping | PASS | 44/44 donors; no missing selected metadata |
| Locked repeated CV | PASS | 20 x 5-fold, donor-stratified, source/fold-contained |
| Identical unadjusted/residualised pipeline | PASS | same folds, scaler, HVG/PCA and C grid |
| GSE135779 additional design blocks | PASS | batch, collection year and complete design |
| Matched restriction controls | PASS | GSE174188 n=89/n=66; GSE135779 n=36 |
| External transfer | PASS | source-only target prediction; paired target inference |
| Learned pooling boundary | PASS | independent-target comparisons retained; biased source-CV headline removed |
| Manuscript numeric validation | PASS | `release_validation.json` |
| New Figures 1-4 | PASS | SVG/PDF/PNG rendered and visually inspected |
| Figure 5 terminology | PASS | n=89 labelled dominant wave 4 |
| Literature/innovation review | PASS | primary records checked through 2026-07-26 |

## Headline evidence supporting the aggressive framing

1. Complete measured design predicts disease at AUC **0.953** in GSE174188 and
   **0.952** in GSE135779.
2. The tested representations recover recorded batch identity at AUC up to
   **0.9997**.
3. Design restriction causes only small matched differences:
   GSE174188 **-0.014 to -0.049** in dominant wave 4 and **-0.026 to -0.044**
   in pure wave 4; GSE135779 **-0.012 to +0.003** after excluding all-case B1.
4. Residualisation can report the opposite conclusion. GSE174188 batch
   residualisation removes **0.230-0.268** AUC. In GSE135779, batch
   residualisation changes almost nothing, while complete-design residualisation
   removes **0.354-0.430** and leaves AUC **0.515-0.526**.
5. Cross-batch and cross-cohort transfer, not internal separability, produce the
   reproducible performance loss. Learned pooling does not repair it.

These results support the title:

> **When design predicts disease: residualisation failure and generalisation gaps
> in donor-level single-cell classifiers**

They also support a narrow first claim:

> To our knowledge, this is the first patient-level single-cell study to quantify,
> across two disease cohorts and on the same donor-matched pipeline, the complete
> discrepancy among design-only predictability, representation-to-design
> predictability, residualisation, design restriction with matched controls, and
> source-only transfer.

## Innovation assessment

**Overall novelty: 7.5/10; editorial distinctiveness: 8/10.**

Patient-level classification, donor-aware folds, learned pooling, external SLE
validation, batch-confounding theory and simple baselines outperforming frozen
foundation models are already occupied by CloudPred, singleDeep, PaSCient, scPhase,
FloREN, recent foundation-model benchmarks and classical batch-confounding work.

The citable increment is the validity chain and its empirical contradiction:

- near-perfect design recognition can coexist with robust within-design disease
  discrimination;
- residualisation damage scales with how strongly the projected covariates predict
  the label;
- design restriction and matched controls can preserve the signal that
  residualisation deletes;
- learned pooling does not turn internal separability into external invariance;
- four executable checks convert this observation into a reporting standard.

FloREN is the closest new collision because it uses both GSE135779 and GSE174188,
but it does not report design-only disease AUC, mixed-design restrictions with
matched donor controls, the residualisation/restriction discrepancy or strict SLE
cross-cohort transfer. The new preprint therefore increases the timeliness of this
paper instead of eliminating its novelty.

## Packaging gate closure

The release-candidate packaging blockers identified in the initial decision are
closed:

1. `CITATION.cff`, `.zenodo.json` and `pyproject.toml` name Hongyu Ying,
   Dandan Yun and Dan Liu.
2. The new manuscript, six main figures, locked source tables, portable donor-level
   inputs, scripts and methods manifests are contained in one release directory.
3. README and DOI-facing metadata describe the design-validity reconstruction.
4. A clean-copy rerun reproduced the locked donor-level outputs.
5. Stable supplementary tables and an initial-submission PDF are included.

The remaining final-release requirements are author declarations, immutable tag/DOI
substitution and one last branch review; they do not block the prerelease branch.

## Recommended release sequence

1. Preserve `v1.0.2` and DOI `10.5281/zenodo.21436893` unchanged.
2. Do not mint `v1.0.3` for the reconstructed paper.
3. Create a `v2.0.0-rc1` directory and branch containing the new manuscript,
   figures, source tables, scripts, verified references, metadata restoration and
   SHA256 manifest.
4. Update creator metadata to the manuscript authors and label FloREN/scContam as
   preprints.
5. Run repository tests, manifest verification and one independent package-open
   check.
6. Publish `v2.0.0` and mint a new version DOI under the existing concept DOI.

## Current disposition

**Scientific reconstruction: GO.**  
**Publish as `v1.0.3`: NO-GO.**  
**Prepare `v2.0.0-rc1`: GO; repository packaging and clean-copy reproduction passed.**
