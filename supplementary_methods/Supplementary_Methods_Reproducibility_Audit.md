# Supplementary Methods: reproducibility and audit record

Version: 2026-06-23  
Project: RheumLens A800 / Mac final manuscript assembly  
Source report archived at: `/Users/lc/Documents/RheumLens/manuscript/source_reports/RheumLens_A800_full_project_integrated_report_20260623.txt`

This document converts the project-level execution history into a manuscript-facing reproducibility supplement. It is intended to support methods, provenance, and auditability. Operational details that do not affect scientific claims, such as transient server failures or SSH connectivity symptoms, are retained only as internal engineering context and are not part of the manuscript claim set.

## 1. Scope of the audit record

The final manuscript uses donor-level, leakage-resistant analyses across three systemic lupus erythematosus (SLE) single-cell cohorts:

- GSE135779 matched-500 benchmark;
- GSE174188 CD4 benchmark and sensitivity analyses;
- GSE285773 CD4 amendment-fold benchmark and source-only transfer analyses.

The audit record covers:

- software environment validation;
- code package integrity;
- data and model provenance;
- donor-disjoint fold construction;
- out-of-fold prediction integrity;
- fixed, learned, and original method evaluation;
- statistical robustness gates;
- superseded or rejected analyses;
- final claim boundaries.

## 2. Computational environments

The project used separate environments to avoid package conflicts between the core benchmark engine, Geneformer, scGPT, and R-based downstream utilities.

| Environment | Python / runtime | Key packages | Validation status |
|---|---:|---|---|
| `rheumlens-core` | Python 3.11 | NumPy, SciPy, scikit-learn, PyTorch 2.5.1 | Unit tests 11/11 passed |
| `rheumlens-geneformer` | Python 3.10 | Geneformer 0.1.0, transformers 4.46.3 | CUDA available; extraction validated |
| `rheumlens-scgpt` | Python 3.10 | scGPT 0.2.4, GeneVocab | compatibility patch applied |
| `rheumlens-r` | R 4.3.3 | coloc 5.2.3 | installed |

The core `src/rheumlens/` implementation was audited against the original v5 execution package. The source tree contained 59 Python files, and the core engine matched the reference package file-by-file at the time of validation.

The release skeleton assembled on Mac was subsequently tested with Python 3.12:

- `pip install -e '.[dev]'`: passed;
- `pytest -q`: 11/11 tests passed.

## 3. Data assets and provenance

### 3.1 Key frozen assets

| Asset | Role |
|---|---|
| `embeddings/geneformer/v2_20260620T0915/GSE135779_matched500.npz` | validated Geneformer V2 embedding for GSE135779 matched-500 |
| `embeddings/scgpt/GSE174188_CD4_v1/GSE174188_CD4_scgpt.npz` | scGPT embedding for GSE174188 CD4 |
| `splits/authoritative_primary/GSE174188_CD4.csv` | authoritative GSE174188 CD4 folds |
| `splits/authoritative_primary/GSE285773.csv` or amendment equivalent | deterministic amendment folds for GSE285773 |
| `GSE174188_donor_covariates.csv` | authoritative donor-level covariate table used for P8.6 V3.1 |
| `configs/query_bank.locked.yaml` | locked FOCUS query bank; retained for exploratory audit only |

### 3.2 GSE174188 gene-symbol repair

The local processed GSE174188 CD4 expression matrix initially contained numeric feature labels. The source h5ad was inspected on the server, and gene symbols were recovered from the h5ad `var.feature_name` field. The local `lognorm.npz` feature names were repaired from this source, while the numeric-label version was retained as a backup.

After repair:

- GSE174188 and GSE285773 expression matrices shared 18,916 genes;
- canonical interferon-stimulated genes including `ISG15`, `MX1`, `IFI44L`, `OAS1`, and `STAT1` were present;
- the mechanism supplement could be computed at donor-level gene-symbol resolution.

## 4. Fold construction and leakage controls

### 4.1 GSE135779 matched-500

GSE135779 used a 5-fold donor-level cross-validation design:

- 44 donors;
- 5 folds;
- each donor appears in the test set exactly once;
- train/test donor sets are disjoint within each fold;
- out-of-fold predictions are donor-unique.

The matched-500 design used 500 cells per donor, giving 22,000 cells total for the validated embedding analyses.

### 4.2 GSE285773 amendment folds

GSE285773 required deterministic amendment folds because authoritative primary folds were unavailable. The amendment design used stratified 5-fold cross-validation:

- 26 donors;
- 16 SLE, 10 control;
- `StratifiedKFold(n_splits=5, shuffle=True, random_state=285773)`;
- 130 fold rows, corresponding to 26 donors × 5 folds;
- train/test donor sets are disjoint;
- each donor is tested exactly once.

The amendment fold SHA256 reported during execution was `90fe8a0c...`.

### 4.3 GSE174188 CD4 folds

GSE174188 CD4 used authoritative donor-level folds after resolving the split file semantics. The split file contains per-donor fold assignment rows, and `role == "test"` identifies each donor's held-out fold. Final analyses used:

- 261 donors;
- 162 SLE, 99 control;
- five donor-disjoint folds;
- fold sizes approximately 53/52/52/52/52 donors;
- donor-level OOF predictions.

## 5. Geneformer provenance and repair

Geneformer required full-chain provenance repair before final use. Two environment defects were identified and corrected:

| Defect | Correction |
|---|---|
| Geneformer package not installed | installed from local source with editable install |
| incompatible transformers version | pinned/downgraded to transformers 4.46.3 |

Validated extraction parameters:

| Parameter | Value |
|---|---|
| checkpoint | Geneformer V2-104M |
| model depth | 12 layers |
| hidden dimension | 768 |
| source commit | `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5` |
| embedding mode | CLS |
| embedding layer | last layer (`-1`) |
| input | 22,000 cells, 44 donors × 500 cells |
| batch size | 8 |
| transformers | 4.46.3 |

Frozen model asset hashes:

| Asset | SHA256 |
|---|---|
| `model.safetensors` | `fff5cba29ddd8792991fa77b4872246fbe548a178cebda3775cdc72b67780e7f` |
| `config.json` | `467d4492f0dd53b4d60afffe20812db484ca1cf9fdbeb6a6e060e93564f70859` |
| token dictionary | `67c445f4385127adfc48dcc072320cd65d6822829bf27dd38070e6e787bc597f` |
| gene median dictionary | `a51c53f6a771d64508dfaf61529df70e394c53bd20856926117ae5d641a24bf5` |
| Ensembl mapping | `0819bcbd869cfa14279449b037eb9ed1d09a91310e77bd1a19d927465030e95c` |

The final validated Geneformer embedding reproduced the expected donor-level signal:

- `geneformer_mean` AUC = 0.9256;
- expected reference approximately 0.92;
- new embedding SHA256 prefix: `87276833...`;
- the earlier poor-quality embedding with unknown provenance, AUC approximately 0.48, was quarantined and excluded from final claims.

## 6. Benchmark method registry

### 6.1 GSE135779 matched-500

The final P5 matched-500 benchmark completed 27/27 methods across reproduce, fixed, original, and learned stages.

| Method | AUC | Stage |
|---|---:|---|
| donor_mean_hvg | 0.8926 | reproduce |
| donor_expression_pca | 0.9146 | reproduce |
| raw_pseudobulk | 0.9725 | reproduce |
| isg_scalar | 0.8678 | reproduce |
| isg_vector | 0.8760 | reproduce |
| scgpt_mean | 0.8457 | reproduce |
| geneformer_mean | 0.9256 | reproduce |
| moments_mean_var@scgpt | 0.9256 | fixed |
| moments_mean_var@geneformer | 0.9421 | fixed |
| quantiles@scgpt | 0.9091 | fixed |
| quantiles@geneformer | 0.9174 | fixed |
| tail_fractions@scgpt | 0.9091 | fixed |
| tail_fractions@geneformer | 0.9614 | fixed |
| kme_multiscale@scgpt | 0.8733 | fixed |
| kme_multiscale@geneformer | 0.8926 | fixed |
| prototype_histogram@scgpt | 0.9091 | fixed |
| prototype_histogram@geneformer | 0.9421 | fixed |
| cckme_u@scgpt | 0.9174 | originals |
| cckme_u@geneformer | 0.9036 | originals |
| uder_meanvar@scgpt | 0.8044 | originals |
| uder_kme@scgpt | 0.8981 | originals |
| focus_lite@scgpt | 0.9669 | originals |
| focus_lite@geneformer | 0.9614 | originals |
| red@scgpt | 0.9559 | originals |
| red@geneformer | 0.9642 | originals |
| attention_mil@scgpt | 0.8733 | learned |
| attention_mil@geneformer | 0.9642 | learned |

Interpretation boundary: raw pseudobulk was the strongest method in this cohort. `focus_lite@scgpt` was the strongest scGPT-derived embedding method but is treated as exploratory unless additional query-specific controls are completed.

### 6.2 GSE285773 CD4

The GSE285773 benchmark completed 11/11 available methods. Key AUCs:

| Method | AUC |
|---|---:|
| donor_expression_pca | 0.9688 |
| raw_pseudobulk | 0.9563 |
| scgpt_mean | 0.9250 |

Geneformer embeddings were unavailable for this cohort.

### 6.3 GSE174188 CD4

GSE174188 CD4 used 261 donors and 380,275 cells. The benchmark completed the major fixed, learned, and original method groups. Key AUCs included:

| Method group | Method | AUC |
|---|---|---:|
| expression | raw_pseudobulk | 0.9879 |
| expression | donor_expression_pca | 0.9857 |
| expression | donor_mean_hvg | 0.9840 |
| embedding | scgpt_mean | 0.9784 |
| embedding | kme_multiscale@scgpt | 0.9779 |
| learned | attention_mil | 0.9699 |
| original | cckme_u_weighted | 0.9698 |
| original | focus_lite | 0.9579 |
| original | red | 0.9393 |
| original | donorclr | 0.9355 |

The formal manuscript emphasizes that direct expression summaries were highly competitive and often strongest.

## 7. Statistical audit gates

### 7.1 P8.1 donor-level bootstrap

P8.1 completed donor-level stratified bootstrap with 10,000 replicates across the available OOF method outputs.

### 7.2 P8.2 formal permutation

P8.2 completed formal 1,000-replicate permutation testing for key primary methods:

| Method | Empirical P | Status |
|---|---:|---|
| scgPT mean | 0.000999 | completed |
| donor expression PCA | 0.000999 | completed |
| kme_multiscale@scgpt | not used for final superiority claim | KME-specific formal permutation interrupted |

Because KME did not show positive incremental discrimination over scGPT mean in later formal gates, the interrupted KME permutation does not support or block any final KME-superiority claim. The manuscript does not claim KME superiority.

### 7.3 P8.3 repeated cross-validation

P8.3 repeated cross-validation used 30 fold-resampling seeds. Key means:

| Method | Mean AUC | SD |
|---|---:|---:|
| scgpt_mean | 0.9773 | 0.0035 |
| donor_expression_pca | 0.9858 | 0.0032 |
| kme_multiscale@scgpt | 0.9747 | 0.0043 |

The repeated-CV results were consistent with P8.2/P6 estimates.

### 7.4 P8.4 cell-count sensitivity

P8.4 used nested, label-corrected sampling manifests and a V5 runner that invoked registered method objects instead of manually collapsing all methods to mean pooling. This corrected an identity-collapse failure from an earlier runner.

At common-support level 50, the final 10-repeat mean AUCs were:

| Method | Mean AUC |
|---|---:|
| scgpt_mean | 0.953 |
| donor_expression_pca | 0.970 |
| kme_multiscale@scgpt | 0.948 |

Across levels 25, 50, 100, 200, and 500 cells/donor, disease discrimination remained high, and KME did not show a stable positive increment over scGPT mean.

### 7.5 P8.5 surrogate nulls

P8.5 was skipped under the pre-specified trigger logic:

- KME did not exceed scGPT mean in P8.3, P8.4, or the original P8.5 baseline;
- original baseline AUCs were scGPT mean 0.9784 and KME 0.9726;
- baseline ΔAUC (KME − scGPT mean) was −0.0057;
- no KME superiority claim required high-order surrogate null testing.

### 7.6 P8.6 covariate sensitivity

P8.6 V3.1 used the authoritative donor covariate file. The final frozen covariate set included:

- numeric: cells per donor, mean UMI per cell, mean genes per cell, mean mitochondrial percentage, age from development stage, number of samples, number of libraries, number of suspensions;
- categorical: sex and processing cohort;
- excluded: label/disease proxies and sparse ethnicity categories.

Final AUCs:

| Method | AUC | 95% bootstrap CI |
|---|---:|---|
| covariates_only | 0.8465 | [0.798, 0.892] |
| scgpt_mean_covres | 0.8683 | [0.819, 0.911] |
| donor_expression_pca_covres | 0.8533 | [0.801, 0.899] |
| donor_mean_hvg_covres | 0.8698 | [0.822, 0.913] |

Unadjusted-to-adjusted AUC changes:

| Method | Unadjusted AUC | Adjusted AUC | ΔAUC |
|---|---:|---:|---:|
| scgpt_mean | 0.9784 | 0.8683 | −0.110 |
| donor_expression_pca | 0.9857 | 0.8533 | −0.132 |
| donor_mean_hvg | 0.9840 | 0.8698 | −0.114 |

Paired bootstrap residual-vs-covariates-only contrasts all included zero:

| Method | ΔAUC vs covariates_only | 95% CI |
|---|---:|---|
| scgpt_mean_covres | +0.022 | [−0.043, +0.085] |
| donor_expression_pca_covres | +0.007 | [−0.062, +0.074] |
| donor_mean_hvg_covres | +0.023 | [−0.042, +0.087] |

Interpretation boundary: residualization is a sensitivity analysis, not causal proof. It may remove disease-related biology that is correlated with measured technical or demographic variables.

### 7.7 P8.7 KME kernel diagnostics

P8.7 confirmed that the KME identity-collapse blocker was resolved:

- bandwidth stable across cell-count levels, approximately 7.55–7.65;
- exact score equality with scGPT mean: zero exact matches at all audited levels;
- Pearson correlation between KME and scGPT mean predictions: approximately 0.94–0.97;
- effective rank: approximately 2.6–2.8 out of 512 diagnostic points.

The final claim is therefore limited: KME generated distinct predictions without identity collapse, but the low effective rank indicates a concentrated spectrum. The analysis does not claim that the kernel is generally well-conditioned or superior.

### 7.8 P8.8 FOCUS controls

FOCUS results existed in P5 and P6, but the manuscript does not make query-specific or focused-cell-importance claims. Earlier “outperform mean pooling” phrasing was withdrawn. P8.8 was therefore closed as:

`SKIPPED_AFTER_CLAIM_RETRACTION_EXPLORATORY_FOCUS_ONLY`

Deferred controls remain available for future work:

- random-query null;
- focused-cell removal;
- random-removal null.

## 8. Superseded and rejected analyses

Several intermediate analyses were explicitly rejected or superseded and are not used for final claims:

| Item | Final status | Reason |
|---|---|---|
| Unknown-provenance Geneformer NPZ | quarantined | poor AUC and unknown source-row provenance |
| P8.4 pre-V5 runner | superseded | method object ignored; KME collapsed to mean-pooling logic |
| P8.4 smoke single-repeat values | not formal | replaced by 10-repeat formal summaries |
| P8.6 V1 | preliminary only | covariate set incomplete |
| P8.6 V2 | rejected | incomplete metadata source |
| P8.6 V3 | superseded | incorrect observed column in bootstrap report |
| P8.6 V3.1 | accepted | corrected final covariate sensitivity analysis |
| P8.8 full FOCUS control suite | deferred/skipped | no final FOCUS specificity claim |

## 9. Source-only transfer and final integration

Mac-side final P9/P10 integration used the locally staged processed data and repaired GSE174188 gene symbols. Source-only transfer was evaluated in both directions between GSE285773 and GSE174188 CD4.

| Direction | Method | AUC |
|---|---|---:|
| GSE285773 → GSE174188 CD4 | scgpt_mean | 0.730 |
| GSE285773 → GSE174188 CD4 | donor_expression_pca | 0.790 |
| GSE285773 → GSE174188 CD4 | donor_mean_hvg | 0.806 |
| GSE174188 CD4 → GSE285773 | scgpt_mean | 0.881 |
| GSE174188 CD4 → GSE285773 | donor_expression_pca | 0.931 |
| GSE174188 CD4 → GSE285773 | donor_mean_hvg | 0.900 |

The transfer analysis supports cross-cohort signal but not clinical deployment, prospective generalization, or causal interpretation.

## 10. Manuscript-facing claim boundary

### Supported claims

The final data support the following bounded claims:

- donor-level SLE signal is strong across the benchmarked cohorts;
- expression/pseudobulk baselines are highly competitive and often strongest;
- scGPT and Geneformer embeddings contain donor-level disease signal;
- Geneformer results require strict provenance control;
- source-only transfer between pediatric and adult CD4 cohorts is detectable;
- measured covariates account for a substantial portion of GSE174188 discrimination;
- residualized embedding/expression scores are not significantly above covariates-only in the final P8.6 V3.1 sensitivity analysis;
- GSE174188 CD4 signal is biologically anchored by a strong interferon-stimulated gene program;
- KME identity collapse was corrected, but KME superiority is not established.

### Unsupported or prohibited claims

The final manuscript should not claim:

- clinical deployment readiness;
- causal disease mechanism;
- foundation-model superiority over expression baselines;
- KME superiority;
- FOCUS query specificity;
- focused cells are disease-causal;
- residualized signal is biologically independent;
- AUC attenuation directly equals “percentage of signal explained”;
- KME is well-conditioned solely because identity collapse was resolved.

## 11. Internal engineering notes

Early execution was affected by shared-storage and SSH stability problems on multiple GPU servers. These failures influenced execution strategy but not final scientific claims. The relevant methodological lessons were:

- avoid high-worker concurrent reads from shared storage;
- stage large matrices to local SSD before repeated resampling;
- prefer registered method execution paths over hand-written method-specific shortcuts;
- maintain fail-closed assertions for donor IDs, fold maps, label sources, and OOF uniqueness.

These engineering notes are useful for reproduction and future workflow design but should remain outside the main manuscript narrative.

