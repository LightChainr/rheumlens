# Methods lock

**Lock date:** 2026-07-26  
**Scope:** design-entanglement manuscript candidate  
**Executable sources:** `scripts/gse135779_metadata_audit.py`,
`scripts/locked_validity_audit.py`,
`scripts/augment_gse135779_design_residualisation.py`,
`scripts/design_entanglement.py`,
`scripts/pure_wave_stratum.py`, `scripts/deconfounded_transfer.py`, and
`../14_learned_pooling_20260726/scripts/`.

This document distinguishes the frozen cell-embedding generation environment from the
local donor-level reanalysis environment. No target-cohort label or target-fitted
statistic enters a source-only transfer model.

## 1. Cohorts and statistical unit

| Cohort | Analysis donors | Cases | Controls | Analysis compartment |
|---|---:|---:|---:|---|
| GSE135779 | 44 | 33 | 11 | childhood cohort |
| GSE174188 | 261 | 162 | 99 | CD4-positive alpha-beta T cells |
| GSE285773 | 26 | 16 | 10 | CD4-positive T cells |

The donor is the unit of inference. All cells from one donor remain together in every
training/test split. Cell-level observations are not treated as independent replicates.

## 2. Frozen Geneformer representation

### 2.1 Checkpoint and assets

- Repository: `ctheodoris/Geneformer`
- Repository revision: `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5`
- Checkpoint: `Geneformer-V2-316M`
- `model.safetensors` SHA256:
  `965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3`
- `config.json` SHA256:
  `2cc4af3442644e84af71814a607b61c958d7b4e5ebde6791ab77b5f534ac6f6e`
- `token_dictionary_gc104M.pkl` SHA256:
  `67c445f4385127adfc48dcc072320cd65d6822829bf27dd38070e6e787bc597f`
- `gene_median_dictionary_gc104M.pkl` SHA256:
  `a51c53f6a771d64508dfaf61529df70e394c53bd20856926117ae5d641a24bf5`
- Archived embedding environment: Python 3.10, PyTorch 2.5.1,
  `transformers` 4.46.3.

The 316M checkpoint asset hashes are preserved in
`extra/00_registry/checkpoints.tsv`. The earlier integrated report's 104M checkpoint
hash must not be cited for the 316M primary representation.

### 2.2 Token construction

The executable extractor is
`06_remote_patientrepbench_scripts/12_run_geneformer_cell_sample_pilot.py`.

1. Ensembl identifiers were version-trimmed. If a feature was supplied only as a gene
   symbol, it was mapped through the archived symbol-to-Ensembl dictionary.
2. Duplicate mapped Ensembl identifiers were removed after first occurrence.
3. For each cell and retained gene, the ranking value was
   `raw_count / total_cell_count * 10,000 / gc104M_gene_median`.
4. Genes were sorted in descending ranking value.
5. The maximum sequence length was 4096. V2 start and end token IDs 2 and 3 occupied
   two positions, leaving at most 4094 ranked genes. Padding token ID 0 was masked.
6. No task label was used during tokenisation or embedding.

### 2.3 Hidden layer and pooling

The checkpoint was loaded with `transformers.AutoModel`. The extractor used
`out.last_hidden_state`, i.e. the final hidden layer (`emb_layer=-1` in the archived
run report), and selected token position 0 as the cell CLS vector. The output dimension
was 1152. No intermediate layer, decoder output, fine-tuning, or task-specific
projection was used.

Cell embeddings were accumulated in float64 and divided by the sampled-cell count;
stored donor means were float32. The resulting representation is therefore a
coordinate-wise donor mean of frozen final-layer cell CLS embeddings.

### 2.4 Cell sampling

- Sampling seed: integer 1.
- GSE135779 primary design audit: cap 500 cells per donor (22,000 cells).
- GSE174188 primary design audit: cap 1,000 cells per donor (227,504 cells available
  under the cap).
- GSE285773 primary cross-cohort reference: cap 1,000 cells per donor (26,000 cells).

Cells were sampled without replacement with NumPy `default_rng(1)` after deterministic
donor ordering. The cap is per donor, not per class or per fold. The cell samples and
frozen donor vectors were created before the donor-level classifier and are reused
unchanged across donor folds.

## 3. Pseudobulk representations

The stored input is donor-by-gene `log1p(CPM)`. For each outer training fold:

1. gene variance was computed on training donors only;
2. the 4,000 highest-variance genes were selected, with stable descending ordering;
3. each selected coordinate was standardised using training-donor mean and standard
   deviation;
4. mean-HVG logistic regression used these standardised values directly;
5. PCA pseudobulk fitted up to 30 components on the same training values and applied
   the fitted transformation to held-out donors.

For cross-cohort transfer, identifiers were stripped and intersected exactly before
source-only HVG selection. Gene means, variances, HVGs, scaler, PCA, classifier and
regularisation were all source-fitted and carried unchanged to the target.

## 4. Internal donor-level evaluation

- Outer design: 20 repeated stratified five-fold donor splits.
- Integer split seeds: 20260801 through 20260820.
- Classifier: `LogisticRegression`, `solver="liblinear"`,
  `class_weight="balanced"`, `max_iter=20,000`.
- Candidate inverse regularisation strengths:
  `10^-4, 10^-3, ..., 10^4`.
- Tie rule: first maximum, hence the smaller `C`.

For each outer fold, representation construction and standardisation were fitted on
the entire outer training set and never used outer test donors. A five-fold
stratified resampling inside that outer training set selected only `C`, conditional on
the outer-training representation. This is an honest outer out-of-fold estimate; it is
not claimed to refit the unsupervised representation separately inside every
inner-`C` split.

Out-of-fold predictions were pooled across all donors within a repeat. The primary
summary is the mean AUC and 2.5th–97.5th percentile of the 20 repeat-level AUCs. That
percentile range is a split-sensitivity interval, not a sampling confidence interval.

## 5. Design metadata

### 5.1 GSE174188

`Processing_Cohort` was read per cell from the distributed CELLxGENE object. Per-donor
wave fractions were computed over the analysed CD4 cells. Dominant wave is the
maximum-fraction wave; a pure-wave donor has maximum fraction at least 0.99.

### 5.2 GSE135779

The childhood-cohort metadata were restored without missing mappings:

- batch, collection year, demographics and clinical fields: original Supplementary
  Table 1b;
- sequencing QC: original Supplementary Table 1c;
- study sample name to analysis donor ID: GEO title/accession mapping.

The executable restoration script writes all 44 rows, cross-tabs and SHA256 provenance
to `results/gse135779_metadata/`. Batch is not inferred from expression.

## 6. Covariate-only models

Within each outer training fold:

- numeric variables were median-imputed and standardised;
- categorical variables were most-frequent-imputed and one-hot encoded with unknown
  held-out levels ignored;
- `C` was selected from the same grid;
- the fitted transform and classifier were applied to held-out donors.

GSE135779 blocks are batch, collection year, sequencing QC, demographics, and all
measured variables. GSE174188 blocks are processing-wave fractions, QC, demographics,
and their union.

## 7. Fold-contained residualisation

Residualisation and the unadjusted comparator now share the same outer splits,
representations, scaling, HVG/PCA construction, `C` grid and classifier.

For each outer fold:

1. the design matrix was fitted on outer training donors;
2. a multivariate ridge regression (`alpha=1`, intercept included) mapped that design
   matrix to every raw representation coordinate using training donors only;
3. training and held-out coordinates were replaced by observed minus predicted values;
4. the identical representation-specific feature construction and classifier pipeline
   was fitted to the training residuals and applied to held-out residuals.

The primary residualisation block in both cohorts was recorded batch. GSE135779
additionally used collection year and the complete measured design block
(batch, collection year, sequencing QC, age, sex, race and ethnicity) in
prespecified sensitivity analyses. The disease label is not included in the
residualiser. Under label-design
collinearity, this procedure deliberately cannot distinguish disease-aligned biology
from design-associated technical signal. Its AUC change is a sensitivity estimand, not
a causal batch-effect estimate.

## 8. Design-restricted and matched analyses

GSE174188 analyses restrict to dominant wave 4 and to pure wave 4. GSE135779 excludes
the all-case B1 batch and analyses B2–B6, in which every retained batch contains both
classes. The full internal pipeline is rerun after restriction.

Each observed subset is compared with 20 donor subsets drawn without replacement from
the full cohort at the identical donor count and case/control count. Each matched draw
uses one prespecified split seed. These controls estimate the performance loss expected
from reduced donor number and class composition; they are not formal equivalence tests.

## 9. Batch predictability

Each representation predicts each observed batch one-versus-rest under repeated
donor-level cross-validation. The number of folds is the smaller of five and the
minority class count. These AUCs quantify retained batch information, not its technical
or biological cause.

## 10. Cross-batch and cross-cohort transfer

Cross-batch analyses train on complete source waves and evaluate untouched held-out
waves. Cross-cohort transfer fits every statistic on the source cohort only. Target
labels are used only after prediction for evaluation.

Target-donor AUC intervals use 5,000 case/control-stratified bootstrap resamples and
the percentile method. They are conditional on one fitted source model and do not
include source-training uncertainty. Paired AUC differences use the same target donors;
paired DeLong tests are adjusted by Benjamini–Hochberg within the prespecified
comparison family.

## 11. Learned pooling

DeepSets, gated-attention MIL and pooling by multi-head attention operate on the same
cap-500 frozen cell embeddings. Cell standardisation, the whitened source PCA32
projection, network weights and hyperparameters are source-fitted. Candidate hidden
widths are 16 and 32; weight decay is 0.01 or 0.1; Adam learning rate is 0.02 for 150
epochs. Five-fold source-donor cross-validation selects the configuration. Three
initialisations are refitted on all source donors and target probabilities are
averaged.

The comparison includes an ordinary donor mean in the identical source PCA32 space so
that any gain is attributable to learned pooling rather than the projection.

The PCA32 projection used during hyperparameter selection was fitted once on all
source-cohort cells before the source-donor folds. It never used target cells or target
labels, so the independent-target comparison remains source-only. However, its
source-CV AUC is conditional on a source-wide unsupervised projection and is not used
as an unbiased estimate of source-cohort generalisation. The manuscript and Figure 6
therefore report independent-target AUC and paired target-donor comparisons rather
than a numerical "source CV minus target" claim.

## 12. Local reanalysis environment

The locked validity audit was executed locally on macOS arm64 with:

- Python 3.13.7
- NumPy 2.4.6
- pandas 3.0.3
- SciPy 1.18.0
- scikit-learn 1.9.0
- PyArrow 21.0.0
- joblib 1.5.3

Input and output hashes are recorded in the audit manifests. Package-version
sensitivity is bounded by deterministic split seeds, explicit tie handling and
rank-based AUC; the final release should still include this exact environment lock.
