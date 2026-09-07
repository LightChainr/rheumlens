# Smoke test, 2026-09-07 (local Mac, before cloud handoff)

The new multi-cohort screen was run against the locked v2 SLE cohort to confirm
it executes end to end and produces numbers consistent with the published
analysis before any cloud compute is spent.

Input: `inputs/design_metadata/sle_gse174188_cd4_donor_covariates.tsv`, derived
from the locked `gse174188_final_donor_covariates.tsv` by adding the log-QC
columns and collapsing `processing_fraction_*` into `batch__processing_wave` -
i.e. the same transformation `locked_validity_audit.py` applies inline. The
cloud run should rebuild this cohort from the h5ad with
`build_donor_level.py` so that every cohort passes through one code path.

Command: `run_design_screen.py --cohorts SLE_GSE174188_CD4 --n-perm 50`

| block | I_D | design AUC (linear) | design AUC (RF) |
|---|---|---|---|
| qc | 0.770 | 0.784 | 0.828 |
| demographic | 0.849 | 0.692 | 0.725 |
| batch | 0.683 | 0.768 | 0.766 |
| **all** | **0.487** | **0.904** | **0.911** |

Disease AUC (pseudobulk, tuned) 0.975; frozen-pipeline AUC 0.963.

Consistency: v2 reported design-only prediction of disease at AUC 0.953 in this
cohort. The screen reaches 0.904 on a design block recomputed from scratch with
different cross-validation settings. Same order, same conclusion, independently
derived - which is the point of the check.

Both permutation p-values sat at the 1/(n+1) floor because only 50 permutations
were run. **This says nothing about whether the two nulls differ.** They can only
separate with the full 1000 permutations, and separating them is the entire
purpose of the design-preserving null. Do not quote these two numbers.

Fixed during the smoke test:
1. Permutation cost - the null originally ran a full nested C search per
   permutation (~135 fits each), which would not finish for 1000 permutations x
   7 contrasts. The null now runs a frozen pipeline (PCA to 50 components,
   C=1.0).
2. A p-value bias - the permutation null was being compared against a *tuned*
   observed AUC. The observed statistic is now recomputed under the identical
   frozen pipeline.
3. `to_numpy()` returns a read-only array under pandas 3.0, breaking median
   imputation.
