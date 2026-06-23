# P9 structured source-only transfer extension

Non-neural structured methods were fit on source donors only and scored on target donors before metric calculation. This extension does not include neural MIL models because their source-target training protocol requires a separate pre-registered transfer definition.

## Key findings

- GSE285773 → GSE174188 CD4: expression baselines remained strongest; FOCUS was the best structured scGPT method but did not exceed donor_mean_hvg or donor_expression_pca.
- GSE174188 CD4 → GSE285773: moments_mean_var@scGPT achieved the highest AUC in this extension, with focus_lite@scGPT comparable to donor_expression_pca. The pediatric target has only 26 donors, so uncertainty remains wide.
- KME did not improve transfer performance in either direction.

## Top methods by direction

| direction                  | method_id              |      auc |   auc_ci_low |   auc_ci_high |   pr_auc |    brier |   n_donors |
|:---------------------------|:-----------------------|---------:|-------------:|--------------:|---------:|---------:|-----------:|
| GSE174188_CD4_to_GSE285773 | moments_mean_var@scgpt | 0.9625   |     0.88125  |      1        | 0.977072 | 0.34616  |         26 |
| GSE174188_CD4_to_GSE285773 | focus_lite@scgpt       | 0.93125  |     0.799844 |      1        | 0.970147 | 0.286793 |         26 |
| GSE174188_CD4_to_GSE285773 | donor_expression_pca   | 0.93125  |     0.8      |      1        | 0.96875  | 0.173172 |         26 |
| GSE285773_to_GSE174188_CD4 | donor_mean_hvg         | 0.806397 |     0.754204 |      0.854285 | 0.892524 | 0.313122 |        261 |
| GSE285773_to_GSE174188_CD4 | donor_expression_pca   | 0.790248 |     0.734939 |      0.841941 | 0.896082 | 0.252029 |        261 |
| GSE285773_to_GSE174188_CD4 | focus_lite@scgpt       | 0.754458 |     0.695408 |      0.809828 | 0.85639  | 0.586893 |        261 |

## Full summary

| direction                  | method_id              |      auc |   auc_ci_low |   auc_ci_high |   pr_auc |    brier |   n_donors | status   |
|:---------------------------|:-----------------------|---------:|-------------:|--------------:|---------:|---------:|-----------:|:---------|
| GSE174188_CD4_to_GSE285773 | moments_mean_var@scgpt | 0.9625   |     0.88125  |      1        | 0.977072 | 0.34616  |         26 | SUCCESS  |
| GSE174188_CD4_to_GSE285773 | focus_lite@scgpt       | 0.93125  |     0.799844 |      1        | 0.970147 | 0.286793 |         26 | SUCCESS  |
| GSE174188_CD4_to_GSE285773 | donor_expression_pca   | 0.93125  |     0.8      |      1        | 0.96875  | 0.173172 |         26 | SUCCESS  |
| GSE174188_CD4_to_GSE285773 | donor_mean_hvg         | 0.9      |     0.75     |      1        | 0.952679 | 0.476457 |         26 | SUCCESS  |
| GSE174188_CD4_to_GSE285773 | scgpt_mean             | 0.88125  |     0.7375   |      0.98125  | 0.935056 | 0.615318 |         26 | SUCCESS  |
| GSE174188_CD4_to_GSE285773 | quantiles@scgpt        | 0.8625   |     0.69375  |      0.9875   | 0.936066 | 0.325896 |         26 | SUCCESS  |
| GSE174188_CD4_to_GSE285773 | red@scgpt              | 0.85     |     0.675    |      0.975    | 0.920191 | 0.384615 |         26 | SUCCESS  |
| GSE174188_CD4_to_GSE285773 | kme_multiscale@scgpt   | 0.640625 |     0.475    |      0.815625 | 0.694928 | 0.384557 |         26 | SUCCESS  |
| GSE174188_CD4_to_GSE285773 | tail_fractions@scgpt   | 0.5      |     0.5      |      0.5      | 0.615385 | 0.384615 |         26 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | donor_mean_hvg         | 0.806397 |     0.754204 |      0.854285 | 0.892524 | 0.313122 |        261 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | donor_expression_pca   | 0.790248 |     0.734939 |      0.841941 | 0.896082 | 0.252029 |        261 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | focus_lite@scgpt       | 0.754458 |     0.695408 |      0.809828 | 0.85639  | 0.586893 |        261 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | scgpt_mean             | 0.730453 |     0.669909 |      0.789065 | 0.844889 | 0.443241 |        261 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | red@scgpt              | 0.719541 |     0.651506 |      0.785386 | 0.741501 | 0.263681 |        261 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | kme_multiscale@scgpt   | 0.61348  |     0.539403 |      0.685561 | 0.713529 | 0.48421  |        261 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | moments_mean_var@scgpt | 0.590971 |     0.518451 |      0.660251 | 0.718336 | 0.469312 |        261 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | quantiles@scgpt        | 0.572827 |     0.498748 |      0.644727 | 0.682854 | 0.475081 |        261 | SUCCESS  |
| GSE285773_to_GSE174188_CD4 | tail_fractions@scgpt   | 0.463493 |     0.386424 |      0.540971 | 0.55923  | 0.553233 |        261 | SUCCESS  |

## Failures

No failures.


## Claim boundary

These are retrospective source-only transfer analyses. They support method-prioritization hypotheses but do not establish clinical deployment, causal biology, or query-specific FOCUS mechanisms. Neural MIL transfer remains a future analysis requiring a separate protocol.
