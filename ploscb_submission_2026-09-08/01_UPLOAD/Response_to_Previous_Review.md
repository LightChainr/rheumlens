# Response to the decision on PCOMPBIOL-D-26-01856

**Previous submission:** *Study design predicts disease and defines the identifiability
boundary for patient-level single-cell classifiers* (rejected after review).
**This submission:** *Recorded metadata predicts the phenotype label in eight of nine
public single-cell cohort comparisons* (new Research Article).

This document lists each comment in the decision letter and in the four reviews, what was
changed or added, and where to find it. Section, figure and table numbers refer to the new
manuscript. Reviewer comments are paraphrased; Reviewers 1 and 4 supplied attachments,
whose points are grouped by topic.

## Editor

| Comment | Response | Location |
|---|---|---|
| The "identifiability boundary" is not theoretically supported; the diagnostic only measures how strongly recorded design variables correlate with the label | Claim withdrawn. The paper no longer states an identifiability result. The primary measure is now the cross-validated metadata-only AUC against the comparison's own permutation null; the residual label variance V_D (formerly 1−R²) is reported beside it as a secondary, explicitly associational summary | Title, Abstract, §2.2, §9.3–9.4 |
| Comparison of adjustment strategies lacks a causal framework | A causal diagram defines what each quantity and adjustment targets; residualisation and restriction are presented as non-equivalent adjustments whose disagreement is informative, not as competing estimators of one effect | §2.1, Figure 1, §4.2, §5.3 |
| "Five required checks" proposed as a general standard from one disease and one cell type | "Required" and "standard" withdrawn. The checks are described as complementary, each with what it shows and does not show. Evidence extended to nine phenotype contrasts from five public datasets (809 donors) across lupus, COVID-19, influenza, cytomegalovirus serostatus and sepsis versus COVID-19 | §3, Table 1, §6 |
| Undefined terms, missing reference numbers, figure errors | All symbols defined at first use; references numbered and verified against Crossref; all figures redrawn; line-numbered PDF | Throughout; References; Figures 1–8, S1–S7 |
| Refine the framework, recalibrate claims, expand validation | See the rows above and the reviewer responses below | — |

## Reviewer 1

| Comment | Response | Location |
|---|---|---|
| The "boundary" is never formally defined; identifiability and estimability are mixed | Boundary language removed. "Estimable / underpowered / not answerable with these data" are defined operationally for the stratified test | §2.4 |
| 1−R² is not an information fraction | Renamed residual label variance, V_D, and described as a linear summary sensitive to the parameterisation of the metadata matrix; cross-fitted and tested against its own null | §2.2, §9.3–9.4 |
| No causal framework; restriction and residualisation target different estimands | Causal diagram added; the two adjustments' assumptions are stated and the comparison is framed as showing that they are not interchangeable | §2.1, Figure 1, §5.3, Table in §6 (check 4) |
| Unmeasured confounding not addressed | Stated as a limitation; the paper measures association with recorded variables only | Abstract (final sentence), §8 |
| Simulation continuous and linear while labels are binary; parameterisation sensitivity | Simulation extended to binary outcomes and six further regimes (nonlinear design effect, heteroscedastic noise, design-by-cell-type interaction, class imbalance with and without threshold offset, different design mechanism in the target cohort); 19,600 simulated cohorts | §4.1–4.3, §9.9, Figure 3, Table S2 |
| Decision-tree thresholds unexplained | Thresholds stated as heuristics; the degeneracy gate is a single structural criterion (no label-mixed collection stratum), and restriction feasibility is a stated rule | §6.1, §9.7, Figure 8 |
| Reference list unnumbered; formulas rendered inconsistently; no line numbers | Numbered references; plain-text formulas; line-numbered PDF | References; PDF |
| Symbols (M₁, ε, ρ, X) undefined; which model produced each AUC; what restriction restricts to | Symbols defined; every fitted pipeline described separately, with a generated hyperparameter table; restriction strata defined | §2.1–2.2, §9.5, §9.14, Table S10 |
| Noise added to X rather than Y | Generating model specified in full | §9.9 |
| Ridge intercept wording | The intercept is fitted and not penalised; wording corrected | §9.13 |
| Author summary: "removed much more signal" | Author summary rewritten | Author summary |

## Reviewer 2

| Comment | Response | Location |
|---|---|---|
| 1. Broad claims, single-disease validation | Nine phenotype contrasts from five datasets and several diseases; claims restricted to association with recorded metadata in public blood datasets. All datasets are blood, which is stated as a limitation | §3, Table 1, §8 |
| 2. Unclear single-cell-specific contribution | The screen operates on donor-level summaries and is stated to apply to any patient-level table. Single-cell-specific analyses are retained: cell-type composition as a representation, pseudobulk versus foundation-model embeddings, and learned pooling | §1, §5.1, §5.4, §9.15–9.16, Table S6 |
| 3. Limited theoretical novelty of 1−R² | No theoretical novelty is claimed for it; it is a secondary summary. The stratified permutation is attributed to Chaibub Neto et al. The contribution is framed as cross-disease empirical evidence | §1, §2.2, §2.3, §7 |
| 4. Checks and thresholds not evaluated as a general standard | Standard claim withdrawn; calibration of both permutation tests measured by simulation; the checks are applied to all nine comparisons and three are traced step by step through the decision tree | §4.4, §6, §6.1, Table S7 |
| 5. Unmeasured factors; persistence after restriction is not biological proof | Both points adopted. A simulation shows the stratified test correctly rejects its null when a recorded variable outside the conditioning set drives both phenotype and expression, so passing it does not certify disease biology | §4.4, §6 (check 2), §8 |
| 6. Reads like a technical report | Restructured into framework (§2), cross-disease screen (§3), simulation interpreting the screen (§4), mechanism in two lupus cohorts (§5), integration (§6), discussion (§7) | Throughout |

## Reviewer 3

| Comment | Response | Location |
|---|---|---|
| 1. Scope does not match general positioning; third lupus cohort used inconsistently | Validation expanded across diseases (see Reviewer 2.1); the role of each lupus cohort is stated | §3, §5, §9.1, Table S8 |
| 2. Abrupt transition to the checks; decision tree never applied step by step; how is the first step determined for GSE174188; why is [residualised, restricted] an attribution range | The tree is traced through three comparisons with the statistic and result behind each node; the routes follow from those results. The "attribution range" is withdrawn: the two adjusted AUCs are not bounds on a biological quantity | §6, §6.1, Figure 8, Table S7, §5.3 |
| 3. Figures 7 and 8 cited only as whole figures; panels undiscussed | Figures reorganised and renumbered; panels are cited where their result is reported | Figures 1–8, S1–S7 |
| 4. Results list AUCs without motivation or synthesis; abstract lists results | Each results section opens with its question and closes with its interpretation; abstract rewritten around the question, the findings and their limit | Abstract, §3–§5 |
| 5. Formatting: unnumbered references, unrendered formulas, wrong figure numbers inside Figures 4 and 5, overlapping annotations | Fixed; figures redrawn from scripts, with a build check for retired labels | References; all figures |

## Reviewer 4

| Comment | Response | Location |
|---|---|---|
| Complete-pipeline label permutation also destroys the design-label association, so it does not test whether the classifier exploits acquisition structure | Adopted as the central methodological change. A collection-stratified permutation (following Chaibub Neto et al.) is run beside the free permutation, with the null of each stated. The two disagree in sepsis versus COVID-19, and the stratified test is uninformative in influenza | §2.3, §3.4, §9.6, Table 1 |
| Representation predicting design does not show the classifier uses it | Adopted; GSE135779 is reported as a counterexample (batch readable from the representation, batch-only phenotype AUC near chance, little change after batch residualisation) | §5.1–5.2, §6 (check 3) |
| I_D measures only one edge of the confounding structure | Stated: metadata predictability is necessary, not sufficient, for the classifier to use metadata. The increment of expression over metadata, AUC(D+X) − AUC(D), is now reported for every comparison | §2.2, §3.5, Table S4 |
| "Required standard" should be "complementary checks" with what each establishes | Adopted | §6 |
| Remaining discrimination after restriction attributed to biology | Removed throughout | §5.3, §5.5, §7 |
| Simulation extensions (nonlinear, heteroscedastic, interaction, imbalance, different target mechanism) | All implemented | §4, §9.9, Table S2 |
| Report metadata-only predictability after restriction and within matched strata | Metadata-only AUC reported by variable group with its own null; restriction and matched donor controls reported together | §3.1–3.2, §5.3, §9.14 |
| Why only random forest and gradient boosting; hyperparameters | Rationale given; every setting generated from the pipeline objects that ran | §9.5, Table S10 |
| Informal writing; sentences such as "Design overlap supports attribution" | Rewritten | Throughout |
| Abstract symbols and GEO accessions unexplained; vague "frozen ... direct workflow" | Abstract contains no undefined symbols; "frozen" reserved for the untrained foundation model | Abstract, §1 |
| Framework schematic | Figure 1 | Figure 1 |
| Missing related work (ProtoCell4P, Liu et al. Cell Systems, scPanel, Wagle et al., HiDDEN) | All cited, with further work on confounding in biomedical machine learning | §1, References |

## Changes beyond the reviews

Before resubmission the analysis code was checked against the text. Three implementation
inconsistencies were found and corrected, and the affected results re-run:

1. The metadata matrix's imputation and one-hot encoding had been fitted on the whole
   cohort. They are now fitted on training donors only. Metadata-only AUCs changed by at
   most 0.012; expression results and all headline counts were unchanged.
2. The calibration simulation had fitted its PCA on all donors, and one arm of the
   residualisation simulation used a differently fitted representation from the other.
   All pipelines now call one shared implementation in which the two arms differ only in
   the adjustment step. Simulation results in §4.4–4.5 and Table S5 are from the
   corrected code.
3. The hyperparameter table described one configuration for pipelines that differ. It is
   now generated from the pipeline objects (Table S10).
