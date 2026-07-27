# Draft cover letter

Dear Editors,

We submit the Research Article **“When study design predicts disease: an
identifiability limit for patient-level single-cell classifiers”** for consideration
in *PLOS Computational Biology*.

Patient-level single-cell studies increasingly report high cross-validated
discrimination after compressing thousands of cells into one patient
representation. The unresolved problem is not simply whether study design is encoded,
but whether disease biology and acquisition effects remain separately estimable when
the design itself predicts the outcome.

We make three linked contributions. First, we define a design-adjusted
label-information scale and construct exact observational equivalence between
biological-only and technical-only mechanisms at perfect design-label collinearity.
Across 11,200 simulations, increasing design-label association depleted estimable
outcome contrast, destabilised biological-effect recovery and created
internal-to-external discrimination gaps. A low-entanglement positive control passed
all five proposed checks, showing what a successful cohort looks like rather than
offering only failure diagnostics.

Second, we establish that the empirical lupus cohorts enter this failure regime. In
GSE135779, complete recorded design predicted disease at AUC 0.952, exceeded 1,000
complete-pipeline label permutations (p=0.000999), and left only 20.3% of centred
label variation after linear adjustment. In GSE174188, CD4-subtype composition
reached AUC 0.899-0.915 but fell to 0.689-0.709 after processing-wave
residualisation. Frozen Geneformer and two expression pseudobulk representations
showed the same contradiction: residualisation removed far more discrimination than
matched restriction to observed design overlap.

Two new robustness results make the interpretation sharper. A training-fold
random-forest residualiser attenuated every representation, so the finding does not
depend on a linear residualiser. Conversely, GSE135779 supplied an internal negative
control: representations identified batch at AUC up to 0.993, but batch alone
predicted disease at 0.499 and batch residualisation changed disease AUC by at most
0.009. Encoding design is therefore exposure, not proof that the disease classifier
used a shortcut.

Third, we convert the evidence into an executable five-part validity standard:
quantify design-label estimability, test the full pipeline under null data, report
representation-to-design predictability, pair residualisation with matched overlap
restriction, and require source-only external transfer. Strict external evaluation
favoured expression pseudobulk over frozen Geneformer in the 261-donor target, and
three learned pooling architectures did not repair the gap.

Confound adjustment and confound leakage have been studied extensively in brain
imaging and clinical prediction. To our knowledge, this is the first study to connect
those failure modes to a quantitative estimability scale, an internal negative
control and strict external transfer in patient-level single-cell classification.
The contribution is not another representation leaderboard; it is a criterion for
deciding when any such leaderboard supports biological attribution.

The work is suited to *PLOS Computational Biology* because it connects a formal
estimability result to an immediately executable protocol for high-dimensional
biological prediction. It explains why leakage-free cross-validation and
fold-contained adjustment can still produce scientifically non-attributable results,
and provides diagnostics that apply beyond lupus or any one representation family.

Reproducibility is a primary deliverable. The MIT-licensed repository and review
archive include `REPRODUCE.md`, two locked conda environments, donor-level inputs and
predictions, complete simulation and null draws, figure source tables,
`SHA256SUMS`, `RELEASE_MANIFEST.json`, and an executable package validator. The
public repository is https://github.com/LightChainr/rheumlens and the permanent
Zenodo record is https://doi.org/10.5281/zenodo.20813922.

The manuscript is original, is not under consideration elsewhere, and has been
approved by all authors. The authors declare no competing interests.

Sincerely,

Hongyu Ying, Dandan Yun and Dan Liu  
Department of Rheumatology and Immunology  
Shanghai Pudong Hospital, Fudan University Pudong Medical Center  
Shanghai 201399, China
