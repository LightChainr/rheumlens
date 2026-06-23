# RheumLens manuscript v0.6 finalization status

Date: 2026-06-23

## Current deliverables

- Working manuscript Markdown: `/Users/lc/Documents/RheumLens/manuscript/RheumLens_manuscript_working_v0.6.md`
- Working manuscript DOCX: `/Users/lc/Documents/RheumLens/manuscript/RheumLens_manuscript_working_v0.6.docx`
- Figures: `/Users/lc/Documents/RheumLens/manuscript/figures_v0.6/`
- Supplementary tables: `/Users/lc/Documents/RheumLens/manuscript/supplementary_tables_v0.6/`
- Open-release skeleton: `/Users/lc/Documents/RheumLens/release/rheumlens_open_repo_v0/`

## Verification completed

- Numeric audit after P8.4 correction: PASS 24/24, CHECK 0/24.
  - `/Users/lc/Documents/RheumLens/manuscript/v0.6_audit/NUMERIC_AUDIT_v0.6_AFTER_FIX.md`
  - `/Users/lc/Documents/RheumLens/manuscript/v0.6_audit/NUMERIC_AUDIT_v0.6_AFTER_FIX.tsv`
- Release repository install/test:
  - Python: `/Users/lc/Documents/RheumLens/.venv-rheumlens-release-py312/bin/python`
  - `pip install -e '.[dev]'`: PASS
  - `pytest -q`: PASS, 11/11 tests
- DOCX rebuild:
  - Rebuilt from current v0.6 Markdown.
  - Structural check after package update: 181 paragraphs, 7 tables, 16 embedded inline figures.
  - DOCX render QA completed: 27 rendered pages, no obvious figure/table clipping or overlap on contact-sheet inspection.
  - P8.4 corrected values confirmed in DOCX: 0.953, 0.970, 0.948.
  - Old single-smoke P8.4 prose value pattern absent.

## 2026-06-23 extension analyses added

Reviewer-facing strengthening analyses were added under:

- `/Users/lc/Documents/RheumLens/manuscript/extension_results_20260623/`

New figures:

- Figure 13: covariate decomposition and processing-cohort encoding sensitivity.
- Figure 14: paired donor-bootstrap ΔAUC versus scGPT mean with Holm correction.
- Figure 15: GSE135779 matched-500 scGPT mean underperformance context.
- Figure 16: P9 structured source-only transfer extension.

New supplementary tables:

- Supplementary Table S15: GSE174188 covariate decomposition.
- Supplementary Table S16: paired ΔAUC / Holm correction.
- Supplementary Table S17: GSE135779 matched-500 method ranking around scGPT mean.
- Supplementary Table S18: P9 structured transfer extension.

Key added findings:

- GSE174188 formal-compatible covariate decomposition: full model AUC 0.848, technical depth/QC AUC 0.776, processing AUC 0.685, demographic AUC 0.653.
- One-hot processing-cohort sensitivity increased covariates-only AUC to 0.930; this is a sensitivity warning, not the formal P8.6 claim.
- Paired bootstrap with Holm correction: raw pseudobulk and focus_lite@scGPT significantly exceeded scGPT mean in GSE135779; GSE174188 expression baselines were numerically close to scGPT mean but not significantly above it after correction.
- P9 structured transfer: GSE285773→GSE174188 remained strongest for expression baselines; GSE174188→GSE285773 had moments_mean_var@scGPT AUC 0.963, with wide uncertainty due to 26-donor target.

## 2026-06-23 package completion update

Three publication-support packages were completed:

1. Data/code availability package:
   - Added `environment.yml`, `Dockerfile`, `Singularity.def`, `CITATION.cff`, `.zenodo.json`.
   - Added `scripts/verify_manifests.sh`.
   - Added `docs/DATA_AND_CODE_AVAILABILITY.md` and `docs/ARCHIVAL_DATA_MANIFEST_TEMPLATE.tsv`.
   - Replaced local-path data availability text with public repository and archival DOI language.
2. Reference and positioning package:
   - Added Related work / contribution framing to the Introduction.
   - Expanded manuscript references.
   - Added `references.bib`.
   - Added `docs/REFERENCE_POSITIONING_CHECKLIST.md`.
3. Figure/table package:
   - Confirmed 16 main figures with PNG and PDF versions.
   - Confirmed 6 embedded main tables and S1-S18 supplementary tables plus mechanism tables.
   - Added `docs/FIGURE_AND_TABLE_INDEX.md`.

Verification after package completion:

- `bash scripts/verify_manifests.sh`: PASS.
- `PYTHON_BIN=/Users/lc/Documents/RheumLens/.venv-rheumlens-release-py312/bin/python bash scripts/reproduce_minimal.sh`: PASS.
- Numeric audit: PASS=24, CHECK=0.
- DOCX render QA: PASS by contact-sheet inspection.

## Scientific status

The current v0.6 draft is ready for substantive manuscript editing and external scientific review. It supports the following bounded claims:

- strong donor-level SLE signal across matched-500 and GSE174188 CD4 analyses;
- source-only transfer signal between pediatric and adult CD4 cohorts;
- expression/pseudobulk baselines are highly competitive and often strongest;
- scGPT/Geneformer embeddings contain useful donor-level signal but do not establish foundation-model superiority;
- measured covariates explain a substantial portion of observed discrimination in GSE174188;
- GSE174188 CD4 signal is biologically anchored by a strong interferon-stimulated gene program;
- KME implementation identity collapse was resolved, while effective-rank diagnostics require a low-rank caveat.

## Next recommended work

1. Human scientific edit of Introduction/Discussion for target journal.
2. Formal citation cleanup with DOI/PMID completion.
3. Fill author, affiliation, correspondence, ethics, funding and conflict-of-interest fields.
4. Mint Zenodo/figshare DOI values and replace all `TBD` DOI placeholders.
5. Fill final SHA256 and file sizes in the archival data manifest after deposition.
