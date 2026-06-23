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
  - Structural check: 146 paragraphs, 7 tables, 12 embedded inline figures.
  - P8.4 corrected values confirmed in DOCX: 0.953, 0.970, 0.948.
  - Old single-smoke P8.4 prose value pattern absent.

## Known limitation

LibreOffice/`soffice` is not installed on the Mac environment used here, so DOCX-to-PNG visual render QA was not performed. Structural DOCX inspection and key numerical text checks passed.

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
3. Optional LibreOffice render QA on a machine with `soffice`.
4. Decide whether to include the lightweight release skeleton as-is or expand it with archival data download scripts.
