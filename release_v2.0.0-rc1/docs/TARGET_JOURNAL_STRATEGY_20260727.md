# 2026 target-journal strategy

**Decision date:** 2026-07-27  
**Scientific line assessed:** design-label entanglement, residualisation failure,
and external generalisation of donor-level single-cell classifiers.

## Executive decision

Use the following submission sequence:

1. **PLOS Computational Biology — Research Article**
2. **GigaScience — Research**
3. **Bioinformatics Advances — Original Article**
4. **Applied Sciences — Special Issue: Research on Computational Biology and
   Bioinformatics**

This sequence is deliberately ambitious without making a structurally incompatible
journal the first target. The manuscript should be written for the PLOS
Computational Biology readership from the outset: the central contribution is a
general validity result for patient-level computational biology, demonstrated with
single-cell foundation-model and pseudobulk representations in SLE.

The public scientific release should remain journal-neutral. Journal-specific
cover letters, title pages and summaries belong under `submission/`, not in the
archived analysis core.

## Decision matrix

| Journal | Scope fit | Editorial increment required | Present fit | Principal risk | Decision |
|---|---:|---|---:|---|---|
| PLOS Computational Biology | 9/10 | Significant computational or biological insight; rigorous evidence | 8/10 | Editors may regard two diseases-cohort audits as too narrow unless the validity standard is made broadly executable | **Primary target** |
| GigaScience | 9/10 | Reproducibility, usability and utility of open data/code/workflows | 8.5/10 | Requires a genuinely self-contained, curator-friendly research object | **First transfer** |
| Bioinformatics Advances | 9/10 | Original bioinformatics contribution or biological insight | 8/10 | Eight-page main-text constraint will require sharp compression | **Second transfer** |
| Applied Sciences special issue | 9/10 | Original computational biology, applied analysis or benchmarking | 9/10 | Lower editorial selectivity and less concentrated computational-biology readership | **Stable fallback** |
| Genome Biology | 8/10 | Field-level genomic or methodological advance | 5.5/10 | Likely desk rejection without multi-disease replication or a new correction method | Stretch only |
| Patterns | 7/10 | Major advance of broad interdisciplinary data-science interest | 5/10 | Current evidence is domain-specific and does not introduce a general algorithm | Do not lead with it |
| Cell Reports Methods | 6/10 | Significant new, robust and reproducible method or tool | 4/10 | The paper audits methods but does not introduce a new method | Do not submit in current form |
| Communications Biology | 6/10 | New biological insight that changes thinking in a specialist field | 5/10 | Primary contribution is methodological validity rather than SLE biology | Not preferred |
| Briefings in Bioinformatics | 2/10 | Reviews, protocols and case studies; not pure original research | 0/10 | Article type is incompatible | **Remove from list** |

## Why PLOS Computational Biology is the lead target

The journal explicitly covers computational work from cells to patient populations,
AI and machine learning, biomarker analysis, and benchmarking. Its Research Article
criteria emphasize originality, field importance, methodological or biological
insight, rigor and substantial evidence. It imposes no fixed limit on manuscript
length, figure count or supporting information at initial submission.

The paper should therefore be submitted as a **Research Article**, not a Methods
Article. The work does not introduce a new classifier. Its advance is the empirical
and operational result that patient-level single-cell representations can separate
disease internally while encoding study design even more strongly, and that
residualisation, restriction and external transfer answer non-equivalent validity
questions.

### PLOS-facing claim

> Patient-level single-cell classifiers can achieve high donor-level discrimination
> while learning a study design that predicts disease almost perfectly. Across two
> independently reconstructed SLE cohorts, restriction preserves signal that
> residualisation can erase, whereas external transfer reveals the reproducible
> generalisation gap. These contradictions define four executable checks for
> validating patient-level single-cell models.

### PLOS submission additions

- Add a 150–200 word non-technical Author Summary.
- Add a short title and 5–8 keywords.
- Retain Figures 1–6 in the main manuscript.
- Put the four-check validity standard in the abstract, final Introduction
  paragraph, Results 2.7 and first Discussion paragraph.
- State that SLE is the empirical system and patient-level single-cell inference is
  the general methodological object.
- Use a cover letter that explains why the contradiction among residualisation,
  restriction and external transfer is a field-level computational insight.
- Submit the locked source tables, donor-fold definitions and release manifest as
  stable supporting files.

## Transfer logic

### GigaScience

GigaScience evaluates open biomedical research by reproducibility, usability and
utility rather than a subjective impact threshold. This project is well aligned
provided the release contains all executable scripts, compact source tables,
environment information, checksums and stable repository identifiers. A submission
here should emphasize the reusable audit protocol and complete research object.

### Bioinformatics Advances

The journal accepts Original Articles in disease bioinformatics, translational
medicine and immunoinformatics. It is a strong scientific fit, but the current
manuscript must be compressed to eight pages, with no more than six main figures
and six tables. The six-figure design already matches this ceiling.

### Applied Sciences special issue

The active special issue **Research on Computational Biology and Bioinformatics**
explicitly welcomes benchmarking studies, single-cell analysis and applied machine
learning. Its current submission deadline is **20 March 2027**. This is the stable
fallback if the more selective computational-biology journals decline the paper.
The Applied Sciences version should emphasize the executable evaluation protocol,
engineering reproducibility and deployment consequences.

## Journals removed from the active route

### Briefings in Bioinformatics

The current author instructions state that Briefings in Bioinformatics is a review
journal and does not publish pure original research. Although benchmarking and
independent tool evaluation are discussed in its aims, the accepted article types
do not match this manuscript. It must not be used as the lead target.

### Cell Reports Methods

The journal's primary criterion is a significant, robust and reproducible new
method or tool. The current paper intentionally does not introduce a new pooling or
correction algorithm. Submission would invite a predictable desk rejection unless
a validated method is added.

## Stable positioning across all targets

Keep the title and central result aggressive. Do not dilute the paper back into a
Geneformer-versus-pseudobulk comparison.

Use:

> **When design predicts disease: residualisation failure and generalisation gaps
> in donor-level single-cell classifiers**

Lead with three increments:

1. design variables alone predict disease in two reused SLE cohorts;
2. residualisation, restriction and external transfer give empirically
   contradictory answers because they target different validity questions;
3. a four-check audit converts this contradiction into an executable standard for
   patient-level single-cell studies.

Confine limitations to one dedicated section. Elsewhere, state scope precisely
without repeatedly apologising for frozen embeddings, mean pooling, SLE or
retrospective public data.

## Official sources checked

- PLOS Computational Biology journal scope and criteria:
  https://journals.plos.org/ploscompbiol/s/journal-information
- PLOS Computational Biology submission guidelines:
  https://journals.plos.org/ploscompbiol/s/submission-guidelines
- PLOS Computational Biology Benchmarking Collection:
  https://collections.plos.org/collection/benchmarking/
- GigaScience instructions and aims:
  https://academic.oup.com/gigascience/pages/instructions_to_authors
- Bioinformatics Advances author guidelines:
  https://academic.oup.com/bioinformaticsadvances/pages/author-guidelines
- Applied Sciences special issue:
  https://www.mdpi.com/journal/applsci/special_issues/SSY60TBE8W
- Genome Biology aims and scope:
  https://link.springer.com/journal/13059/aims-and-scope
- Communications Biology aims and scope:
  https://www.nature.com/commsbio/aims
- Cell Reports Methods aims and scope:
  https://www.sciencedirect.com/journal/cell-reports-methods
- Patterns aims and scope:
  https://www.sciencedirect.com/journal/patterns
- Briefings in Bioinformatics author instructions:
  https://academic.oup.com/bib/pages/msprep_submission
