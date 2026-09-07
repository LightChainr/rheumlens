# v2.0.0-rc1 readiness

**Audit date:** 2026-07-27  
**Candidate:** `release/v2.0.0-rc1`

## Decision

The research-object package is published for review on the dedicated
`release/v2.0.0-rc1` repository branch. It is not yet ready for a final public Zenodo
version DOI because author declarations and the immutable release tag remain open.

## Completed

- The locked manuscript is copied without numeric changes.
- Figures 1–4 use the locked two-cohort validity figures.
- Figures 5–6 use the corrected generalisation and learned-pooling figures.
- Restored GSE135779 metadata, repeated OOF predictions, transfer predictions,
  learned-pooling predictions and source tables are included.
- The methods lock, literature audit, technical audit, verified references, release
  gate and 2026 journal strategy are included.
- Creator metadata names Hongyu Ying, Dandan Yun and Dan Liu.
- Exact analysis and Geneformer environment specifications are included.
- Donor-level pseudobulk, primary Geneformer embeddings, labels and compact design
  metadata are included for portable reruns of the core audit and strict transfer.
- Release-facing scripts use package-relative paths with environment-variable
  overrides; no machine-specific absolute paths remain.
- A package-level validator checks required files, figure triplets, manuscript hash,
  stale claims, metadata donor count, creator metadata and symlinks.
- The candidate passes validation with no errors.
- A clean-copy rerun reproduced the restored GSE135779 metadata, strict source-only
  transfer, all 540 locked repeat metrics, summary tables and decompressed OOF
  predictions.
- A 15-sheet stable supplementary workbook was generated and every worksheet was
  rendered for visual inspection.
- A 16-page initial-submission PDF containing the complete manuscript, Author Summary,
  six main figures and full legends was rendered and visually inspected.

## Remaining before public v2.0.0

1. Confirm author contributions, ORCIDs, funding, ethics and competing interests.
2. Replace future-tense repository language with the final Git tag and Zenodo DOI.
3. Review the dedicated branch diff, then publish a
   prerelease before minting the final version DOI.

## Release boundary

Legacy drafts, atlas figures and scripts containing superseded claims remain outside
`release/v2.0.0-rc1`. They are preserved in the working archive but are not part of
the candidate manifest.
