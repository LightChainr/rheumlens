# v2.0.0-rc1 readiness

**Audit date:** 2026-07-27  
**Candidate:** `release/v2.0.0-rc1`

## Decision

The local research-object package is ready to enter a dedicated repository branch.
It is not yet ready for a public Zenodo DOI.

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
- A package-level validator checks required files, figure triplets, manuscript hash,
  stale claims, metadata donor count, creator metadata and symlinks.
- The candidate passes validation with no errors.

## Remaining before public v2.0.0

1. Parameterise scripts that retain original project-relative paths.
2. Run the package from a clean directory against documented upstream inputs.
3. Produce stable supplementary tables and a submission PDF.
4. Confirm author contributions, ORCIDs, funding, ethics and competing interests.
5. Replace future-tense repository language with the final Git tag and Zenodo DOI.
6. Commit the candidate on a dedicated branch, review the diff, then publish a
   prerelease before minting the final version DOI.

## Release boundary

Legacy drafts, atlas figures and scripts containing superseded claims remain outside
`release/v2.0.0-rc1`. They are preserved in the working archive but are not part of
the candidate manifest.
