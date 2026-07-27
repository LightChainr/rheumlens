# Independent clean-copy reproduction

**Run date:** 2026-07-27  
**Candidate:** `v2.0.0-rc1`  
**Disposition:** PASS

The candidate was copied to a new directory outside the repository, excluding
transient Python bytecode. No path in the original project workspace was available
through the release scripts. The clean copy then restored GSE135779 metadata,
recomputed strict source-only transfer, reran the complete locked two-cohort audit,
added the prespecified GSE135779 design blocks and executed both scientific and
package validators.

## Commands

```bash
python3 scripts/design_validity/gse135779_metadata_audit.py
python3 scripts/strict_source_only_transfer/recompute_strict_source_only_transfer.py
N_JOBS=6 python3 scripts/design_validity/locked_validity_audit.py
python3 scripts/design_validity/augment_gse135779_design_residualisation.py
python3 scripts/design_validity/validate_design_validity_release.py
python3 scripts/validate_rc_package.py .
```

## Reproduction results

The following regenerated files were byte-identical to the locked candidate:

- GSE135779 restored donor metadata and label cross-tabs;
- strict source-only transfer metrics;
- paired DeLong comparisons;
- strict transfer donor predictions;
- feature-harmonisation record;
- locked audit summary;
- all 540 repeat-level locked metrics;
- GSE135779 matched-restriction estimates;
- batch-predictability estimates and summary.

The compressed OOF file differed only in the gzip container timestamp. After
decompression, its 5,645,546 bytes were identical, with SHA256:

```text
235fa66e1abbc51c7886e117099962d285c6d9a264b035b44ea0ee1bdf9b0164
```

The regenerated scientific validation record was semantically identical after
removing the expected absolute location of the clean-copy manuscript. Its locked
summary hash was:

```text
8e6120bc13dafdb852d96f054b0f0f3ee05daa0a4feb1cd8744065e602b88b62
```

Both validators returned `status: passed`.

## Interpretation

This is a clean-copy deterministic reproduction of the packaged donor-level
analyses. It does not regenerate the upstream frozen cell embeddings from raw
single-cell matrices, and it does not retrain the exploratory learned-pooling models,
whose cell-level training inputs are intentionally not duplicated in the repository.
Those upstream boundaries are separately locked by checkpoint, tokenizer,
environment, sampling and input-output provenance records.
