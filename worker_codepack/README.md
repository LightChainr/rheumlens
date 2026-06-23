# RheumLens Worker Code Pack v2

This package replaces the Markdown-only 20260621 reference pack. It contains an executable, reviewed P8.2 runner and strict validation utilities.

## Status

```text
Package status: CODE_READY_FOR_HOST_GATES
Formal 1000-rep status: NOT AUTHORIZED UNTIL HOST GATES PASS
```

## Important changes

- Uses the exact 1000-rep seed table recovered from accepted `scgpt_mean` output.
- Parent loads one locally staged NPZ, then Linux `fork` shares it read-only.
- Workers never read/write shared storage.
- The coordinator is the sole checkpoint/error writer.
- Local checkpoints are periodically published as compact durable checkpoints.
- Resume identity includes runner, project code, config, folds, seed table, input manifest and staged input hashes.
- Strict output requires exactly rep IDs 0..999, matching seeds, finite AUCs and no errors.
- No environment-variable dump and no automatic selection of arbitrary result CSVs.

## Required order

1. Finish migration and compare source/destination key manifests.
2. Promote incoming directory only if `/autodl-fs/data/rheumlens` does not already exist.
3. Rebuild the packaged `rheumlens-core` environment locally; require compileall, 11/11 pytest and 9/9 smoke.
4. Generate the frozen input manifest with `runner/make_input_manifest.py`.
5. Audit migrated scgpt_mean with `runner/audit_scgpt_mean.py`.
6. Run package tests.
7. Run 1-vs-4 worker numerical gate, real TERM/resume gate, bad-hash resume rejection, then 16-worker resource smoke.
8. Supervisor may authorize 24/32-worker formal only after evidence review.

## Do not do

- Do not use the old `perm_parallel_kme.py` or any 48-spawn partial.
- Do not regenerate P4–P7 assets.
- Do not launch PCA and KME concurrently.
- Do not treat this package's presence as permission to start formal.
