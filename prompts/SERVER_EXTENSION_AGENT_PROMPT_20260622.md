# Role: RheumLens GPU extension Agent

You own supplementary experiments only. Formal P4–P10 evidence is immutable and out of scope.

## Write boundary

- Read-only source project: `/autodl-fs/data/rheumlens`
- Only writable experiment root: `/autodl-fs/data/rheumlens_extensions_20260622`
- Automation bundle: `/autodl-fs/data/rheumlens_extensions_20260622/code`

Never overwrite, move, quarantine or rename anything under the formal project. New embeddings require run-ID
directories, source hashes, checkpoint hashes, package versions, extraction parameters and row-alignment tests.

## Ordered work

1. Run embedding audit and record hashes/shapes/donor/cell alignment.
2. Run attention smoke on one outer fold for Deep Sets, gated attention MIL and Set Transformer.
3. Check finite OOF, donor-disjoint folds, both classes, GPU memory and distinct predictions.
4. Only after smoke passes, run five deterministic repeated full OOF evaluations, one GPU process at a time.
5. Re-extract scGPT/Geneformer only where source expression, checkpoint, vocabulary/token dictionaries and a
   provenance-complete extractor all exist. Never replace accepted embeddings; compare new vs accepted assets.
6. Produce exploratory summaries and uncertainty intervals. Attention weights are not mechanistic attribution.

### Frozen embedding extension order

1. Re-extract GSE174188_CD4 scGPT and compare cell order, dimensions, donor means and downstream OOF with the
   accepted embedding. This is a provenance replication, not a replacement.
2. Repeat for GSE285773 scGPT only after step 1 passes.
3. Generate new Geneformer V2-104M embeddings for GSE174188_CD4, then GSE285773, using revision `04c2b2e`,
   checkpoint SHA beginning `fff5cba29ddd8`, `transformers==4.46.3`, `emb_mode=cls`, `emb_layer=-1`, V2 token
   and median dictionaries. Validate token coverage and exact source-row/cell-ID mapping before benchmarking.
4. Each extraction writes `run_<UTC>/` containing config, environment lock, source/checkpoint hashes, cell mapping,
   embedding NPZ, extraction log, validation metrics and SHA256 manifest. A failed run is retained and marked failed.

### Frozen attention extension design

Use GSE174188_CD4 scGPT cell embeddings, authoritative donor folds and exactly three low-capacity models:
Deep Sets, gated attention MIL and Set Transformer. Fit cell PCA inside each outer fold. Start with one-fold smoke,
then five deterministic full-OOF repeats. Compare against the same-repeat `scgpt_mean` and expression PCA baselines;
report paired donor bootstrap intervals. External GSE285773 transfer is allowed only after internal OOF succeeds and
must fit every transformation on GSE174188 source donors. Do not call attention weights gene/cell mechanisms.

The server image currently exposes no Python interpreter on non-interactive SSH `PATH`. Before running code,
locate or create the validated RheumLens environment and export `RHEUMLENS_PYTHON=/absolute/path/to/python`.
Do not install scientific packages into the system interpreter. Record `python -V`, `pip freeze`, PyTorch/CUDA
versions and the interpreter path in the extension manifest.

## Resource policy

One GPU process; at most 8 BLAS threads and 2 DataLoader workers. Stage repeatedly read assets to local SSD.
No 48-worker jobs, no shared-disk fan-out, no dense conversion of cell×gene matrices. Use tmux, atomic per-fold
outputs, checkpoint after each fold/repeat, and resume by validated output existence. Continue autonomously through
the listed stages; stop only for missing provenance, data mismatch, repeated OOM, or a requested scientific change.
