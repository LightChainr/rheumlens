# Supervisor deployment instructions

## Package decision

```text
Code package: READY_FOR_HOST_GATES
Formal execution: NOT AUTHORIZED BY PACKAGE ALONE
```

The previous Markdown-only pack is superseded. Do not copy code blocks from it.

## After migration promotion

1. Place this directory at:

```text
/autodl-fs/data/rheumlens/supervisor_control/rheumlens_worker_codepack_v2
```

2. Verify `MANIFEST_SHA256.txt`.
3. Build `rheumlens-core` on local disk using `env/environment-core-p82.yml` (adjust only the environment prefix, not package versions).
4. Require project compileall, 11/11 pytest and 9/9 project smoke.
5. Generate the input manifest:

```bash
python runner/make_input_manifest.py \
  --root /autodl-fs/data/rheumlens \
  --output /autodl-fs/data/rheumlens/results/P8_2_formal_permutation_runs/manifests/GSE174188_CD4_P8_2_input_manifest.json
```

6. Audit `scgpt_mean` against the packaged frozen seed table.
7. Stage only the input required by the current method.
8. Run `runner/gate_controller.py` for PCA. It must report `PASS_ALL_GATES`.
9. Review resource evidence, then authorize one 24/32-worker formal job using `scripts/run_p82_job.sh`.
10. Run `scripts/validate_p82_job.sh`; Supervisor, not Worker, records final ACCEPTED.
11. Repeat gates/formal/validation for KME only after PCA acceptance.

## Frozen scientific values

```text
base_seed: 20260619
seed table: recovered verbatim from accepted scgpt_mean 1000-rep output
donor_expression_pca observed AUC: 0.9856590597331338
kme_multiscale@scgpt observed AUC: 0.977865070457663
folds: splits/authoritative_primary/GSE174188_CD4.csv
cohort: GSE174188_CD4
estimand: matched_500_cells_per_donor
```

Do not regenerate the seed table, change the methods, or infer observed AUC from a newly run model.

## Safety properties

- One parent loads the staged dataset before Linux fork.
- Workers perform computation only; the coordinator is the sole writer.
- Local checkpoint every 5 completed reps; durable checkpoint every 10.
- Resume refuses any runner/project/config/folds/seed/input hash change.
- Formal output requires exactly 1000 unique finite rows with exact frozen seeds.
- No environment dump, arbitrary CSV selection, silent NaN or automatic formal launch.
