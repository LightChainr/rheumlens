# Changes from rheumlens_worker_codepack_20260621.zip

- Added executable `perm_parallel_v3.py`; removed dependence on an unavailable candidate runner.
- Recovered and packaged the exact accepted scgpt_mean 1000-rep seed table.
- Corrected the actual scgpt permutation path to include `GSE174188_CD4`.
- Replaced arbitrary CSV discovery with fixed result schemas and paths.
- Replaced string-based p-value acceptance with strict recomputation.
- Added exact rep/seed/finite/row assertions.
- Added parent-load + Linux fork architecture and sole-writer coordinator.
- Added local/durable checkpoint layers and hash-guarded resume identity.
- Added real SIGTERM recovery and bad-hash resume rejection gates.
- Removed unfiltered environment dumps.
- Added safe promotion refusal when a canonical project root already exists.
- Added pinned CPU environment and local integration tests.
- Corrected manifest generation to exclude the manifest itself.
