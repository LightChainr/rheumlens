# Local validation report

Date: 2026-06-22

## Passed

- Python compileall for runner and tests.
- Bash syntax checks for all shell scripts.
- Package frozen seed table contains 1000 unique rep IDs/seeds recovered from accepted scgpt_mean output.
- Toy-data 1-worker versus 4-worker AUCs match exactly for identical rep seeds.
- Cooperative stop/checkpoint/resume matches uninterrupted reference exactly.
- Real SIGTERM/checkpoint/resume matches uninterrupted reference exactly.
- Resume after config content change is rejected.
- Gate controller reaches `PASS_ALL_GATES` on toy data.
- Static scan found no spawn context, full environment dump, `np.save('.tmp')` rename pattern, or silent pytest masking in executable package code.

## Still required on target host

- Migration source/destination key manifest comparison.
- Rebuild pinned rheumlens-core and obtain compileall + 11/11 pytest + 9/9 smoke.
- Full-data scgpt_mean audit.
- Full-data 1/4/8/16-worker gates and resource measurements.
- Supervisor authorization before formal 1000-rep execution.

Local tests validate runner mechanics; they do not constitute scientific acceptance of target-host outputs.
