# REF-001 execution and recovery runbook

Runtime 2.1 keeps its execution evidence in `work/runtime-v2-fastpause`. Before using an existing calculation tree, inspect its actual processes, control mode, file identities, current repository revision and recorded capabilities. The [reference archive](../ref001-progress-20260912/README.md) gives the latest independently verified inventory; historical summaries are not process-state authority.

## Pause, continue and save

- `4/hold` pauses the owned process threads while retaining RAM, memory allocation and Job handles. `6/resume` continues the same processes and restores dispatch capacity. RAM pause is not a disk save; power loss or closing the coordinator can destroy that state.
- `10/pause` stops dispatch after an accepted optimization segment or the current stage reaches its boundary.
- `8/interrupt` requires `SAVE`, stops only the owned Windows Job tree, waits for all descendants, then performs cold copy, flush, hash and checkpoint probing. Frequency recovery preserves named RWF/INT/D2E/CHK files; optimization preserves CHK/input/log and retains the original RWF. Unwritten work may be lost.
- `9/shutdown` exits after the current boundary. A RAM-paused stage is resumed to reach that boundary; shutdown is not an immediate close operation.

Thread control uses owned handles and reverses only the supervisor's suspend increment. Process-name termination, stale-PID control and unreviewed native freeze interfaces are excluded. The full startup-affinity compatibility path is retained for Gaussian 16W A.03. The 72-hour limit counts active elapsed time and defaults to RAM hold on expiry.

## Native capability gates

Optimization segments use an input-declared `%KJob L103 2` boundary followed by cold verification and `Opt=Restart`. Log-line detection followed by a timed process kill is not a checkpoint boundary. Acceptance requires the expected KJob termination, actual optimization progress and agreement with a baseline calculation.

Analytic frequency recovery copies an immutable cold snapshot into a new attempt and uses `#p Restart`; `Freq=Restart` is the numerical-frequency route. Missing or unaccepted RWF recovery enters `needs_review`. A reviewed `replay` decision can rerun frequency while retaining the completed optimization.

Every capability must match the runtime source identity, Gaussian launcher/link/DLL fingerprints, R/U method and evidence hashes. CHK/FCHK parse success does not certify all optimization or RWF restart histories. The current archived receipt accepts RAM pause for RPBE1PBE and UPBE1PBE; link-103 segmentation and analytic RWF restart remain unaccepted.

## Existing calculations

Keep legacy sources, inputs, identities, archives and original calculation trees. Cold import requires exclusion locks and verified file identities. OLD-018's timeout checkpoint and OLD-051's failed stability attempt require explicit review; deleting a review marker must not convert a retained checkpoint into a fresh calculation. Unique checkpoints are never used for destructive capability experiments.

Completed stages are reused by their recorded reference and binary identities. The current snapshot has 90 collected candidates, two partially collected cases, two cases needing review and 179 not-started cases. Collection does not establish scientific acceptance. All failed attempts and recovery evidence remain available.

## Resource and scientific contract

Use measured physical topology rather than SMT thread count: the aggregate budgets are 14 cores for at least 16 physical cores, 10 for 12, 6 for 8, 4 for 6, 2 for 4, and 1 for one or two. At most two jobs run concurrently. Actual REF basis size, open-shell status and free resources determine 1–14 cores per new attempt or segment; allocations are not changed inside an active stage.

Each job uses 24 GB input memory and a 32 GiB Job-tree limit while retaining system memory and disk headroom. Global and legacy-directory locks enforce exclusion. Model, basis, charge/spin, VeryTight convergence, SuperFineGrid, fixed-geometry purpose and the two stability-repair limit remain unchanged. SCF fluctuations and log size are not completion percentages.

## Validation and deployment

Run `test_runtime.py`, `test_fast_pause.py` and the frozen recovery tests before adopting implementation changes. `python -X utf8 resume.py native-acceptance` performs real Gaussian tests in isolated water/NO copies while holding the same exclusion locks. It exercises baseline, RAM pause, segmented optimization and analytic-frequency restart, retaining logs and comparing energy, geometry and frequency. A path that did not actually interrupt or failed is not accepted.

Simulation and Windows Python process-tree tests do not replace Gaussian native acceptance. Preserve failed evidence, correct the implementation and validate again without weakening scientific criteria. Migrate only from verified records, keep original directories, and record the resulting implementation/configuration identities, capability evidence and remaining limits.
