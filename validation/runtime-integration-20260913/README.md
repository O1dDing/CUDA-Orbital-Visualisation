# Runtime 2.2 integration and validation — 2026-09-13

The reference runner now shares one process-control implementation, one state store, and compatible local entry points. The program PRE11 / PR #3 remains separate from this validation-branch update.

## Verified implementation changes

- The pause transaction verifies live Job membership and thread creation identity, owns exactly one suspend increment, checks its deadline and cancellation inside enumeration, requires coverage of every computational process, and rolls back partial holds. Failed rollback keeps the owned handles visible to Job cleanup.
- The compute launcher starts detached with explicit standard-file handles. Windows may still create a console host for descendants. The exact System32/conhost.exe path is classified as OS infrastructure and remains responsive; unknown console hosts fail closed. Resource accounting and termination still cover the full owned Job.
- During native Gaussian link transitions, QueryFullProcessImageNameW can temporarily fail while the process is still listed. Such a process remains required, unresolved coverage. A later successful identification or confirmed exit is required before pause acknowledgement; the original scan deadline is unchanged. The accepted native report includes the transient observations.
- JSON readers and replacement share a Windows named mutex. Neither lock naming nor optional-read checks open/stat the mutable file outside that mutex. Python 3.12's Windows realpath and slow stat implementations can use share mode zero, so resolving the final file or testing its existence before locking still races with replacement. The fix resolves the parent directory for lock naming and reads optional files directly under the lock. Permanent replacement errors preserve the previous document and propagate. See [CPython 3.12.10 implementation](https://github.com/python/cpython/blob/v3.12.10/Modules/posixmodule.c#L4518).
- Control commands have serialized sequence numbers and a sticky interrupt sequence. Explicit release after an active-time limit grants the next bounded interval. Job cleanup, callback failure, and guard termination confirm zero owned processes.
- Migration validates source binding, input and binary identity, original receipts, and SHA-256 values before publishing metadata. Original producer identity and file paths remain visible. Native capability receipts from 2.1 are never imported as new acceptance.

## Local validation

The final controlled regression passed **102 tests**. These include 20 separately bounded early-pause trials, parent/child pause and resume with an independent test process continuing to run, cancellation, callback cleanup, repeated active-time limits, concurrent commands, 3,000 independent reads racing 600 replacements, migration integrity, and compatible entry activation/rollback. The separate frozen recovery suite passed **7 tests**. The watchdog self-test deliberately hung its own worker and child; both exited after approximately one second.

Fresh native acceptance ran isolated distorted water and NO copies. Their Opt and Freq attempts each acknowledged RAM pause and resumed the same attempt. All four holds measured zero tree CPU increment while held. Pause acknowledgement was approximately **0.297–0.313 seconds**. Final energy, pair distances and frequencies matched their uninterrupted baselines exactly at the recorded precision. RPBE1PBE and UPBE1PBE RAM pause passed for the runtime identity in the receipt. This is not scientific acceptance of formal candidates.

L103 Opt segmentation and analytic RWF restart remain unaccepted and gated. The new runtime does not silently substitute a full frequency replay.

## Historical failures and remaining limits

The original early-pause series had **4 failures in 20 trials**, with roughly 30–53 second delays. Seventy controlled baseline trials later did not reproduce that long delay. Source inspection establishes the old missing inner deadline/cancellation checks and console-host inclusion, but the precise Windows call responsible for those four historical stalls remains unproven. Passing newer trials does not erase those failures.

The five-second scan limit is cooperative between API calls. It cannot preempt a blocked kernel call. The test harness adds an independently owned outer Job and a separate deadline; it does not touch Codex, terminals, or unrelated processes.

Intermediate experiments are retained. The first new native harness only considered the final Freq hold and could miss an Opt pause failure; that receipt is explicitly not accepted. The harness now requires both stages and empty failure histories. A subsequent stricter native trial exposed the transient image-query failure. The final report records successful resolution before acknowledgement.

## Preserved REF-001 state

The current read-only audit verified **337 collected stages** and all **273 frozen inputs**. Migration yields **103 candidate_collected, 1 partial, 2 needs_review, 167 not_started**. Scientific passes remain **0**. No formal candidate queue was started by this integration.

OLD-018's saved checkpoint retains SHA-256:
`12c50a01ca0d2b37b1f1234df94ea7d0b2c62553ef605a1108212238de08087e`.
OLD-018 / OLD-051 review holds and original recovery material remain in place.

The prior data PRE described 90 collected cases. The new delta archives only newly observed bytes and reuses existing assets by digest. Its manifest records 4,308 source files, 335 newly unique files, and 15 new archive parts. The source data, legacy results, old runtime, failed trials and new run directory remain distinguishable.

## Entry points, publication and evidence

Resume/START.cmd and FastPause-2.1/validation/runtime_v2/START.cmd use the active bundle and shared configuration under Resume/runtime. The work directory is identity-specific and selected by runtime_directory. Every replaced entry has a byte-verified backup and rollback receipt. Activation requires no Gaussian/coordinator and accepted R/U RAM receipts.

- [Integration summary](integration-summary.json)
- [Accepted native report](native-report.json) and [capability receipt](native-capabilities.json)
- [Current stage audit](current-ref001-stage-audit.json)
- [Migration receipt](migration-complete.json)
- [Final controlled regression log](release-regression.log)
- [Original long-delay failures](historical-early-pause-failures.log)
- [Data PRE](https://github.com/O1dDing/CUDA-Orbital-Visualisation/releases/tag/cov-ref-paused-20260906)
- [Program PRE11](https://github.com/O1dDing/CUDA-Orbital-Visualisation/releases/tag/v0.3.0-pre.11)
- [Program PR #3](https://github.com/O1dDing/CUDA-Orbital-Visualisation/pull/3)

The release bundle manifest records the exact validation commit, CI run, file hashes and upload receipts. PR #3 stays draft and unmerged; no auto-merge is enabled.
