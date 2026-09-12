# Validation records

## Current software validation

[COMMON-ROOT-008E](common-root-008e-20260912/README.md) records the COV source revision `3752f8a`: independent numerical checks, scoped orbital and diagram checks, dedicated controls, the fixed 40-case precheck, the new 273-case software round and ordinary-viewer interaction. The two builds each passed 39 component tests. A single unexplained blank startup remains recorded separately from the successful interactions.

The [pre.11 release notes](../docs/releases/v0.3.0-pre.11.md) describe the distributable programs, evidence and source identities. Software-regression passes apply to the declared checks; they do not establish physical electronic states or full scientific acceptance.

## REF-001 calculation data

The [2026-09-12 reference snapshot](ref001-progress-20260912/README.md) contains 90 collected candidates, two partially collected cases, two cases needing review and 179 not-started cases. It verifies 294 main Runtime v2 stage records plus 10 native-capability probe stage records, with legacy imports identified separately. Complete scientific acceptance remains zero.

The detailed file manifest resolves new data, earlier archive members and identical-content aliases by SHA-256. Large calculation files are retained in the [REF-001 data archive](https://github.com/O1dDing/CUDA-Orbital-Visualisation/releases/tag/cov-ref-paused-20260906); verified asset identities are published with the snapshot. Inputs, logs, checkpoints, timeout/failure records and recovery material remain available.

The [2026-09-09 snapshot](ref001-progress-20260909/README.md) is a historical baseline with 17 collected candidates and 37 completed legacy stages. Its raw receipts and uploaded assets retain their original bytes. The later snapshot supersedes its aggregate status, not its evidence.

## Runtime and recovery

[Runtime 2.1](runtime_v2/README.md), its [pause/recovery controls](runtime_v2/FAST_PAUSE.md) and the [execution runbook](runtime_v2/WORK_PROMPT.md) describe the current execution layer. Runtime evidence is stored separately in `work/runtime-v2-fastpause`; old calculation directories are imported only after identity and lock checks.

The current native-capability evidence accepts RAM pause for its recorded RPBE1PBE/UPBE1PBE runtime and Gaussian installation. Optimization link-103 segmentation and analytic RWF restart remain unaccepted. Simulation and CI results do not replace installation-specific native evidence.

`paused-20260906/` remains an immutable legacy source/input snapshot. Runtime resource limits use detected physical cores, preserve aggregate CPU/memory/disk budgets and keep the established model, basis, state, convergence and grid requirements.

## Earlier software rounds and scientific scope

The [COMMON-ROOT-008C record](common-root-008-20260909/README.md) preserves an earlier completed software round. Later ordinary-viewer checks exposed the details-window defect corrected in the E revision. Other dated directories retain their original failures, contracts and results.

Scientific requirements and implementation details are described in [the validation specification](../docs/FCHK_VALIDATION_AGENT.md) and [native validation](../docs/native-validation.md). The remaining physical work comprises the original 273-case reference/state/geometry adjudication and at least 50 distinct external molecules. Candidate collection, normal process termination, display success and automated software tests are separate evidence categories.
