# COV Helper 2.2.1 desktop deployment

The desktop Helper now opens from a single `Resume` deployment. It reads current
status without importing, dispatching or resuming calculations. Manual calculation
controls are disabled at every startup and cannot be restored from saved settings.

The former standalone directory has been retired after its complete frozen
calculation block was copied to `Resume/recovery/frozen-reference-20260910` and
verified against every original file. Existing Resume calculations, checkpoints,
progress records, archived inputs and recovery dependencies remain in place.

## Delivered behavior

- `Resume/START.cmd` opens the visible `COV-Helper.pyw` desktop entry using the
  configured existing Python 3.12 windowed interpreter.
- All status, checking, manual dispatch and pause/recovery controls use the single
  active Runtime 2.2 bundle and `runtime/config.json`.
- Startup, refresh and close issue no calculation or control commands. A separate
  read-only deployment check verifies the 273 frozen inputs and resource budget.
- The accepted computational runtime identity remains
  `84fe2a76288c49870928b3809af258a0535cf42ad693939b910e1deff5d90235`.
  Scientific policies and the previous native-capability gates are unchanged.
- The old Helper package, configuration, entries, program backups and superseded
  source bundles are retired through the Windows Recycle Bin. Local deployment
  receipts retain exact original paths, Recycle Bin locations and content hashes.

## Non-computing verification

- 120 controlled local regression tests passed in 31.609 seconds, with zero
  remaining owned processes. The added Helper tests use fake dispatch and a hidden
  Tk window; they launch no Gaussian calculations.
- After final cleanup classification changes, 21 focused Helper/deployment/entry
  tests passed. They cover read-only startup and close, explicit per-window enable,
  malformed selections, existing-coordinator admission, protected cleanup paths,
  exact frozen migration, destination conflicts and unclassified old results.
- A disposable file and directory verified actual Recycle Bin preservation and
  receipt lookup. Explorer can expose directory metadata shortly after a successful
  shell operation; receipt lookup is bounded and matches the exact original path.
- The deployment verification checks all protected file paths, sizes, original
  timestamps and file IDs, plus content hashes for metadata and all migrated files.
  The complete local inventory is retained with the deployment evidence.

See `verification.json` and `cleanup-summary.json` for measured deployment results.
Local receipts are stored under `Resume/runtime/deployments`; they retain recovery
locations without publishing machine account identifiers.

[Desktop operation and recovery instructions](../runtime_v2/HELPER.md) describe the
new entry, explicit calculation controls and restoration procedure. Candidate
collection still does not constitute scientific acceptance. No formal calculation
is part of this deployment.
