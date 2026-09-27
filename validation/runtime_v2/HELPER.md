# COV Helper 3.0

The desktop Helper opens in read-only mode. It shows the current coordinator lease,
control state, candidate counts, runtime identity, and calculation paths. Opening,
refreshing, checking the deployment, and closing the panel never dispatch work.
Runtime 3.0 adds atomic whole-tier resource admission and creation-time Job ownership. Scientific inputs and pause/recovery gates are preserved; native acceptance is bound to the new runtime identity.

## Deployment layout

`Resume` is the only active deployment directory:

- `START.cmd` launches the configured existing Python 3.12 windowed interpreter.
- `COV-Helper.pyw` loads `helper_ui.py` from `runtime/active.json`.
- `runtime/config.json` is the single calculation configuration.
- `runtime/launcher.json` records the interpreter, Helper version and read-only startup policy.
- `runtime/releases/<identity>` contains the verified active source bundle and frozen dependencies.
- `data`, `work`, `assets`, `OLD-018-SALVAGE` and `tools` retain original calculations,
  archives, checkpoints, identities, progress, and frozen recovery dependencies.
- `recovery/frozen-reference-20260910` preserves the complete frozen calculation
  block migrated from the retired standalone deployment, including its original bytes.
- `runtime/deployments` records migration, cleanup, new entries and Recycle Bin locations.

There is no second active runtime or compatibility directory. References inside
preserved historical records remain provenance; they are not used as active launch paths.

## Desktop use

Open `Resume/START.cmd`. The status screen is idle unless a separately started
coordinator already holds the calculation lease. The screen reports that lease
instead of treating old status text as proof of a live calculation.

`刷新状态` and `检查部署（不计算）` are available without enabling calculation controls.
Manual controls must be enabled explicitly in each new window. This permission is
kept only in memory and is never restored from settings. A new batch defaults to one
case and one worker slot; explicit case IDs may narrow the selection. The runtime
continues to enforce its admission locks, physical core budget, checkpoint identities,
and scientific capability gates.

RAM hold, resume, stage-boundary pause, stop/cold-save, and stage-boundary shutdown
operate only through the shared runtime commands. No control is written when no
coordinator owns the work directory. Closing the panel does not terminate a running
coordinator. RAM hold is not a disk checkpoint; wait for cold-save evidence before shutdown.

## Cleanup and recovery

`helper_deployment.py plan` pins the retired program files, verifies replacement
frozen dependencies, and rejects unclassified additions. `apply` requires idle
runtime/Helper leases, verifies all planned identities again, copies the frozen
calculation block without overwriting different files, and verifies every copied byte.
Only then is the explicitly named old directory retired.

Cleanup rejects reparse points and targets outside the two named directories.
It excludes all protected calculation roots. Retired content is sent to the
[Windows Recycle Bin](https://learn.microsoft.com/en-us/windows/win32/api/shellapi/nf-shellapi-shfileoperationw);
the receipt records its original and actual Recycle Bin paths and verifies the retained bytes.

To recover a retired file, close the Helper and ensure no coordinator or Gaussian is
running, consult its deployment receipt, and restore the recorded original location
from the Recycle Bin. Preserve any newer file before restoring an older version.
The migrated calculation block and existing Resume calculation directories do not
need restoration: they remain in Resume.

The focused regression suite uses fake command dispatch and a hidden Tk window.
It verifies that startup and close execute only a status read, that every mutating
operation is disabled by default, and that malformed selections and duplicate
coordinators cannot dispatch work. It launches no Gaussian calculation.
