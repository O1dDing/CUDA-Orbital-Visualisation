# REF-001 calculation archive — 2026-09-12

This snapshot preserves existing PBE0-D3BJ reference calculations, stage identities and recovery material. It does not start calculations or assign scientific acceptance. The 273-case reference set remains distinct from the COV software-regression corpus and the planned external-molecule validation.

| Case state | Count | Meaning |
|---|---:|---|
| Candidate collected | 90 | The recorded calculation workflow has collected its required stages; physical adjudication remains pending. |
| Partially collected | 2 | OLD-093 and OLD-094 have completed optimization records; frequency is next. |
| Needs review | 2 | OLD-018 retains its legacy optimization timeout; OLD-051 retains a failed stability attempt after optimization and frequency. |
| Not started | 179 | OLD-095 through OLD-273. |
| Scientific acceptance | 0 | Candidate collection is not proof of a ground state, physical geometry or complete scientific agreement. |

The aggregate is reconstructed from all case and stage files. A coordinator status snapshot describes its most recent selection and can contain historical process flags; it is not the aggregate case count.

## Identity checks

All 273 frozen reference inputs, basis files and reference identities were checked. The archive verifies 294 main Runtime v2 stage records and 10 native-capability probe stage records against checkpoint/FCHK/log hashes, Gaussian completion records, recorded input hashes, recorded formchk success and the actual FCHK method, basis dimension, atoms, charge and spin identity. No recorded-stage identity checks failed.

The 37 completed legacy stage records are also verified and retained. Imported v2 optimization/frequency records can reference the same earlier computation, so legacy and imported record counts must not be added together as independent new calculations. Native-capability probes are listed separately from reference cases.

See [case status](case-status.json), [stage verification](stage-verification.json), [runtime capability evidence](native-capabilities-verification.json) and the [file manifest](manifest.json).

## Archive organization

The inventory covers 4,880 files and 20,562,025,672 source bytes. Content already present in verified earlier archives is referenced by asset name, member path, byte count and SHA-256. Identical files retain all original paths through `same_content_as` entries. The 2,090 new unique files total 15,464,659,896 uncompressed bytes and are stored in 34 ZIP volumes totaling 4,423,489,633 bytes.

Large files are published in the existing [REF-001 data archive](https://github.com/O1dDing/CUDA-Orbital-Visualisation/releases/tag/cov-ref-paused-20260906). The `cov-ref001-progress-20260912-001-of-034.zip` through `034-of-034.zip` volumes preserve calculation inputs, logs, CHK/FCHK, available RWF and scratch state, failed attempts, native-acceptance evidence and the exact runtime Python sources. The original `cov-ref001-progress-20260912-metadata.zip` preserves the metadata snapshot taken before upload. The companion `cov-ref001-progress-20260912-published-metadata.zip` supplies the complete public index, detailed manifests and [upload verification](upload-verification.json) for the 34 data volumes and original metadata snapshot. The original manifest retains its historical pre-upload state; the upload receipt records completed publication. The final metadata archive has its own [published asset identity](published-metadata.json).

The original 23 data volumes and the 2026-09-09 increment remain available with their original bytes. Source files, failed evidence, checkpoints and cold copies are retained. Interpreter caches and ephemeral OS lock files are excluded from restore material; live locks must not be recreated from an archive.

## Recovery state

OLD-018 retains the legacy 24-hour timeout and its checkpoint with SHA-256 `12c50a01ca0d2b37b1f1234df94ea7d0b2c62553ef605a1108212238de08087e`. Its optimization is incomplete. The timeout is not classified as an SCF convergence failure. OLD-051 has verified optimization and frequency records but its failed stability attempt requires review. OLD-093 and OLD-094 retain their completed optimization checkpoints.

The existing native-capability receipt accepts RAM pause for RPBE1PBE and UPBE1PBE on its recorded binary/runtime identity. Optimization link-103 segmentation and analytic RWF restart are not accepted. A RAM pause is not a disk checkpoint, and a parsed CHK/FCHK does not by itself certify every restart route. These capability limits remain enforced; this archival pass runs no Gaussian or formchk process.

Restore into a new isolated directory. Resolve each manifest entry from its new archive, `existing_archive` locator or `same_content_as` entry, then verify its byte count and SHA-256. Preserve the relative tree and retain original metadata as provenance. Relocate only the working restore copy, revalidate the frozen reference/runtime/binary identities and use an accepted recovery route. Existing unique calculation directories must not be overwritten by restoration.

## Software validation

The independent [COMMON-ROOT-008E software record](../common-root-008e-20260912/README.md) covers the COV executable and its numerical/UI/export checks. Its 273 passing software cases do not increase the scientific acceptance count in this reference archive.
