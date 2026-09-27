# PR #3 plan — FCHK-first Gaussian interface

PR #3 starts from the merged PR #2 baseline and changes the input/data-analysis layer while preserving the GPU-resident orbital rendering architecture.

## Primary interface

- FCHK becomes the preferred Gaussian wavefunction input.
- Molden remains a supported compatibility/exchange format for Gaussian and non-Gaussian producers.
- `.chk` may be accepted through `formchk` when Gaussian tools are locally available.
- `.log/.out` is an optional enrichment source for producer-reported symmetry, population, NBO/Wiberg-style results and calculation metadata.

## FCHK parser targets

Map the following into the existing `Wavefunction` model:

- atomic numbers and Cartesian coordinates
- shell types and shell-to-atom mapping
- primitive exponents and contraction coefficients
- Cartesian vs pure-spherical shell conventions
- alpha/beta electron counts
- alpha/beta MO energies and coefficients
- available density matrices and related machine-readable metadata

The parser must remain strict about dimensions and must never silently reinterpret a shell convention to make a file fit.

## Native symmetry engine

COV does **not** take a runtime dependency on libmsym or copy its implementation. libmsym and the older molecular-symmetry literature are useful design references, but PR #3 builds a new symmetry layer around COV's own wavefunction model and basis conventions.

The symmetry stack is deliberately split into verifiable layers:

1. generate candidate axes and planes from centred molecular geometry
2. validate every symmetry operation by an explicit same-element atom permutation and a recorded mapping error
3. classify point-group families from the validated operation set rather than molecule-specific templates
4. preserve producer-reported Gaussian symmetry as higher-precedence metadata
5. build AO/MO representations on top of the same validated operations in a later layer, with Cartesian and real-spherical basis conventions handled explicitly
6. assign a conventional irrep only when the representation evidence is sufficient; otherwise retain an unlabeled symmetry fingerprint rather than guessing

The first native layer covers linear groups and the common axial/polyhedral families needed by the regression set, including D3h, D5h, Td and Oh. Improper rotations are represented explicitly so S_n-only cases are not silently reduced to C_n.

For MO symmetry, the intended design is independent of chemistry heuristics: transform the actual AO basis under a validated operation, evaluate the MO/subspace representation using the AO overlap metric, and classify the resulting characters. Degenerate subspaces are treated as subspaces rather than forcing arbitrary individual members to carry one-dimensional labels.

## Bond/connectivity analysis

PR #2 demonstrated that geometry-only bond perception is insufficient for coordination and multicentre systems. PR #3 uses a provenance-aware hierarchy:

1. producer-reported bond/population information when explicitly available
2. wavefunction-derived analysis from FCHK/Molden data when scientifically supportable
3. conservative geometry fallback
4. unresolved/unknown rather than forced connectivity

The current wavefunction-derived path recovers an AO overlap matrix from a complete canonical MO block when possible, reconstructs density matrices when safe, and computes pairwise Mayer-type indices. Pairwise bond indices must not be treated as proof of 3c2e/3c4e; multicentre classification requires separate MO-subspace/participation analysis.

## Gaussian log/out enrichment

A sibling Gaussian `.log/.out` may enrich an FCHK with producer metadata. The enrichment path accepts the final `Orbital symmetries:` block and point-group report, deliberately excluding `Initial guess orbital symmetries:`. Orbital-count mismatches are rejected rather than shifted or padded.

## Regression cases

- H2
- Cp- Cartesian and spherical references
- H3+ — 3c2e
- XeF2 — 3c4e
- [TiF6]2- — Oh
- [ZnCl4]2- — Td

Equivalent Gaussian FCHK and Molden exports should converge on the same internal orbital semantics where their source data are equivalent. Real Gaussian 09/16 FCHK files are preferred over synthetic fixtures whenever a redistributable public regression source is available; synthetic fixtures remain useful for narrowly targeted parser invariants.

## Deferred item 11

The larger validation issue explicitly deferred during PR #2 remains a first-class PR #3 workstream. Its exact corrective implementation must be recovered from the stored PR #2/manual-validation context before coding; the FCHK migration by itself must not be presented as solving it.

## Scientific constraints

- no molecule-specific chemistry hard-coding
- no `occupied -> bonding` or `virtual -> antibonding` shortcuts
- no strict SALC claims without sufficient information
- inferred annotations retain source, confidence and heuristic flags
- producer metadata outranks derived metadata
- use `unavailable` / `unclassified` when evidence is insufficient

## Initial acceptance gates

- restricted and unrestricted FCHK parser tests
- shell-convention coverage through the angular momenta represented by the project model
- FCHK-vs-Molden regression for shared Gaussian calculations
- native geometry symmetry regressions for Dinfh, D3h, D5h, Td and Oh
- Windows + Ubuntu CPU CI
- Windows CUDA 12.8 / `sm_120` viewer build remains healthy
- no regression of PR #2 UI, exports, localisation or rendering
