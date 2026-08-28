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

## Bond/connectivity analysis

PR #2 demonstrated that geometry-only bond perception is insufficient for coordination and multicentre systems. PR #3 should use a provenance-aware hierarchy:

1. producer-reported bond/population information when explicitly available
2. wavefunction-derived analysis from FCHK/Molden data when scientifically supportable
3. conservative geometry fallback
4. unresolved/unknown rather than forced connectivity

A Mayer-type pathway may be built from AO overlap and density information. Pairwise bond indices must not be treated as proof of 3c2e/3c4e; multicentre classification requires separate MO-subspace/participation analysis.

## Regression cases

- H2
- Cp- Cartesian and spherical references
- H3+ — 3c2e
- XeF2 — 3c4e
- [TiF6]2- — Oh
- [ZnCl4]2- — Td

Equivalent Gaussian FCHK and Molden exports should converge on the same internal orbital semantics where their source data are equivalent.

## Deferred item 11

The larger validation issue explicitly deferred during PR #2 remains a first-class PR #3 workstream. Its exact corrective implementation must be recovered from the stored PR #2/manual-validation context before coding; the FCHK migration by itself must not be presented as solving it.

## Scientific constraints

- no molecule-specific chemistry hard-coding
- no `occupied -> bonding` or `virtual -> antibonding` shortcuts
- no strict SALC claims without sufficient information
- inferred annotations retain source, confidence and heuristic flags
- use `unavailable` / `unclassified` when evidence is insufficient

## Initial acceptance gates

- restricted and unrestricted FCHK parser tests
- shell-convention coverage through the angular momenta represented by the project model
- FCHK-vs-Molden regression for shared Gaussian calculations
- Windows + Ubuntu CPU CI
- Windows CUDA 12.8 / `sm_120` viewer build remains healthy
- no regression of PR #2 UI, exports, localisation or rendering
