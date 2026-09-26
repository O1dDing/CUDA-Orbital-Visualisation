# Focused canonical MO / NAO composition

This extension uses the existing COV central MO diagram to define the orbitals
of interest. Its valence, compact, delocalised-pi and multicentre settings remain
authoritative. Selecting a three-dimensional inspection orbital does not silently
add a level to the compact diagram. The NAO composition view uses the diagram's
actual member and opposite-spin counterpart indices, including unclassified and
delta-labelled members. A chemical annotation is not a prerequisite for a
numerically verified decomposition.

## Input and preservation

Open the original Gaussian canonical FCHK, then attach an NBO report, ARCHIVE,
AONBO and NBOMO as in `nbo-initial-integration.md`. The three additional explicit
paths are NAOMO (NAO rows / canonical MO columns), AONAO (AO rows / NAO columns),
and NAONBO (NAO rows / NBO columns). AONAO and NAOMO are required together;
NAONBO supplies an additional independent composition check. No filename suffix
is interpreted as a matrix type. File contents, axes, size and spin must agree.

The tested producer requests `NAOMO=W51 AONAO=W52 NAONBO=W53` in addition to
`AONBO=W37 NBOMO=W49 SAO=W50 ARCHIVE`. On Gaussian 16W A.03 the installed NBO7
interface is invoked through `Pop=NBO6Read`. One observed interface path wrote
NAONBO to its default FILE.30 despite the requested number; the campaign retains
that output and uses a complete GenNBO reanalysis bundle with explicit numbered
matrices. Do not concatenate matrices or NBO identities from different analyses.

GenNBO can reuse an existing archive when it includes the required data:
LCAOMO for NAOMO, density and overlap for the basis/density checks, and Fock for
E2 analysis. A reanalysis is a new analysis version of the archived wavefunction.
It is not a new SCF, and it is not silently merged into the original NBO report.

Canonical coefficients, ordering, energies, spin and source indices remain
unchanged. Localised orbital checkpoint options must not overwrite the canonical
identity. NBO diagonal Fock values are not canonical eigenvalues; E2 is a separate
perturbative interaction estimate, not a total bond or reaction energy.

## Numerical interpretation

For compatible AO representations, A is AONAO, S is AO overlap, C is the
canonical coefficient matrix and T is NAOMO. The import verifies A^T S A,
A T = C, T = A^T S C, and the available NAONBO/AONBO/NBOMO composite relations.
An atom or specified shell weight is the sum of squared NAOMO coefficients in
that group. Raw nonorthogonal AO coefficient squares are not percentages.

Complete NAO rows retain literal atom, angular and core/valence/Rydberg labels.
Principal shell labels are parsed when present rather than inferred from the
element or a guessed valence configuration. Tiny and zero entries, unknown labels,
core and Rydberg entries remain in the full numerical export. No omitted subset
is renormalised. A stored null column is not a physical zero-weight MO. Unsupported
rank-deficient representations fail explicitly rather than receiving fabricated
normalisation.

Electron contribution is enabled only when the spin-appropriate occupation
matrix closes against the actual density. Closed-shell total occupations and
unrestricted alpha/beta occupations are distinct. The initial CH3 paired FCHK has
only alpha coefficient columns; its archive nevertheless contains both spins.
Its alpha NAOMO electron contribution uses occupation one, not an inferred shared
occupation two. ECP-replaced core electrons are separate from the explicit
density and from printed NPA totals that restore the core.

Source compatibility has two separate layers:

- Coefficients, ordered atoms, basis, overlap, density and transformations are
  verified numerically. Whole-column phase reversal may be aligned for comparison;
  raw matrices and the phase record remain available. Arbitrary MO permutation or
  rotation does not count as a direct column match.
- The importer cannot infer the Gaussian job/step history from numerical equality.
  `producer_step_identity_status` therefore remains `producer_step_unverified`;
  the producer manifest supplies job, checkpoint, reanalysis and hash evidence.
  The legacy `strict_same_source` status name is retained for compatibility and
  must not be interpreted as proof of the temporal production chain.

The full-spin CH3 fixture uses the initial SCF FCHK. Its alpha and beta columns
were independently compared with the analysis archive using column phase and the
overlap metric. It is explicitly **not** a FCHK from the later NBO checkpoint.
The later checkpoint's missing-beta FCHK is retained as a required limited-input
fixture. No beta FCHK coefficients were fabricated or silently copied into it.

## Capability and missing-value fields

`NboNaoValidation.available` means that the spin's archive transform, labels and
density checks are usable. It can be true for `verified_archive_only` while
`direct_fchk_coefficients` is false. `NboMoDecomposition.available` is stronger:
the actual canonical column must exist and pass direct coefficient association.
An archive-only beta transform never creates a beta canonical MO.

Omitting all new NAO files leaves the original NBO functions available and reports
missing NAO transforms. Supplying NAOMO without AONAO explicitly requests an
incomplete extension and rejects that association. Malformed matrices, wrong spin,
wrong headings, inconsistent coefficients and ambiguous analysis segments are
rejected. An absent, omitted or unsupported quantity is not silently zero.

## View and export

The composition section appears with the existing central MO diagram. Choose a
member MO and the atoms or specified shells relevant to the current question.
Valence shells are the initial category; core and Rydberg categories and explicit
shell selections remain available. Thus a metal shell marked Rydberg by the
producer can still be deliberately inspected. The diagram uses actual composition
nodes and lines, with no universal percentage or top-N edge cutoff and no
remainder/other node. Unshown relationships are retained in the data and are not
declared chemically irrelevant.

Atomic full weights and the subtotal of selected shells are distinct. Explicit
ligand groups retain their atom membership rather than being inferred from a
particular test molecule. The graph carries the selected canonical identity and
its original energy/spin. NAO composition relations and E2 interactions have
different meanings; composition lines must not be read as donor-acceptor energy
arrows. Full NBOMO and NPA/NBO/Wiberg/E2 source records remain accessible in the
NBO report and exports.

Multiple saved, disjoint ligand groups can be shown together. The optional
ligand s/p/d grouping sums each group's actual angular components and retains
ungrouped atoms as individual n/l shells; it never counts a grouped atom twice.
Draft atom selections do not become a ligand until the user adds the group.
The explicit "Inspect this MO in 3D" action changes the inspected canonical
orbital without changing the central diagram's member set.

The original export bundle is extended with `.focus.json`, `.focus.csv`, `.focus.groups.csv`,
`.focus.svg` and `.focus.png`. JSON and CSV retain all canonical decompositions,
including those outside the central view and unavailable identities. View metadata
records the central member set, focused MO, selected groups/shells, source evidence
and composition nodes/edges. The central snapshot and focus selection are shared
across the exported numerical and visual records. The full `.nbo.json` retains
the matrices, provenance and residuals.

## Machine-readable acceptance

The existing local `COV_VALIDATION 1` text plan / JSONL protocol drives real UI
events. New input targets are `nbo.naomo`, `nbo.aonao`, `nbo.naonbo`; composition
controls use the `nbo.focus.*` namespace. Session completion is not a scientific
verdict. An external checker compares source files, raw vendor matrices, literal
NAO report labels, every canonical decomposition, grouped values, graph selection,
export identities and selected actual GPU textures.

The probe remains compatible with its original arguments. Append three explicit
paths for the extended mode:

```text
cov_nbo_probe canonical.fchk analysis.log FILE.47 FILE.37 FILE.49 new-output-directory FILE.51 FILE.52 FILE.53
```

Python entry points are `tests/nbo_focus_suite.py` and
`tests/nbo_focus_negative.py`. They consume a `cases.json` manifest and a build
directory; `--native` on the suite adds real UI checks. Independent reference
libraries (NumPy, IOData, GBasis) are provided separately via the configured
reference-dependency path; they are not replaced by the production parser.
Every invocation uses a fresh evidence directory and never starts chemistry.

The focused real cases are H2O, CH3, CO and Ni(CO)4. The original ten cases remain
the compatibility regression set (twelve distinct molecules in their union).
Fixed geometries and modest basis sets validate the recorded interface and
representation. They do not establish optimised minima, true ground states,
quantitative reaction energetics or all Gaussian/NBO/rank-deficient variants.
Failing data, reports and attempted producer recovery are retained. A general fix
is followed by a targeted retest and affected regressions; final evidence belongs
to one frozen binary/source version. No case is dropped or numerical threshold
relaxed merely to obtain a passing report.
