# Initial NBO integration

This implementation attaches an NBO report and its matrix sidecars to an open
Gaussian FCHK wavefunction. Canonical orbitals and NBOs keep separate identities,
selections and GPU evaluators. The NBO report is not a replacement canonical
wavefunction, and NBO diagonal Fock values are not canonical eigenvalues.

## Open an analysis

1. Open the FCHK formatted from the **same Gaussian calculation step that ran NBO**.
   A previous SCF checkpoint may describe rotated canonical orbitals, including
   rotations within degenerate spaces, and is not interchangeable column by column.
2. Scroll to **NBO analysis**. Enter the NBO output (`analysis.log` or standalone
   NBO output), ARCHIVE `.47`, AONBO matrix, and NBOMO matrix paths explicitly.
3. If the output contains several analyses, select its zero-based analysis segment.
   The default `-1` requires a unique analysis and rejects ambiguity.
4. Attach the data. A successful source association enables **NBO orbitals**.
   Select an NBO to render it; return to **Canonical MOs** to restore its selection.

The matrix file numbers are producer choices. The validation inputs request
`AONBO=W37 NBOMO=W49 SAO=W50 ARCHIVE PLOT` and therefore use `FILE.37` and
`FILE.49`; these numbers are not inferred from a filename or the output stem.
The report alone can be inspected and exported, but does not establish renderable
AO coefficients. Printed CMO percentages are incomplete summaries, never silently
promoted to a full transformation.

## Data meaning and source association

- NPA charge, core/valence/Rydberg/total populations and printed spin density;
  NAO occupations, printed diagonal Fock or spin-density columns; NBO labels,
  occupation, local components and optional diagonal Fock values remain distinct.
- Wiberg entries retain ordered atom identities and spin. E2 interactions retain
  donor/acceptor NBO identities, spin, kcal/mol units, printed threshold, energy
  gap and Fock coupling. Missing sections and threshold-omitted entries are not zero.
- NPA for an ECP calculation can include the replaced core. The printed value,
  derived ECP core count and derived explicit-electron population are separate.
- Association checks ordered atoms/effective nuclear charges/geometry, primitive
  basis, complete available canonical coefficient columns, unique AO row mapping,
  AO overlap, electron density and trace. AONBO occupations and NBOMO algebra must
  match. Dimension equality alone does not establish common provenance.
- Some Gaussian checkpoints contain producer total/spin densities and only alpha
  MO coefficients even when the NBO archive contains both spins. If alpha columns
  uniquely establish the AO mapping and the independent producer density checks
  pass for each spin, NBO beta rendering can be density-supported. Direct beta
  FCHK coefficient verification remains unavailable; beta canonical columns in
  that NBO archive are labelled as archive provenance. The canonical FCHK dataset
  is not filled with guessed beta columns, and its NBOMO contribution view only
  maps directly verified canonical columns.

Every parsed row records source path, line, analysis segment and literal source
text. Complete matrix data and AO mapping evidence are included in JSON.

## Export

The NBO export bundle includes `.nbo.json`, `.npa.csv`, `.nao.csv`, `.nbo.csv`,
`.wiberg.csv`, `.e2.csv`, `.e2-sections.csv`, `.view.json`, `.view.svg` and
`.view.png`. Empty optional CSV fields and JSON nulls mean unavailable values.
The image is an occupation view in electrons; it does not put NBO Fock diagonals
on a canonical energy diagram. The view snapshot identifies both the NBO report
selection and the actual displayed orbital set. Canonical-to-NBO contributions
use raw squared NBOMO coefficients, with a display threshold, shown weight and
remaining weight. They are not renormalized to hide omitted entries.

## Machine-readable validation

The local protocol is `COV_VALIDATION 1`: a text plan plus JSON/JSONL evidence,
not a network service. Only the validation executable accepts the plan options.
Actions inject ordinary UI input, using recorded semantic targets to navigate
the real scroll panes and controls.

```text
COV_VALIDATION 1
scene 0.92 0.7 0.25 2.2 0.03 80
seek "panel.nbo"
text "nbo.path" "C:/case/analysis.log"
text "nbo.archive47" "C:/case/FILE.47"
text "nbo.aonbo" "C:/case/FILE.37"
text "nbo.nbomo" "C:/case/FILE.49"
click "nbo.attach"
click "nbo.set.nbo"
click "nbo.orbital.0"
volume "first-nbo" "0"
capture "first-nbo"
export-name "nbo-bundle"
click "nbo.export"
```

```text
cov_validation.exe C:/case/paired-canonical.fchk --validation-plan C:/case/plan.txt --validation-output C:/new-evidence --validation-background
```

`nbo.orbital.<index>` uses the zero-based list index from the exported dataset,
while the producer NBO ID remains one-based and scoped by spin/analysis.
`volume` reads the actual GPU texture and records up to 8192 deterministic sample
indices. `volume_full` additionally saves the complete float32 texture. Captures
contain a framebuffer and UI state/targets; frames record rendered/applied set,
dataset, spin, source index and generation. `session.json` reports command
completion only. Scientific and software acceptance belongs to an external checker.

`tests/nbo_integration_suite.py` compares real producer reports and matrices with
the production parser, exercises native UI/export/switching and compares selected
GPU textures with independent IOData/GBasis values. `tests/nbo_negative.py`
checks rejected and report-only inputs. Both use new output directories and keep
failures. Reference libraries are supplied separately through `COV_REFERENCE_DEPS`.

## Scope

The initial ten fixed-geometry single-point cases cover ordinary bonds, lone pairs,
delocalization, an anion, an open-shell radical, a donor/acceptor adduct,
three-centre bonding and an ECP metal example. These checks establish the recorded
software/data paths and selected sampled textures only. They do not establish
optimized geometries, physical ground states, quantitative E2 chemical conclusions,
every voxel or universal compatibility with all Gaussian/NBO versions. The broader
reference validation program retains its separate status.
