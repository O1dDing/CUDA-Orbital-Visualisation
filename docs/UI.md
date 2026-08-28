# UI and localisation

The interface keeps the molecular-orbital viewport dominant and treats controls as a compact scientific inspector rather than a desktop form.

## COV Slate

`COV Slate` uses a neutral graphite/slate hierarchy so red/blue orbital phase colours remain visually dominant, with a restrained cool-blue interaction accent and explicit success/error states. The panel overlays the viewport rather than shrinking it.

## Locales and fonts

Runtime locales: English, 简体中文, 日本語 and Français. English is the fallback. User-visible strings belong in the localisation table; stable `##` identifiers ensure a language switch does not change orbital selection or control state. No font binaries are redistributed.

## Molecular representation

The default molecular skeleton is an **enhanced ball-and-stick** presentation: atoms and bonds are intentionally larger than the initial MVP so the molecular framework reads immediately, but the MO surface remains the subject.

Common elements use conventional chemical colours (H white with a dark visual outline, C graphite, N blue, O red, halogens green where defined, P orange, S yellow). The palette should continue towards CPK/Jmol conventions rather than decorative application-specific colours.

`Stick + delocalisation` is a second, orbital-friendly view. Dashed/segmented bonds come from a conservative geometric ring heuristic only. They are a visual aid, not a calculated bond-order assignment. Controls remain available for atom size, bond size, molecule opacity, orbital opacity and hydrogen visibility.

## Orbital browser semantics

Human-facing grouped labels such as `17-a` / `17-b` may represent an energy-degenerate compatible set. They never replace source identity: raw one-based Molden MO numbers and internal zero-based indices remain available in tooltips and machine-readable exports.

Default degeneracy tolerance is `1e-5 Ha`. If meaningful producer symmetry labels disagree, coincident printed energies are not automatically collapsed into one group. High virtual levels are hidden only through non-destructive filters.

## Energy units

Hartree remains the source of truth. UI presentation may be switched consistently between Ha, eV, J/mol, kJ/mol, cal/mol and kcal/mol. Changing units does not launch a CUDA calculation.

## Valence-focused central MO diagram

The automatic MO diagram is deliberately a **central molecular-orbital energy diagram**, not an AO-interaction diagram and not a claimed strict SALC reconstruction.

It is designed around information that ordinary Molden files actually provide:

- MO energy;
- occupation;
- spin;
- producer symmetry label where present;
- raw MO numbering;
- energy-tolerance degeneracy grouping.

### Selection

The default `ValenceCentral` plan is non-destructive. It selects:

1. occupied orbitals classified as valence rather than deep core;
2. low-lying frontier virtual orbitals inside the configured virtual-energy window;
3. the current selected non-core MO if it would otherwise fall outside the compact view;
4. complete degenerate sets when a diagram boundary would otherwise split one.

A capacity derived from the diagram UI control prevents hundreds of levels from being plotted at once. The export records the number shown and the number hidden. Hidden orbitals remain present in the parsed wavefunction and in machine-readable metadata.

### Layout

Energy is quantitative on the vertical axis: lower levels are below higher levels. Degenerate orbitals retain the same y coordinate and expand horizontally around the diagram centre. For example, a twofold set is drawn left/right; threefold sets use left/centre/right; larger sets continue symmetrically.

Labels are allowed to move slightly to prevent collisions, with leader lines back to the quantitative level. The physical level line itself must never be displaced to fake an energy separation.

Each plotted level prioritizes:

- grouped display label;
- electron occupancy;
- producer symmetry;
- energy in the selected unit;
- optional orbital-family / bonding-class annotation only when supported.

Raw MO number and internal index remain secondary/tooltip/export information.

### Orbital-family and bonding annotations

The project separates presentation from scientific provenance.

`σ / π / δ / φ` family labels are currently derived only when producer text explicitly supports such a mapping (for example `sigma` / `pi` in a supplied symmetry/label string). A familiar molecule or a familiar energy pattern is not enough to assign a family.

`bonding / nonbonding / antibonding` is stricter still. Energy and occupation alone are never used to assign bonding character. If producer text does not explicitly support a classification, the value remains `unclassified`.

Machine metadata records the annotation source, confidence, and heuristic flag. This allows human-friendly labels where defensible without concealing uncertainty from downstream analysis.

### Export

`Export diagram + metadata` writes PNG, SVG, JSON and CSV together.

PNG/SVG use a white/light report-friendly central layout. JSON/CSV retain:

- exact raw MO number and internal index;
- grouped display label;
- Hartree source energy and converted display energy;
- occupation and spin;
- producer symmetry;
- degeneracy size;
- core/valence/virtual browsing region;
- whether the MO is included in the compact diagram;
- visibility/selection state;
- optional orbital family and bonding class;
- classification source, confidence and heuristic status.

The current feature **does not claim universal AO→MO reconstruction or strict SALC generation** from ordinary Molden files.

## GPU interaction rules

- changing isovalue: rendering only;
- hover/search/filter/unit/degeneracy-display changes: no CUDA recompute if selection is unchanged;
- selecting a different MO: one frame-end CUDA evaluation;
- changing grid resolution: CUDA recompute;
- language changes: no scientific-state mutation.
