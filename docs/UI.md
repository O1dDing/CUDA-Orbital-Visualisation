# UI and localisation

The interface keeps the molecular-orbital viewport dominant and treats controls as a compact scientific inspector rather than a desktop form.

## COV Slate

`COV Slate` uses a neutral graphite/slate hierarchy so red/blue orbital phase colours remain visually dominant, with a restrained cool-blue interaction accent and explicit success/error states. The control panel and molecular scene occupy separate rectangles. Molecular geometry and orbital surfaces share the same scene viewport, projection and depth buffer, including after window resizing and display scaling.

## Locales and fonts

Runtime locales: English, 简体中文, 日本語 and Français. English is the fallback. User-visible strings belong in the localisation table; stable `##` identifiers ensure a language switch does not change orbital selection or control state. No font binaries are redistributed.

## Molecular representation

The default molecular skeleton is an **enhanced ball-and-stick** presentation: atoms and bonds are intentionally larger than the initial MVP so the molecular framework reads immediately, but the MO surface remains the subject.

Common elements use conventional chemical colours (H white with a dark visual outline, C graphite, N blue, O red, halogens green where defined, P orange, S yellow). The palette should continue towards CPK/Jmol conventions rather than decorative application-specific colours.

`Stick + delocalisation` is a second view with less atom clutter. Interaction styles distinguish ordinary bonds, coordination, multicentre and cage support, hydrogen bonds, noncovalent contacts and ionic interactions according to the available interaction graph. A separate geometric ring heuristic adds delocalisation strokes. These visual styles do not by themselves establish a bond order or a physical assignment. Controls remain available for atom size, bond size, molecule opacity, orbital opacity and hydrogen visibility.

## Orbital browser semantics

Human-facing grouped labels such as `17-a` / `17-b` may represent an energy-degenerate compatible set. They never replace source identity: raw one-based source MO numbers and internal zero-based indices remain available in tooltips and machine-readable exports.

Default degeneracy tolerance is `1e-5 Ha`. If meaningful producer symmetry labels disagree, coincident printed energies are not automatically collapsed into one group. High virtual levels are hidden only through non-destructive filters.

## Energy units

Hartree remains the source of truth. UI presentation may be switched consistently between Ha, eV, J/mol, kJ/mol, cal/mol and kcal/mol. Changing units does not launch a CUDA calculation.

## Valence-focused central MO diagram

The automatic MO diagram is deliberately a **central molecular-orbital energy diagram**, not an AO-interaction diagram and not a claimed strict SALC reconstruction.

It combines parsed wavefunction data with explicitly identified producer metadata and derived model results. The basic orbital identity uses:

- MO energy;
- occupation;
- spin;
- producer symmetry label where present;
- raw MO numbering;
- energy-tolerance degeneracy grouping.

### Selection

The `ValenceCentral` plan is non-destructive. Compact mode first considers an eligible haptic metal-pi family; otherwise an available ligand-field environment uses the valence-central view. Other inputs may use a supported multicentre or delocalised-pi active space. The displayed and exported mode describes the diagram actually built, including a fallback when a requested active space is unavailable.

The valence-central plan selects:

1. occupied orbitals classified as valence rather than deep core;
2. low-lying frontier virtual orbitals inside the configured virtual-energy window;
3. complete compatible groups when a diagram boundary would otherwise split one.

The compact diagram stays stable when the inspected MO changes. An MO outside the compact view can still be selected in the browser and inspected in the scene and details; selection does not insert, remove or reprioritise compact rows. Selection may affect the neighbourhood in the expanded view.

A capacity derived from the diagram UI control prevents hundreds of levels from being plotted at once. The export records the number shown and the number hidden. Hidden orbitals remain present in the parsed wavefunction and in machine-readable metadata.

### Layout

The linear energy axis preserves quantitative energy spacing. The nonlinear focus mode expands crowded regions with a labelled transform: its pixel gaps are not proportional to energy differences. Exact source member energies remain unchanged. A grouped row can use a representative mean energy and reports its spread; grouping is not a claim that all member energies are exactly equal. Levels are placed at the coordinates defined by the chosen transform, with members spread horizontally.

Crowded display footprints receive separate horizontal lanes. The interactive canvas can scroll horizontally while the control panel stays in place. Export labels can also move vertically to avoid collisions, with leader lines back to their levels; the level lines retain the energy-axis coordinates.

Each plotted level prioritizes:

- grouped display label;
- electron occupancy;
- producer symmetry;
- energy in the selected unit;
- optional orbital-family / bonding-class annotation only when supported.

Raw MO number and internal index remain secondary/tooltip/export information.

### Orbital-family and bonding annotations

The detail panel identifies where an annotation comes from: producer text, a computed model or a candidate explanation. Derived `sigma / pi / delta / phi` and bonding classifications can use the available orbital-chemistry decomposition; familiar molecule names, energies and occupations alone do not establish those assignments. Missing evidence remains unavailable or unclassified. Source and support fields are retained in machine metadata.

### Details window and narrow panels

The details window keeps its title and close controls within the current viewport when the main window shrinks. It can be moved, resized, scrolled, closed and reopened. These layout changes retain the inspected source member and its scientific values. Long coloured chemistry fields wrap before the sidebar edge so a continuation does not become a narrow vertical column.

### Whole-orbital and local symmetry

Whole-orbital symmetry, local-centre explanations and candidate labels have separate scopes. Each explanation retains its target MO members, alpha/beta identities, point-group context and source. A local label does not replace an unavailable whole-orbital label. Gaussian's detected point group and reported Abelian subgroup remain producer context; a missing or invalid printed field is not turned into a point-group assignment.

The local metric projection presents three different percentages: centre coverage, label purity within the assessed local shell, and the labelled component of the complete target. A high local-shell purity can coexist with a small contribution to the full orbital. These quantities describe the recorded projection and local frame; they are not probabilities that a physical symmetry assignment is correct. If no single label is resolved, the available angular decomposition and residual remain visible.

### Metal-ligand details

The panel labels group analysis separately from the selected member identity. Its metal populations, ligand populations, channel fractions and overlap have independent availability checks:

- **Not applicable:** the source does not support this metal-ligand analysis, for example an organic molecule without a metal centre. Numeric rows for that scope are omitted.
- **Unavailable:** the scope could apply, but required member populations, first-shell pairs or channel denominators are missing or unresolved. The affected quantity is shown as unavailable.
- **Determined:** an actual computed zero remains a zero. Missing values are never presented as computed zeros.

### Energy gaps and topology

Crystal-field gaps and pi-partner splittings have different types and retain their participating levels, local environment, spin and source evidence. Their displayed support is a heuristic score, not a probability or an independent chemical-state determination. Catalogue ligand priors and the orbital evidence are shown separately, including a contradiction between them.

A delocalised-pi family and a topology assignment are separate claims. Topology details retain the contributing channels and paths through the molecular skeleton. A spiro assignment requires the recorded skeleton-ring witness; multiple orientation channels alone are insufficient. These graph and local-model descriptions do not establish a chemical mechanism or a physical ground state.

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

The current feature **does not claim universal AO→MO reconstruction or strict SALC generation** from an arbitrary parsed wavefunction.

## GPU interaction rules

- changing isovalue: rendering only;
- hover/search/filter/unit/degeneracy-display changes: no CUDA recompute if selection is unchanged;
- selecting a different MO: one frame-end CUDA evaluation;
- changing grid resolution: CUDA recompute;
- language changes: no scientific-state mutation.
