# UI and localisation

The interface keeps the molecular-orbital viewport dominant and treats controls as a compact scientific inspector rather than a desktop form.

## COV Slate

`COV Slate` uses a neutral graphite/slate hierarchy so red/blue orbital phase colours remain visually dominant, with a restrained cool-blue interaction accent and explicit success/error states. The panel overlays the viewport rather than shrinking it.

Open-source UI references remain Dear ImGui, ImThemes, dear-imgui-styles and the Blender Human Interface Guidelines. No third-party theme source, icons or font binaries are copied into the repository.

## Locales and fonts

Runtime locales:

| Locale | Display name |
| --- | --- |
| `en` | English |
| `zh-Hans` | 简体中文 |
| `ja` | 日本語 |
| `fr` | Français |

English is the fallback. User-visible strings must remain in the localisation table; widget state uses stable `##` identifiers so a language switch does not change orbital selection or control state.

The repository does not redistribute fonts. Windows prefers Segoe UI, Microsoft YaHei and Yu Gothic/Meiryo; Linux/macOS use suitable system/Noto fallbacks. The font atlas contains only glyphs required by bundled UI strings.

## Molecular representation

The default molecular skeleton is an **enhanced ball-and-stick** presentation: atoms and bonds are intentionally larger than the initial MVP so the molecular framework reads immediately, but the MO surface remains the subject.

Common elements use conventional chemical colours (H white with a dark visual outline, C graphite, N blue, O red, halogens green where defined, P orange, S yellow). The palette should continue towards CPK/Jmol conventions rather than decorative application-specific colours.

`Stick + delocalisation` is a second, orbital-friendly view. Dashed/segmented bonds come from a conservative geometric ring heuristic only. They are a visual aid, not a calculated bond-order assignment.

Controls remain available for atom size, bond size, molecule opacity, orbital opacity and hydrogen visibility.

## Orbital browser semantics

Human-facing grouped labels such as `17-a` / `17-b` may represent an energy-degenerate compatible set. They never replace source identity: raw one-based Molden MO numbers and internal zero-based indices remain available in tooltips and machine-readable exports.

Default degeneracy tolerance is `1e-5 Ha`. If meaningful producer symmetry labels disagree, coincident printed energies are not automatically collapsed into one group.

High virtual levels are hidden only through non-destructive filters. `Auto · reasonable`, `All`, `Occupied`, `Virtual`, `Core` and `Valence` are browsing modes; they do not mutate wavefunction data.

## Energy units

The scientific source of truth remains Hartree. UI presentation may be switched consistently between Ha, eV, J/mol, kJ/mol, cal/mol and kcal/mol. Changing presentation units does not launch a CUDA orbital calculation.

## MO diagram presentation

The interactive energy-level diagram and exported report diagrams serve different visual needs. The GUI can remain compact/dark; PNG/SVG exports are light-background, print-friendly figures.

### Textbook diatomic mode

For a small homonuclear diatomic such as H₂, the exported figure uses the familiar textbook structure: atomic-orbital levels on the sides, molecular-orbital levels in the centre, interaction lines, an energy axis and direct electron occupancy.

H₂ therefore shows side `H 1s` levels, central `σ1s` and `σ*1s`, AO→MO interaction lines, an energy axis and electron arrows. The side AO position is explicitly qualitative unless isolated-atom AO energies are actually available; central MO energies come from the parsed wavefunction.

### Complex systems and SALC honesty

A complex calculation with useful producer symmetry labels can be shown as a **symmetry-grouped MO diagram** with separate symmetry lanes, electron occupancy and collision-managed labels.

A strict universal SALC diagram is a stronger scientific claim. Ordinary Molden data generally do not provide the point-group operations and basis-function transformation matrices required to construct arbitrary SALCs rigorously. When that information is not available, the UI/export must say **`SALC not confidently available`** rather than invent a SALC decomposition.

The supported hierarchy is therefore:

1. textbook/simple MO diagram where scientifically defensible;
2. symmetry-grouped MO diagram when producer symmetry information is sufficient;
3. explicit SALC-unavailable fallback otherwise.

## Diagram export and machine readability

`Export diagram + metadata` writes PNG, SVG, JSON and CSV together. PNG/SVG prioritize human readability: white/light background, adequate lane spacing, direct electron occupancy and collision-managed labels. JSON/CSV retain exact raw MO number, internal index, grouped display label, Hartree energy, converted display energy, occupation, spin, symmetry, degeneracy size, browsing-region heuristic, visibility and selection.

The diagram generator may move **labels** to avoid collisions; it must not move the underlying physical energy level to fake separation. Near-degenerate levels may be given a small horizontal offset for click/readability while their source energies remain unchanged.

## GPU interaction rules

The UI must not turn browsing into GPU work:

- changing isovalue: rendering only;
- hover/search/filter/unit/degeneracy-display changes: no CUDA recompute if selection is unchanged;
- selecting a different MO: one final frame-end CUDA evaluation;
- changing grid resolution: CUDA recompute;
- language changes: no scientific-state mutation.

These rules preserve the project’s GPU-resident architecture while keeping the UI responsive and machine-readable state deterministic.
