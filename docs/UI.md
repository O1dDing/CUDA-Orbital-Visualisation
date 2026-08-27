# UI, orbital semantics and localisation

The UI is deliberately **human-first, machine-readable underneath**. The molecular-orbital viewport remains the dominant visual surface; controls expose scientific state without turning the application into a raw debug panel.

## COV Slate

`COV Slate` is the project's own Dear ImGui theme. It uses a neutral graphite/slate hierarchy so the red/blue orbital phases remain visually dominant, a single cool-blue interaction accent, restrained rounding/borders, and green/red only for healthy/error state.

The design was synthesised after reviewing open-source UI conventions rather than copying a theme implementation or asset:

- Dear ImGui — <https://github.com/ocornut/imgui>
- ImThemes — <https://github.com/Patitotective/ImThemes>
- dear-imgui-styles — <https://github.com/GraphicsProgramming/dear-imgui-styles>
- Blender Human Interface Guidelines — <https://developer.blender.org/docs/features/interface/human_interface_guidelines/>

No third-party theme source, icons or font binaries are copied into this repository.

## Locales and stable IDs

Runtime locales:

| Locale | Display name |
| --- | --- |
| `en` | English |
| `zh-Hans` | 简体中文 |
| `ja` | 日本語 |
| `fr` | Français |

User-facing strings are strongly enumerated in `include/cov/ui.hpp` and stored as four-language rows in `src/ui/ui.cpp`. Internal ImGui identities use stable `##...` IDs, so changing language does not reset MO selection, filter state, sliders or path input.

MSVC compiles the application/UI sources as UTF-8.

## Font policy

The repository does **not** redistribute fonts. The application uses system fonts at runtime:

- Windows Latin: Segoe UI / Arial
- Simplified Chinese: Microsoft YaHei / SimHei fallback
- Japanese: Yu Gothic / Meiryo / MS Gothic fallback
- Linux/macOS: Noto/system fallbacks where available

The font atlas is built from the glyphs actually present in the four bundled locales, rather than adding the entire CJK Unicode range. This keeps the atlas compact while supporting instant language switching.

## Information architecture

The scrollable overlay inspector now follows the user's workflow rather than implementation modules:

```text
┌──────────────────────────────────┐
│ CUDA Orbital Visualisation       │ Language
│ ● status                         │
├──────────────────────────────────┤
│ FILE                             │
│ [Open File…] current file        │
│ manual path / recent files       │
├──────────────────────────────────┤
│ WAVEFUNCTION                     │
│ atoms / shells / basis / MOs     │
├──────────────────────────────────┤
│ ORBITAL BROWSER                  │
│ HOMO-1 · HOMO · LUMO · LUMO+1   │
│ search / filter / units          │
│ degeneracy tolerance             │
│ scrollable MO table              │
├──────────────────────────────────┤
│ ENERGY-LEVEL DIAGRAM             │
│ occupied ↑↓ / virtual / selected │
│ clickable levels                 │
│ [export diagram + metadata]      │
├──────────────────────────────────┤
│ RENDERING                        │
│ molecule style / sizes / opacity │
│ isovalue / grid                  │
├──────────────────────────────────┤
│ PERFORMANCE                      │
│ CUDA device / kernel time        │
│ ● GPU resident                   │
└──────────────────────────────────┘
```

The panel overlays rather than reallocates the viewport framebuffer.

## Orbital numbering and degeneracy

Raw Molden/user numbering is one-based; the internal C++ index remains zero-based.

Degenerate sets are detected from **adjacent same-spin orbital energies** within a configurable tolerance. The default is `1e-5 Ha`. A two-member set beginning at raw MO17 is displayed as:

```text
17-a    raw MO 17
17-b    raw MO 18
```

Larger sets continue `17-c`, `17-d`, etc. The grouped label is presentation metadata only: raw MO numbers and internal indices remain available in tooltips and JSON/CSV export.

Energy-only degeneracy detection cannot prove group-theoretical degeneracy. The tolerance is therefore visible and adjustable; near-degenerate levels outside it are not silently grouped.

## Frontier navigation and open-shell behaviour

HOMO/LUMO helpers use occupation rather than assuming every occupied orbital has occupation 2. If beta-spin orbitals are present, quick navigation resolves the frontier within the selected spin set. The browser retains spin and occupation columns/metadata so an unrestricted calculation is not presented as a single closed-shell ladder.

## Orbital filters

Filters never delete or rewrite wavefunction data.

- **All** — raw complete MO list
- **Occupied** — occupation above the configured threshold
- **Virtual** — unoccupied levels
- **Core** — occupied + below the current core-energy heuristic
- **Valence** — occupied + above the core-energy heuristic
- **Auto · reasonable** — all occupied levels plus virtual levels within a configurable window above LUMO

The current Core/Valence split is a UI heuristic, not a population/localisation analysis. The high-virtual filter is also purely presentational. Machine metadata records whether each raw orbital was visible under the active filter.

## Energy units

Hartree is the immutable source unit. UI and export display conversion supports:

- Ha
- eV
- J/mol
- kJ/mol
- cal/mol
- kcal/mol

The browser and energy diagram share one unit selection. JSON always keeps `energy_hartree` alongside the selected display-unit value.

## Energy-level diagram

The diagram is drawn with `ImDrawList`, avoiding a plotting-framework dependency. It displays the selected neighbourhood rather than compressing hundreds of levels into one unreadable panel.

Visual semantics:

- occupied ordinary level: light neutral
- virtual: muted slate
- HOMO: green accent
- LUMO: amber accent
- selected: blue accent
- electron occupancy: up/down arrows
- degenerate members: small horizontal offsets so each remains clickable

Clicking a level requests the same orbital selection used by the browser. Hovering shows raw MO, grouped label, energy, occupation, spin, symmetry and degeneracy size. Hovering never launches CUDA.

## SALC / symmetry grouping boundary

A generic Molden file does not guarantee enough information to derive a strict general SALC basis: robust SALC construction needs the molecular point group, symmetry operations and the transformation of the relevant AO/basis functions under those operations.

Therefore the application uses three explicit diagram modes:

1. **Simple ordering** — suitable for small/simple systems.
2. **Symmetry-grouped** — complex system where producer-supplied MO symmetry labels have sufficient coverage; these labels can form diagram lanes.
3. **SALC not confidently available** — complex system without enough trustworthy symmetry information. The program falls back to energy/occupation/degeneracy ordering and states why.

The third state is intentional. It is preferable to a visually impressive but chemically fabricated SALC diagram.

## Automatic MO-diagram export

The export bundle writes:

```text
<calculation>.mo.png
<calculation>.mo.svg
<calculation>.mo.json
<calculation>.mo.csv
```

PNG and SVG display orbital levels and electron occupancy. SVG preserves full text labels; PNG uses a deliberately dependency-free compact bitmap annotation layer. JSON/CSV provide exact machine-readable orbital state, including raw number, grouped label, Hartree/display energy, occupation, spin, symmetry, degeneracy size, Core/Valence/Virtual region, visible/filter state and selected state.

## Molecular representation

### Medium ball + thin stick (default)

The default molecular overlay deliberately increases atom readability compared with the original tiny-point renderer while keeping bonds thin so the orbital surface remains the primary visual object. Controls expose atom scale, bond width, molecular opacity, orbital opacity and hydrogen visibility.

### Stick + delocalisation

Molden does not generally carry bond order. Dashed bonds therefore use a **conservative geometry heuristic**, not an asserted Lewis/resonance structure. A candidate must belong to a compact 5–7 member ring of common π-capable elements and the normalized ring bond lengths must be both short enough and near-equal. If those conditions are not met, the program draws a normal solid stick.

This boundary must remain visible in the UI/documentation: dashed = “geometry suggests a delocalised ring”, not “bond order has been proven”.

## CUDA interaction rules

1. Browser hover does not compute an orbital.
2. Energy-diagram hover does not compute an orbital.
3. Search, filtering, unit changes and degeneracy-label changes do not compute an orbital.
4. An actual selected-MO change requests recomputation.
5. Selection requests are coalesced at frame end so only the latest selection in that frame is evaluated.
6. Isovalue changes remain shader-only and do not recompute the CUDA field.
7. Grid-resolution changes resize/re-register the texture and recompute the field.
8. The interactive volumetric path remains GPU-resident; the MO-diagram export does not read back the orbital volume.

## Rules for future UI work

1. Human-readable labels must never destroy raw scientific identifiers.
2. Heuristics must be named as heuristics and remain non-destructive.
3. Do not claim strict SALC/group-theoretical results unless the required symmetry transformations are actually available and validated.
4. Keep raw MO number, energy, occupation, spin and symmetry available to both humans and machine export.
5. Do not commit proprietary/system fonts.
6. Keep the 3D viewport dominant at ordinary desktop resolutions.
7. A green CI check is not a substitute for Windows/CUDA/OpenGL visual validation.
