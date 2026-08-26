# UI and localisation

The first UI pass keeps the molecular-orbital viewport as the dominant surface and treats the controls as a compact scientific inspector rather than a desktop form.

## COV Slate

`COV Slate` is the project's own Dear ImGui theme. It uses:

- a neutral slate/graphite hierarchy so orbital phase colours remain visually dominant;
- one cool-blue interaction accent;
- green for healthy GPU-resident state and red only for errors;
- modest 6–10 px rounding, restrained borders and medium-density spacing;
- grouped cards for File, Wavefunction, Orbital, Rendering and Performance;
- an overlay inspector so the 3D viewport remains full-window;
- stable `##` ImGui identifiers so changing language does not change widget identity.

The design was synthesised after reviewing open-source tooling conventions rather than copying a theme implementation or asset.

### Open-source references

- Dear ImGui — <https://github.com/ocornut/imgui>
- ImThemes theme browser/editor — <https://github.com/Patitotective/ImThemes>
- dear-imgui-styles collection — <https://github.com/GraphicsProgramming/dear-imgui-styles>
- Blender Human Interface Guidelines — <https://developer.blender.org/docs/features/interface/human_interface_guidelines/>

ImThemes and dear-imgui-styles were useful mainly for comparing control density, contrast and rounding. Blender's UI guidance reinforced the decision to keep a dark, low-distraction viewport-first layout. No third-party theme source, icons or font binaries are copied into this repository.

## Locales

The runtime UI currently ships four locales:

| Locale | Display name |
| --- | --- |
| `en` | English |
| `zh-Hans` | 简体中文 |
| `ja` | 日本語 |
| `fr` | Français |

English remains the fallback language. Translation keys are strongly enumerated in `include/cov/ui.hpp`, while the UTF-8 string tables live in `src/ui/ui.cpp`.

The program is compiled with `/utf-8` for MSVC C++ sources so the literals have deterministic encoding on Windows.

## Font policy

The repository deliberately does **not** redistribute font binaries.

At startup, the application builds an ImGui font atlas from locally installed system fonts. On Windows it prefers:

- Segoe UI for the Latin UI;
- Microsoft YaHei for Simplified Chinese;
- Yu Gothic, Meiryo or MS Gothic for Japanese.

Linux and macOS have Noto/system-font fallbacks. Because this project currently pins Dear ImGui 1.90.9, CJK glyph ranges are built explicitly from the strings actually used by the Chinese and Japanese localisations with `ImFontGlyphRangesBuilder`. This avoids allocating a huge all-CJK atlas while still allowing instant runtime language switching.

The selected font chain is shown in the Performance card so missing CJK fallback fonts are visible during testing.

## Layout

The control panel is intentionally narrow and scrollable:

```text
┌──────────────────────────────┐
│ CUDA Orbital Visualisation   │  Language
│ GPU-first ...                │
│ ● Ready                      │
├──────────────────────────────┤
│ FILE                         │
│ [ Molden path .... ] [Load]  │
├──────────────────────────────┤
│ WAVEFUNCTION                 │
│ Atoms                 10/100 │
│ Shells                   ... │
│ Basis functions          ... │
├──────────────────────────────┤
│ ORBITAL                      │
│ [ MO slider ]                │
│ Energy / occupation / spin   │
├──────────────────────────────┤
│ RENDERING                    │
│ Isovalue / grid              │
│ [Recompute] [Reset camera]   │
├──────────────────────────────┤
│ PERFORMANCE                  │
│ RTX 5090 / kernel time       │
│ ● GPU resident               │
└──────────────────────────────┘
```

The panel overlays the viewport instead of shrinking it. The viewport therefore retains the full framebuffer for CUDA/OpenGL rendering.

## Rules for future UI work

1. Scientific state must be explicit: MO number, energy, occupation, spin, symmetry and grid resolution must never be hidden behind decorative UI.
2. Changing isovalue stays a rendering-only interaction; it must not trigger CUDA grid recomputation.
3. Language changes must not alter internal widget IDs or scientific state.
4. Do not commit proprietary/system fonts; detect them at runtime.
5. Prefer text and restrained colour over decorative icons until an icon set has a clear functional purpose and compatible licence.
6. Keep controls compact enough that a 1440×900 window leaves most of the viewport unobstructed.
