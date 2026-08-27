# CUDA Orbital Visualisation

**Wavefunction in, GPU-resident orbital visualisation out — no CPU volumetric bottleneck.**

CUDA Orbital Visualisation is an experimental GPU-first molecular-orbital viewer for NVIDIA CUDA GPUs, developed first for **RTX 5090 / Blackwell (`sm_120`)**. The current target is small and medium molecular systems of **up to 100 atoms**.

Molden geometry/basis/MO data are parsed on the CPU, uploaded to CUDA, evaluated directly on a 3D grid, written into an OpenGL 3D texture through CUDA/OpenGL interop, and ray-marched on the GPU. The interactive volumetric path deliberately avoids full-field GPU → CPU readback.

## Current scope

- [x] C++20 core
- [x] CUDA 12.8+ build path
- [x] Blackwell `sm_120` cubin + PTX by default
- [x] Streaming/two-pass Molden parser
- [x] `[Atoms]` in Angstrom or atomic units
- [x] `s`, `p`, `sp`, `d`, `f`, `g` shells
- [x] Cartesian `6D / 10F / 15G`
- [x] Real spherical `5D / 7F / 9G`
- [x] Strict MO coefficient/basis-count consistency check
- [x] CUDA orbital grid evaluator
- [x] Direct CUDA surface writes into an OpenGL 3D texture
- [x] GPU-resident `+iso / -iso` ray-marched visualisation
- [x] 64³ / 128³ / 256³ / 512³ grids
- [x] Compact viewport-first `COV Slate` Dear ImGui interface
- [x] Runtime localisation: English / 简体中文 / 日本語 / Français
- [x] System-font CJK fallback without redistributing font binaries
- [x] Windows native Molden Open File dialog + drag/drop + manual path
- [x] Session-local recent files
- [x] Orbital browser with HOMO/LUMO quick navigation and `ImGuiListClipper`
- [x] Configurable degeneracy grouping (`17-a`, `17-b`, …) while retaining raw MO numbering
- [x] Non-destructive Core / Valence / Virtual / high-virtual filters
- [x] Interactive energy-level diagram with electron occupancy glyphs
- [x] Ha / eV / J mol⁻¹ / kJ mol⁻¹ / cal mol⁻¹ / kcal mol⁻¹ display
- [x] Medium-ball + thin-stick molecular representation
- [x] Alternate stick + conservative dashed-delocalisation representation
- [x] MO diagram export: PNG + SVG + JSON + CSV
- [x] CPU CI for parser, localisation, grouping/filtering, molecule heuristics and diagram export
- [ ] Gaussian FCHK
- [ ] CHK → `formchk`
- [ ] Gaussian `.log/.out`
- [ ] Multi-MO fused CUDA evaluator
- [ ] GPU Marching Cubes for mesh export
- [ ] Cube export
- [ ] Cp⁻ numerical regression against trusted Gaussian/Multiwfn reference grids

> **Scientific-status note:** this is still an experimental viewer. Passing UI/parser tests is not equivalent to validating wavefunction values. Producer-specific Molden phase/order conventions, especially spherical shells, still require reference-grid regression before a scientific release.

## Architecture invariant

```text
Molden parser (CPU)
    ↓
basis + selected MO → GPU
    ↓
CUDA evaluates ψ(x,y,z)
    ↓
OpenGL 3D texture stays GPU-resident
    ↓
GPU raymarch finds +iso and -iso
    ↓
display
```

Changing the isovalue, orbital-browser filters, energy units, labels or energy-diagram hover state does **not** recompute the CUDA scalar field. Only an actual MO selection or grid change requests new CUDA work. MO-selection requests are coalesced at frame end so rapid UI interactions evaluate only the latest selection in that frame.

## Human-friendly orbital workflow

The default workflow is intentionally closer to a scientific viewer than a debug panel:

```text
Open wavefunction
      ↓
inspect basis / shell convention
      ↓
filter or search orbital browser
      ↓
HOMO / LUMO / degeneracy-aware selection
      ↕
interactive energy-level diagram
      ↓
GPU orbital viewport
      ↓
optional MO-diagram + metadata export
```

### Degenerate orbitals

Adjacent same-spin MOs within a configurable energy tolerance are presented as one human-readable set. For example, raw Molden **MO17** and **MO18** can appear as **17-a** and **17-b** when they satisfy the tolerance. This never rewrites the wavefunction: tooltips, JSON/CSV export and internal state retain the original one-based MO number and zero-based internal index.

The default tolerance is `1e-5 Ha`; it is intentionally user-adjustable because numerical output conventions differ between producers.

### High virtual orbitals

`Auto · reasonable` is a **view filter**, not deletion. It keeps all occupied orbitals and virtual orbitals up to a configurable energy window above the LUMO. `All`, `Occupied`, `Virtual`, `Core` and `Valence` views are also available.

`Core` versus `Valence` is currently a transparent UI heuristic based on occupation and an energy cutoff; it is not claimed to be a universal chemical partition. Raw MOs are always preserved.

## Molecular representation

The default overlay uses **medium atoms + thin sticks** so the molecular frame is readable without hiding the orbital surface. Controls are provided for atom size, bond width, molecule opacity, orbital opacity and hydrogen visibility.

An alternate **Stick + delocalisation** view uses dashed bonds only when a conservative geometry heuristic detects a compact, near-equal 5–7 member ring composed of common π-capable elements. Molden does not generally provide bond orders, so this is deliberately labelled as a heuristic and does not claim a resonance assignment when the evidence is insufficient.

## Energy-level and MO diagrams

The interactive energy-level diagram is synchronized with the orbital browser and 3D selection. It displays occupancy using electron arrows and emphasizes HOMO, LUMO and the selected level. When sufficiently complete symmetry labels are supplied by the producer, complex systems can be organized into symmetry lanes.

For complex systems without sufficient point-group/symmetry information, the program reports **“SALC not confidently available”** and falls back to reliable energy/occupation/degeneracy ordering. A general strict SALC derivation requires point-group operations and AO transformation information that a generic Molden file does not guarantee; this program does not fabricate it.

The export action writes four siblings next to the wavefunction:

```text
calculation.mo.png
calculation.mo.svg
calculation.mo.json
calculation.mo.csv
```

PNG/SVG are human-facing diagrams; JSON/CSV preserve machine-readable raw MO number, grouped label, Hartree energy, selected display-unit value, occupation, spin, symmetry, degeneracy size, region and filtering state.

## Energy units

Orbital energies can be displayed consistently as:

- Hartree (Ha)
- electronvolt (eV)
- J/mol
- kJ/mol
- cal/mol
- kcal/mol

Conversions use one Hartree as the internal source of truth; unit switching never modifies parsed energies.

## UI and localisation

The control surface is a dark inspector over the molecular viewport. All user-facing controls added by the orbital workflow are translated at the same time in:

- English
- 简体中文
- 日本語
- Français

No font files are committed. The application builds an ImGui atlas from local system fonts and only requests the CJK glyphs needed by bundled strings. On Windows it prefers Segoe UI + Microsoft YaHei + Yu Gothic/Meiryo. See [`docs/UI.md`](docs/UI.md) for UI semantics and heuristic boundaries.

## Requirements

### GPU build

- NVIDIA GPU with CUDA support
- **CUDA Toolkit 12.8+**
- CMake 3.28+
- C++20 compiler
- OpenGL 2.1+ compatibility context
- Git/network access for the pinned GLFW and Dear ImGui source archives

Default CUDA target:

```text
sm_120 native cubin
compute_120 PTX
```

Override `CMAKE_CUDA_ARCHITECTURES` for another CUDA GPU.

## Build

### Windows / Visual Studio 2022

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --parallel
.\build\Release\cov.exe .\examples\h2.molden
```

### Linux

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/cov ./examples/h2.molden
```

### CPU/core tests without CUDA

```bash
cmake -S . -B build-cpu -DCOV_ENABLE_CUDA=OFF -DCOV_BUILD_TESTS=ON
cmake --build build-cpu --parallel
ctest --test-dir build-cpu --output-on-failure
```

The CPU suite covers the parser plus degeneracy labels, HOMO/LUMO logic, energy conversions, high-virtual filtering, MO-diagram export/PNG structure, conservative delocalised-ring detection and four-language UI/font-atlas strings.

## Usage

Start with:

```bash
cov calculation.molden
```

or launch the application and use **Open File…**, manual path entry, or drag a Molden wavefunction onto the viewport. The native Windows picker intentionally advertises Molden-like extensions only; a generic MDL `.mol` file is not presented as supported input.

Core interactions:

- **Left mouse drag:** orbit camera
- **Mouse wheel:** zoom
- **Orbital browser / energy diagram:** select MO
- **HOMO-1 / HOMO / LUMO / LUMO+1:** fast frontier navigation
- **Isovalue:** immediate shader-only threshold change
- **Grid:** 64³, 128³, 256³ or 512³; grid changes recompute CUDA
- **Molecule style:** medium ball + stick or stick + dashed-delocalisation heuristic
- **Energy unit:** global display-unit selection
- **Language:** English / 简体中文 / 日本語 / Français at runtime

Positive and negative orbital phases are rendered simultaneously with different colours.

## Molden parser rules

The parser intentionally fails hard on inconsistent basis data:

```text
n_basis_from_shells == n_coefficients_per_MO
```

A coefficient index inconsistent with the derived shell convention is an error, not something to guess or silently repair. Large Molden files are read in two streaming passes rather than loading the entire text into one giant string.

## CUDA evaluator

For each grid point:

\[
\psi_i(\mathbf r)=\sum_\mu C_{\mu i}\chi_\mu(\mathbf r)
\]

The kernel does not allocate an `Ngrid × Nbasis` AO matrix. The current evaluator computes one selected MO per pass; the next performance milestone is fused multi-MO evaluation so AO work can be reused across 4/8/16 MOs.

## First regression target: cyclopentadienyl anion

Cp⁻ (`C5H5−`) at `wB97XD/aug-cc-pVTZ` remains the first serious numerical reference:

- spherical `5D/7F`: 345 basis functions
- Cartesian `6D/10F`: 400 basis functions
- 18 occupied MOs
- Cartesian MO16 ≈ −0.236093 Ha, `a2''` π
- Cartesian MO17 ≈ MO18 ≈ −0.068048 Ha, `e1''` π

The final scientific gate is pointwise comparison against trusted Gaussian/Multiwfn values on identical grid points with maximum/RMS error reporting.

## Repository layout

```text
include/cov/              public C++ interfaces
src/parser/               Molden parser
src/cuda/                 CUDA orbital evaluator
src/render/               OpenGL raymarch + molecule display analysis
src/orbitals/             orbital semantics, filtering and MO diagram export
src/ui/                   theme, localisation, browser and energy-level UI
src/platform/             native file-dialog integration
src/main.cpp              GLFW + Dear ImGui application shell
tests/                    CPU parser/orbital/UI/export smoke tests
docs/UI.md                UI and semantic design notes
examples/                 tiny checked-in examples only
.github/workflows/        CPU/core CI
```

## Near-term roadmap

1. Complete Windows/RTX 5090 end-to-end validation for this UI branch.
2. Validate Cartesian s/p/d/f/g values against reference cubes.
3. Validate spherical 5D/7F/9G ordering, phase and normalisation.
4. Add Cp⁻ GPU-vs-reference maximum/RMS regression.
5. Add fused multi-MO evaluation.
6. Add FCHK, then CHK → `formchk` integration.
7. Add GPU Marching Cubes only for explicit mesh export.
8. Add Cube export.

## Licence

Apache License 2.0.
