# CUDA Orbital Visualisation

**Wavefunction in, GPU-resident orbital visualisation out — no CPU volumetric bottleneck.**

CUDA Orbital Visualisation is an experimental GPU-first molecular-orbital viewer for NVIDIA CUDA GPUs, developed first for **RTX 5090 / Blackwell (`sm_120`)**. The initial target is small and medium molecular systems of **up to 100 atoms**.

The viewer parses Molden wavefunctions on the CPU, uploads basis data to CUDA, evaluates a selected molecular orbital on a 3D grid, writes the scalar field **directly into an OpenGL 3D texture through CUDA/OpenGL interop**, and ray-marches the positive and negative isosurfaces on the GPU.

## Current scope

- [x] C++20 core
- [x] CUDA 12.8+ build path
- [x] Blackwell `sm_120` cubin + PTX by default
- [x] Streaming/two-pass Molden parser
- [x] Cartesian and real-spherical `s/p/d/f/g` basis support
- [x] Strict MO coefficient/basis-count consistency checks
- [x] CUDA orbital grid evaluator
- [x] GPU-resident `+iso / -iso` ray-marched visualisation
- [x] 64³ / 128³ / 256³ / 512³ grids
- [x] Viewport-first Dear ImGui UI
- [x] English / 简体中文 / 日本語 / Français
- [x] Native Windows Open File workflow + drag/drop/manual path
- [x] Orbital browser, search, HOMO/LUMO navigation and filters
- [x] Degeneracy-aware human labels such as `17-a`, `17-b` while retaining raw MO numbers
- [x] Ha / eV / J mol⁻¹ / kJ mol⁻¹ / cal mol⁻¹ / kcal mol⁻¹ display
- [x] Enhanced ball-and-stick defaults plus stick/delocalisation display mode
- [x] Interactive energy-level diagram with electron occupancy
- [x] Print-friendly MO diagram export: PNG + SVG + JSON + CSV
- [x] Textbook-style homonuclear-diatomic AO → MO presentation
- [x] Symmetry-grouped complex-system fallback without fabricating strict SALCs
- [ ] Universal strict SALC construction
- [ ] Gaussian FCHK
- [ ] CHK → `formchk`
- [ ] Gaussian `.log/.out`
- [ ] Multi-MO fused CUDA evaluator
- [ ] GPU Marching Cubes for mesh export
- [ ] Cube export
- [ ] Cp⁻ numerical regression against Gaussian/Multiwfn reference grids

> **Scientific-status note:** this is an early project, not yet a validated production quantum-chemistry package. Producer-specific Molden phase/order conventions and the spherical-harmonic path still require reference-grid regression before scientific release.

## Architecture

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

Changing the isovalue does **not** recompute the CUDA grid or rebuild a mesh.

## Human-facing orbital workflow

The UI is designed around the sequence:

```text
Open calculation
→ inspect wavefunction
→ browse orbitals
→ select MO
→ inspect energy / occupation / spin / symmetry
→ view orbital
→ generate/export an MO diagram
```

The orbital browser keeps raw one-based Molden numbering for machine/scientific traceability while optionally presenting compact degeneracy groups. With the default `1e-5 Ha` tolerance, a compatible two-level degenerate set may be presented as `17-a` and `17-b`; the raw source MO numbers remain available in tooltips and exported metadata.

High virtual orbitals are never deleted. The default `Auto · reasonable` view is a non-destructive filter; `All` restores the complete parsed orbital list. `Core`, `Valence`, and high-virtual hiding are human-facing heuristics rather than new quantum-chemical assignments.

## Molecular display

The default skeleton is deliberately stronger than the first MVP: larger atoms and thicker bonds make the molecule immediately legible while retaining orbital surfaces as the visual subject. A second `Stick + delocalisation` mode uses a conservative geometric ring heuristic for dashed bonds. Those dashed bonds are a visual aid and **do not claim an ab-initio bond order**.

The renderer uses conventional chemical element colours for common atoms (e.g. H white, C graphite, N blue, O red, halogens green where defined, P orange, S yellow) and will continue to expand its CPK/Jmol-style palette.

## MO diagrams and SALC boundary

There are deliberately different presentation levels rather than one diagram claiming to solve every symmetry problem.

### Textbook diatomic presentation

For small homonuclear diatomics, the exporter can render a textbook-style AO → MO diagram: atomic-orbital levels on the sides, molecular-orbital levels in the centre, interaction lines, an energy axis, and electron arrows. For H₂ this gives the familiar `1s + 1s → σ1s / σ*1s` presentation.

The side AO level is explicitly **qualitative** unless an actual isolated-atom AO reference energy is available. Molecular-orbital energies shown next to the central levels are quantitative values from the parsed wavefunction.

### Complex systems

If a complex calculation carries sufficiently useful producer symmetry labels, the program may render a **symmetry-grouped MO diagram**. If the data are insufficient for a rigorous SALC derivation, the program explicitly reports **`SALC not confidently available`** and falls back to reliable energy/occupation/degeneracy/symmetry information.

A strict universal SALC builder requires point-group operations plus basis-function transformation information that ordinary Molden data do not generally provide. The project therefore does not fabricate SALCs from energies alone.

### Export

`Export diagram + metadata` writes:

```text
calculation.mo.png
calculation.mo.svg
calculation.mo.json
calculation.mo.csv
```

PNG/SVG are light-background, report/textbook-oriented views. JSON/CSV retain machine-readable raw MO number, internal index, grouped label, Hartree energy, selected display-unit energy, occupation, spin, symmetry, degeneracy size, orbital-region heuristic, visibility and selection state.

## Requirements

### GPU build

- NVIDIA GPU with CUDA support
- **CUDA Toolkit 12.8+**
- CMake 3.28+
- C++20 compiler
- OpenGL 2.1+ compatibility context
- Git, because GLFW and Dear ImGui are fetched by CMake

Default CUDA target:

```text
sm_120 native cubin
compute_120 PTX
```

For another CUDA GPU, override `CMAKE_CUDA_ARCHITECTURES`.

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
cmake -S . -B build -DCOV_ENABLE_CUDA=OFF -DCOV_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Molden parser rules

The parser intentionally fails hard on inconsistent data. After expanding shells it requires:

```text
n_basis_from_shells == n_coefficients_per_MO
```

An MO coefficient index larger than the derived basis count is treated as a shell-convention/parser error, not guessed or silently repaired.

## CUDA evaluator

For each grid point, the CUDA kernel evaluates

\[
\psi_i(\mathbf r)=\sum_\mu C_{\mu i}\chi_\mu(\mathbf r)
\]

directly. It deliberately does **not** allocate an `Ngrid × Nbasis` AO matrix.

The next performance milestone is a fused multi-MO kernel so expensive Gaussian/AO work can be reused across multiple orbitals.

## First regression target: cyclopentadienyl anion

Cp⁻ (`C5H5−`) at `wB97XD/aug-cc-pVTZ` remains the first serious scientific regression target:

- spherical `5D/7F`: 345 basis functions
- Cartesian `6D/10F`: 400 basis functions
- 18 occupied MOs
- Cartesian MO16 ≈ −0.236093 Ha, `a2''` π
- Cartesian MO17 ≈ MO18 ≈ −0.068048 Ha, `e1''` π

The numerical regression gate will compare CUDA values against a trusted Gaussian/Multiwfn grid at identical points and report maximum/RMS error.

## Repository layout

```text
include/cov/              public C++ interfaces
src/parser/               Molden parser
src/orbitals/             orbital semantics + diagram generation
src/cuda/                 CUDA orbital evaluator
src/render/               volume and molecule rendering
src/ui/                   UI, localisation and orbital browser
src/platform/             native platform integration
tests/                    CPU/core regression tests
examples/                 tiny checked-in examples only
.github/workflows/        CPU/core CI
```

## Licence

Apache License 2.0.
