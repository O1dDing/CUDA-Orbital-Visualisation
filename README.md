# CUDA Orbital Visualisation

**Wavefunction in, GPU-resident orbital visualisation out — no CPU volumetric bottleneck.**

CUDA Orbital Visualisation is an experimental GPU-first molecular-orbital viewer for NVIDIA CUDA GPUs, developed first for **RTX 5090 / Blackwell (`sm_120`)**. The initial target is small and medium molecular systems of **up to 100 atoms**.

The current MVP parses Molden wavefunctions on the CPU, uploads basis data to CUDA, evaluates a selected molecular orbital on a 3D grid, writes the scalar field **directly into an OpenGL 3D texture through CUDA/OpenGL interop**, and ray-marches the positive and negative isosurfaces on the GPU.

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
- [x] Basic atom/bond overlay
- [x] 64³ / 128³ / 256³ / 512³ grids
- [x] Interactive MO, isovalue, camera and grid controls
- [x] CPU-only parser CI
- [ ] Gaussian FCHK
- [ ] CHK → `formchk`
- [ ] Gaussian `.log/.out`
- [ ] Multi-MO fused CUDA evaluator
- [ ] GPU Marching Cubes for mesh export
- [ ] Cube export
- [ ] Cp⁻ numerical regression against Gaussian/Multiwfn reference grids

> **Scientific-status note:** this is an early MVP, not yet a validated production quantum-chemistry package. The spherical-harmonic path and producer-specific Molden phase/order conventions must pass reference-grid regression tests before scientific release.

## Why this architecture?

Conventional orbital workflows often become:

```text
wavefunction
    ↓
CPU grid generation
    ↓
CPU mesh extraction
    ↓
GPU display
```

or, in some GPU-enabled viewers:

```text
GPU scalar field
    ↓
GPU → CPU readback
    ↓
CPU/JS Marching Cubes
    ↓
GPU display
```

This project instead aims for:

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

## Requirements

### GPU build

- NVIDIA GPU with CUDA support
- **CUDA Toolkit 12.8+**
- CMake 3.28+
- C++20 compiler
- OpenGL 2.1+ compatibility context
- Git, because GLFW and Dear ImGui are fetched by CMake

The default CUDA target is RTX 5090 / Blackwell:

```text
sm_120 native cubin
compute_120 PTX
```

For a different CUDA GPU, override `CMAKE_CUDA_ARCHITECTURES`.

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

### Parser-only build without CUDA

Useful for CI or parser development:

```bash
cmake -S . -B build -DCOV_ENABLE_CUDA=OFF -DCOV_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Usage

Start with a Molden file:

```bash
cov calculation.molden
```

or launch the application and either:

- enter a `.molden` path in the UI; or
- drag a Molden file onto the window.

Controls:

- **Left mouse drag:** orbit
- **Mouse wheel:** zoom
- **MO slider:** select orbital
- **Isovalue:** changes the GPU raymarch threshold immediately
- **Grid:** 64³, 128³, 256³ or 512³; changing this recomputes the CUDA scalar field

Positive and negative orbital phases are rendered simultaneously with different colours.

## Molden parser rules

The parser intentionally fails hard on inconsistent data.

After expanding the parsed shells, it requires:

```text
n_basis_from_shells == n_coefficients_per_MO
```

An MO coefficient index larger than the derived basis count is treated as a shell-convention/parser error, not guessed or silently repaired.

For large Molden files the parser uses two streaming passes rather than `read()`-ing the whole file into one giant string. The first pass parses geometry/basis metadata; the second parses MO blocks.

## CUDA evaluator

For each grid point, the CUDA kernel evaluates

\[
\psi_i(\mathbf r)=\sum_\mu C_{\mu i}\chi_\mu(\mathbf r)
\]

directly. It deliberately does **not** allocate an `Ngrid × Nbasis` AO matrix.

The current MVP evaluates one selected MO per kernel pass. The next performance milestone is a fused multi-MO kernel so expensive Gaussian/AO work is reused across 4/8/16 orbitals.

## First regression target: cyclopentadienyl anion

The first serious validation case is Cp⁻ (`C5H5−`) at `wB97XD/aug-cc-pVTZ`.

Reference expectations from the development handoff:

- spherical `5D/7F`: 345 basis functions
- Cartesian `6D/10F`: 400 basis functions
- 18 occupied MOs
- Cartesian MO16 ≈ −0.236093 Ha, `a2''` π
- Cartesian MO17 ≈ MO18 ≈ −0.068048 Ha, `e1''` π

The regression gate will compare CUDA values against a trusted Gaussian/Multiwfn grid at identical points and report maximum/RMS error.

## Repository layout

```text
include/cov/              public C++ interfaces
src/parser/               Molden parser
src/cuda/                 CUDA orbital evaluator
src/render/               OpenGL loader + volume raymarch renderer
src/main.cpp              GLFW + Dear ImGui application
tests/                    CPU parser tests
examples/                 tiny checked-in examples only
.github/workflows/        CPU parser CI
```

## Near-term roadmap

1. Validate Cartesian s/p/d/f/g against reference cubes.
2. Validate spherical 5D/7F/9G ordering, phase and normalisation.
3. Add Cp⁻ regression data and GPU-vs-reference error report.
4. Replace generic spherical-harmonic math with hard-coded low-`l` polynomial fast paths.
5. Add fused multi-MO evaluation.
6. Add FCHK.
7. Add CHK → `formchk`.
8. Add GPU Marching Cubes only for mesh export.
9. Add Cube export.

## Licence

Apache License 2.0.
