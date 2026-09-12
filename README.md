# CUDA Orbital Visualisation

**Wavefunction in, GPU-resident orbital visualisation out — no CPU volumetric bottleneck.**

CUDA Orbital Visualisation is an experimental GPU-first molecular-orbital viewer for NVIDIA CUDA GPUs, developed first for **RTX 5090 / Blackwell (`sm_120`)**. The initial target is small and medium molecular systems of **up to 100 atoms**.

The viewer parses Gaussian FCHK/FCH and Molden wavefunctions on the CPU, uploads basis data to CUDA, evaluates a selected molecular orbital on a 3D grid, writes the scalar field **directly into an OpenGL 3D texture through CUDA/OpenGL interop**, and ray-marches the positive and negative isosurfaces on the GPU.

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
- [x] Degeneracy-aware labels such as `17-a`, `17-b` while retaining raw MO numbers
- [x] Ha / eV / J mol⁻¹ / kJ mol⁻¹ / cal mol⁻¹ / kcal mol⁻¹ display
- [x] Enhanced ball-and-stick defaults plus stick/delocalisation display mode
- [x] Interactive energy-level diagram with electron occupancy
- [x] **Valence-focused central MO diagram** export: PNG + SVG + JSON + CSV
- [x] Non-destructive valence/core/high-virtual filtering for readable diagrams
- [x] Machine-readable annotation provenance for orbital-family/bonding labels
- [ ] Universal strict SALC construction
- [x] Gaussian FCHK/FCH parsing, source identities and independent numerical regression
- [x] Local CHK → `formchk` integration (requires an installed converter; configure `COV_FORMCHK` when needed)
- [x] Gaussian `.log/.out` companion enrichment with explicit source and missing/invalid metadata states
- [ ] Multi-MO fused CUDA evaluator
- [ ] GPU Marching Cubes for mesh export
- [ ] Cube export
- [x] All-MO samples and complete frontier textures compared with independent IOData/GBasis references under frozen input/grid identities

> **Scientific-status note:** this remains an experimental viewer. The [validation records](validation/README.md) distinguish implemented features, completed numerical/UI subsets and pending physical or full product acceptance. Local symmetry explanations, ligand-field heuristics and topology models do not by themselves establish a physical electronic state or chemical mechanism.

## Architecture

```text
FCHK/Molden parser (CPU)
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

```text
Open calculation
→ inspect wavefunction
→ browse orbitals
→ select MO
→ inspect energy / occupation / spin / symmetry
→ view orbital
→ generate/export a valence MO diagram
```

The orbital browser keeps one-based source MO numbers for the parsed FCHK/FCH or Molden orbital sequence, together with each member's spin and zero-based internal index. With the default `1e-5 Ha` tolerance, a compatible two-level set may be presented as `17-a` and `17-b`. This grouping is a presentation rule; exact member energies and identities remain available in orbital details and exported metadata.

High virtual orbitals are never deleted. The default `Auto · reasonable` view and the diagram's valence selection are non-destructive filters. `All` restores the complete parsed orbital list. `Core`, `Valence`, and high-virtual hiding are human-facing selection heuristics rather than new quantum-chemical assignments.

## Molecular display

The default skeleton uses larger atoms and thicker bonds while retaining orbital surfaces as the visual subject. A second `Stick + delocalisation` mode reduces atom clutter. The renderer distinguishes ordinary bonds, coordination, multicentre and cage support, hydrogen bonds, noncovalent contacts and ionic interactions according to the available interaction graph. A separate geometric ring heuristic can add delocalisation strokes. These styles show the recorded model and evidence; a dashed line is not by itself a calculated bond order or a physical assignment.

The renderer uses conventional chemical element colours for common atoms (H white, C graphite, N blue, O red, common halogens green where defined, P orange, S yellow) and will continue to expand its CPK/Jmol-style palette.

## Valence-focused central MO diagram

The automatic diagram combines the parsed FCHK/FCH or Molden wavefunction with explicitly identified producer metadata and derived local models. It displays the central molecular orbitals; it does **not** claim a universal reconstruction of fragment AOs or strict SALCs.

The default diagram:

- shows only the central molecule's molecular orbitals;
- uses a vertical energy axis, low energy at the bottom and high energy at the top;
- separates crowded levels horizontally while preserving the declared energy transform;
- distinguishes exact member energies from a grouped row's representative mean and energy spread;
- shows electron occupancy directly on each level;
- distinguishes producer, derived local, candidate and unavailable symmetry information;
- focuses on occupied valence orbitals plus a bounded frontier-virtual window;
- hides deep core and very high virtual orbitals non-destructively;
- keeps the complete raw orbital set in JSON/CSV metadata.

### Valence selection

`ValenceCentral` selection keeps occupied valence orbitals and a bounded set of frontier virtual orbitals. A diagram capacity is derived from the UI window control, and compatible groups are preserved at selection boundaries. Compact mode can instead select a supported delocalised-pi family or multicentre active space; its recorded mode describes the diagram actually built. Selecting an MO for three-dimensional inspection does not insert it into or reorder the compact diagram.

The exported selection summary reports how many orbitals were shown, how many remained hidden, and how many occupied-valence/frontier-virtual levels were included.

### Orbital type / bonding annotations

The program distinguishes **direct data** from **derived presentation**.

Source orbital identity, energy and spin remain separate from chemical interpretation. Occupation and symmetry fields retain their producer, derived or unavailable status. Energy-tolerance grouping and compact selection are constructed presentation rules, not producer observations.

Orbital family labels such as `σ`, `π`, `δ` or `φ` can come from supported producer text or the computed orbital-chemistry decomposition. Derived labels retain their source, participating atoms or subspace, and support information. Familiar molecule names, energies and occupations alone do not establish these assignments.

`bonding / nonbonding / antibonding` can likewise come from producer annotations or an available computed decomposition. The selected MO's pair contribution is based on overlap population; a displayed Mayer index describes the total density at the atom-pair level. They are different quantities. Missing or unresolved evidence remains unavailable or unclassified, and model support scores are not probabilities of physical correctness.

Every JSON/CSV orbital row therefore includes provenance fields such as `family_source`, `bonding_class_source`, confidence, and a heuristic flag.

### Export

`Export diagram + metadata` writes:

```text
calculation.mo.png
calculation.mo.svg
calculation.mo.json
calculation.mo.csv
```

PNG/SVG use a light report-friendly central energy layout. JSON/CSV retain raw MO number, internal index, grouped label, Hartree energy, selected display-unit energy, occupation, spin, symmetry, degeneracy size, region/filter state, whether the orbital was included in the diagram, optional orbital-family/bonding annotation, annotation source/confidence, visibility and selection state.

The project does **not** claim universal AO→MO textbook reconstruction or universal strict SALC generation from ordinary Molden files.

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

## Historical regression fixture: cyclopentadienyl anion

The original Cp⁻ (`C5H5−`) fixture at `wB97XD/aug-cc-pVTZ` motivated the initial numerical checks. Its historical reference configuration was:

- spherical `5D/7F`: 345 basis functions
- Cartesian `6D/10F`: 400 basis functions
- 18 occupied MOs
- Cartesian MO16 ≈ −0.236093 Ha, `a2''` π
- Cartesian MO17 ≈ MO18 ≈ −0.068048 Ha, `e1''` π

The current numerical gate compares actual CUDA values with independently evaluated references at identical recorded points, reports phase-aware shape and amplitude errors, and preserves failures. See the dated [validation records](validation/README.md) for the tested inputs, program identities, tolerances and remaining acceptance scope.

## Repository layout

```text
include/cov/              public C++ interfaces
src/parser/               FCHK/Molden parsers, CHK conversion and Gaussian companion metadata
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
