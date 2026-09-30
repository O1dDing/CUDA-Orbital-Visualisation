# Build Chemical Orbital Visualiser

[简体中文](BUILD.zh-CN.md)

These instructions are for building COV from source. The ordinary Windows download does not require the CUDA Toolkit, CMake or a C++ compiler. The supplied Windows program is built for RTX 50-series GPUs; a suitable source build is needed for another NVIDIA GPU architecture.

## Requirements

- NVIDIA GPU with CUDA support and a compatible display driver
- CUDA Toolkit 12.8 or newer
- CMake 3.28 or newer
- C++20 compiler
- OpenGL 2.1 compatibility context
- Git, used by CMake to fetch GLFW and Dear ImGui

The default CUDA architecture is `sm_120` with `compute_120` PTX. Set `CMAKE_CUDA_ARCHITECTURES` to the architecture of the GPU you are building for when using another target.

## Choose a version

Build the `v0.3.0` tag for the stable viewer, or `v0.4.0-pre.1` for the NBO preview. The commands below run from the source folder for the version you chose.

## Windows / Visual Studio 2022

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --parallel
.\build\Release\cov.exe .\examples\h2.molden
```

## Linux

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/cov ./examples/h2.molden
```

## Core tests without CUDA

```bash
cmake -S . -B build -DCOV_ENABLE_CUDA=OFF -DCOV_BUILD_TESTS=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Inputs and parser

COV reads Gaussian FCHK/FCH and Molden wavefunctions. Local CHK conversion requires a separately installed `formchk`; set `COV_FORMCHK` to its path if automatic discovery fails. A Gaussian `.log` or `.out` file can add information to the corresponding wavefunction, but does not replace the FCHK/FCH input.

The current input limit is 100 atoms. Cartesian and real-spherical `s/p/d/f/g` basis functions are supported. For a Molden file, the expanded basis count must match the coefficients of each molecular orbital. An out-of-range coefficient index is treated as an input or shell-convention error; the parser does not guess a repair.

## Rendering architecture

The CPU parses the input and sends the basis and selected orbital coefficients to the GPU. CUDA evaluates the orbital value on a three-dimensional grid and writes the result into an OpenGL 3D texture. The renderer draws positive and negative isosurfaces from that texture. Changing the isovalue changes the display without recalculating the orbital grid.

For each grid point, the evaluator computes

\[
\psi_i(\mathbf r)=\sum_\mu C_{\mu i}\chi_\mu(\mathbf r).
\]

It does not allocate a full grid-points-by-basis-functions matrix. Available grid resolutions are 64³, 128³, 256³ and 512³.

For Gaussian/NBO file preparation, see the [one-job template](NBO_ONE_JOB.md).
