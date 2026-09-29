# Chemical Orbital Visualiser (COV)

**English** · [简体中文](README.zh-CN.md) · [日本語](README.ja.md) · [Français](README.fr.md)

Explore orbital energies, occupations and the connections between atomic and molecular orbitals.

COV displays orbital energies and occupations in interactive energy-level diagrams, with figures and data available for export. It reads Gaussian FCHK/FCH and Molden files.

The v0.4 preview adds NBO analysis. With the calculation's NBO output, you can see how atomic and localised orbitals contribute to molecular orbitals, follow their connections in the diagram, and inspect charge and bonding information.

Select an orbital to view its shape in 3D. Rendering uses NVIDIA CUDA, and the interface is available in English, Simplified Chinese, Japanese and French.

## What you can do

- **Read energy-level diagrams.** See orbital energies and electron occupations, change energy units, and export a diagram as PNG or SVG.
- **Follow orbital connections — v0.4 preview.** See how atomic and localised orbitals contribute to molecular orbitals. Select a connection to inspect its contribution.
- **Inspect charge and bonding — v0.4 preview.** Read NPA charges, spin populations, Wiberg bond indices and donor–acceptor interactions when the NBO files contain them.
- **Export figures and data.** Save energy diagrams as PNG or SVG, and the corresponding data as CSV or JSON.
- **Browse orbitals.** Jump to HOMO or LUMO, search the orbital list, and show occupied, virtual, core or valence orbitals.
- **View orbitals in 3D.** Select an orbital, adjust its isosurface, and rotate or zoom the molecular view.

## Download

| Version | Includes | Windows download |
|---|---|---|
| [Stable v0.3.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0) | Orbital energies, occupations, energy-level diagrams and 3D views | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.3.0/CUDA-Orbital-Visualisation-v0.3.0-Windows-sm120.zip) |
| [Preview v0.4.0-pre.1](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0-pre.1) | Adds NBO analysis, orbital composition and connections | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.1/CUDA-Orbital-Visualisation-v0.4.0-pre.1-Windows-sm120.zip) |

Existing downloads still show the earlier product name.

The current Windows downloads are built for NVIDIA GeForce RTX 50-series GPUs. Other NVIDIA architectures need a suitable build; see [Building from source](docs/BUILD.md). Install a compatible NVIDIA driver. The CUDA Toolkit is only needed when building from source.

## Get started

1. Download and extract a Windows ZIP, then run `cov.exe`.
2. Open a Gaussian FCHK/FCH or compatible Molden file, or drag it into the window.
3. Browse the energies and occupations in the orbital list and energy diagram. Select a level to view its orbital, or export the diagram.
4. In the v0.4 preview, open the calculation folder to load its wavefunction and NBO files together. If the folder contains several calculations, choose the one to open.

Convert a Gaussian `.chk` file to `.fchk` with Gaussian's `formchk` before opening it. For NBO features, use files for the same geometry and electronic state as the wavefunction. COV reads the results; it does not run Gaussian or NBO.

## Input and requirements

- **Wavefunctions:** Gaussian `.fchk` / `.fch`, or compatible `.molden` / `.mol` / `.input` files. Gaussian `.chk` needs conversion first.
- **NBO files — v0.4 preview:** the report, archive and orbital matrices provide different parts of the analysis. See [Using NBO results](docs/AOMO_NBO.md) for the files each view needs.
- **Molecule size:** up to 100 atoms per input.
- **Graphics:** an NVIDIA GPU, a compatible driver and OpenGL 2.1 or newer. Use a build made for your GPU architecture.

COV reads existing calculation results. Automatic orbital and bonding labels can be wrong; use the original output when interpreting them.

## Documentation

- [Using COV](docs/UI.md)
- [Using NBO results](docs/AOMO_NBO.md) — v0.4 preview
- [Building from source](docs/BUILD.md)
- [Stable release notes](docs/releases/v0.3.0.md) · [Preview release notes](docs/releases/v0.4.0-pre.1.md)

## Licence

[Apache License 2.0](LICENSE). Bundled library licences are listed in [Third-party notices](THIRD_PARTY_NOTICES.md).
