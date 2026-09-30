# Chemical Orbital Visualiser (COV)

**English** · [简体中文](README.zh-CN.md) · [日本語](README.ja.md) · [Français](README.fr.md)

Explore orbital energies, occupations and the connections between atomic and molecular orbitals.

COV displays orbital energies and occupations in interactive energy-level diagrams, with figures and data available for export. It reads Gaussian FCHK/FCH and Molden files.

The v0.4 preview adds NBO analysis. With the calculation's NBO output, you can see how atomic and localised orbitals contribute to molecular orbitals, follow their connections in the diagram, and inspect charge and bonding information.

Select an orbital to view its shape in 3D. Rendering uses an NVIDIA GPU, and the interface is available in English, Simplified Chinese, Japanese and French.

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

The current Windows packages support RTX 50-series GPUs; packages for other architectures are not yet available.

## Get started

1. Download and extract a Windows ZIP. For v0.3.0, run `cov.exe`. For v0.4.0-pre.1, run `program/cov.exe`, or double-click an `Open-*.cmd` launcher to open an example.
2. Open a Gaussian FCHK/FCH or compatible Molden file, or drag it into the window.
3. Browse the energies and occupations in the orbital list and energy diagram. Select a level to view its orbital, or export the diagram.
4. In the v0.4 preview, open the calculation folder to load its wavefunction and NBO files together. If the folder contains several calculations, choose the one to open.

With `formchk` installed, COV can convert Gaussian CHK files to FCHK. FCHK/FCH and Molden files provide MO energies, occupations and shapes. NBO orbital shapes and composition analysis also need the matching report, `.47` archive and orbital matrices.

## Input and requirements

- **Wavefunctions:** Gaussian `.fchk` / `.fch`, or compatible `.molden` / `.mol` / `.input` files. Gaussian `.chk` can be opened with an installed `formchk`.
- **NBO files — v0.4 preview:** [Prepare calculation files](docs/NBO_ONE_JOB.md) explains how to produce them; [Using NBO results](docs/AOMO_NBO.md) lists the files for each view.
- **Molecule size:** up to 100 atoms per input.
- **Graphics:** an NVIDIA GPU, a compatible driver and OpenGL 2.1 or newer. Use a build made for your GPU architecture.

## Documentation

- [Using COV](docs/UI.md)
- [Using NBO results](docs/AOMO_NBO.md) — v0.4 preview
- [Prepare calculation files](docs/NBO_ONE_JOB.md) — includes the source-tree job template
- [Building from source](docs/BUILD.md)
- [Stable release notes](docs/releases/v0.3.0.md) · [Preview release notes](docs/releases/v0.4.0-pre.1.md) · [Older v0.3 previews](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0-pre-archive)

[Report a problem or suggest a feature](https://github.com/O1dDing/Chemical-Orbital-Visualiser/issues/new/choose).

## Licence

[Apache License 2.0](LICENSE). Bundled library licences are listed in [Third-party notices](THIRD_PARTY_NOTICES.md).
