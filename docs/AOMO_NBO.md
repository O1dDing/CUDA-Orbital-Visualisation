# Explore AO–MO and NBO relationships

[简体中文](AOMO_NBO.zh-CN.md)

Chemical Orbital Visualiser (COV) brings orbital energies, occupations and composition into one view. In the v0.4.0 preview, you can open an existing Gaussian FCHK and its corresponding NBO results, then move from an MO, atom or bond to related orbitals and contributions. Inspect a selected orbital in 3D when you want to see its shape.

## Which files do I need?

COV reads results that already exist; it does not run Gaussian or NBO. Drop the files together, drop their directory, or enter a path. If the directory contains several calculations, choose the one you want. The advanced path controls help when the files are stored in different places.

| What you want to see | Files to provide |
| --- | --- |
| MO energies, occupations, energy diagram, Gaussian AO–MO coefficient links and 3D MO | Gaussian FCHK/FCH; Molden supports basic MO views. |
| Values printed in the NBO report | FCHK and an NBO report containing those sections. |
| NAO shapes, NPA charge and Wiberg colouring | Matching FCHK, report, `.47` archive and AONAO; colouring also needs matching density. |
| MO composition in NAOs | Matching FCHK, report, archive, AONAO and NAOMO; NAONBO provides another transform. |
| NBO shapes and canonical MO–NBO links | Matching FCHK, report, archive and AONBO; MO–NBO links also need NBOMO. |
| NHO shapes and components | Base files, AONAO, AONHO and NAONHO; NHO–NBO links also use NHONBO and AONBO. |
| NLMO shapes, main NBO and tail | Base files, AONBO, AONLMO, NBONLMO and matching report entries; NLMO–MO links also use NLMOMO. |
| PNAO shapes | Base files, AONAO and AOPNAO. |
| E(2) donor–acceptor orbitals | Matching report, archive Fock data and AONBO. |

Base files mean the matching FCHK, NBO report and `.47` archive. The FCHK and NBO data must describe the same calculation step. The NBO report and orbital matrices must come from the same analysis; use a GenNBO reanalysis report and its matrices together. The archive needs the overlap, density, canonical MO and Fock data used by the selected view. Spin colouring also needs complete spin-resolved density and report data.

Matrix file numbers are chosen when the NBO calculation is produced; they are not fixed by COV. For example, a calculation requesting `AONBO=W37 NBOMO=W49 NAOMO=W51 AONAO=W52 NAONBO=W53` would write files with those numbers if the producer follows the requests. Check the matrix headings and your calculation output rather than assuming that `FILE.37` or another numbered file has a particular role. This five-matrix example covers MO composition in NAOs and MO–NBO links. For the full file set, see the [one-job template](NBO_ONE_JOB.md).

If you start with a CHK file, use an installed `formchk` to create an FCHK, or let COV call the installed converter. Use that step's canonical checkpoint. A checkpoint whose orbitals were replaced by SaveNBOs or SaveNLMOs is not the canonical MO input. If NBO files are missing, the ordinary MO view still works; an analysis that needs missing data will be unavailable.

## Follow an orbital through the diagram

Choose an MO in the browser or energy diagram. The diagram focuses on valence and nearby unoccupied levels; **All** in the browser shows the full imported orbital list. Click an MO to highlight its available links to Gaussian atomic orbitals (AO), natural atomic orbitals (NAO) and atom groups. Click an AO, NAO or link to inspect a contribution or find related MOs. You can select several terms and show their signed partial sum or overlay them in 3D. Folding groups tidies the diagram without removing their members.

NAOs form an orthogonal representation, so squared NAO coefficients can describe weights in that representation. Squared coefficients of the original, generally nonorthogonal Gaussian AOs are not atomic populations. When NAO or SALC energies are available, they can use the same numerical axis as the MOs. Their values are operator expectation values in the molecular environment; the central values are canonical MO energies. The illustrative side layout arranges orbitals without using their energies. A group of orbitals is a SALC only when the available data supports that symmetry meaning.

A local basis may contain fewer orbitals than the Gaussian AO basis. The diagram then shows the available contributions and the uncovered part, without rescaling the contributions to 100%.

**Overview** shows a compact view. **Research analysis** adds detail, and **Full basis** includes core and Rydberg orbitals. These presets change what is shown.

## Start from the molecule

Click an atom or bond in the 3D view to open its related values and orbital links. With the report and matching matrices available, atoms can be coloured by NPA charge or spin population; bonds can show Wiberg indices, bonding and antibonding orbitals, and coordination or multicentre relationships. A Wiberg index is a continuous value, not an integer bond order.

NHO views can show directional orbital lobes and their angular components. Select an E(2) interaction to see its donor and acceptor orbitals together; the E(2) value is a perturbation estimate, not a bond or reaction energy. NLMO views can separate the full orbital, its main NBO component and the remaining tail when the needed data is available.

A small contribution may be hidden by the current isosurface threshold. Use **Fit component** to adjust the display threshold or **Reveal bonds** to make the molecule easier to see. Keep the same threshold when comparing the sizes of different orbitals.

## Export

**Export whole diagram** saves `.aomo.svg`, `.aomo.png`, `.aomo.json` and `.aomo.csv`. SVG/PNG are static figures; JSON/CSV retain selections, groups and values. Interactive 3D exploration remains in COV.
