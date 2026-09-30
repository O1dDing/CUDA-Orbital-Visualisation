# Prepare one calculation for COV

[简体中文](NBO_ONE_JOB.zh-CN.md)

The [one-job template](../examples/nbo-one-job/cov_one_job.py) runs Gaussian/NBO, `formchk` and GenNBO in sequence, then puts the wavefunction and NBO files in a folder that COV can open. GenNBO reanalyses the archived wavefunction without another SCF calculation.

This template runs from the COV source tree on Windows with Python 3.11 or newer, Gaussian 16W A.03 and NBO7 i8. Gaussian and NBO must already be installed, with `gaunbo6.bat` configured for that NBO installation. A.03 uses `Pop=NBO6Read` to call NBO7. Other Gaussian revisions use their own interface setup; this template checks for A.03.

The existing v0.3.0 and v0.4.0-pre.1 download packages do not include this template. Its process supervisor is [tests/validation_process.py](../tests/validation_process.py), included in the source tree.

## Submit a calculation

Copy [water.json](../examples/nbo-one-job/water.json) and set the molecule, geometry, method, basis, charge and multiplicity. The example is a fixed-geometry PBE1PBE / 6-31G(d) water single point, charge 0, multiplicity 1. `basis_ecp_tail` accepts the explicit basis/ECP block when needed.

Run from the source-tree root, using your installation paths and a new output directory:

```powershell
python examples/nbo-one-job/cov_one_job.py `
  --recipe examples/nbo-one-job/water.json `
  --output work/water-job `
  --gaussian-bin 'C:/Gaussian/Gaussian 16 W' `
  --nbo-bin 'C:/NBO/bin'
```

This previews the input and checks the paths without creating files or starting a calculation. Add `--run` to the same command to submit it. The output directory must not already exist. The template keeps the supplied geometry and uses one single-point step; optimization and frequency calculations need their own inputs.

The template permits up to three Gaussian threads. The process tree is limited to three cores and 20 GiB; timeouts are 900 seconds for Gaussian, 300 for `formchk` and 600 for GenNBO. `--memory-limit-gib` and `--gaussian-timeout` adjust the documented limits; `--help` lists the options.

## Open the result

Open `work/water-job/cov-package` in the NBO preview, or drag its `drop.covnbopkg` into COV. The folder contains `canonical.fchk`, one `analysis.nbo` report, the `.47` archive and the requested orbital matrices.

The canonical checkpoint and original Gaussian archive are preserved. The report, archive and matrices in `cov-package` come from the same GenNBO analysis. Run records stay outside this input folder. Failed jobs retain their output; use a new directory for another attempt.

## Files for each view

The FCHK and NBO data must describe the same calculation step. The archive supplies overlap, density, canonical MO and Fock data; the NBO report supplies the orbital identities and printed values. Matrix numbers below are the template's settings; COV reads their headings.

| View | Matrix data written by the template |
| --- | --- |
| MO energies, occupations, Gaussian AO–MO coefficient links and 3D MO | `canonical.fchk` |
| NAO shapes, NPA charge and Wiberg colouring | AONAO (W52), report and matching density |
| MO composition in NAOs | AONAO (W52), NAOMO (W51); NAONBO (W53) provides another transform |
| NBO shapes and MO–NBO links | AONBO (W37), NBOMO (W49) |
| NHO shapes and components | AONHO (W54), NAONHO (W57), AONAO (W52); NHONBO (W59) links to NBOs |
| NLMO shapes, main NBO and tail | AONLMO (W55), NBONLMO (W60), AONBO (W37); NLMOMO (W61) links to MOs |
| PNAO shapes | AOPNAO (W56), AONAO (W52) |
| E(2) donor–acceptor orbitals | AONBO (W37), report and archive Fock data |
| NAO/fragment energies and SALC data | AONAO (W52), matching density and Fock data; available symmetry depends on the geometry and selected fragment |

The template also writes SAO, NAONLMO and AOMO. Missing data affect their corresponding views. Open-shell spin displays need the actual spin-resolved data.

Matrix and archive options are described in the [NBO 7 manual](https://nbo.chem.wisc.edu/nboman.pdf). COV's views and exports are described in [Using NBO results](AOMO_NBO.md).
