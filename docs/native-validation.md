# Native single-case validation

The `test/fchk-validation-native` branch starts at PR #3 commit
`bdc3ece61b3847a782873bfd2654521c9808efb5`. It adds observation and a local
test driver. It does not fix the scientific failures retained in the original
validation ledger. Normal builds use `COV_ENABLE_VALIDATION=OFF` and produce
`cov.exe`; ON builds produce `cov_validation.exe`.

## Build and run

Use the existing CUDA 12.8 / Visual Studio 2022 / ImGui 1.90.9 / OpenGL 2.1
configuration. Configure `COV_ENABLE_CUDA=ON`, `COV_BUILD_TESTS=ON`, and
`COV_ENABLE_VALIDATION=ON`. An independent OFF build is used by the checker.
No ImGui upgrade, server, or ImGui Test Engine dependency is introduced.

```
cov_validation.exe input.fch --validation-plan case.plan --validation-output new-directory
```

The first plan line is `COV_VALIDATION 1`. Supported commands are:

| Command | Meaning |
| --- | --- |
| `scene opacity yaw pitch distance iso resolution` | Set scene parameters only |
| `seek "semantic.id"` | Reach a drawn target through actual wheel input |
| `click "semantic.id"` | Send mouse position, down and up to the observed hit rectangle |
| `hover "semantic.id"` | Move over a drawn target |
| `text "browser.search" "query"` | Click, Ctrl+A, type, Enter through ImGui input |
| `key "Home/Down/Up/Enter/Escape"` | Keyboard navigation, including dropdowns |
| `capture "name"` | Capture the completed visible viewer + ImGui back buffer |
| `volume "name" "zero-based-MO"` | Read the texture just consumed by the renderer |
| `wait "label"` | Wait four drawn frames |

The plan owns input while active. GLFW cursor events are cleared and the
ordered native event batch is processed in the same frame, with event
trickling disabled. The real production controls remain the sole writers of
selection/filter/expansion/export state. A missing or clipped target fails;
there is no direct-selection fallback. An `executed` action records that input
was delivered; the external checker separately verifies its outcome.

The scene is drawn before deferred selection/evaluation. Frames therefore
record requested, UI-drawn, rendered and subsequently applied MO indices,
volume and diagram generations, and kernel/recompute information. A new
applied selection after drawing is a normal transition, not automatically a
scene/UI mismatch. Readback occurs before the later CUDA update, after the
previous CUDA resource unmap. It uses the actual renderer texture, not a
separate evaluator. Sample indices and grid interpolation conventions are
saved. Captures use the final GL back buffer after ImGui rendering.

`tests/native_validation_case.py` runs one supplied molecule with explicit
input/log/build/OFF-build/audit-helper/reference-dependency/output arguments.
It imports the retained independent audit helpers without modifying their
ledgers. Pillow is installed in the build's `python-deps` directory only.
It checks direct coefficients/energies, independent overlap/density/Mayer,
symmetry operations, all-MO sampled actual textures, individually selectable
compact members, expansion/restoration, units/languages and actual exports.
The output includes a review page, raw evidence and phase-specific timings.
No file in an existing output directory is overwritten.

## Evidence and limits

- All-MO means every MO is selected in the actual browser, then sampled at
  up to 8192 deterministic texture indices (duplicates removed). It does not
  mean every voxel or an integral over all space.
- A reference-only odd-|m| phase experiment can explain the known failure;
  it never changes a production failure into a pass.
- The independent checker owns the verdict. Native collector success is
  distinct from scientific success. The checker returns failure for known
  scientific failures instead of counting expected failures as passing cases.
- Scientific scope is the supplied wavefunction and calculation. A fixed
  geometry single point is not evidence of an optimized stable ground state.
- Image review is a separate measured stage; machine-only timings exclude
  development, compilation, and subsequent review/reporting.
- ON/OFF science data and original regressions are compared. This first
  implementation is not a complete certification of non-interference across
  every UI timing, resolution, degenerate group and molecule.
- Full tooltip clipping, all four-language terminology, complete general
  active-space selection, all point groups and all proposed inspection APIs
  still require additional coverage. Unsupported/insufficient checks stay
  explicit. There is no automatic expansion to the 273-case queue or new
  Gaussian calculations from this single-case command.

References: the installed Gaussian `doc/formchk.txt`,
[GBasis evaluation documentation](https://gbasis.qcdevs.org/tutorial/Evaluations_basis_and_potential.html),
and the [D2h character table](https://www.staff.ncl.ac.uk/j.p.goss/symmetry/D2h.html).
