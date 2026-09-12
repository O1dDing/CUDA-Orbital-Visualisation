# Native single-case validation

The `test/fchk-validation-native` branch originally added observation and a local
test driver at PR #3 commit `bdc3ece61b3847a782873bfd2654521c9808efb5`.
It now also contains the subsequent parser, numerical, scope, navigation,
diagram and UI fixes recorded in the [dated validation history](../validation/README.md).
Each pass belongs to its recorded build and scope; it does not certify a physical
electronic state or all product behavior. Normal builds use
`COV_ENABLE_VALIDATION=OFF` and produce `cov.exe`; ON builds produce
`cov_validation.exe`.

## Build and run

Use the existing CUDA 12.8 / Visual Studio 2022 / ImGui 1.90.9 / OpenGL 2.1
configuration. Configure `COV_ENABLE_CUDA=ON`, `COV_BUILD_TESTS=ON`, and
`COV_ENABLE_VALIDATION=ON`. An independent OFF build is used by the checker.
No ImGui upgrade, server, or ImGui Test Engine dependency is introduced.

Commit the intended source, then reconfigure and rebuild the validation target.
CMake embeds the current Git commit and appends `-dirty` for modified tracked
files. Committing after a build does not update that binary. Before freezing a
round, run a short native plan and verify that its `identity.json.git_commit`
matches the intended source commit and that the executable SHA-256 matches the
program inventory. Preserve the source, inputs, plans, criteria and build receipt.
A mismatch invalidates the round even if its actions succeed; rebuild and create
a new round instead of editing the old manifest or identity output.

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
| `window width height` | Resize the actual GLFW window; defaults remain 2100x1250 |
| `drag "semantic.id" dx dy` | Move to the observed target, press, move by logical pixels, then release |
| `wheel "semantic.id" amount` | Move to the observed target and send a vertical wheel event |

Window controls accept integral client sizes from 640x360 to 7680x4320 and
verify the actual framebuffer against the request. Pointer controls feed
normal ImGui events; they never write the camera or selection directly.
Their delivered events and resulting camera/GL viewport are recorded so the
external reviewer can check scene/panel separation and actual resize behavior.

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
- Tooltip and details-window defects are adjudicated individually in the dated
  validation records, including actual resize, scroll, close and language checks.
  Those checks do not cover every possible term, point group, active space or
  display environment. Unsupported or insufficient checks remain explicit.
  This single-case command does not automatically start the 273-case queue or
  new Gaussian calculations.

References: the installed Gaussian `doc/formchk.txt`,
[GBasis evaluation documentation](https://gbasis.qcdevs.org/tutorial/Evaluations_basis_and_potential.html),
and the [D2h character table](https://www.staff.ncl.ac.uk/j.p.goss/symmetry/D2h.html).

## Corpus collection without scientific analysis

`tests/native_collection_batch.py` recursively inventories `.fch` and `.fchk`,
then uses a thread pool to launch independent COV processes. Its raw default is
four workers; a frozen campaign's resource policy overrides that default. The
current COV campaign uses two workers, each limited to two physical cores and
24 GiB, inside a four-core/64-GiB parent limit. Every case has an isolated
current directory, input/log snapshot,
production dump, native plan, exports, captures, action trace and process IDs.
The new `--validation-background` flag creates a hidden window/context and
still draws the real production scene and ImGui into the same back buffer.
Framebuffer size is required to remain 2100x1250; captures are losslessly
converted and checked. See the [GLFW offscreen-context documentation](https://www.glfw.org/docs/latest/context_guide.html#context_offscreen).

Collection includes every alpha/beta MO texture and the current compact
members, including spin counterparts reached through the actual browser.
Only collection integrity is checked at this stage. The scientific and image
review fields remain deferred, even for a completely collected case.
`progress.json` is updated atomically, failures retain their evidence, and
`COMPLETED.json` is written only after all cases return and the index is saved.
The final summary distinguishes complete collection, collection with gaps,
and runtime errors. Completion does not imply scientific correctness.
