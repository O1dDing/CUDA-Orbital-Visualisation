# Native single-case checks

The native driver runs a text plan through COV's ordinary UI controls and records
the displayed scene, selections, exports and GPU texture samples. Historical
run reports are retained in the private development archive. Current runtime
and recovery instructions are in [validation/README.md](../validation/README.md).
Normal builds use
`COV_ENABLE_VALIDATION=OFF` and produce `cov.exe`; ON builds produce
`cov_validation.exe`.

## Build and run

Use the existing CUDA 12.8 / Visual Studio 2022 / ImGui 1.90.9 / OpenGL 2.1
configuration. Configure `COV_ENABLE_CUDA=ON`, `COV_BUILD_TESTS=ON`, and
`COV_ENABLE_VALIDATION=ON`. An independent OFF build is used by the checker.

CMake embeds the Git commit in the executable and adds `-dirty` when tracked
files have local changes. Rebuild after changing the source. The run's
`identity.json` records the executable's source version.

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
The output includes a review page, captured data and timings for each stage.
No file in an existing output directory is overwritten.

## Reading the output

- All-MO means every MO is selected in the actual browser, then sampled at
  up to 8192 distinct texture positions.
- The collector records input, actions and renderer output. The separate
  checker compares those results with the supplied reference data.
- Image review has its own timing. Program timings exclude builds and review.
- The single-case command uses the supplied files and does not start Gaussian.

References: the installed Gaussian `doc/formchk.txt`,
[GBasis evaluation documentation](https://gbasis.qcdevs.org/tutorial/Evaluations_basis_and_potential.html),
and the [D2h character table](https://www.staff.ncl.ac.uk/j.p.goss/symmetry/D2h.html).

## Batch collection

`tests/native_collection_batch.py` recursively inventories `.fch` and `.fchk`,
then uses a thread pool to launch independent COV processes. Its raw default is
four workers; set the worker and resource limits for the current run explicitly.
Every case has an isolated
current directory, input/log snapshot,
production dump, native plan, exports, captures, action trace and process IDs.
The `--validation-background` flag creates a hidden window/context and
still draws the real production scene and ImGui into the same back buffer.
Framebuffer size is required to remain 2100x1250; captures are losslessly
converted and checked. See the [GLFW offscreen-context documentation](https://www.glfw.org/docs/latest/context_guide.html#context_offscreen).

Collection includes every alpha/beta MO texture and the current compact
members, including spin counterparts reached through the actual browser.
The numerical and image checks run separately from collection.
`progress.json` records progress, failures retain their output, and
`COMPLETED.json` is written only after all cases return and the index is saved.
The final summary distinguishes complete collection, collection with gaps,
and runtime errors.
