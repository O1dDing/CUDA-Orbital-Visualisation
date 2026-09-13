# COMMON-ROOT-008E software validation

The measured ordinary viewer and native-validation programs are built from `3752f8afe7b1e5687f9960fb7f15cd7c591a3e79`. The revision corrects unreachable orbital details after a main-window resize and narrow vertical wrapping of colored chemistry text. The [pre.11 release notes](../../docs/releases/v0.3.0-pre.11.md) describe the public distribution.

| Scope | Verified result | Record |
|---|---|---|
| Dedicated controls | 39 component tests in each build, eight mathematical probes, viewport checks, 20 language/availability/spin frames and existing g-shell fixtures | [Complete controls](complete-local-controls.json) |
| Fixed 40-case precheck | 39,123 native actions, zero failed actions; numerical, view, scope and selected visual checks passed | [Fixed 40](fixed40-adjudication.json) |
| Existing pure/Cartesian g inputs | 2,038 actions, 500 sampled MOs and 12 complete frontier textures | [g-shell controls](g-controls-adjudication.json) |
| New 273-case software round | 187,792 native actions, 31,871 sampled MOs, 1,214 complete frontier textures and 2,113 export bundles; recorded checks passed | [Full 273](full273-adjudication.json) |
| Representative images | Direct inspection of the specified current-round and fixed-precheck frames | [Full-round images](full273-visual-review.json), [precheck images](fixed40-visual-review.json) |
| Ordinary viewer | 13 recorded frames covering manual load, member selection, details movement/resize/scroll/close/reopen and English/Chinese text | [Interaction](ordinary-interaction-review.json) |
| Distribution directory | Startup, manual wavefunction load, rendered scene and normal exit | [Directory check](ordinary-delivery-review.json) |

All 273 cases remain in the corpus. Independent reference values were reused only under matching mathematical identities; current COV outputs were collected and compared again. Cached pass verdicts were not used.

## Corrected behavior

- The details title and close control remain reachable when the main window shrinks. Valid manually chosen positions remain stable.
- Colored sidebar segments wrap onto a usable line before their remaining width becomes a narrow vertical column.
- Whole-orbital and local symmetry, selected-member and group values, missing/inapplicable/measured-zero states, typed energy gaps, ligand priors and skeleton-ring witnesses retain their declared scopes.

## Limits and retained observations

The first ordinary startup showed a blank client area. The same executable subsequently started normally four consecutive times, including the full interaction sequence and the distribution-directory check. The cause remains unknown; no source fix or definitive attribution is claimed.

Representative visual checks do not certify every pixel, glyph or environment. The 11 small-viewport checks address scene placement and shape; the 640×360 scrolling sidebar does not prove that all browser content is visible simultaneously. See the [viewport record](viewport-visual-review.json).

The numerical and local-model results do not establish physical electronic states, geometries, symmetry assignments or chemical mechanisms. Complete scientific acceptance remains zero. Current physical-reference progress is reported separately in the [REF-001 snapshot](../ref001-progress-20260912/README.md), including the accepted and unaccepted runtime capabilities.

## Build identity and history

The D precheck retained an executable configured before its source commit, so all 40 cases were invalidated by the embedded-commit gate despite successful native actions. Its [failure record](prior-D-invalid-identity.json) is preserved. A clean source commit, reconfiguration, rebuild and actual identity preflight produced the E validation executable. The ordinary and audit executable bytes were unchanged.

The earlier [C-round result](prior-C-full273-adjudication.json) belongs to its own program identity. Later ordinary-viewer evidence exposed the resize defect; it is not used as this revision's pass verdict.

[Program hashes](programs.json), [build identity](build-identity-preflight.json), [workflow gate](workflow-gate.json), [frozen plan](candidate-plan-frozen.json) and [artifact index](artifact-index.json) preserve the exact records. Frozen JSON bytes are not normalized by Git. The [local package receipt](delivery-receipt.json) records the earlier local distribution; the public prerelease includes separate distribution and checksum manifests.

Full raw inputs, textures and action traces retain their recorded source locations. The public evidence bundle contains detailed numerical/view/scope reports and the representative captured images; it is not a duplicate of the entire raw acquisition directory.
