# R16 NBO integration acceptance

The acceptance documents record the scoped final R16 replay. `source-capture.json` binds the 50 previously uncommitted files to the accepted 3,435-file source snapshot. The release commit adds publication documentation and CMake project-version metadata and the explicit parser dependency required for portable static-library linking.

The optional release assets retain the exact build-input source, raw replay evidence, native-validation executable and reference inputs. Local paths and process IDs inside raw evidence describe the machine on which the tests ran; binaries and inputs are identified by hashes. Failed original runs remain distinguishable from corrected passing runs. See [release notes](../../docs/releases/v0.4.0-pre.1.md).
