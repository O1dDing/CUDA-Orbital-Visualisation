# Windows CI portability

The first CI run of d2a9b7d passed all Linux runtime tests. Windows ran 130 tests and reported one admission error: the existing parent/child pause test requested two cores, while the hosted runner exposes two physical cores and four SMT logical processors. The new conservative small-host profile correctly reserves one physical core for the system, leaving a one-core managed budget.

The test now requests one core. Its assertions still verify that both owned processes stop updating during RAM pause, the foreign process keeps updating, both owned processes resume and finish, and the entire owned tree is empty on return. The test exercises lifecycle isolation rather than parallel throughput. The 16/8 host profiles and production runtime files are unchanged, so the d9f48dc70e1e native Gaussian capability identity remains valid.

The raw failure log and the targeted passing local rerun are retained under the local delivery output as `ci-windows-d2a9b7d.log` and `test-ci-small-host`. Remote CI runs the complete suite again on the corrected source. No Gaussian or formal reference calculation was restarted for this test-only correction.
