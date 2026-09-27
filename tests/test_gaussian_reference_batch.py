from pathlib import Path
import tempfile
import unittest

from gaussian_reference_batch import identity_gate, needs_relaxation, stage_diagnostics


class ReferenceChainTests(unittest.TestCase):
    def test_initial_instability_requires_reoptimization_even_if_finally_stable(self):
        with tempfile.TemporaryDirectory() as raw:
            log = Path(raw) / "job.log"
            log.write_text("The wavefunction has an RHF -> UHF instability.\n"
                           "The wavefunction is stable under the perturbations considered.\n"
                           "Normal termination of Gaussian\n")
            diagnostics = stage_diagnostics(log)
            before = {"fields": {"method": "RPBE1PBE", "Total Energy": -1.0}}
            after = {"fields": dict(before["fields"]), "diagnostics": diagnostics}
            self.assertTrue(diagnostics["stability_reported"])
            self.assertTrue(needs_relaxation(before, after))

    def test_same_stable_state_does_not_start_an_extra_optimization(self):
        before = {"fields": {"method": "RPBE1PBE", "Total Energy": -1.0}}
        after = {"fields": dict(before["fields"]), "diagnostics": {"instability_encountered": False}}
        self.assertFalse(needs_relaxation(before, after))
        after["fields"]["method"] = "UPBE1PBE"
        self.assertTrue(needs_relaxation(before, after))

    def test_ecp_and_g_conventions_are_checked_against_planned_basis(self):
        reference = {"explicit_electrons": 15, "alpha_electrons": 8, "beta_electrons": 7,
                     "ecp_core_electrons_by_atomic_number": {"75": 60}, "expected_g_shell_type": -4}
        fields = {"Atomic numbers": [75], "Charge": 0, "Multiplicity": 2,
                  "Number of electrons": 15, "Number of alpha electrons": 8,
                  "Number of beta electrons": 7, "Nuclear charges": [15.0], "Shell types": [0, -4]}
        original = {"fchk_identity": dict(fields)}
        self.assertTrue(identity_gate(reference, original, fields)["consistent"])
        fields["Shell types"] = [0, 4]
        self.assertIn("g-shell representation", identity_gate(reference, original, fields)["mismatches"])
        fields["Nuclear charges"] = [75.0]
        self.assertIn("ECP nuclear charges", identity_gate(reference, original, fields)["mismatches"])


if __name__ == "__main__":
    unittest.main()
