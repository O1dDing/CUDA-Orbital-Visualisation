import unittest

from gaussian_reference_inputs import make_followup_input, make_initial_input, reference_purpose


class ReferenceInputTests(unittest.TestCase):
    def setUp(self):
        self.case = {"purpose": {"kind": "scan_frame", "reference_geometry": "fixed"},
                     "original_input": {"route": "#p PBE1PBE/def2SVP SP"},
                     "fchk_identity": {"Number of atoms": 2, "Atomic numbers": [1, 1],
                         "Current cartesian coordinates": [0, 0, 0, 0, 0, 3],
                         "Charge": 0, "Multiplicity": 1}}
        self.reference = {"case_id": "test", "purpose": reference_purpose(self.case),
                          "gaussian_basis_keyword": "Gen", "shell_flags": "5D 7F 9G"}

    def test_scan_stays_fixed_and_contains_force_evidence(self):
        result = make_initial_input(self.reference, self.case, "H 0\nS 1 1.0\n1 1\n****")
        self.assertIn(" Force ", result)
        self.assertNotIn("Opt=", result)
        self.assertIn("Integral=SuperFineGrid", result)
        self.assertNotIn("Guess=Read", result)
        repaired = make_followup_input(self.reference, "UPBE1PBE", "relax_after_stability")
        self.assertIn("UPBE1PBE/ChkBasis Force", repaired)
        self.assertNotIn("Opt=", repaired)

    def test_changed_unrestricted_state_is_retained_and_refrequenced(self):
        self.reference["purpose"] = {"reference_geometry": "optimize"}
        result = make_followup_input(self.reference, "UPBE1PBE", "relax_after_stability")
        self.assertIn("UPBE1PBE/ChkBasis Opt=", result)
        self.assertIn(" Freq ", result)
        self.assertIn("Guess=Read", result)

    def test_original_optimize_takes_precedence_over_name_classification(self):
        self.case["purpose"] = {"kind": "specified_geometry", "reference_geometry": "fixed"}
        self.case["original_input"]["route"] = "#p PBE1PBE/def2SVP Opt=(Tight,CalcFC) Freq"
        self.assertEqual(reference_purpose(self.case)["reference_geometry"], "optimize")

    def test_atom_has_no_internal_vibrations(self):
        self.case["fchk_identity"]["Number of atoms"] = 1
        purpose = reference_purpose(self.case)
        self.assertEqual(purpose["reference_geometry"], "single_atom")
        self.assertIn("no_internal_nuclear_coordinates", purpose["vibration_applicability"])


if __name__ == "__main__":
    unittest.main()
