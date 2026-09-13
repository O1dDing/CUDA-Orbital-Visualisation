import argparse
import json
from pathlib import Path
import tempfile
import unittest

from reproduction_compare_batch import require_barrier


class ComparisonBarrierTests(unittest.TestCase):
    def test_missing_partial_and_wrong_version_batches_cannot_start_analysis(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            args = argparse.Namespace(output=output)
            identity = {"case_set_identity": "frozen_original_cases"}
            self.assertFalse(require_barrier(args, identity))
            (output/"reproduction").mkdir()
            path = output/"reproduction/batch.json"
            for terminal, count in ((False, 273), (True, 272)):
                path.write_text(json.dumps(dict(identity, all_terminal=terminal, case_count=count)))
                self.assertFalse(require_barrier(args, identity))
            path.write_text(json.dumps({"case_set_identity": "other_cases", "all_terminal": True, "case_count": 273}))
            with self.assertRaisesRegex(ValueError, "case sets differ"):
                require_barrier(args, identity)
            path.write_text(json.dumps(dict(identity, all_terminal=True, case_count=273)))
            self.assertTrue(require_barrier(args, identity))


if __name__ == "__main__":
    unittest.main()
