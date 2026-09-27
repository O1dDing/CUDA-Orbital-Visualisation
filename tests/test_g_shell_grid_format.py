import unittest

from g_shell_fixture_batch import grid_spec


class GaussianGridSyntaxTests(unittest.TestCase):
    def test_integer_valued_coordinates_are_serialized_as_real_tokens(self):
        text = grid_spec({'grid': {'origin_bohr': [-6, -6, -6], 'step_bohr': .375, 'shape': [33, 33, 33]}})
        rows = [row.split() for row in text.splitlines()]
        self.assertEqual([int(row[0]) for row in rows], [-1, -33, 33, 33])
        for row in rows:
            self.assertEqual(len(row), 4)
            self.assertTrue(all('.' in token for token in row[1:]))
        self.assertEqual([float(v) for v in rows[0][1:]], [-6., -6., -6.])
        self.assertEqual([[float(v) for v in row[1:]] for row in rows[1:]],
                         [[.375, 0., 0.], [0., .375, 0.], [0., 0., .375]])


if __name__ == '__main__':
    unittest.main()
