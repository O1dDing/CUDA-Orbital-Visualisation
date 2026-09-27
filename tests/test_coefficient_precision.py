from decimal import Decimal, localcontext
import struct
import unittest

from coefficient_precision import measure_coefficient_rounding, np


class CoefficientPrecisionTests(unittest.TestCase):
    def test_nonorthogonal_cancellation_against_decimal_polynomial_integral(self):
        # This analytic two-function integral is computed independently with
        # Decimal arithmetic and struct's binary32 conversion, not Cholesky.
        c = np.array([[100000.003], [-99999.997]])
        s = np.array([[1.0, .99999999], [.99999999, 1.0]])
        untouched = c.copy()
        measured = measure_coefficient_rounding(c, s)
        with localcontext() as context:
            context.prec = 60
            a, b, overlap = map(Decimal.from_float, (float(c[0, 0]), float(c[1, 0]), float(s[0, 1])))
            qa = Decimal.from_float(struct.unpack('<f', struct.pack('<f', float(a)))[0])
            qb = Decimal.from_float(struct.unpack('<f', struct.pack('<f', float(b)))[0])
            da, db = qa-a, qb-b
            norm_squared = a*a + b*b + 2*a*b*overlap
            error_squared = da*da + db*db + 2*da*db*overlap
            expected = float((error_squared/norm_squared).sqrt())
        self.assertAlmostEqual(float(measured['relative_errors'][0])/expected, 1.0, places=8)
        self.assertTrue(np.array_equal(c, untouched))
        self.assertGreater(expected, 1e-5)

    def test_exact_binary32_values_and_rectangular_matrix(self):
        c = np.array([[1.0], [.5], [-2.0]])
        measured = measure_coefficient_rounding(c, np.diag([2., 3., 4.]))
        self.assertEqual(float(measured['error_norms'][0]), 0.0)
        self.assertAlmostEqual(float(measured['original_norms'][0]), np.sqrt(18.75))
        self.assertEqual(measured['relative_errors'].shape, (1,))

    def test_missing_or_invalid_metric_evidence_is_rejected(self):
        fixtures = [(np.zeros((2, 1)), np.eye(2)),
                    (np.ones((2, 1)), np.array([[1., 1.], [1., 1.]])),
                    (np.ones((2, 1)), np.array([[1., .1], [.2, 1.]])),
                    (np.ones((2, 1)), np.eye(3)),
                    (np.array([[np.nan]]), np.eye(1))]
        for c, s in fixtures:
            with self.subTest(c=c.tolist(), s=s.tolist()), self.assertRaises((ValueError, np.linalg.LinAlgError)):
                measure_coefficient_rounding(c, s)


if __name__ == '__main__':
    unittest.main()
