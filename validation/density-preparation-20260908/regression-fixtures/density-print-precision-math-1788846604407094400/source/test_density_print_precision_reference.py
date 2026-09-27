from decimal import Decimal,localcontext
import unittest

from density_print_precision_reference import (np,printed_value_interval,outer_sum_envelope,
                                                compare_matrix_to_envelope,gamma)

class PrintPrecisionControls(unittest.TestCase):
    def test_fortran_final_place_and_declared_rounding(self):
        for token,quantum in [('1.23456789D-03','1E-11'),('-2.300e2','0.1'),('.125','0.001'),('0.00000000E+00','1E-8')]:
            nearest=printed_value_interval(token,rounding_mode='nearest')
            whole=printed_value_interval(token,rounding_mode='one-last-place')
            self.assertEqual(Decimal(nearest['decimal_quantum']),Decimal(quantum))
            self.assertEqual(Decimal(nearest['exact_decimal_radius']),Decimal(quantum)/2)
            self.assertGreaterEqual(whole['radius'],nearest['radius'])
        for value in ['NaN','Inf','1.0junk','1.0 2.0','1e9999']:
            with self.assertRaises(ValueError):printed_value_interval(value,rounding_mode='nearest')
        with self.assertRaises(ValueError):printed_value_interval('1.0',rounding_mode='unknown')

    def test_decimal_nearest_endpoints_are_inside_float_intervals(self):
        for text in ['1.2345678901234567','0.00000000E+00','-3.1415926535897932','1E-999']:
            record=printed_value_interval(text,rounding_mode='nearest')
            center=Decimal(text);r=Decimal(record['exact_decimal_radius'])
            with localcontext() as ctx:
                ctx.prec=1100
                for endpoint in [center-r,center+r]:
                    self.assertLessEqual(abs(endpoint-Decimal.from_float(record['value'])),Decimal.from_float(record['radius']))
            self.assertGreater(record['radius'],0)

    def test_exact_decimal_perturbations_and_signed_cancellation(self):
        random=np.random.default_rng(20260908)
        c=random.normal(size=(4,17));r=np.full_like(c,5e-7)
        w=np.array([1.,-1.,.25,-.5,2.,0.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,1.])
        rw=np.full_like(w,1e-8)
        envelope=outer_sum_envelope(c,w,r,rw)
        for _ in range(40):
            exact_c=c+r*random.uniform(-1,1,c.shape)
            exact_w=w+rw*random.uniform(-1,1,w.shape)
            # Decimal arithmetic gives an independent non-BLAS multiplication path.
            with localcontext() as ctx:
                ctx.prec=100
                matrix=[[float(sum((Decimal.from_float(float(exact_c[i,k]))*
                                   Decimal.from_float(float(exact_w[k]))*
                                   Decimal.from_float(float(exact_c[j,k])) for k in range(len(w))),Decimal(0)))
                         for j in range(4)] for i in range(4)]
            self.assertTrue(compare_matrix_to_envelope(matrix,envelope)['entrywise_compatible'])

    def test_strict_same_input_check_does_not_inherit_print_allowance(self):
        c=np.array([[.512345678,-.314159265],[.271828182,.712345678]])
        weights=np.array([2.,1.])
        strict=outer_sum_envelope(c,weights)
        rounded_source=outer_sum_envelope(c,weights,np.full_like(c,.0005))
        damaged=strict['nominal'].copy();damaged[0,1]+=1e-4;damaged[1,0]+=1e-4
        self.assertFalse(compare_matrix_to_envelope(damaged,strict)['entrywise_compatible'])
        self.assertTrue(compare_matrix_to_envelope(damaged,rounded_source)['entrywise_compatible'])

    def test_orbital_phase_reorder_and_ao_permutation_invariants(self):
        c=np.array([[.12345,-.33333,.9],[.1,.99,-.22222]])
        w=np.array([1.,-.5,2.]);r=np.full_like(c,1e-6)
        base=outer_sum_envelope(c,w,r)
        order=[2,0,1];sign=np.array([-1.,1.,-1.])
        changed=outer_sum_envelope(c[:,order]*sign,w[order],r[:,order])
        self.assertTrue(compare_matrix_to_envelope(changed['nominal'],base)['entrywise_compatible'])
        np.testing.assert_allclose(base['source_radius'],changed['source_radius'],rtol=3e-14,atol=0)
        permuted=outer_sum_envelope(c[::-1],w,r[::-1])
        np.testing.assert_allclose(permuted['nominal'],base['nominal'][::-1,::-1],rtol=3e-15,atol=0)

    def test_zero_matrix_has_distinct_zero_and_unknown_source_uncertainty(self):
        c=np.array([[.5,.5],[.75,.75]])
        exact=outer_sum_envelope(c,[1.,-1.])
        uncertain=outer_sum_envelope(c,[1.,-1.],np.full_like(c,5e-9))
        self.assertTrue(np.all(exact['nominal']==0))
        self.assertTrue(np.all(uncertain['source_radius']>exact['source_radius']))
        self.assertTrue(compare_matrix_to_envelope(np.zeros((2,2)),exact)['entrywise_compatible'])
        self.assertFalse(compare_matrix_to_envelope(np.eye(2)*1e-8,exact)['entrywise_compatible'])

    def test_dimension_nonfinite_negative_radius_and_overflow_rejected(self):
        for c,w,r in [([[1.,2.]],[1.],None),([[float('nan')]],[1.],None),([[1.]],[1.],[[-1.]])]:
            with self.assertRaises(ValueError):outer_sum_envelope(c,w,r)
        with self.assertRaises((ValueError,FloatingPointError)):outer_sum_envelope([[1e300]],[1e300])
        with self.assertRaises(ValueError):gamma(True)
        with self.assertRaises(ValueError):compare_matrix_to_envelope([[float('inf')]],outer_sum_envelope([[1.]],[1.]))

if __name__=='__main__':unittest.main(verbosity=2)
