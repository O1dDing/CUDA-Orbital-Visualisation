"""Exercise false-pass hazards before reviewing actual COV density output."""
from copy import deepcopy
from pathlib import Path
import json
import unittest
from review_density_source_controls import check_control, unpack

ROOT = Path(r'F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\regression-fixtures\density-source-v1')
MANIFEST = json.loads((ROOT/'manifest.json').read_text(encoding='utf-8'))
SPECS = {r['control_id']:r for r in MANIFEST['controls']}

def declared_good(identifier):
    spec = SPECS[identifier]
    if spec['expected_total_status'].startswith(('invalid-input','parser-rejects-')):
        return dict(parse_status='rejected',error='The supplied control is invalid under its declared ordinary-state purpose.')
    observed = dict(parse_status='returned',nbasis=3,bond_orders=[],bond_order_provenance=0)
    for field in ('total','spin'):
        values = spec.get(field+'_density_lower_packed_decimal')
        observed[field+'_density_packed'] = [float(v) for v in values] if values else []
        observed[field+'_density_provenance'] = (1 if identifier=='DEN-C19' else 2) if values else 0
        observed[field+'_density_status'] = 'available' if values else 'missing-input'
        observed[field+'_density_reason'] = 'Explicit expected control definition' if values else 'Input does not define this complete density.'
        observed[field+'_density_model'] = 'declared shared determinant' if identifier in ('DEN-C02','DEN-C03') else 'declared control state'
        observed['production_reconstructed_'+field] = dict(returned=False,error='Missing occupied subspace') if identifier in ('DEN-C18','DEN-C19') else dict(returned=True,packed=observed[field+'_density_packed'])
    if observed['total_density_packed'] and observed['spin_density_packed']:
        s = [[float(x) for x in row] for row in MANIFEST['overlap_decimal']]
        p = unpack(observed['total_density_packed'],3)
        q = unpack(observed['spin_density_packed'],3)
        ps = [[sum(p[i][k]*s[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        qs = [[sum(q[i][k]*s[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(i+1,3):
                value = ps[i][j]*ps[j][i] + qs[i][j]*qs[j][i]
                if abs(value) >= .02: observed['bond_orders'].append(dict(atoms=[i,j],mayer=value,provenance=2))
        observed['bond_order_provenance'] = 2
    return observed

class ReviewerTests(unittest.TestCase):
    def review(self, ident, observed=None):
        return check_control(SPECS[ident],observed or declared_good(ident),MANIFEST)

    def test_all_declared_outcomes_and_zero_are_usable(self):
        for ident in SPECS:
            with self.subTest(ident=ident):
                result=self.review(ident)
                self.assertTrue(result['control_requirements_satisfied'],result['issues'])
                self.assertFalse(result['formal_case_pass'])
        self.assertEqual(declared_good('DEN-C09')['spin_density_packed'],[0.0]*6)

    def test_unknown_occupation_cannot_be_derived_zero(self):
        for field in ('total','spin'):
            bad=declared_good('DEN-C01');bad[field+'_density_packed']=[0.0]*6;bad[field+'_density_provenance']=2
            self.assertFalse(self.review('DEN-C01',bad)['control_requirements_satisfied'])
            bad=declared_good('DEN-C01');bad[field+'_density_provenance']=2
            self.assertFalse(self.review('DEN-C01',bad)['control_requirements_satisfied'])

    def test_ambiguous_spin_has_no_full_mayer_even_if_no_pairs(self):
        bad=declared_good('DEN-C04');bad['bond_order_provenance']=2
        result=self.review('DEN-C04',bad)
        self.assertIn('DENSITY-MAYER-INCOMPLETE',[i['code'] for i in result['issues']])
        bad=declared_good('DEN-C04');bad['spin_density_packed']=[0.0]*6;bad['spin_density_provenance']=2
        self.assertFalse(self.review('DEN-C04',bad)['control_requirements_satisfied'])

    def test_equal_counts_nonzero_spin_and_wrong_direction_detected(self):
        self.assertGreater(max(abs(x) for x in declared_good('DEN-C06')['spin_density_packed']),.1)
        bad=declared_good('DEN-C06');bad['spin_density_packed']=[0.0]*6
        self.assertFalse(self.review('DEN-C06',bad)['control_requirements_satisfied'])
        bad=declared_good('DEN-C02');bad['spin_density_packed']=declared_good('DEN-C16')['spin_density_packed']
        self.assertFalse(self.review('DEN-C02',bad)['control_requirements_satisfied'])

    def test_invalid_json_numbers_dimensions_and_over_rejection_fail(self):
        for value in (None,float('nan'),float('inf'),True):
            bad=declared_good('DEN-C06');bad['total_density_packed'][0]=value
            self.assertFalse(self.review('DEN-C06',bad)['control_requirements_satisfied'])
        bad=declared_good('DEN-C17');bad['spin_density_packed'].pop()
        self.assertFalse(self.review('DEN-C17',bad)['control_requirements_satisfied'])
        bad=dict(parse_status='rejected',error='Counts missing')
        self.assertFalse(self.review('DEN-C20',bad)['control_requirements_satisfied'])

    def test_producer_density_and_incomplete_reconstruction_are_distinct(self):
        for field in ('total','spin'):
            bad=declared_good('DEN-C19');bad[field+'_density_provenance']=2
            self.assertFalse(self.review('DEN-C19',bad)['control_requirements_satisfied'])
            bad=declared_good('DEN-C19');bad['production_reconstructed_'+field]=dict(returned=True,packed=[0.0]*6)
            self.assertFalse(self.review('DEN-C19',bad)['control_requirements_satisfied'])

    def test_shared_model_must_be_declared_or_explicitly_unresolved(self):
        bad=declared_good('DEN-C02');bad.pop('spin_density_model')
        self.assertFalse(self.review('DEN-C02',bad)['control_requirements_satisfied'])
        unresolved=declared_good('DEN-C02')
        unresolved.update(spin_density_packed=[],spin_density_provenance=0,spin_density_status='missing-input',
                          spin_density_reason='No shared determinant model was declared',bond_orders=[],bond_order_provenance=0)
        self.assertTrue(self.review('DEN-C02',unresolved)['control_requirements_satisfied'])

if __name__ == '__main__': unittest.main(verbosity=2)
