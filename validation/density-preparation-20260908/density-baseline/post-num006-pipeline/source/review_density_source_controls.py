"""Independent finite-AO control review; only runs after the entire control batch.

The pre-existing manifest supplies matrices, interpretation requirements and the
fixed tolerance. This checker never infers an unknown spin partition from COV.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

def unpack(values, n):
    if not isinstance(values, list) or len(values) != n*(n+1)//2:
        raise ValueError('Missing or dimensionally invalid packed density')
    if any(isinstance(v, bool) or not isinstance(v, (float,int)) or not math.isfinite(v) for v in values):
        raise ValueError('Nonfinite or nonnumeric density entry')
    matrix = [[0.0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1):
            matrix[i][j] = matrix[j][i] = float(values[k]); k += 1
    return matrix

def matrix_check(values, expected, n, tolerance):
    try:
        unpack(values, n)
    except ValueError as error:
        return dict(passed=False, reason=str(error))
    difference = max(abs(v-float(e)) for v,e in zip(values,expected))
    return dict(passed=difference <= tolerance, maximum_absolute_error=difference,
                fixed_absolute_tolerance=tolerance)

def check_control(spec, observed, manifest):
    """Return bounded diagnostics; a mathematical control never passes a molecule."""
    n = len(manifest['atomic_numbers'])
    tolerance = manifest['fixed_matrix_absolute_tolerance']
    assert tolerance > 0 and math.isfinite(tolerance)
    issues = []; checks = {}
    rejected = observed.get('parse_status') == 'rejected'
    def fail(code, detail):
        issues.append(dict(code=code, detail=detail))
    # Malformed coefficient blocks and explicitly invalid ordinary-state
    # occupations may be rejected. Missing metadata must preserve orbital use.
    expected_p = spec['expected_total_status']
    if rejected:
        permitted = expected_p.startswith(('invalid-input','parser-rejects-'))
        checks['parser'] = dict(passed=permitted, status='rejected', reason=observed.get('error'))
        if not permitted: fail('DENSITY-OVER-REJECTION', 'Usable orbitals or producer density were rejected as a whole.')
        if not observed.get('error'): fail('DENSITY-REASON-MISSING', 'Parser rejection lacks an explanatory error.')
        return dict(control_id=spec['control_id'], checks=checks, issues=issues,
                    control_requirements_satisfied=not issues, formal_case_pass=False)
    if observed.get('parse_status') != 'returned':
        fail('DENSITY-EVIDENCE-MISSING','Neither successful parse nor retained parser rejection.')
        return dict(control_id=spec['control_id'], checks=checks, issues=issues,
                    control_requirements_satisfied=False, formal_case_pass=False)
    if expected_p.startswith('parser-rejects-'):
        fail('DENSITY-COEFFICIENT-DIMENSION','Malformed coefficient block was accepted.')
    if observed.get('nbasis') != n:
        fail('DENSITY-AO-DIMENSION','Production AO dimension differs from the independent input.')
    available = {}
    for field in ('total','spin'):
        status = spec[f'expected_{field}_status']
        values = observed.get(f'{field}_density_packed')
        provenance = observed.get(f'{field}_density_provenance')
        expected = spec.get(f'{field}_density_lower_packed_decimal')
        reason = observed.get(f'{field}_density_reason')
        model = observed.get(f'{field}_density_model')
        observed_status = observed.get(f'{field}_density_status')
        conditional = status.startswith(('available-under-','same-interpretation-as-'))
        if expected is not None:
            # A stricter explicitly unknown-model answer is allowed only by the
            # already frozen conditional shared-determinant controls C02/C03.
            if conditional and values == [] and provenance == 0 and reason and observed_status == 'missing-input':
                checks[field] = dict(passed=True, status='explicitly-unresolved-model', reason=reason)
                available[field] = False
                continue
            result = matrix_check(values, expected, n, tolerance)
            checks[field] = result
            available[field] = result['passed'] and provenance in (1,2)
            if not result['passed']: fail('DENSITY-'+field.upper()+'-MATRIX', result)
            if provenance not in (1,2): fail('DENSITY-'+field.upper()+'-PROVENANCE','A known matrix lacks a valid source.')
            if conditional and not model:
                fail('DENSITY-SPIN-MODEL-EVIDENCE','The shared-integer determinant interpretation is not exposed by the production result.')
            if spec['control_id'] == 'DEN-C19' and provenance != 1:
                fail('DENSITY-PRODUCER-OVERWRITE','Complete stored density must retain producer provenance.')
        else:
            # Unavailable means both no fabricated matrix and no Derived tag.
            absent = values == [] and provenance == 0
            checks[field] = dict(passed=absent, expected_status=status,
                                 observed_status=observed_status, provenance=provenance,
                                 matrix_entries=len(values) if isinstance(values,list) else None)
            available[field] = False
            if not absent:
                fail('DENSITY-'+field.upper()+'-FALSE-AVAILABILITY','Unknown/invalid density is represented by entries or a usable provenance.')
            if not reason or observed_status not in ('missing-input','invalid-input'):
                fail('DENSITY-'+field.upper()+'-STATUS-EVIDENCE','No explicit availability status and explanation for unknown/invalid density.')
    # C19 is deliberately producer-complete but MO-incomplete. Raw
    # reconstruction must not pretend the missing occupied direction is zero.
    if spec['control_id'] in ('DEN-C18','DEN-C19'):
        for field in ('total','spin'):
            raw = observed.get('production_reconstructed_'+field,{})
            incomplete_reported = raw.get('returned') is False or raw.get('packed') == []
            checks['raw_reconstructed_'+field] = dict(passed=incomplete_reported)
            if not incomplete_reported:
                fail('DENSITY-TRUNCATED-OCCUPIED-SPACE','Raw '+field+' reconstruction returned a complete-looking matrix from incomplete occupied space.')
    # A total/spin Mayer value needs both P and Q. A P-only measure would
    # need a separately named definition; the legacy field is the full Mayer.
    if not all(available.values()):
        blocked = observed.get('bond_order_provenance') == 0 and observed.get('bond_orders') == []
        checks['full_mayer_availability'] = dict(passed=blocked)
        if not blocked:
            fail('DENSITY-MAYER-INCOMPLETE','Full total/spin Mayer analysis is marked usable without both known densities.')
    else:
        p = unpack([float(v) for v in spec['total_density_lower_packed_decimal']],n)
        q = unpack([float(v) for v in spec['spin_density_lower_packed_decimal']],n)
        s = [[float(v) for v in row] for row in manifest['overlap_decimal']]
        mul = lambda a,b: [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        ps = mul(p,s); qs = mul(q,s)
        expected = {(i,j):ps[i][j]*ps[j][i]+qs[i][j]*qs[j][i] for i in range(n) for j in range(i+1,n)}
        # This pre-existing production display floor determines records, not
        # the independently fixed matrix tolerance or a physical bond cutoff.
        floor = 0.02
        required = {pair:value for pair,value in expected.items() if abs(value) >= floor}
        bonds = observed.get('bond_orders')
        seen = {}
        if not isinstance(bonds,list): bonds = []; fail('DENSITY-MAYER-EVIDENCE','Missing actual pair records.')
        for bond in bonds:
            pair = tuple(bond.get('atoms',[])); value = bond.get('mayer')
            if pair in seen or pair not in expected or not isinstance(value,(int,float)) or isinstance(value,bool) or not math.isfinite(value):
                fail('DENSITY-MAYER-RECORD','Duplicate, invalid or unexpected atom-pair record.'); continue
            seen[pair] = value
            if abs(value-expected[pair]) > tolerance:
                fail('DENSITY-MAYER-VALUE',dict(atoms=pair,observed=value,expected=expected[pair]))
        if set(seen) != set(required):
            fail('DENSITY-MAYER-COVERAGE',dict(required=list(required),observed=list(seen)))
        if observed.get('bond_order_provenance') not in (1,2):
            fail('DENSITY-MAYER-PROVENANCE','Valid full Mayer result is not available.')
        checks['full_mayer'] = dict(expected_records=len(required), actual_records=len(seen), record_floor=floor)
    return dict(control_id=spec['control_id'], checks=checks, issues=issues,
                control_requirements_satisfied=not issues, formal_case_pass=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    read = lambda p: json.loads(p.read_text(encoding='utf-8'))
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    collection = read(args.baseline/'collection-complete.json')
    manifest = read(args.manifest)
    assert collection['all_terminal'] and collection['terminal_controls'] == manifest['control_count'] == 20
    assert len(collection['records']) == 20
    assert sha(args.manifest) == read(args.baseline/'identity.json')['controls_manifest_sha256']
    if args.output.exists(): raise FileExistsError('Keep prior control adjudication; use a new review version')
    args.output.mkdir(parents=True)
    (args.output/'reviewer.py').write_bytes(Path(__file__).read_bytes())
    records = []; actual = {}; receipts = {r['control_id']:r for r in collection['records']}
    assert len(receipts) == 20
    for spec in manifest['controls']:
        ident = spec['control_id']; receipt = receipts[ident]
        observed_path = args.baseline/'cases'/ident/'observed.json'
        if receipt['status'] != 'collected' or not observed_path.exists():
            records.append(dict(control_id=ident,control_requirements_satisfied=False,formal_case_pass=False,
                issues=[dict(code='DENSITY-COLLECTION-GAP',detail=receipt)])); continue
        assert sha(observed_path) == receipt['observed_sha256']
        assert sha(args.manifest.parent/spec['file']) == receipt['input_sha256'] == spec['sha256']
        observed = read(observed_path); actual[ident] = observed
        result = check_control(spec,observed,manifest)
        result['observed_sha256'] = sha(observed_path)
        records.append(result)
    relations = []
    for spec in manifest['controls']:
        relation = spec.get('metamorphic_relation')
        if not relation: continue
        for field in ('total','spin'):
            try:
                left = actual[spec['control_id']][field+'_density_packed']
                right = actual[relation['control']][field+'_density_packed']
                unpack(left,len(manifest['atomic_numbers'])); unpack(right,len(manifest['atomic_numbers']))
                sign = -1 if relation[field] == 'negated' else 1
                error = max(abs(a-sign*b) for a,b in zip(left,right))
                result = dict(status='measured', passed=error <= manifest['fixed_matrix_absolute_tolerance'],maximum_absolute_error=error)
            except (KeyError,ValueError,TypeError) as error:
                result = dict(status='evidence-missing',passed=False,reason=str(error))
            relations.append(dict(control_id=spec['control_id'],reference_control=relation['control'],field=field,**result))
    summary = dict(all_control_inputs_terminal=True, controls=20,
        controls_satisfying_recorded_requirements=sum(r['control_requirements_satisfied'] for r in records),
        controls_with_failures_or_evidence_gaps=sum(not r['control_requirements_satisfied'] for r in records),
        formal_cases_added=0,formal_case_passes=0,baseline=str(args.baseline),
        manifest_sha256=sha(args.manifest),records=records,metamorphic_relations=relations,
        scientific_scope='Finite-AO parser/density/Mayer diagnostic controls only; no molecule, physical-reference or rendering acceptance.')
    (args.output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('records','metamorphic_relations')},ensure_ascii=False))

if __name__ == '__main__': main()
