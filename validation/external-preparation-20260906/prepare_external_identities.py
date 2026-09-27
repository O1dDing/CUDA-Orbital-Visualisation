"""Read-only public identity lookup; does not freeze cases or start Gaussian."""
from pathlib import Path
import collections
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request

workspace = Path(__file__).resolve().parent.parent
root = workspace / 'outputs/cov-complete-validation-20260906/external-preparation-v1'
spec_path = root / 'candidates.json'
spec = json.loads(spec_path.read_text(encoding='utf-8'))
source_sha = hashlib.sha256(spec_path.read_bytes()).hexdigest()
campaign = json.loads((root.parent / 'campaign.json').read_text(encoding='utf-8'))
symbols = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn '
           'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr '
           'Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn '
           'Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()

def atomic(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)

def formula_counts(text):
    parts = re.findall(r'([A-Z][a-z]?)([0-9]*)', text)
    if ''.join(a+b for a,b in parts) != text:
        raise ValueError('Unsupported formula syntax: ' + text)
    result = collections.Counter()
    for element, number in parts:
        result[element] += int(number or 1)
    return dict(sorted(result.items()))

old_formulas = collections.defaultdict(list)
for case in campaign['cases']:
    counts = collections.Counter(symbols[z] for z in case['fchk_identity']['Atomic numbers'])
    old_formulas[tuple(sorted(counts.items()))].append(case['case_id'])

records = []
for candidate in spec['candidates']:
    ident = candidate['candidate_id']
    out = root / 'identities' / ident
    out.mkdir(parents=True, exist_ok=True)
    status_path = out / 'identity.json'
    if status_path.exists():
        cached = json.loads(status_path.read_text(encoding='utf-8'))
        if cached['candidate_spec_sha256'] != source_sha:
            raise RuntimeError('Candidate list changed; preserve this preparation and create a new version')
        records.append(cached)
        continue
    expected = formula_counts(candidate['expected_formula'])
    record = {**candidate, 'candidate_spec_sha256': source_sha,
              'original_cases_with_same_element_counts': old_formulas.get(tuple(expected.items()), []),
              'lookup_epoch': time.time(), 'lookup_status': 'source-coordinate-work-pending',
              'identity_candidates': [], 'attempts': [], 'formal_case': False,
              'dedup_scope': 'formula screen and PubChem identity only; molecular graph and metal connectivity still require review'}
    if not candidate.get('skip_name_lookup'):
        url = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' +
               urllib.parse.quote(candidate['query'], safe='') +
               '/property/MolecularFormula,ConnectivitySMILES,InChIKey,IUPACName,Charge/JSON')
        record['query_url'] = url
        for attempt in range(1, 3):
            started = time.time()
            evidence = {'attempt': attempt, 'started_epoch': started}
            try:
                request = urllib.request.Request(url, headers={'User-Agent': 'COV-validation-identity-preparation/1.0', 'Accept': 'application/json'})
                with urllib.request.urlopen(request, timeout=20) as response:
                    payload = response.read()
                    evidence['http_status'] = response.status
                path = out / ('response-' + str(attempt) + '.json')
                path.write_bytes(payload)
                evidence.update(path=str(path), sha256=hashlib.sha256(payload).hexdigest(), bytes=len(payload))
                props = json.loads(payload)['PropertyTable']['Properties']
                for prop in props:
                    prop['formula_matches_request'] = formula_counts(prop['MolecularFormula']) == expected
                    prop['charge_matches_request'] = prop.get('Charge') == candidate['expected_charge']
                    prop['compound_page'] = 'https://pubchem.ncbi.nlm.nih.gov/compound/' + str(prop['CID'])
                record['identity_candidates'] = props
                matches = [p for p in props if p['formula_matches_request'] and p['charge_matches_request']]
                record['lookup_status'] = 'unique-metadata-match' if len(matches) == 1 else 'ambiguous-or-mismatched-metadata'
                evidence['elapsed_seconds'] = time.time() - started
                record['attempts'].append(evidence)
                break
            except Exception as exc:
                evidence.update(error=type(exc).__name__ + ': ' + str(exc), elapsed_seconds=time.time()-started)
                if isinstance(exc, urllib.error.HTTPError):
                    evidence['http_status'] = exc.code
                    payload = exc.read()
                    path = out / ('response-' + str(attempt) + '.error')
                    path.write_bytes(payload)
                    evidence.update(path=str(path), sha256=hashlib.sha256(payload).hexdigest(), bytes=len(payload))
                record['attempts'].append(evidence)
                record['lookup_status'] = 'lookup-failed-retained'
                if isinstance(exc, urllib.error.HTTPError) and exc.code in (400, 404):
                    break
                if attempt < 2:
                    time.sleep(2)
            finally:
                time.sleep(max(0, 0.6 - (time.time() - started)))
    atomic(status_path, record)
    records.append(record)
    atomic(root / 'progress.json', {'candidate_spec_sha256': source_sha, 'terminal_candidates': len(records),
           'expected_candidates': len(spec['candidates']), 'counts': dict(collections.Counter(r['lookup_status'] for r in records)),
           'formal_external_cases': 0, 'gaussian_started': False, 'REF-001': 'paused-by-user'})
    if len(records) % 10 == 0:
        print(json.dumps({'terminal_candidates': len(records), 'last': ident}), flush=True)

connectivity = collections.defaultdict(list)
for record in records:
    if record['lookup_status'] == 'unique-metadata-match':
        match = next(p for p in record['identity_candidates'] if p['formula_matches_request'] and p['charge_matches_request'])
        connectivity[match['InChIKey'].split('-')[0]].append(record['candidate_id'])
atomic(root / 'identity-summary.json', {'candidate_spec_sha256': source_sha,
       'all_metadata_lookups_terminal': len(records) == len(spec['candidates']), 'candidate_count': len(records),
       'counts': dict(collections.Counter(r['lookup_status'] for r in records)),
       'same_connectivity_key_candidates': [ids for ids in connectivity.values() if len(ids)>1],
       'formula_disjoint_from_original_count': sum(not r['original_cases_with_same_element_counts'] for r in records),
       'formal_external_cases': 0, 'gaussian_started': False, 'REF-001': 'paused-by-user',
       'limitations': ['Metadata is not validated geometry, electronic state, stability, or COV evidence.',
                      'Formula screening cannot distinguish isomers; identical formula does not establish duplication.',
                      'Disconnected or standardized metal SMILES cannot establish a coordination shell.',
                      'All failed and ambiguous identities remain in this preparation; no member has entered the formal set.'],
       'records': records})
print(json.dumps({'metadata_terminal': len(records), 'formal_external_cases': 0}), flush=True)
