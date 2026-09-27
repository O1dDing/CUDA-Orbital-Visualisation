"""Archive PubChem computed conformers for reviewed identities, never reference geometries."""
from pathlib import Path
import collections
import ctypes
import hashlib
import json
import math
import re
import sys
import time
import urllib.error
import urllib.request

workspace = Path(__file__).resolve().parent.parent
source = workspace / 'outputs/cov-complete-validation-20260906/external-preparation-v1/identity-review-v2/summary.json'
prior = workspace / 'outputs/cov-complete-validation-20260906/external-coordinate-preparation-v1'
root = workspace / 'outputs/cov-complete-validation-20260906/external-coordinate-preparation-v2'
root.mkdir(exist_ok=True)
sys.path.insert(0, r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json, physical_core_masks
kernel = ctypes.WinDLL('kernel32', use_last_error=True)
kernel.GetCurrentProcess.restype = ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(), physical_core_masks(12)[8])
symbols = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn '
           'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr '
           'Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn').split()
identity_sha = hashlib.sha256(source.read_bytes()).hexdigest()
spec = dict(identity_review=str(source), identity_review_sha256=identity_sha,
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            purpose='Computed starting conformers only; not experimental coordinates or validated quantum references',
            cpu_core_slot=8, maximum_concurrent_requests=1, request_timeout_seconds=20,
            maximum_response_bytes=4*1024*1024, minimum_request_interval_seconds=0.7,
            automatic_retries=0, formal_external_cases=0, gaussian_started=False,
            independently_managed_gaussian='F:/CalChem/COV/Resume; not controlled by this preparation',
            source_documentation='https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest')
spec['revision_reason'] = 'Use reviewed matching_cids for original and identifier-followup records alike. Retain v1 responses and interrupted collector.'
spec_path = root / 'preparation.json'
if spec_path.exists():
    if json.loads(spec_path.read_text(encoding='utf-8')) != spec:
        raise RuntimeError('Preparation identity changed; preserve it and use a new version')
else:
    atomic_json(spec_path, spec)
records = []
for candidate in json.loads(source.read_text(encoding='utf-8'))['records']:
    ident = candidate['candidate_id']
    out = root / ident
    out.mkdir(exist_ok=True)
    receipt = out / 'receipt.json'
    if receipt.exists():
        row = json.loads(receipt.read_text(encoding='utf-8'))
        if row['identity_review_sha256'] != identity_sha:
            raise RuntimeError('Saved conformer identity differs from the reviewed candidate list')
        records.append(row)
        continue
    prior_receipt = prior / ident / 'receipt.json'
    if prior_receipt.exists():
        row = json.loads(prior_receipt.read_text(encoding='utf-8'))
        if row['identity_review_sha256'] != identity_sha:
            raise RuntimeError('Prior conformer identity differs from the reviewed candidate list')
        for file_key, sha_key in [('raw_file','raw_sha256'),('normalized_file','normalized_sha256'),('error_file','error_sha256')]:
            if row.get(file_key) and hashlib.sha256(Path(row[file_key]).read_bytes()).hexdigest() != row[sha_key]:
                raise RuntimeError('Prior evidence bytes changed')
        row['retained_prior_receipt'] = str(prior_receipt)
        row['retained_prior_receipt_sha256'] = hashlib.sha256(prior_receipt.read_bytes()).hexdigest()
        atomic_json(receipt, row)
        records.append(row)
        continue
    row = dict(candidate_id=ident, query=candidate['query'], family=candidate['family'],
               expected_formula=candidate['expected_formula'], expected_charge=candidate['expected_charge'],
               identity_review_sha256=identity_sha, started_epoch=time.time(),
               input_geometry_accepted=False, formal_case=False, gaussian_started=False)
    if candidate['lookup_status'] != 'unique-metadata-match':
        row.update(status='identity-unresolved-retained', lookup_status=candidate['lookup_status'])
    else:
        matches = [item for item in candidate['identity_candidates'] if item['CID'] in candidate['matching_cids']]
        assert len(matches) == 1
        match = matches[0]
        cid = match['CID']
        row.update(cid=cid, inchikey=match['InChIKey'],
                   disconnected_identity=match.get('connectivity_representation_disconnected', False),
                   original_cases_with_same_element_counts=candidate['original_cases_with_same_element_counts'])
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON?record_type=3d'
        row['url'] = url
        started = time.monotonic()
        try:
            req = urllib.request.Request(url, headers={'User-Agent':'COV-validation-conformer-preparation/1.0', 'Accept':'application/json'})
            with urllib.request.urlopen(req, timeout=20) as response:
                raw = response.read(4*1024*1024+1)
                row['http_status'] = response.status
            if len(raw)>4*1024*1024:
                raise ValueError('Conformer exceeds frozen response size bound')
            raw_path = out / 'pubchem-computed-3d.json'
            raw_path.write_bytes(raw)
            row.update(raw_file=str(raw_path), raw_sha256=hashlib.sha256(raw).hexdigest(), raw_bytes=len(raw))
            compounds = json.loads(raw)['PC_Compounds']
            if len(compounds) != 1:
                raise ValueError('Expected one explicitly identified compound')
            compound = compounds[0]
            if compound['id']['id']['cid'] != cid:
                raise ValueError('Compound id differs from requested identity')
            atoms = compound['atoms']
            aids = atoms['aid']
            elements = atoms['element']
            expected = collections.Counter({el:int(n or 1) for el,n in re.findall(r'([A-Z][a-z]?)([0-9]*)',candidate['expected_formula'])})
            observed = collections.Counter(symbols[z] for z in elements)
            if observed != expected or len(aids)!=len(elements) or len(set(aids))!=len(aids):
                raise ValueError('Explicit atom identities/counts differ from expected formula')
            if compound.get('charge', 0) != candidate['expected_charge']:
                raise ValueError('3D compound charge differs from requested identity')
            conformer_records = []
            for coord in compound.get('coords', []):
                if 2 not in coord.get('type', []):
                    continue
                caids = coord['aid']
                if len(caids)!=len(aids) or set(caids)!=set(aids):
                    raise ValueError('3D coordinates omit or duplicate atom identities')
                for conf in coord.get('conformers', []):
                    axes = [conf[key] for key in ('x','y','z')]
                    if any(len(axis)!=len(caids) for axis in axes):
                        raise ValueError('Coordinate arrays have inconsistent lengths')
                    by_aid = {aid:[axis[i] for axis in axes] for i,aid in enumerate(caids)}
                    xyz = [by_aid[aid] for aid in aids]
                    if any(not math.isfinite(v) for p in xyz for v in p):
                        raise ValueError('Nonfinite coordinate')
                    minimum = min((math.dist(xyz[i],xyz[j]) for i in range(len(xyz)) for j in range(i)),default=None)
                    if minimum is not None and minimum <= 0:
                        raise ValueError('Coincident atoms in supplied conformer')
                    conformer_records.append(dict(atom_ids=aids, atomic_numbers=elements,
                        coordinates_angstrom=xyz, minimum_separation_angstrom=minimum,
                        pubchem_coord_type=coord.get('type'), source_conformer_id=conf.get('id')))
            if not conformer_records:
                raise ValueError('No explicit 3D conformer in response')
            normalized = out / 'starting-conformers.json'
            atomic_json(normalized, dict(source_raw_sha256=row['raw_sha256'], cid=cid,
                provenance='PubChem computed 3D conformer; not an experimental reference',
                conformers=conformer_records, formal_case=False,
                independent_geometry_and_electronic_state_review='pending'))
            row.update(status='computed-starting-conformer-retrieved', atom_count=len(aids),
                       conformer_count=len(conformer_records), normalized_file=str(normalized),
                       normalized_sha256=hashlib.sha256(normalized.read_bytes()).hexdigest())
        except Exception as error:
            row.update(status='conformer-unavailable-or-rejected-retained', error=f'{type(error).__name__}: {error}')
            if isinstance(error, urllib.error.HTTPError):
                row['http_status'] = error.code
                raw = error.read(4*1024*1024)
                error_path = out / 'response.error'
                error_path.write_bytes(raw)
                row.update(error_file=str(error_path), error_sha256=hashlib.sha256(raw).hexdigest())
        finally:
            time.sleep(max(0, .7-(time.monotonic()-started)))
    row['finished_epoch'] = time.time()
    atomic_json(receipt, row)
    records.append(row)
    progress = dict(terminal_candidates=len(records), expected_candidates=67,
                    counts=dict(collections.Counter(record['status'] for record in records)),
                    formal_external_cases=0, gaussian_started=False)
    atomic_json(root/'progress.json', progress)
    if len(records)%10 == 0:
        print(json.dumps(progress), flush=True)
atomic_json(root/'summary.json', {**spec, **progress, 'all_terminal':len(records)==67, 'records':records})
print(json.dumps(progress),flush=True)
