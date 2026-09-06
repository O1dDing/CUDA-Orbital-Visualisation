"""Review all 67 public identity lookups, preserving their original receipts."""
from pathlib import Path
from collections import Counter, defaultdict
import copy
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request

root = Path(__file__).resolve().parent.parent / 'outputs/cov-complete-validation-20260906/external-preparation-v1'
destination = root / 'identity-review-v2'
destination.mkdir(exist_ok=True)
if (destination / 'summary.json').exists():
    raise FileExistsError('Preserve completed identity review')
initial = json.loads((root / 'identity-summary.json').read_text(encoding='utf-8'))
assert initial['all_metadata_lookups_terminal']

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path, data): path.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')

def formula(text):
    # PubChem puts charge AFTER composition as +, - or +2/-2. This is not a
    # general chemical-formula parser and does not interpret parentheses/salts.
    charge_match = re.search(r'([+-])([0-9]*)$', text)
    suffix_charge = 0
    composition = text
    if charge_match:
        suffix_charge = (1 if charge_match[1]=='+' else -1) * int(charge_match[2] or 1)
        composition = text[:charge_match.start()]
    parts = re.findall(r'([A-Z][a-z]?)([0-9]*)', composition)
    if not parts or ''.join(a+b for a,b in parts)!=composition:
        raise ValueError('Unsupported PubChem formula: '+text)
    counts = Counter()
    for a,b in parts: counts[a] += int(b or 1)
    return dict(sorted(counts.items())), suffix_charge

followups = {
    'PREP-045': ('inchikey', 'WTFNSXYULBQCQV-UHFFFAOYSA-N', 'https://cccbdb.nist.gov/exp2x.asp?casno=2143580&charge=0'),
    'PREP-058': ('cid', '11963622', 'https://pubchem.ncbi.nlm.nih.gov/compound/11963622'),
}
records=[]
for original in initial['records']:
    record = copy.deepcopy(original)
    record['original_lookup_status'] = original['lookup_status']
    record['review_evidence'] = []
    props=[]
    for attempt in original['attempts']:
        if attempt.get('http_status')==200 and attempt.get('path'):
            path=Path(attempt['path'])
            assert sha(path)==attempt['sha256']
            props=json.loads(path.read_bytes())['PropertyTable']['Properties']
            record['review_evidence']=[{'path':str(path),'sha256':sha(path),'source':'original-query'}]
    if record['candidate_id'] in followups:
        namespace, value, source = followups[record['candidate_id']]
        url=('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'+namespace+'/'+urllib.parse.quote(value,safe='')+
             '/property/MolecularFormula,ConnectivitySMILES,InChIKey,IUPACName,Charge/JSON')
        path=destination/(record['candidate_id']+'-followup.json')
        receipt_path=destination/(record['candidate_id']+'-followup-receipt.json')
        if not path.exists():
            try:
                with urllib.request.urlopen(url, timeout=20) as response: payload=response.read()
                path.write_bytes(payload)
                save(receipt_path, {'url':url,'source_for_identifier':source,'sha256':sha(path),'epoch':time.time()})
            except Exception as exc:
                record['followup_error']=type(exc).__name__+': '+str(exc)
        if path.exists():
            receipt=json.loads(receipt_path.read_text(encoding='utf-8'))
            assert sha(path)==receipt['sha256']
            props=json.loads(path.read_bytes())['PropertyTable']['Properties']
            record['review_evidence'].append({'path':str(path),'sha256':sha(path),'url':url,'identifier_source':source})
    expected, _ = formula(record['expected_formula'])
    for prop in props:
        counts, suffix_charge = formula(prop['MolecularFormula'])
        prop['formula_matches_request']=counts==expected
        prop['charge_matches_request']=prop.get('Charge')==record['expected_charge']
        prop['formula_charge_consistent']=suffix_charge==prop.get('Charge')
        prop['compound_page']='https://pubchem.ncbi.nlm.nih.gov/compound/'+str(prop['CID'])
        prop['connectivity_representation_disconnected']='.' in prop.get('ConnectivitySMILES','')
    matches=[p for p in props if p['formula_matches_request'] and p['charge_matches_request'] and p['formula_charge_consistent']]
    identities={p['InChIKey'] for p in matches}
    record['identity_candidates']=props
    record['matching_cids']=[p['CID'] for p in matches]
    record['lookup_status']=('unique-metadata-match' if len(identities)==1 else 'ambiguous-or-mismatched-metadata') if props else original['lookup_status']
    record['geometry_accepted']=False
    record['formal_case']=False
    records.append(record)
keys=defaultdict(list)
for r in records:
    if r['lookup_status']=='unique-metadata-match':
        p=next(p for p in r['identity_candidates'] if p['CID'] in r['matching_cids'])
        keys[p['InChIKey'].split('-')[0]].append(r['candidate_id'])
summary={'candidate_spec_sha256':initial['candidate_spec_sha256'],'original_lookup_summary_sha256':sha(root/'identity-summary.json'),
         'review_script_sha256':sha(Path(__file__)),'candidate_count':len(records),'all_metadata_reviews_terminal':True,
         'counts':dict(Counter(r['lookup_status'] for r in records)),
         'same_connectivity_key_candidates':[v for v in keys.values() if len(v)>1],
         'formula_disjoint_from_original_count':initial['formula_disjoint_from_original_count'],
         'formal_external_cases':0,'gaussian_started':False,'REF-001':'paused-by-user',
         'corrections':['v1 formula reader rejected PubChem charge suffixes; v2 reuses the preserved response bytes and checks charge explicitly.',
                        'Methylperoxy and nitroprusside use separately sourced identifier followups; original failed/mismatched name responses remain intact.'],
         'limitations':initial['limitations'],'records':records}
save(destination/'summary.json',summary)
print(json.dumps({k:v for k,v in summary.items() if k!='records'},ensure_ascii=False,indent=2))
