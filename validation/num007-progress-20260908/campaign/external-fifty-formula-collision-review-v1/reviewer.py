"""Resolve the remaining structural identity collisions, keeping every source."""
from collections import Counter
from pathlib import Path
import ctypes
import hashlib
import json
import math
import re
import sys
import time
base=Path(r'F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906')
sys.path.insert(0,r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json,physical_core_masks
k=ctypes.WinDLL('kernel32',use_last_error=True);k.GetCurrentProcess.restype=ctypes.c_void_p
k.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert k.SetProcessAffinityMask(k.GetCurrentProcess(),physical_core_masks(12)[8])
read=lambda p:json.loads(p.read_text(encoding='utf-8'));sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
out=base/'external-fifty-formula-collision-review-v1';out.mkdir(exist_ok=False)
cutoffs=[1.65,1.70,1.75,1.80]
atomic_json(out/'protocol.json',dict(created_epoch=time.time(),cutoffs_angstrom=cutoffs,
    declaration_before_measurements=True,scope='Remaining C/H and C/N/O/H constitutional-identity collisions only; distance graphs are checked across every declared cutoff.',
    preceding_identity_screen_sha256=sha(base/'external-fifty-identity-screen-v1.json'),formal_cases_added=0,gaussian_started=False))
(out/'reviewer.py').write_bytes(Path(__file__).read_bytes())
def graph(elements,xyz,cutoff):
    ids=[i for i,z in enumerate(elements) if z!=1]
    edges=[(a,b) for a,i in enumerate(ids) for b,j in enumerate(ids) if a<b and math.dist(xyz[i],xyz[j])<cutoff]
    degrees=[sum(i in e for e in edges) for i in range(len(ids))]
    return dict(source_atom_indices=ids,edges=[list(e) for e in edges],sorted_degrees=sorted(degrees),
        oxygen_nitrogen_edges=sum({elements[ids[a]],elements[ids[b]]}=={7,8} for a,b in edges))
def swept(elements,xyz):return [dict(cutoff_angstrom=c,**graph(elements,xyz,c)) for c in cutoffs]
records={r['candidate_id']:r for r in read(base/'external-coordinate-preparation-v2/summary.json')['records']}
pair=[]
for ident in ('PREP-022','PREP-023'):
    row=records[ident];raw_path=Path(row['raw_file']);normalized=Path(row['normalized_file'])
    assert sha(raw_path)==row['raw_sha256'] and sha(normalized)==row['normalized_sha256']
    raw=read(raw_path)['PC_Compounds'][0];z=raw['atoms']['element'];aids=raw['atoms']['aid']
    xyz=read(normalized)['conformers'][0]['coordinates_angstrom'];heavy=[a for a,n in zip(aids,z) if n!=1];index={a:i for i,a in enumerate(heavy)}
    declared=sorted(tuple(sorted((index[a],index[b]))) for a,b in zip(raw['bonds']['aid1'],raw['bonds']['aid2']) if a in index and b in index)
    results=swept(z,xyz);agrees=all(sorted(map(tuple,r['edges']))==declared for r in results)
    assert agrees
    pair.append(dict(candidate_id=ident,raw_sha256=row['raw_sha256'],normalized_sha256=row['normalized_sha256'],
        inchikey=row['inchikey'],source_declared_heavy_edges=declared,coordinate_sweeps=results,
        declared_and_coordinate_edges_agree_at_all_cutoffs=agrees))
pair_distinct=all(a['oxygen_nitrogen_edges']!=b['oxygen_nitrogen_edges'] for a,b in zip(pair[0]['coordinate_sweeps'],pair[1]['coordinate_sweeps']))
atomic_json(out/'oxazole-isoxazole.json',dict(candidates=pair,constitutional_graph_nonisomorphism_established=pair_distinct,
    reason='An element-preserving graph isomorphism must preserve the number of oxygen-nitrogen edges. The two source graphs have zero and one respectively.'))
folder=base/'benchmark-geometry-preparation-v1/primary-records/PREP-043';receipt=read(folder/'receipt.json')
assert sha(folder/'directory.json')==receipt['directory_sha256']
author=next(r for r in receipt['files'] if r['name']=='coord');raw=(folder/'coord').read_bytes()
assert sha(folder/'coord')==author['sha256']
assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==author['git_blob_sha1']
assert next(r for r in read(folder/'directory.json') if r['name']=='coord')['sha']==author['git_blob_sha1']
z=[];xyz=[];inside=False
for line in raw.decode('ascii').splitlines():
    if line.strip()=='$coord':inside=True;continue
    if inside and line.lstrip().startswith('$'):break
    if inside and line.strip():
        fields=line.split();assert len(fields)==4
        z.append({'c':6,'h':1}[fields[3].lower()]);xyz.append([float(v)*.529177210903 for v in fields[:3]])
benzyl=swept(z,xyz)
case=next(r for r in read(base/'campaign.json')['cases'] if r['case_id']=='OLD-097')
source=Path(case['original_fchk']);assert sha(source)==case['original_fchk_sha256'];lines=source.read_text(encoding='ascii').splitlines()
def field(label,convert):
    pos=[i for i,l in enumerate(lines) if l[:40].strip()==label];assert len(pos)==1
    i=pos[0];n=int(re.search(r'N=\s*(\d+)',lines[i]).group(1));tokens=[]
    for line in lines[i+1:]:
        if len(tokens)>=n:break
        tokens.extend(line.split())
    assert len(tokens)==n
    return [convert(v.replace('D','E')) for v in tokens]
old_z=field('Atomic numbers',int);old_flat=field('Current cartesian coordinates',float)
assert old_z==case['fchk_identity']['Atomic numbers'] and old_flat==case['fchk_identity']['Current cartesian coordinates']
assert Counter(z)==Counter(old_z)==Counter({6:7,1:7})
old_xyz=[[v*.529177210903 for v in old_flat[i:i+3]] for i in range(0,len(old_flat),3)];tropylium=swept(old_z,old_xyz)
stable=lambda rows:all(r['edges']==rows[0]['edges'] for r in rows)
benzyl_distinct=stable(benzyl) and stable(tropylium) and all(a['sorted_degrees']!=b['sorted_degrees'] for a,b in zip(benzyl,tropylium))
atomic_json(out/'benzyl-tropylium.json',dict(candidate_id='PREP-043',original_case_id='OLD-097',
    author_coord_sha256=author['sha256'],author_coord_git_blob_sha1=author['git_blob_sha1'],author_commit=receipt['commit'],
    original_fchk_sha256=case['original_fchk_sha256'],benzyl_sweeps=benzyl,tropylium_sweeps=tropylium,
    constitutional_graph_nonisomorphism_established=benzyl_distinct,
    reason='The carbon degree sequences differ at every declared cutoff: a six-membered ring with a pendant carbon versus a seven-membered ring. Charge differences are not used to count new structures.'))
prior=read(base/'external-fifty-identity-screen-v1.json');assert prior['external_same_formula_groups']==[['PREP-022','PREP-023']]
for r in prior['records']:
    if r['candidate_id']=='PREP-043':
        assert r['original_same_formula_cases']==['OLD-097'];r['nonidentity_to_originals_established']=benzyl_distinct
        r['constitutional_distinction']=dict(path=str(out/'benzyl-tropylium.json'),sha256=sha(out/'benzyl-tropylium.json'))
summary=dict(created_epoch=time.time(),prepared_candidates=50,original_cases=273,
    nonidentity_to_originals_established=sum(r['nonidentity_to_originals_established'] for r in prior['records']),
    all_fifty_pairwise_constitutionally_distinct=pair_distinct,remaining_unresolved_identity_collisions=0 if pair_distinct and benzyl_distinct else 1,
    records=prior['records'],new_evidence={p.name:sha(p) for p in out.glob('*.json')},
    formal_external_cases=0,accepted_quality_references=0,gaussian_started=False,
    limitation='Identity screening only: none of these results certify equilibrium geometry, electronic state, physical reference quality, or COV correctness.')
atomic_json(out/'summary.json',summary)
print(json.dumps({k:v for k,v in summary.items() if k not in ('records','new_evidence')},ensure_ascii=False),flush=True)
