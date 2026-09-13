"""Formula/charge screen of available local XYZ data; no chemical identity assignment."""
from pathlib import Path
import collections
import ctypes
import hashlib
import json
import math
import re
import sys
import time

sys.path.insert(0,r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import physical_core_masks
kernel=ctypes.WinDLL('kernel32',use_last_error=True);kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
workspace=Path(__file__).resolve().parent.parent
source=Path(r'F:\CalChem\ORCA\Executor\datasets\Geometries\GMTKN55')
out=workspace/'outputs/cov-complete-validation-20260906/benchmark-geometry-preparation-v1'
out.mkdir(exist_ok=True)
path=out/'local-formula-screen.json'
if path.exists():raise FileExistsError('Retain previous source screen')
summary=json.loads((out.parent/'external-coordinate-preparation-v2/summary.json').read_text(encoding='utf-8'))
targets=summary['records']
formula=lambda s:tuple(sorted((a,int(n or 1)) for a,n in re.findall(r'([A-Z][a-z]?)([0-9]*)',s)))
index=collections.defaultdict(list);errors=[];files=list(source.rglob('*.xyz'));started=time.time()
for file in files:
    try:
        raw=file.read_bytes();rows=raw.decode('utf-8-sig').splitlines();n=int(rows[0]);charge,multiplicity=map(int,rows[1].split())
        if len(rows[2:])!=n and any(line.strip() for line in rows[2+n:]):raise ValueError('Unexpected coordinate records')
        coords=[line.split() for line in rows[2:2+n]]
        if len(coords)!=n or any(len(r)!=4 for r in coords):raise ValueError('Incomplete XYZ coordinates')
        for row in coords:
            if any(not math.isfinite(float(v.replace('D','E'))) for v in row[1:]):raise ValueError('Nonfinite XYZ coordinate')
        counts=tuple(sorted(collections.Counter(row[0] for row in coords).items()))
        record=dict(file=str(file),relative_path=file.relative_to(source).as_posix(),sha256=hashlib.sha256(raw).hexdigest(),
                    atom_count=n,charge=charge,multiplicity=multiplicity,formula_counts=dict(counts),
                    original_source_coordinates_verified=False)
        index[(counts,charge)].append(record)
    except Exception as error:errors.append(dict(file=str(file),error=f'{type(error).__name__}: {error}'))
records=[]
for candidate in targets:
    hits=index.get((formula(candidate['expected_formula']),candidate['expected_charge']),[])
    records.append(dict(candidate_id=candidate['candidate_id'],query=candidate['query'],expected_formula=candidate['expected_formula'],
        expected_charge=candidate['expected_charge'],pubchem_status=candidate['status'],formula_charge_matches=hits,
        geometry_identity_confirmed=False,formal_case=False))
result=dict(source=str(source),source_provenance='Local package metadata has empty source/citations; verify against primary distribution before accepting',
            xyz_files_seen=len(files),xyz_files_parsed=len(files)-len(errors),errors=errors,records=records,
            matched_candidates=sum(bool(r['formula_charge_matches']) for r in records),
            matched_candidates_without_pubchem_3d=sum(bool(r['formula_charge_matches']) and r['pubchem_status']!='computed-starting-conformer-retrieved' for r in records),
            formal_external_cases=0,gaussian_started=False,seconds=time.time()-started)
path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['records','errors']},ensure_ascii=False))
print(json.dumps([dict(candidate_id=r['candidate_id'],query=r['query'],matches=[x['relative_path'] for x in r['formula_charge_matches']]) for r in records if r['formula_charge_matches'] and r['pubchem_status']!='computed-starting-conformer-retrieved'],ensure_ascii=False))
