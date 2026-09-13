"""Extract literal Jmol XYZ strings as data, without executing page JavaScript."""
from collections import Counter
from decimal import Decimal
from pathlib import Path
import ctypes
import hashlib
import json
import math
import re
import sys
import time

workspace=Path(r'F:\Codex\2026-09-05\branch-15')
sys.path.insert(0,r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json,physical_core_masks
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
root=workspace/'outputs/cov-complete-validation-20260906/nist-followed-geometries-v1'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
collection=read(root/'collection-complete.json')
assert collection['all_terminal'] and collection['terminal_sources']==2
out=root/'starting-geometries-v1'
if out.exists():raise FileExistsError('Preserve the extracted coordinate source')
out.mkdir();(out/'extractor.py').write_bytes(Path(__file__).read_bytes())
specs={'PREP-007':dict(elements={'Cl':1,'F':3},name='chlorine trifluoride',multiplicity=1,state='1 A 1',point_group='C2v',key='JOHWNGGYGAVMGU-UHFFFAOYSA-N'),
       'PREP-045':dict(elements={'C':1,'H':3,'O':2},name='methylperoxy radical',multiplicity=2,state='2 A "',point_group='Cs',key='WTFNSXYULBQCQV-UHFFFAOYSA-N')}
records=[]
quoted=re.compile(r'"(?:\\.|[^"\\])*"')
for source in collection['records']:
    ident=source['candidate_id'];spec=specs[ident]
    assert source['status']=='geometry-page-retrieved-awaiting-review'
    raw_path=root/ident/'geometry.html';raw=raw_path.read_text(encoding='utf-8')
    assert sha(raw_path)==source['requests'][-1]['sha256']
    data=read(root/ident/'tables.json')
    table=next(t for t in data['tables'] if t['index']==1)
    assert spec['key'] in table['preceding_text'] and spec['state'] in table['preceding_text']
    assert 'geometry from B3LYP/6-31G*' in table['preceding_text']
    matches=re.findall(r'\bvar\s+xyzmodel\s*=\s*(.*?);',raw,re.DOTALL)
    assert len(matches)==1
    expression=matches[0]
    remainder=quoted.sub('',expression)
    assert re.fullmatch(r'[+\s]*',remainder),'Only literal string concatenation is accepted'
    payload=''.join(json.loads(token) for token in quoted.findall(expression))
    lines=payload.splitlines();n=int(lines[0]);assert n==sum(spec['elements'].values())
    assert [line.strip() for line in lines[n+2:] if line.strip()]==['END']
    atom_lines=lines[2:n+2];atoms=[]
    for line in atom_lines:
        fields=line.split();assert len(fields)==4 and re.fullmatch(r'[A-Z][a-z]?',fields[0])
        values=[float(x) for x in fields[1:]];assert all(math.isfinite(v) for v in values)
        atoms.append(dict(element=fields[0],literal_coordinates=fields[1:],xyz_angstrom=values))
    assert dict(Counter(a['element'] for a in atoms))==spec['elements']
    tabular=table['rows'][2:];assert len(tabular)==n
    maximum=Decimal(0)
    for index,(atom,row) in enumerate(zip(atoms,tabular),1):
        assert row[0]==atom['element']+str(index)
        for a,b in zip(atom['literal_coordinates'],row[1:4]):
            difference=abs(Decimal(a)-Decimal(b));quantum=Decimal(1).scaleb(Decimal(b).as_tuple().exponent)
            assert difference<=quantum/2,'Visualizer and displayed internal-coordinate table disagree beyond their printed precision'
            maximum=max(maximum,difference)
    directory=out/ident;directory.mkdir()
    xyz=directory/'starting.xyz'
    xyz.write_text('\n'.join([str(n),f'{ident}; literal NIST visualization coordinates in angstrom; B3LYP/6-31G* starting geometry',*atom_lines])+'\n',encoding='ascii',newline='\n')
    record=dict(candidate_id=ident,name=spec['name'],charge=0,multiplicity=spec['multiplicity'],state=spec['state'],
        source_reported_point_group=spec['point_group'],inchikey=spec['key'],element_counts=spec['elements'],atom_count=n,
        source_url=source['followed_link']['url'],session_identity_url=source['initial_url'],raw_source_sha256=sha(raw_path),
        coordinate_method='B3LYP/6-31G*',model_chemistry_page='G3B3',units='angstrom',
        source='Public NIST page literal xyzmodel payload, checked against the same page internal-coordinate table',
        extraction='Decoded JSON-compatible quoted string literals and concatenation only; no JavaScript execution or coordinate fitting.',
        actual_literal_coordinates=atoms,maximum_difference_from_four_decimal_table_angstrom=str(maximum),
        xyz_file=str(xyz.relative_to(out)),xyz_sha256=sha(xyz),status='literal-starting-geometry-prepared',
        minimum_frequency_and_stability_not_verified=True,formal_case=False,accepted_quality_reference=False,
        gaussian_started=False,cov_executed=False)
    atomic_json(directory/'source-and-state.json',record);records.append(record)
atomic_json(out/'summary.json',dict(created_epoch=time.time(),starting_geometries_prepared=2,records=records,
    formal_external_cases=0,accepted_quality_references=0,gaussian_started=False,cov_executed=False))
print(json.dumps(dict(output=str(out),starting_geometries=2,formal_external_cases=0)))
