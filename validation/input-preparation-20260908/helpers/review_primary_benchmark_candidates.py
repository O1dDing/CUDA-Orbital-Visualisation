"""Compare archived author geometries with local copies; preserve charge/spin sources."""
from pathlib import Path
import hashlib
import json
import math
import re

root=Path(__file__).resolve().parent.parent/'outputs/cov-complete-validation-20260906/benchmark-geometry-preparation-v1'
collection=json.loads((root/'primary-collection.json').read_text(encoding='utf-8'))
assert collection['terminal_count']==collection['selected_count']==4
out=root/'review-v2'
if out.exists():raise FileExistsError('Retain reviewed source evidence')
out.mkdir()
bohr_to_angstrom=0.529177210903
records=[]
for source in collection['records']:
    record=dict(candidate_id=source['candidate_id'],query=source['query'],primary_commit=source['commit'],
                formal_case=False,accepted_reference_geometry=False)
    try:
        assert source['status']=='primary-records-collected'
        folder=root/'primary-records'/source['candidate_id']
        for item in source['files']:
            assert hashlib.sha256((folder/item['name']).read_bytes()).hexdigest()==item['sha256']
        rows=(folder/'coord').read_text().splitlines();atoms=[];bohr=[];in_coords=False;state=[]
        for row in rows:
            if row.startswith('$coord'):
                if row.strip()!='$coord':raise ValueError('Explicitly handle any alternate coordinate units')
                in_coords=True;continue
            if row.startswith('$'):
                in_coords=False
                if row.startswith('$eht'):state.append(row)
                continue
            if in_coords and row.strip():
                x,y,z,symbol=row.split();atoms.append(symbol.capitalize());bohr.append([float(v) for v in (x,y,z)])
        xyz=[[v*bohr_to_angstrom for v in p] for p in bohr]
        local=Path(source['local_candidate']['file']);assert hashlib.sha256(local.read_bytes()).hexdigest()==source['local_candidate']['sha256']
        local_rows=local.read_text().splitlines();n=int(local_rows[0]);charge,mult=map(int,local_rows[1].split())
        local_atoms=[];local_xyz=[]
        for row in local_rows[2:2+n]:
            symbol,x,y,z=row.split();local_atoms.append(symbol);local_xyz.append(list(map(float,[x,y,z])))
        if atoms!=local_atoms:raise ValueError('Local atom order/elements differ from primary geometry')
        difference=max(abs(a-b) for p,q in zip(xyz,local_xyz) for a,b in zip(p,q))
        distance_difference=max(abs(math.dist(xyz[i],xyz[j])-math.dist(local_xyz[i],local_xyz[j])) for i in range(n) for j in range(i))
        coordinate_extent=max(abs(v) for p in bohr for v in p)
        state_records=[]
        for key in ['.CHRG','.UHF']:
            if (folder/key).exists():state_records.append(dict(field=key,value=int((folder/key).read_text().strip()),source_file=key))
        for line in state:
            state_records.append(dict(field='$eht',raw_record=line))
        known_charge=next((r['value'] for r in state_records if r['field']=='.CHRG'),None)
        known_unpaired=next((r['value'] for r in state_records if r['field']=='.UHF'),None)
        if known_charge is not None and known_charge!=charge:raise ValueError('Local charge differs from explicit author state')
        if known_unpaired is not None and known_unpaired+1!=mult:raise ValueError('Local multiplicity differs from explicit author unpaired-electron setting')
        for line in state:
            m=re.search(r'charge=(-?\d+)\s+unpaired=(\d+)',line)
            if m and (int(m[1])!=charge or int(m[2])+1!=mult):raise ValueError('Primary state records disagree')
        prepared=dict(candidate_id=source['candidate_id'],symbols=atoms,coordinates_bohr=bohr,
                      coordinates_angstrom=xyz,bohr_to_angstrom=bohr_to_angstrom,
                      coordinate_source=f"https://github.com/grimme-lab/GMTKN55/blob/{source['commit']}/{source['source_directory']}/coord",
                      coordinate_provenance='Author-distributed computational benchmark geometry; not an experimental coordinate set',
                      electronic_state_records=state_records,
                      local_proposed_charge=charge,local_proposed_multiplicity=mult,
                      electronic_state_provenance='explicit-author-records' if state_records else 'source-state-fields-absent; proposed state requires independent confirmation',
                      formal_case=False,gaussian_started=False,accepted_reference_geometry=False)
        target=out/(source['candidate_id']+'-starting-geometry.json')
        target.write_text(json.dumps(prepared,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
        record.update(status='author-coordinate-identity-established',atom_count=n,
                      maximum_local_coordinate_difference_angstrom=difference,
                      primary_maximum_coordinate_bohr=coordinate_extent,
                      maximum_corresponding_pair_distance_difference_angstrom=distance_difference,
                      comparison_scope='Same atom ordering, coordinate frames may differ. Pair distances are rotation/translation invariant; no fitted scaling or threshold adjustment.',
                      electronic_state_source=prepared['electronic_state_provenance'],
                      source_geometry=prepared['coordinate_source'],
                      normalized_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    except Exception as error:
        record.update(status='source-review-failed-retained',error=f'{type(error).__name__}: {error}')
    records.append(record)
summary=dict(all_4_records_reviewed=True,records=records,formal_external_cases=0,gaussian_started=False,
             prior_review_correction='v1 raw coordinate difference did not establish a common coordinate frame. Retain its raw values; v2 adds corresponding pair-distance differences without declaring a scientific reference pass.',
             accepted_quality_references=0,attribution='Goerigk, Hansen, Bauer, Ehrlich, Najibi, Grimme and coauthors, PCCP 2017, 19, 32184-32215, DOI 10.1039/C7CP04913G; dataset CC-BY-4.0',
             source_repository_commit=collection['commit'])
(out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summary,ensure_ascii=False),flush=True)
