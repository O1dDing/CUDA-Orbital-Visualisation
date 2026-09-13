"""Review a completed conformer preparation without running COV or Gaussian."""
from pathlib import Path
import collections
import hashlib
import json
import math

workspace = Path(__file__).resolve().parent.parent
root = workspace/'outputs/cov-complete-validation-20260906/external-coordinate-preparation-v2'
summary = json.loads((root/'summary.json').read_text(encoding='utf-8'))
assert summary['all_terminal'] and summary['terminal_candidates']==summary['expected_candidates']==67
out = root/'review-v1'
if out.exists(): raise FileExistsError('Retain the completed source review')
out.mkdir()
records=[]
for row in summary['records']:
    record={key:row.get(key) for key in ['candidate_id','query','family','status','cid','http_status','error']}
    record.update(input_geometry_accepted=False, formal_case=False)
    for file_key,sha_key in [('raw_file','raw_sha256'),('normalized_file','normalized_sha256'),('error_file','error_sha256')]:
        if row.get(file_key):
            assert hashlib.sha256(Path(row[file_key]).read_bytes()).hexdigest()==row[sha_key]
    if row['status']=='computed-starting-conformer-retrieved':
        raw=json.loads(Path(row['raw_file']).read_text())['PC_Compounds'][0]
        normal=json.loads(Path(row['normalized_file']).read_text())
        aids=raw['atoms']['aid']; aid_to_pos={aid:i for i,aid in enumerate(aids)}
        edges=raw.get('bonds',{})
        adjacency={aid:set() for aid in aids}; lengths=[]
        conf=normal['conformers'][0]; xyz=conf['coordinates_angstrom']
        for a,b,order in zip(edges.get('aid1',[]),edges.get('aid2',[]),edges.get('order',[])):
            assert a in adjacency and b in adjacency and a!=b
            adjacency[a].add(b);adjacency[b].add(a)
            lengths.append(dict(aid1=a,aid2=b,order=order,distance_angstrom=math.dist(xyz[aid_to_pos[a]],xyz[aid_to_pos[b]])))
        unseen=set(aids);components=[]
        while unseen:
            todo=[next(iter(unseen))]; reached=set()
            while todo:
                aid=todo.pop()
                if aid in reached:continue
                reached.add(aid);todo.extend(adjacency[aid]-reached)
            unseen-=reached;components.append(sorted(reached))
        properties=[]
        for prop in raw.get('props',[]):
            urn=prop.get('urn',{})
            if urn.get('label') in ['Molecular Formula','Charge','InChIKey','SMILES','Conformer','Conformer Model','Conformer RMSD']:
                properties.append(prop)
        record.update(atom_count=len(aids),conformer_count=len(normal['conformers']),
                      formula_charge_identity_checked=True,raw_and_normalized_sha256_verified=True,
                      graph_components=components,component_count=len(components),
                      stored_bond_lengths=lengths, pubchem_properties=properties,
                      original_cases_with_same_element_counts=row['original_cases_with_same_element_counts'],
                      further_checks=['Confirm molecular connectivity and requested stereochemistry',
                                      'Confirm electronic state and environment',
                                      'Gaussian optimization, frequency and applicable stability/reference checks',
                                      'Formal external identity freeze and full COV evidence only after original loop acceptance'])
    records.append(record)
result=dict(source_summary_sha256=hashlib.sha256((root/'summary.json').read_bytes()).hexdigest(),
            all_67_preparation_records_reviewed=True,
            retrieved_conformers=sum(r['status']=='computed-starting-conformer-retrieved' for r in records),
            retrieved_single_component=sum(r.get('component_count')==1 for r in records),
            retrieved_formula_disjoint_from_original=sum(r['status']=='computed-starting-conformer-retrieved' and not r.get('original_cases_with_same_element_counts') for r in records),
            unavailable_retrieval_http_counts=dict(collections.Counter(str(r['http_status']) for r in records if r['status']=='conformer-unavailable-or-rejected-retained')),
            raw_identity_integrity_failures=0, formal_external_cases=0,accepted_reference_geometries=0,
            gaussian_started=False,records=records)
(out/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k!='records'},ensure_ascii=False))
print(json.dumps([dict(candidate_id=r['candidate_id'],query=r['query'],status=r['status'],http_status=r.get('http_status'),error=r.get('error')) for r in records if r['status']!='computed-starting-conformer-retrieved'],ensure_ascii=False))
