"""Review preserved native export frames only after the whole collection ends.

Checks identity, members, raw values, selection, draw-call semantics and selected
PNG strokes. This is not a substitute for complete screenshot/chemistry review.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import sys
import time
import xml.etree.ElementTree as ET

from validation_process import atomic_json


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def load(path): return json.loads(Path(path).read_text(encoding='utf-8'))


def review_case(root,record,Image):
    terminal=load(root/'cases'/record['case_id']/'terminal.json')
    evidence=root/terminal['evidence_directory']
    manifest_file=root/terminal['evidence_manifest']
    assert sha(manifest_file)==terminal['evidence_manifest_sha256']
    artifacts={Path(k).as_posix():v for k,v in load(manifest_file).items()}
    consumed={}
    def consume(relative):
        relative=Path(relative).as_posix()
        path=evidence/relative
        identity=sha(path)
        assert identity==artifacts[relative]['sha256'],relative+' changed'
        consumed[relative]=identity
        return path
    native=evidence/'native'
    expected=load(consume('expected-evidence.json'))
    production=load(consume('production.json'))
    events=[json.loads(line) for line in consume('native/events.jsonl').read_text(encoding='utf-8').splitlines()]
    traces={x['frame']:x['data'] for x in events if x['kind']=='export.frame-trace'}
    exports=[x for x in events if x['kind']=='export.actual']
    failures=[]; bundles=[]
    def check(ok,code,detail):
        if not ok: failures.append({'check':code,'detail':detail})
    names=expected.get('export_names',[])
    check(bool(names),'VIEW-COVERAGE','No frozen export-state coverage specified')
    check(len(names)==len(set(names))==len(exports),'VIEW-EXPORT-COUNT',
          {'expected':len(names),'measured':len(exports)})
    for name in names:
        begin=len(failures)
        matches=[x for x in exports if Path(x['data']['base']).name==name]
        check(len(matches)==1,'VIEW-EXPORT-EVENT',name)
        if len(matches)!=1: continue
        event=matches[0]; capture_trace=traces.get(event['frame'],[])
        snapshot_records=[x['data'] for x in capture_trace if x['kind']=='diagram.snapshot']
        check(len(snapshot_records)==1,'VIEW-FRAME-SNAPSHOT',name)
        if len(snapshot_records)!=1: continue
        drawn=snapshot_records[0]
        paths={ext:consume('native/'+name+'.mo.'+ext) for ext in ('json','csv','svg','png')}
        data=load(paths['json']); context=data.get('view_snapshot') or {}
        csv_rows=list(csv.DictReader(paths['csv'].open(encoding='utf-8',newline='')))
        svg=ET.parse(paths['svg']).getroot()
        with Image.open(paths['png']) as image:
            png_id=image.info.get('cov.view.snapshot')
            image.load(); pixels=image.convert('RGB')
        check(event['data']['success'] and context.get('origin')=='interactive-canvas' and
              drawn['id']==event['data'].get('snapshot_id')==context.get('id')==
              svg.attrib.get('data-view-snapshot')==png_id,'VIEW-IDENTITY',name)
        selected=context.get('inspected_orbital_index')
        check(selected==drawn['inspected_orbital_index']==event['data'].get('selected_index'),
              'VIEW-INSPECTION',name)
        metadata={m['index']:m for m in data['orbitals']}
        csv_metadata={int(m['index']):m for m in csv_rows}
        check(len(metadata)==len(data['orbitals']) and set(metadata)==set(csv_metadata),
              'VIEW-RAW-METADATA-COVERAGE',name)
        for i,mo in metadata.items():
            source=production['orbitals'][i]; item=csv_metadata.get(i,{})
            check(abs(mo['energy_hartree']-source['energy_hartree'])<=1e-12 and
                  abs(mo['occupation']-source['occupation'])<=1e-7,
                  'VIEW-RAW-VALUES',{'bundle':name,'mo':i})
            check(mo['selected']==(i==selected) and item.get('selected')==str(int(i==selected)) and
                  item.get('view_snapshot_id')==context.get('id') and
                  item.get('energy_unit')==data['energy_unit'] and
                  math.isclose(float(item.get('energy_hartree','nan')),mo['energy_hartree'],abs_tol=1e-12),
                  'VIEW-CSV-JSON-SELECTION',{'bundle':name,'mo':i})
            check((mo['spin'].lower().startswith('b'))==(source['spin']==1),
                  'VIEW-RAW-SPIN',{'bundle':name,'mo':i})
        rows=data['diagram_rows']
        svg_rows=[node for node in svg if node.attrib.get('class')=='mo-level']
        check(len(rows)==data['diagram_row_count']==drawn['row_count']==len(svg_rows),
              'VIEW-ROW-COUNT',name)
        native_members=[x['data'] for x in capture_trace if x['kind']=='draw.level']
        expected_native=[]
        selected_png=[]
        for row_index,row in enumerate(rows):
            members=row.get('members',[])
            check(bool(members),'VIEW-MEMBERS-MISSING',{'bundle':name,'row':row_index})
            marks=[] if row_index>=len(svg_rows) else list(svg_rows[row_index].iter())
            marks=[m for m in marks if m.attrib.get('class')=='mo-member']
            check(len(marks)==len(members),'VIEW-SVG-MEMBER-COUNT',{'bundle':name,'row':row_index})
            energies=[]
            for member_index,member in enumerate(members):
                indices=[member['orbital_index']]+([member['spin_counterpart']]
                    if member['spin_counterpart'] is not None else [])
                expected_selected=selected in indices
                expected_used=selected if expected_selected else indices[0]
                check(member['selected']==expected_selected and member['inspected_orbital_index']==expected_used,
                      'VIEW-MEMBER-SELECTION',{'bundle':name,'row':row_index,'member':member_index})
                for i in indices:
                    check(i in metadata and metadata[i]['diagram_row_index']==row_index and
                          metadata[i]['included_in_diagram'],'VIEW-MEMBER-RAW-LINK',{'bundle':name,'mo':i})
                    energies.append(production['orbitals'][i]['energy_hartree'])
                expected_native.append((expected_used,row['representative_metadata_index'],
                    row['layout_energy_hartree'],member['alpha_arrows'],member['beta_arrows']))
                if member_index>=len(marks): continue
                attrs=marks[member_index].attrib
                check(int(attrs['data-orbital-index'])==indices[0] and
                      int(attrs['data-inspected-orbital-index'])==expected_used and
                      (attrs['data-selected']=='true')==expected_selected,
                      'VIEW-SVG-SELECTION',{'bundle':name,'row':row_index,'member':member_index})
                if expected_selected:
                    # Inspect actual decoded PNG pixels at the selected stroke,
                    # away from the two electron-arrow shafts.
                    x0,x1,y=float(attrs['x1']),float(attrs['x2']),float(attrs['y1'])
                    probes=[(round(x),round(y)+dy) for x in (x0+2,x1-2) for dy in (-1,0,1)]
                    visible=any(0<=x<pixels.width and 0<=yy<pixels.height and
                                pixels.getpixel((x,yy))==(40,109,224) for x,yy in probes)
                    check(visible,'VIEW-PNG-SELECTED-STROKE',{'bundle':name,'row':row_index,'probes':probes})
                    selected_png.append(visible)
            check(bool(energies) and math.isclose(min(energies),row.get('all_members_energy_min_hartree',math.nan),abs_tol=1e-12)
                  and math.isclose(max(energies),row.get('all_members_energy_max_hartree',math.nan),abs_tol=1e-12),
                  'VIEW-ROW-ENERGY-RANGE',{'bundle':name,'row':row_index})
        check(len(expected_native)==len(native_members),'VIEW-FRAME-MEMBER-COUNT',name)
        for wanted,actual in zip(expected_native,native_members):
            check(wanted[0]==actual['used_internal_mo'] and wanted[1]==actual['level_metadata_internal_mo'] and
                  math.isclose(wanted[2],actual['layout_energy_hartree'],abs_tol=1e-12) and
                  wanted[3]==actual['alpha_arrows'] and wanted[4]==actual['beta_arrows'],
                  'VIEW-DRAW-CALL-CONSISTENCY',{'bundle':name,'wanted':wanted,'actual':actual})
        bundles.append({'name':name,'status':'pass' if len(failures)==begin else 'fail',
            'snapshot_id':context.get('id'),'export_frame':event['frame'],'selected_mo':selected,
            'rows':len(rows),'members':len(expected_native),'selected_png_checks':len(selected_png)})
    return {'case_id':record['case_id'],'status':'view_subset_pass' if not failures else 'view_subset_fail',
        'bundles':bundles,'failures':failures,'consumed_artifacts':consumed,'formal_case_pass':False}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--round',type=Path,required=True)
    parser.add_argument('--review-id',default='view-review-v1')
    args=parser.parse_args(); root=args.round.resolve()
    if Path(args.review_id).name!=args.review_id: raise ValueError('review-id must be a directory name')
    manifest=load(root/'manifest.json'); barrier=load(root/'collection-complete.json')
    assert barrier['all_terminal'] and barrier['terminal_cases']==manifest['case_count']
    assert barrier['round_identity']==manifest['round_identity']
    out=root/args.review_id
    if out.exists(): raise FileExistsError('Existing review evidence is never overwritten')
    for relative,identity in manifest['image_dependency_artifacts'].items():
        assert sha(root/'image-deps'/relative)==identity
    sys.path.insert(0,str(root/'image-deps'))
    from PIL import Image
    out.mkdir(); start=time.perf_counter(); results=[]
    for record in manifest['cases']:
        try: result=review_case(root,record,Image)
        except Exception as error:
            result={'case_id':record['case_id'],'status':'view_subset_error',
                    'error':repr(error),'formal_case_pass':False}
        atomic_json(out/(record['case_id']+'.json'),result); results.append(result)
    result={'round_identity':manifest['round_identity'],'case_count':len(results),
        'all_reviewed':len(results)==manifest['case_count'],
        'status_counts':{s:sum(x['status']==s for x in results) for s in sorted({x['status'] for x in results})},
        'bundle_count':sum(len(x.get('bundles',[])) for x in results),
        'failure_count':sum(len(x.get('failures',[])) for x in results),
        'wall_seconds':time.perf_counter()-start,'reviewer_sha256':sha(Path(__file__)),
        'collection_barrier_sha256':sha(root/'collection-complete.json'),'formal_case_passes':0,
        'limitations':['Draw-call semantics do not prove all visible text/glyph positions or clipping.',
                       'PNG stroke probes do not certify all output pixels.',
                       'Chemical/symmetry interpretation and paused physical references remain separate.'],
        'case_reports':{p.name:sha(p) for p in out.glob('OLD-*.json')}}
    atomic_json(out/'summary.json',result)
    print(json.dumps(result),flush=True)
    return 0 if result['status_counts']=={'view_subset_pass':manifest['case_count']} else 2


if __name__=='__main__': raise SystemExit(main())
