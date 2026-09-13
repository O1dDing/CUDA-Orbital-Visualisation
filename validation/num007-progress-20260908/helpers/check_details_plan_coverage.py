"""Compare new controls against every preserved full-round plan before execution."""
from pathlib import Path
import ctypes
import hashlib
import importlib.util
import json
import os
import sys
import time
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
repo=Path(r'F:\Dev\cov-native-validation-20260905')
sys.path[:0]=[str(repo/'tests'),r'F:\Dev\cov-validation-20260905\reference-deps']
from validation_process import atomic_json,physical_core_masks
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
base=Path(r'F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906')
old=base/'NUM-FIX-006-original'
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result
previous=module('old_plan',old/'runner-source/native_collection_batch.py')
current=module('new_plan',repo/'tests/native_collection_batch.py')
from review_export_snapshot_round import review_details_frame
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
records=[]
for case in read(old/'manifest.json')['cases']:
    terminal=read(old/'cases'/case['case_id']/'terminal.json')
    data=read(old/terminal['evidence_directory']/'production.json')
    before,a=previous.make_plan(data);after,b=current.make_plan(data)
    assert a['mos']==b['mos'] and a['primary_members']==b['primary_members'] and a['spin_counterparts']==b['spin_counterparts']
    assert a['export_names']==b['export_names'] and set(a['captures']).issubset(b['captures'])
    assert [l for l in before.splitlines() if l.startswith('volume ')]==[l for l in after.splitlines() if l.startswith('volume ')]
    assert b['primary_members'] and b['orbital_details_member']==b['primary_members'][0]
    assert len(b['captures'])==len(a['captures'])+4
    records.append({'case_id':case['case_id'],'mos':b['mos'],'primary_members':len(b['primary_members']),
        'exports':len(b['export_names']),'retained_captures':len(a['captures']),'added_details_captures':4})
state=dict(scene_matches_applied=True,rendered_mo=4,applied_mo=4,drawn_ui_mo=4,requested_mo=4)
target=dict(id='diagram.details.close',visible=True,rect=[20,20,70,50],clip_rect=[0,0,150,150])
fixture={'state':state,'targets':[target]}
assert not review_details_frame(fixture,4,target['id'],(100,100))
assert review_details_frame(fixture,3,target['id'],(100,100))==['selected-member-mismatch']
assert review_details_frame(fixture,4,None,(100,100))==['window-remains-open']
assert not review_details_frame({'state':state,'targets':[]},4,None,(100,100))
for altered,expected in [({**target,'visible':False},'required-control-not-visible'),
                         ({**target,'rect':[20,20,120,50]},'required-control-outside-frame'),
                         ({**target,'rect':[20,20,float('nan'),50]},'nonfinite-control-bounds'),
                         ({**target,'clip_rect':[0,0,150,45]},'required-control-partially-clipped')]:
    assert review_details_frame({'state':state,'targets':[altered]},4,target['id'],(100,100))==[expected]
out=base/('details-plan-preflight-'+str(time.time_ns()));out.mkdir()
atomic_json(out/'summary.json',{'all_original_cases_checked':len(records)==273,'case_count':len(records),
    'retained_all_mos':sum(r['mos'] for r in records),'retained_primary_targets':sum(r['primary_members'] for r in records),
    'retained_export_bundles':sum(r['exports'] for r in records),'new_details_captures':4*len(records),
    'control_checker_positive_and_negative_checks':8,'formal_case_passes':0,
    'source_hashes':{n:hashlib.sha256((repo/'tests'/n).read_bytes()).hexdigest() for n in ('native_collection_batch.py','review_export_snapshot_round.py')},
    'records':records})
print(json.dumps({'evidence':str(out),'cases':len(records),'retained_primary_targets':sum(r['primary_members'] for r in records),
    'new_details_captures':4*len(records),'checker_checks':8}),flush=True)
