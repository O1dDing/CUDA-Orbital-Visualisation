"""Aggregate terminal numerical and drawn-view evidence, retaining all failures."""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json
import time

def read(path): return json.loads(path.read_text(encoding='utf-8'))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

parser=argparse.ArgumentParser()
parser.add_argument('--round',type=Path,required=True)
parser.add_argument('--view-review-id',default='view-review-v1')
args=parser.parse_args()
root=args.round.resolve()
manifest=read(root/'manifest.json')
collection=read(root/'collection-complete.json')
numeric=read(root/'review-v1/unified-review/summary.json')
view_root=root/args.view_review_id
view=read(view_root/'summary.json')
identity=manifest['round_identity']
expected={c['case_id'] for c in manifest['cases']}
assert len(expected)==manifest['case_count']
assert collection['all_terminal'] and view['all_reviewed']
assert all(d['round_identity']==identity for d in (collection,numeric,view))
assert {c['case_id'] for c in collection['cases']}==expected
assert {c['case_id'] for c in numeric['cases']}==expected
assert {Path(p).stem for p in view['case_reports']}==expected
assert view['collection_barrier_sha256']==sha(root/'collection-complete.json')
destination=root/'unified-view-numeric-review'
if destination.exists(): raise FileExistsError('Preserve earlier unified reports')
failures=[]; reports=[]; bundles=Counter(); probes=0
for name,digest in view['case_reports'].items():
    path=view_root/name
    assert sha(path)==digest, str(path)
    result=read(path)
    assert result['case_id']==path.stem
    reports.append({'case_id':result['case_id'],'status':result['status'],
                    'report':str(path),'sha256':digest})
    if result['status']=='view_subset_error':
        failures.append({'case_id':result['case_id'],'check':'VIEW-REVIEW-EXECUTION',
                         'detail':result.get('error'),'evidence':str(path)})
    for f in result.get('failures',[]):
        failures.append({'case_id':result['case_id'],**f,'evidence':str(path)})
    for bundle in result.get('bundles',[]):
        bundles[bundle['status']]+=1
        probes+=bundle.get('selected_png_checks',0)
result={'round_identity':identity,'production_commit':manifest['git_commit'],
    'created_epoch':time.time(),'summary_script_sha256':sha(Path(__file__)),
    'collection_counts':collection['collection_counts'],
    'numeric_subset_counts':numeric['numeric_subset_counts'],
    'view_subset_counts':view['status_counts'],'bundle_counts':dict(bundles),
    'selected_png_strokes_checked':probes,'case_count':len(expected),
    'view_failures':failures,'view_failures_by_check':dict(Counter(f['check'] for f in failures)),
    'view_case_reports':reports,'numeric_failures':numeric['failures'],
    'numeric_texture_failures':numeric['texture_failures'],
    'sampled_mos':numeric['sampled_mos'],
    'complete_frontier_textures':numeric['complete_frontier_textures'],
    'timing':{**numeric['timing'],'view_review_wall_seconds':view['wall_seconds'],
        'OLD-015':next(c for c in numeric['cases'] if c['case_id']=='OLD-015')},
    'consumed_summaries':{str(p):sha(p) for p in (root/'collection-complete.json',
        root/'review-v1/unified-review/summary.json',view_root/'summary.json')},
    'formal_case_passes':0,'external_molecule_count':0,'REF-001':'paused-by-user',
    'limitations':view['limitations']+numeric['remaining']}
destination.mkdir()
(destination/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:result[k] for k in ('case_count','numeric_subset_counts','view_subset_counts',
    'bundle_counts','selected_png_strokes_checked','view_failures_by_check','timing')},ensure_ascii=False))
