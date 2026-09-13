"""Finish the authorized post-collection stages without restarting Gaussian.

Waits for the actual collection process tree to finish before launching any
review. Each child reviewer owns the existing aggregate Windows Job limits.
"""
from pathlib import Path
import argparse
import hashlib
import json
import msvcrt
import os
import subprocess
import sys
import time

parser=argparse.ArgumentParser()
parser.add_argument('--round',type=Path,required=True)
args=parser.parse_args(); root=args.round.resolve()
work=Path(__file__).resolve().parent
def read(p): return json.loads(p.read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def atomic(p,d):
    q=p.with_name(p.name+'.tmp');q.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');q.replace(p)
lock=(root/'post-collection-pipeline.lock').open('a+b');lock.seek(0)
if lock.read(1)==b'':lock.write(b'0');lock.flush()
lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
manifest=read(root/'manifest.json'); identity=manifest['round_identity']
destination=root/('post-collection-pipeline-'+str(time.time_ns()));destination.mkdir()
steps=[]; start=time.time()
atomic(destination/'identity.json',{'round_identity':identity,'created_epoch':start,
    'pipeline_sha256':sha(Path(__file__)),'gaussian_started':False})
def status(stage):
    atomic(root/'pipeline-progress.json',{'round_identity':identity,'stage':stage,
        'started_epoch':start,'updated_epoch':time.time(),'steps':steps,'formal_case_passes':0,
        'REF-001':'paused-by-user'})
def run(name,argv):
    status(name); log=destination/(name+'.log'); before=time.time()
    with log.open('wb') as output:
        done=subprocess.run([sys.executable,*map(str,argv)],stdout=output,stderr=subprocess.STDOUT)
    steps.append({'stage':name,'exit_code':done.returncode,'wall_seconds':time.time()-before,
                  'log':str(log),'log_sha256':sha(log)})
    status(name+'-terminal')
    return done.returncode
status('waiting-for-complete-collection-process-tree')
while True:
    barrier=root/'collection-complete.json'
    supervisors=sorted(root.glob('supervisor-*'))
    receipt=supervisors[-1]/'process.json' if supervisors else None
    if barrier.exists() and receipt and receipt.exists():
        data=read(barrier)
        if data.get('round_identity')!=identity: raise RuntimeError('Collection identity changed')
        if data.get('all_terminal') and data.get('terminal_cases')==manifest['case_count']: break
    time.sleep(5)
run('numerical-review',[work/'run_numeric_review_batch.py','--round',root])
if (root/'review-v1/review-complete.json').exists():
    if not (root/'review-v1/unified-review').exists():
        run('numerical-summary',[work/'summarize_numeric_review.py','--round',root,'--review-id','review-v1'])
    if (root/'review-v1/unified-review').exists() and not (root/'review-v1/unified-review/extrema.json').exists():
        run('numerical-extrema',[work/'numeric_round_extrema.py','--round',root])
view=root/'view-review-v1/summary.json'
if not view.exists():
    review_index=1
    while (root/f'view-review-v{review_index}').exists(): review_index+=1
    name=f'view-review-v{review_index}'
    run('view-snapshot-review',[work/'run_view_snapshot_review.py','--round',root,'--review-id',name])
    view=root/name/'summary.json'
result={'round_identity':identity,'wall_seconds':time.time()-start,'steps':steps,
    'numerical_summary':str(root/'review-v1/unified-review/summary.json'),
    'view_summary':str(view),'view_summary_exists':view.exists(),
    'REF-001':'paused-by-user','formal_case_passes':0,
    'next':'Unified root-cause analysis; no algorithm or input has been changed by this pipeline.'}
atomic(destination/'pipeline-complete.json',result);status('reviews-terminal-awaiting-root-cause-analysis')
print(json.dumps(result),flush=True)
