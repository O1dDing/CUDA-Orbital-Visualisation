"""Run the prepared real-production baseline only after NUM006 reviews terminate."""
from pathlib import Path
import ctypes
import hashlib
import json
import msvcrt
import os
import shutil
import sys
import time

workspace=Path(r'F:\Codex\2026-09-05\branch-15')
controls=workspace/'outputs/cov-complete-validation-20260906/regression-fixtures/density-source-v1'
baseline=controls/'baseline-production-21036b2-v3'
round_root=workspace/'outputs/general-fixes-20260906/NUM-FIX-006-original'
root=baseline/'post-num006-pipeline'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sys.path.insert(0,str(baseline/'source/helpers'))
from validation_process import atomic_json, physical_core_masks, run_tree
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])

if len(sys.argv)==2 and sys.argv[1]=='prepare':
    if root.exists():raise FileExistsError('Retain this frozen queue instead of replacing it')
    (root/'source').mkdir(parents=True)
    for name in ['queue_density_source_baseline.py','review_density_source_controls.py']:
        shutil.copy2(workspace/'work'/name,root/'source'/name)
    tested=workspace/'outputs/cov-complete-validation-20260906/regression-fixtures/density-reviewer-preflight-1788847493973793400'
    assert read(tested/'process.json')['exit_code']==0
    assert sha(root/'source/review_density_source_controls.py')==sha(tested/'source/review_density_source_controls.py')
    atomic_json(root/'identity.json',dict(created_epoch=time.time(),baseline_identity_sha256=sha(baseline/'identity.json'),
        source_hashes={p.name:sha(p) for p in (root/'source').iterdir()},
        preceding_round_identity=read(round_root/'manifest.json')['round_identity'],
        tested_reviewer=str(tested),resource_policy=dict(own_cores=4,own_memory_gib=64,core_slot_offset=8,
        independent_gaussian_reserved_cores=8,independent_gaussian_reserved_memory_gib=48),
        gaussian_controlled=False,algorithm_changed=False,formal_cases_added=0))
    print(json.dumps(dict(status='queue-prepared-not-started',root=str(root))),flush=True)
    raise SystemExit(0)

identity=read(root/'identity.json')
assert sha(Path(__file__))==identity['source_hashes']['queue_density_source_baseline.py']
assert sha(baseline/'identity.json')==identity['baseline_identity_sha256']
for name,digest in identity['source_hashes'].items():assert sha(root/'source'/name)==digest
lock=(root/'queue.lock').open('a+b');lock.seek(0)
if lock.read(1)==b'':lock.write(b'0');lock.flush()
lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
if (root/'started.json').exists():raise FileExistsError('Do not repeat this started baseline queue')
started=time.time();atomic_json(root/'started.json',dict(started_epoch=started,pid=os.getpid()))
atomic_json(root/'progress.json',dict(stage='waiting-for-NUM006-full-collection-and-reviews',updated_epoch=time.time(),gaussian_controlled=False))
while True:
    barrier=round_root/'collection-complete.json'
    if barrier.exists():
        collection=read(barrier);pipeline=read(round_root/'pipeline-progress.json')
        assert collection['round_identity']==identity['preceding_round_identity']
        if collection.get('all_terminal') and collection.get('terminal_cases')==273 and pipeline.get('stage')=='reviews-terminal-awaiting-root-cause-analysis':break
    time.sleep(5)
mask=sum(physical_core_masks(12)[8:12]);env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
steps=[]
for name,command,timeout in [
    ('production-baseline',[sys.executable,'-X','utf8',baseline/'source/run_density_source_baseline.py','run'],3600),
    ('independent-control-review',[sys.executable,'-X','utf8',root/'source/review_density_source_controls.py',
      '--baseline',baseline,'--manifest',controls/'manifest.json','--output',baseline/'review-v1'],180)]:
    atomic_json(root/'progress.json',dict(stage=name,updated_epoch=time.time(),steps=steps,gaussian_controlled=False))
    directory=root/name;directory.mkdir()
    result=run_tree(command,directory,env,mask,64,timeout)
    atomic_json(directory/'process.json',result)
    steps.append(dict(stage=name,**result))
    if result['exit_code'] or result['timed_out']:
        atomic_json(root/'progress.json',dict(stage='retained-execution-error-needs-review',updated_epoch=time.time(),steps=steps))
        raise SystemExit(result['exit_code'] or 1)
atomic_json(root/'progress.json',dict(stage='all-density-controls-reviewed-awaiting-common-root-fix',updated_epoch=time.time(),steps=steps,
    formal_case_passes=0,formal_cases_added=0,gaussian_controlled=False))
print(json.dumps(dict(status='density-baseline-and-review-terminal',baseline=str(baseline),steps=len(steps))),flush=True)
