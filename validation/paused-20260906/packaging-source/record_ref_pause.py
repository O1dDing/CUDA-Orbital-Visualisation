"""Record an explicit user pause without changing the frozen REF runner/state."""
from pathlib import Path
from collections import Counter
import ctypes as ct
from ctypes import wintypes as wt
import hashlib
import json
import os
import sys
import time

ROOT = Path(r'F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906')
JOBS = Path(r'F:\Dev\cov-cycle-20260906\jobs')
PAUSE = ROOT / 'runtime-incidents/ref-user-pause-20260906'
sys.path.insert(0, r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json

before = json.loads((PAUSE / 'before-stop.json').read_text())
stopped = json.loads((PAUSE / 'coordinator-stop.json').read_text())
assert stopped['terminate_ok'] and stopped['wait_after_1000ms'] == 0
k = ct.WinDLL('kernel32', use_last_error=True)
k.OpenProcess.argtypes = [wt.DWORD, wt.BOOL, wt.DWORD]
k.OpenProcess.restype = wt.HANDLE
k.CloseHandle.argtypes = [wt.HANDLE]
k.WaitForSingleObject.argtypes = [wt.HANDLE, wt.DWORD]
processes = []
for pid in [stopped['pid'], 39504] + [s['state']['attempts'][-1]['pid'] for s in before['running_stages']]:
    handle = k.OpenProcess(0x100000 | 0x1000, False, pid)
    if handle:
        wait = k.WaitForSingleObject(handle, 0)
        k.CloseHandle(handle)
        assert wait == 0, (pid, wait)
        processes.append({'pid': pid, 'process_signaled': True})
    else:
        error = ct.get_last_error()
        assert error == 87, (pid, error)
        processes.append({'pid': pid, 'absent_error': error})
# Verify the actual kernel reservation was released, not merely a stale PID.
import msvcrt
with (JOBS / 'supervisor.lock').open('r+b') as lock:
    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
    msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)

manifest = json.loads((ROOT / 'references/REF-001/reference-candidates.json').read_text())
records = []
for reference in manifest['candidates']:
    case_id = reference['case_id']
    directory = JOBS / 'REF-001' / case_id
    raw_result_path = directory / 'result.json'
    raw_result = json.loads(raw_result_path.read_text()) if raw_result_path.exists() else None
    if raw_result is None:
        state = 'not_started'
    elif raw_result['status'] == 'candidate_collected':
        state = 'candidate_collected'
    elif raw_result['status'] == 'running':
        state = 'paused_interrupted'
    else:
        state = raw_result['status']
    stages = []
    for path in sorted(directory.glob('*/stage.json')):
        raw = json.loads(path.read_text())
        stages.append({'stage': path.parent.name, 'last_runner_status': raw['status'],
                       'current_execution_status': 'paused_interrupted' if raw['status'] == 'running' else raw['status'],
                       'attempts': len(raw['attempts']), 'raw_state_sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    records.append({'case_id': case_id, 'reference_identity': reference['identity'],
                    'state': state, 'stages': stages, 'scientific_pass': False})
counts = dict(Counter(r['state'] for r in records))
assert len(records) == 273 and counts == {'candidate_collected': 13, 'paused_interrupted': 2, 'not_started': 258}, counts
value = {'schema_version': 2, 'stage': 'paused_by_user_for_incremental_work',
         'updated_epoch': time.time(), 'runner_identity': before['progress']['runner_identity'],
         'reference_set_identity': manifest['reference_set_identity'],
         'coordinator_pid': None, 'gaussian_jobs_running': 0,
         'counts': counts, 'collection_barrier_reached': False,
         'scientifically_passed': 0, 'all_required_evidence_complete': False,
         'exclusive_reservation_available': True, 'process_exit_checks': processes,
         'cases': {r['case_id']: r['state'] for r in records},
         'records': records,
         'raw_runner_state_policy': 'Raw stage/result files preserve the last coordinator write. Running there is historical, not a live process. Frozen resume marks interrupted attempts and starts new attempts; complete stages are hash-checked and reused.',
         'partial_checkpoint_policy': 'Preserve incomplete checkpoint/scratch bytes, never promote them to a complete reference. Default recovery starts the interrupted stage again from its verified prior complete checkpoint, or a fresh initial guess.',
         'execution_priority': 'Proceed with all established common-root COV fixes independent of REF-001; physical reference conclusions remain pending.',
         'pause_evidence': 'runtime-incidents/ref-user-pause-20260906'}
atomic_json(PAUSE / 'pause-summary.json', value)
atomic_json(ROOT / 'references/REF-001/progress.json', value)
atomic_json(ROOT / 'progress.json', value)

inventory = []
for scope, directory in [('reference_jobs', JOBS / 'REF-001'), ('reproduction_jobs', JOBS / 'R001'),
                         ('stability_jobs', JOBS / 'S-001'), ('campaign', ROOT),
                         ('original_collection', ROOT.parent / 'fchk-collection-20260905')]:
    bins = Counter()
    counts_by_type = Counter()
    for current, folders, files in os.walk(directory):
        for name in files:
            path = Path(current) / name
            relative = path.relative_to(directory)
            kind = 'scratch' if 'scratch' in relative.parts else path.suffix.lower() or '(none)'
            bins[kind] += path.stat().st_size
            counts_by_type[kind] += 1
    inventory.append({'scope': scope, 'directory': str(directory), 'bytes_by_type': dict(bins),
                      'files_by_type': dict(counts_by_type), 'total_bytes': sum(bins.values())})
atomic_json(PAUSE / 'transfer-inventory.json', inventory)
print(json.dumps({'counts': counts, 'reservation_available': True, 'inventories': inventory}, ensure_ascii=False))
