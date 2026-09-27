from pathlib import Path
import hashlib
import json
import os
import shutil
import sys
import time

workspace=Path(r'F:\Codex\2026-09-05\branch-15')
sys.path.insert(0,r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json,physical_core_masks,run_tree
out=workspace/'outputs/cov-complete-validation-20260906/regression-fixtures'/('density-reviewer-preflight-'+str(time.time_ns()))
(out/'source').mkdir(parents=True)
names=['review_density_source_controls.py','test_density_source_reviewer.py','run_density_reviewer_tests.py']
for name in names:shutil.copy2(workspace/'work'/name,out/'source'/name)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
atomic_json(out/'identity.json',dict(source_hashes={n:sha(out/'source'/n) for n in names},
    controls_manifest_sha256=sha(workspace/'outputs/cov-complete-validation-20260906/regression-fixtures/density-source-v1/manifest.json'),
    cpu_core_slot=8,memory_gib=1,scope='Reviewer false-pass prevention tests; no COV execution or current-round result analysis'))
result=run_tree([sys.executable,'-X','utf8',out/'source/test_density_source_reviewer.py'],out,
    dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'),physical_core_masks(12)[8],1,120)
atomic_json(out/'process.json',result)
print(json.dumps(dict(output=str(out),exit_code=result['exit_code'],wall_seconds=result['wall_seconds'],console_log=result['console_log'])),flush=True)
raise SystemExit(result['exit_code'])
