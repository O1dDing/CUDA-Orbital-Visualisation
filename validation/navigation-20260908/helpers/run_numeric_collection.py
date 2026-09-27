from pathlib import Path
import argparse
import json
import os
import sys
import time

parser=argparse.ArgumentParser()
parser.add_argument('--kind',choices=['pilot','original','g-controls','diagnostic'],required=True)
parser.add_argument('--round-id',default='NUM-FIX-002')
args=parser.parse_args()
root=Path(r'F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906')/(args.round_id+'-'+args.kind)
sys.path.insert(0,str(root/'runner-source'))
from validation_process import atomic_json,physical_core_masks,run_tree
quota=json.loads((root/'manifest.json').read_text(encoding='utf-8'))['resources']
cores=quota['aggregate_cpu_cores'];offset=quota.get('core_slot_offset',0)
assert 1<=cores<=4 and offset==8 and 1<=quota['aggregate_memory_gib']<=64, 'Reserve the independent Gaussian allocation'
directory=root/('supervisor-'+str(time.time_ns()))
directory.mkdir()
env=dict(os.environ,COV_CPU_THREADS=str(quota['threads_per_worker']),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
         MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
command=[Path(r'F:\Dev\Python312\python.exe'),'-u',root/'runner-source'/'general_fix_validation_round.py',
         '--round',root,'--reference-deps',Path(r'F:\Dev\cov-validation-20260905\reference-deps')]
# The containing Job limits the whole coordinator and all its nested case Jobs.
# A 7-day emergency ceiling is not used as a polling wait; progress is atomic.
result=run_tree(command,directory,env,sum(physical_core_masks(offset+cores)[offset:]),quota['aggregate_memory_gib'],7*86400)
atomic_json(directory/'process.json',result)
print(json.dumps({'exit':result['exit_code'],'seconds':result['wall_seconds'],'evidence':str(directory)}),flush=True)
raise SystemExit(result['exit_code'])
