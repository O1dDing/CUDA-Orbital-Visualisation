from pathlib import Path
import argparse
import json
import os
import sys
import time
parser=argparse.ArgumentParser()
parser.add_argument('--round',type=Path,required=True)
parser.add_argument('--review-id',default='view-review-v1')
args=parser.parse_args()
root=args.round.resolve()
sys.path.insert(0,str(root/'runner-source'))
from validation_process import atomic_json,physical_core_masks,run_tree
out=root/('view-review-process-'+str(time.time_ns())); out.mkdir()
env=dict(os.environ,COV_CPU_THREADS='3',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
result=run_tree([sys.executable,root/'runner-source'/'review_export_snapshot_round.py',
    '--round',root,'--review-id',args.review_id],out,env,sum(physical_core_masks(3)),24,1800)
atomic_json(out/'process.json',result)
summary=root/args.review_id/'summary.json'
print(json.dumps({'exit_code':result['exit_code'],'wall_seconds':result['wall_seconds'],
    'process_directory':str(out),'summary':str(summary)},ensure_ascii=False),flush=True)
if summary.exists():
    data=json.loads(summary.read_text(encoding='utf-8'))
    print(json.dumps({k:data[k] for k in ('case_count','status_counts','bundle_count','failure_count','wall_seconds')},ensure_ascii=False),flush=True)
raise SystemExit(result['exit_code'])
