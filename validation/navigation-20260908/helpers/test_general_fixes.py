from pathlib import Path
import argparse
import json
import os
import sys
import time
parser=argparse.ArgumentParser()
parser.add_argument('--off',action='store_true')
parser.add_argument('--match')
settings=parser.parse_args()
ROOT=Path(r'F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906')
OUT=ROOT/('component-tests-'+('off-' if settings.off else 'on-')+str(time.time_ns()))
OUT.mkdir()
sys.path.insert(0,r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json,physical_core_masks,run_tree
env=dict(os.environ,COV_CPU_THREADS='2',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
         MSBUILDDISABLENODEREUSE='1')
args=[Path(r'C:\Program Files\CMake\bin\ctest.exe'),'--test-dir',Path(r'F:\Dev\cov-general-fixes-20260906-build'+('-off' if settings.off else '')),
      '-C','Release','--parallel','2','--output-on-failure','--output-junit',OUT/'ctest.xml']
if settings.match:args.extend(['-R',settings.match])
result=run_tree(args,OUT,env,sum(physical_core_masks(12)[8:12]),64,1800)
atomic_json(OUT/'process.json',result)
print(json.dumps({'exit_code':result['exit_code'],'wall_seconds':result['wall_seconds'],'output':str(OUT)}),flush=True)
raise SystemExit(result['exit_code'])
