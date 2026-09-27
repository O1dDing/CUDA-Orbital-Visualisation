"""Build in an isolated directory under the actual Windows job resource limits."""
from pathlib import Path
import argparse
import json
import os
import sys
import time

parser=argparse.ArgumentParser()
parser.add_argument('--off',action='store_true')
parser.add_argument('--cores',type=int,default=4)
parser.add_argument('--memory-gib',type=int,default=64)
args=parser.parse_args()
assert 1<=args.cores<=4 and 1<=args.memory_gib<=64, 'Reserve 8 cores/48 GiB for the independent reference workflow'
REPO=Path(r'F:\Dev\cov-native-validation-20260905')
BUILD=Path(r'F:\Dev\cov-general-fixes-20260906-build'+('-off' if args.off else ''))
OUT=Path(r'F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906')
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(REPO/'tests'))
from validation_process import atomic_json, physical_core_masks, run_tree
cmake=Path(r'C:\Program Files\CMake\bin\cmake.exe')
env=dict(os.environ,COV_CPU_THREADS=str(args.cores),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
         MSBUILDDISABLENODEREUSE='1')
mask=sum(physical_core_masks(12)[8:8+args.cores])
helpers=[Path(r'F:\Dev\VS2022BuildTools\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64\vctip.exe'),
         Path(r'F:\Dev\VS2022BuildTools\VC\Tools\MSVC\14.44.35207\bin\Hostx86\x64\vctip.exe'),
         Path(r'F:\Dev\VS2022BuildTools\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64\mspdbsrv.exe')]
configure=[cmake,'-S',REPO,'-B',BUILD,'-G','Visual Studio 17 2022','-A','x64',
    '-DCMAKE_GENERATOR_INSTANCE=F:/Dev/VS2022BuildTools','-DCOV_ENABLE_CUDA=ON',
    '-DCOV_ENABLE_VALIDATION='+('OFF' if args.off else 'ON'),'-DCOV_BUILD_TESTS=ON',
    '-DFETCHCONTENT_SOURCE_DIR_GLFW=F:/CalChem/COV/DEV/build/_deps/glfw-src',
    '-DFETCHCONTENT_SOURCE_DIR_IMGUI=F:/CalChem/COV/DEV/build/_deps/imgui-src']
run=OUT/('build-preflight-'+('off-' if args.off else 'on-')+str(time.time_ns()))
run.mkdir()
for name,argv in [('configure',configure),('build',[cmake,'--build',BUILD,'--config','Release','--parallel',str(args.cores),'--','/nodeReuse:false'])]:
    directory=run/name
    directory.mkdir()
    result=run_tree(argv,directory,env,mask,args.memory_gib,3600,finite_driver_helpers=helpers)
    atomic_json(directory/'process.json',result)
    print(json.dumps({'stage':name,'exit':result['exit_code'],'wall_seconds':result['wall_seconds'],
                      'log':result.get('console_log',str(directory/'launcher.log'))}),flush=True)
    if result['exit_code'] or result['timed_out']:
        raise SystemExit(result['exit_code'] or 1)
print(json.dumps({'build_complete':True,'directory':str(BUILD)}),flush=True)
