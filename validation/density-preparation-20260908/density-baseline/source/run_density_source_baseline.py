"""Prepare, then execute a production-library density probe after the full-round barrier."""
from pathlib import Path
import argparse
import ctypes
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

parser=argparse.ArgumentParser()
parser.add_argument('action',choices=['prepare','run'])
args=parser.parse_args()
workspace=Path(r'F:\Codex\2026-09-05\branch-15')
repo=Path(r'F:\Dev\cov-native-validation-20260905')
build=Path(r'F:\Dev\cov-general-fixes-20260906-build')
controls=workspace/'outputs/cov-complete-validation-20260906/regression-fixtures/density-source-v1'
root=controls/'baseline-production-21036b2-v3'
round_root=workspace/'outputs/general-fixes-20260906/NUM-FIX-006-original'
sys.path.insert(0,str(repo/'tests' if args.action=='prepare' else root/'source/helpers'))
from validation_process import atomic_json, physical_core_masks, run_tree
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
git=r'F:\Dev\Git\cmd\git.exe'
scientific_commit=read(controls/'manifest.json')['scientific_source_commit_to_probe']

if args.action=='prepare':
    if root.exists():raise FileExistsError('The prepared library/header snapshot is immutable')
    changed=subprocess.check_output([git,'diff','--name-only',scientific_commit,'--','src','include','CMakeLists.txt'],cwd=repo,text=True)
    assert not changed,'Scientific source differs from the original control baseline'
    (root/'source/density-probe-project').mkdir(parents=True)
    shutil.copytree(repo/'include/cov',root/'source/include/cov')
    for name in ['density_source_probe.cpp','density-probe-project/CMakeLists.txt','run_density_source_baseline.py']:
        shutil.copy2(workspace/'work'/name,root/'source'/name)
    (root/'source/helpers').mkdir()
    shutil.copy2(repo/'tests/validation_process.py',root/'source/helpers/validation_process.py')
    (root/'libraries').mkdir()
    original_libraries=[build/'Release/cov_parser.lib',Path(r'F:\Dev\CUDA\v12.8\lib\x64\cudart_static.lib')]
    library_records=[]
    for source in original_libraries:
        before=source.stat();destination=root/'libraries'/source.name
        shutil.copy2(source,destination)
        assert sha(source)==sha(destination) and before.st_mtime_ns==source.stat().st_mtime_ns
        library_records.append(dict(source=str(source),source_modified_ns=before.st_mtime_ns,
                                    file=str(destination.relative_to(root)),sha256=sha(destination)))
    for item in read(controls/'manifest.json')['controls']:
        assert sha(controls/item['file'])==item['sha256']
    records={str(p.relative_to(root)).replace('\\','/'):sha(p) for p in sorted(root.rglob('*')) if p.is_file()}
    atomic_json(root/'identity.json',dict(scientific_source_commit=scientific_commit,
        preparation_repo_commit=subprocess.check_output([git,'rev-parse','HEAD'],cwd=repo,text=True).strip(),
        preparation_epoch=time.time(),controls_manifest_sha256=sha(controls/'manifest.json'),
        files=records,libraries=library_records,
        library_provenance='Existing ON build used in NUM-FIX-006 preflight; only diagnostic wrapper will be compiled.',
        build_preflight='outputs/general-fixes-20260906/navigation-preflight-1788840607516376700',
        production_density_algorithms_changed=False,probe_executed=False,
        resource_policy=dict(cores=4,memory_gib=64,core_slot_offset=8,
                             independent_gaussian_reserved_cores=8,independent_gaussian_reserved_memory_gib=48)))
    print(json.dumps(dict(status='prepared-not-compiled-or-executed',root=str(root),files=len(records))),flush=True)
    raise SystemExit(0)

identity=read(root/'identity.json')
assert sha(Path(__file__))==identity['files']['source/run_density_source_baseline.py'],'Use the frozen baseline runner'
assert sha(controls/'manifest.json')==identity['controls_manifest_sha256']
for name,digest in identity['files'].items():assert sha(root/name)==digest,name
collection=read(round_root/'collection-complete.json')
pipeline=read(round_root/'pipeline-progress.json')
assert collection['all_terminal'] and collection['terminal_cases']==273
assert pipeline['stage']=='reviews-terminal-awaiting-root-cause-analysis'
assert (round_root/'review-v1/unified-review/summary.json').exists()
assert (round_root/'view-review-v1/summary.json').exists()
if (root/'execution-started.json').exists():raise FileExistsError('Retain the completed or interrupted baseline; do not silently rerun')
atomic_json(root/'execution-started.json',dict(started_epoch=time.time(),round_identity=collection['round_identity'],
    preceding_pipeline=pipeline,gaussian_controlled=False))
cmake=Path(r'C:\Program Files\CMake\bin\cmake.exe')
env=dict(os.environ,COV_CPU_THREADS='2',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MSBUILDDISABLENODEREUSE='1')
mask=sum(physical_core_masks(12)[8:12])
helpers=[Path(r'F:\Dev\VS2022BuildTools\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64\vctip.exe'),
         Path(r'F:\Dev\VS2022BuildTools\VC\Tools\MSVC\14.44.35207\bin\Hostx86\x64\vctip.exe'),
         Path(r'F:\Dev\VS2022BuildTools\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64\mspdbsrv.exe')]
configure=[cmake,'-S',root/'source/density-probe-project','-B',root/'build','-G','Visual Studio 17 2022','-A','x64',
    '-DCMAKE_GENERATOR_INSTANCE=F:/Dev/VS2022BuildTools',
    '-DCOV_PRODUCTION_INCLUDE='+str(root/'source/include'),
    '-DCOV_PRODUCTION_LIBRARY='+str(root/'libraries/cov_parser.lib'),
    '-DCOV_CUDART_LIBRARY='+str(root/'libraries/cudart_static.lib')]
for name,argv in [('configure',configure),('build',[cmake,'--build',root/'build','--config','Release','--parallel','2','--','/nodeReuse:false'])]:
    directory=root/name if name=='configure' else root/'build-process';directory.mkdir(exist_ok=True)
    result=run_tree(argv,directory,env,mask,64,600,finite_driver_helpers=helpers)
    atomic_json(directory/'process.json',result)
    if result['exit_code'] or result['timed_out']:raise SystemExit(result['exit_code'] or 1)
exe=root/'build/Release/cov_density_source_probe.exe'
spec=read(controls/'manifest.json');records=[];start=time.time()
for item in spec['controls']:
    directory=root/'cases'/item['control_id'];directory.mkdir(parents=True)
    record=dict(control_id=item['control_id'],input_sha256=item['sha256'],exe_sha256=sha(exe),started_epoch=time.time())
    try:
        assert sha(controls/item['file'])==item['sha256']
        result=run_tree([exe,controls/item['file']],directory,env,mask,64,90)
        atomic_json(directory/'process.json',result)
        raw=Path(result['console_log']).read_text(encoding='utf-8')
        record.update(status='collected',exit_code=result['exit_code'],timed_out=result['timed_out'],
                      wall_seconds=result['wall_seconds'],log_sha256=sha(Path(result['console_log'])))
        if result['timed_out']:record['status']='timeout-retained'
        else:
            observed=json.loads(raw)
            atomic_json(directory/'observed.json',observed)
            record['observed_sha256']=sha(directory/'observed.json')
    except Exception as error:
        record.update(status='collection-error-retained',error=f'{type(error).__name__}: {error}')
    record['finished_epoch']=time.time();records.append(record)
    atomic_json(directory/'receipt.json',record)
    atomic_json(root/'progress.json',dict(expected_controls=20,terminal_controls=len(records),records=records,
                                        scientific_review='pending complete control batch'))
atomic_json(root/'collection-complete.json',dict(all_terminal=True,expected_controls=20,terminal_controls=len(records),
    wall_seconds=time.time()-start,records=records,scientific_review='pending',formal_cases_added=0))
print(json.dumps(dict(status='all-control-inputs-terminal',root=str(root),controls=len(records))),flush=True)
