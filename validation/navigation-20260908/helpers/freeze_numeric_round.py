"""Freeze reproducible COV programs, actual source inputs and collection rules."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import time

REPO=Path(r'F:\Dev\cov-native-validation-20260905')
CAMPAIGN=Path(r'F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906')
BASE=Path(r'F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906')
parser=argparse.ArgumentParser()
parser.add_argument('--kind',choices=['pilot','original','diagnostic'],required=True)
parser.add_argument('--case-id',help='One explicitly identified diagnostic case; never a formal replacement')
parser.add_argument('--round-id',default='NUM-FIX-002')
parser.add_argument('--production-commit',default='a9bf30a18a16c16e876abec133691b911696da3f')
args=parser.parse_args()
root=BASE/(args.round_id+'-'+args.kind)
if root.exists(): raise FileExistsError('Existing frozen round is never replaced')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1<<20),b''):h.update(block)
    return h.hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def git(*args):return subprocess.check_output([r'F:\Dev\Git\cmd\git.exe',*args],cwd=REPO).decode().strip()
commit=git('rev-parse',args.production_commit)
runner_commit=git('rev-parse','HEAD')
assert not git('status','--porcelain')
assert all(p.startswith('tests/') and p.endswith('.py') for p in git('diff','--name-only',commit,runner_commit).splitlines())
campaign=read(CAMPAIGN/'campaign.json')
selected=campaign['cases']
if args.kind=='pilot':selected=[x for x in selected if x['case_id']==campaign['pilot_case_id']]
if args.kind=='diagnostic':
    assert args.case_id, 'An explicit case is required for a diagnostic snapshot'
    selected=[x for x in selected if x['case_id']==args.case_id]
else:assert args.case_id is None, 'Formal collection and the designated pilot cannot be replaced by a subset'
assert len(selected)==(273 if args.kind=='original' else 1)
root.mkdir()
(root/'programs').mkdir()
programs={}
for source,name in [
    (Path(r'F:\Dev\cov-general-fixes-20260906-build\Release\cov_validation.exe'),'cov_validation.exe'),
    (Path(r'F:\Dev\cov-general-fixes-20260906-build\Release\cov_scientific_audit_dump.exe'),'cov_scientific_audit_dump.exe'),
    (Path(r'F:\Dev\cov-general-fixes-20260906-build-off\Release\cov.exe'),'cov.exe'),
    (Path(r'F:\Dev\cov-general-fixes-20260906-build-off\Release\cov_scientific_audit_dump.exe'),'cov_scientific_audit_dump_off.exe')]:
    shutil.copy2(source,root/'programs'/name);programs[name]=sha(root/'programs'/name)
(root/'runner-source').mkdir()
runners={}
for name in ('general_fix_validation_round.py','native_collection_batch.py','validation_process.py',
             'review_numeric_fix_case.py','grid_reference_compare.py','independent_grid_reference.py',
             'review_export_snapshot_round.py'):
    shutil.copy2(REPO/'tests'/name,root/'runner-source'/name);runners[name]=sha(root/'runner-source'/name)
image_deps=Path(r'F:\Dev\cov-native-validation-20260905-build\python-deps')
image_artifacts={}
for directory in ('PIL','pillow-12.3.0.dist-info'):
    for source in (image_deps/directory).rglob('*'):
        if not source.is_file() or '__pycache__' in source.parts: continue
        relative=source.relative_to(image_deps)
        destination=root/'image-deps'/relative
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,destination);image_artifacts[str(relative)]=sha(destination)
cases=[]
input_hashes={}
for old in selected:
    archived=Path(old['original_fchk'])
    assert sha(archived)==old['original_fchk_sha256']
    identity=read(archived.parent/'source-identity.json')
    target=root/'inputs'/'.chk'/identity['relative_path']
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(archived,target)
    input_hashes[str(target.relative_to(root))]=sha(target)
    logs=[]
    for item in identity['logs']:
        original=Path(item['original'])
        try:relative=original.relative_to(Path(r'F:\CalChem\Gaussian\Files\.log'))
        except ValueError:relative=Path(identity['relative_path']).parent/original.name
        log=root/'inputs'/'.log'/relative
        log.parent.mkdir(parents=True,exist_ok=True)
        frozen=archived.parent/item['copy']
        if log.exists():assert sha(log)==sha(frozen)
        else:shutil.copy2(frozen,log)
        logs.append(str(log));input_hashes[str(log.relative_to(root))]=sha(log)
    stat=target.stat()
    cases.append({'case_id':old['case_id'],'input':str(target),'relative_path':identity['relative_path'],
        'source_size':stat.st_size,'source_mtime_ns':stat.st_mtime_ns,'log_candidates':logs,
        'sha256':old['original_fchk_sha256'],'original_snapshot':str(archived)})
manifest={'schema':2,'name':args.round_id+'-'+args.kind,'kind':args.kind,'git_commit':commit,'runner_source_commit':runner_commit,
    'case_count':len(cases),'formal_original_case_count':273,'cases':cases,
    'programs':programs,'runner_sources':runners,'input_artifacts':input_hashes,
    'image_dependency_artifacts':image_artifacts,'image_dependency_version':'Pillow 12.3.0',
    'resources':{'workers':2 if args.kind=='original' else 1,'threads_per_worker':2,
                 'memory_gib_per_worker':24,'aggregate_cpu_cores':4,'aggregate_memory_gib':64,
                 'core_slot_offset':8,'reserved_independent_reference_cpu_cores':8,
                 'reserved_independent_reference_memory_gib':48},
    'native_timeout_seconds':1800,'created_epoch':time.time(),
    'collection_policy':'Collect all cases to terminal states, retaining failures, before any scientific comparison or algorithm change',
    'thresholds':{'all_mo_nrms_max':1e-4,'all_mo_abs_cosine_min':1-1e-7,
       'all_mo_relative_peak_error_max':1e-3,'same_grid_minimum_unique_samples':8000,
       'basis_overlap_max_absolute_error':1e-8,'mo_orthonormality_max_error':5e-6,
       'electron_trace_error':5e-6,'raw_mo_coefficients_max_absolute_error':1e-12,
       'mapped_mo_coefficients_max_absolute_error':1e-12,'energy_error_hartree':1e-10},
    'scope':['actual all-MO textures and frame identities','real browser/member/spin/control paths',
             'actual PNG/SVG/CSV/JSON export','complete frontier textures at 64 and 128',
             'independent AO integrals and all spin MO metrics','ON/OFF scientific-data consistency'],
    'not_a_final_acceptance':['Physical reference acceptance remains pending; independent reference workflow has its own stage-boundary pause',
                            'remaining known product issues still open',
                            'full invariant and visual/chemical adjudication still required'],
    'formal_case_passes':0,'gaussian_calculations_started':0}
manifest['round_identity']=hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
(root/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'round':str(root),'cases':len(cases),'commit':commit,'round_identity':manifest['round_identity']},ensure_ascii=False),flush=True)
