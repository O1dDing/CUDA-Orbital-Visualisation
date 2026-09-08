"""Anti-false-pass checks and diagnostic use of already terminal native evidence."""
from pathlib import Path
import copy
import ctypes
import hashlib
import json
import os
import sys
import time
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
repo=Path(r'F:\Dev\cov-native-validation-20260905')
sys.path[:0]=[str(repo/'tests'),r'F:\Dev\cov-validation-20260905\reference-deps']
from validation_process import atomic_json,physical_core_masks
kernel=ctypes.WinDLL('kernel32',use_last_error=True);kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
import numpy as np
from iodata import load_one
from gbasis.integrals.overlap import overlap_integral
from gbasis.wrappers import from_iodata
from review_numeric_fix_case import internal_basis_map
from review_density_matrices import same_input_density,compare_same_input,review_density
base=Path(r'F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906')
out=base/('density-matrix-checker-preflight-'+str(time.time_ns()));out.mkdir()
tests=[]
def positive(name,condition):
    assert condition,name;tests.append(name)
c=np.array([[1.,.25,-.7],[.3,.6,.1]])
nominal,radius,_=same_input_density(c,np.array([1.,-1.,0.]))
positive('signed-spin-weights',compare_same_input(nominal,nominal,radius)['pass_same_parsed_input'])
wrong=nominal.copy();wrong[0,1]+=1e-9
positive('small-corruption-is-rejected',not compare_same_input(wrong,nominal,radius)['pass_same_parsed_input'])
perm=np.array([2,0,1]);phase=np.array([-1,1,-1])
b,br,_=same_input_density(c[:,perm]*phase,np.array([1.,-1.,0.])[perm])
positive('same-spin-permutation-and-phase-invariant',compare_same_input(b,nominal,radius+br)['pass_same_parsed_input'])
z,zr,_=same_input_density(c,np.zeros(3));positive('known-zero',np.array_equal(z,np.zeros((2,2))))
for invalid in (np.full((2,2),np.nan),np.ones((3,3))):
    try:compare_same_input(invalid,z,zr)
    except ValueError:tests.append('invalid-matrix-rejected')
    else:raise AssertionError('Invalid matrix accepted')
pilot=base/'NUM-FIX-007-pilot';read=lambda p:json.loads(p.read_text(encoding='utf-8'))
assert read(pilot/'collection-complete.json')['all_terminal']
record=read(pilot/'manifest.json')['cases'][0];terminal=read(pilot/'cases/OLD-001/terminal.json')
evidence=pilot/terminal['evidence_directory'];production=read(evidence/'production.json')
events_path=evidence/'native/events.jsonl'
events=[json.loads(l) for l in events_path.read_text(encoding='utf-8').splitlines()]
density=next(e['data'] for e in events if e['kind']=='input.density_evidence')
mol=load_one(record['input']);s=overlap_integral(tuple(from_iodata(mol)));indices,scales=internal_basis_map(mol)
s=s[np.ix_(indices,indices)]*scales[:,None]*scales[None,:];actual_s=np.asarray(production['overlap']).reshape(s.shape)
def run(name,data):
    folder=out/name;folder.mkdir()
    checks=review_density(record['input'],mol,indices,scales,s,actual_s,data,5e-6,folder)
    atomic_json(folder/'checks.json',checks);return checks
actual=run('measured-pilot',density)
positive('preserved-native-pilot-density-all-checks-pass',all(c['status']=='pass' for c in actual))
bad=copy.deepcopy(density);bad['actual']['total']['packed'][0]+=1e-6
positive('wrong-producer-matrix-rejected',any(c['check']=='DENSITY-PRODUCER-AO-TRANSFORM-total' and c['status']=='fail' for c in run('corrupt-producer-negative',bad)))
bad=copy.deepcopy(density);bad['production_mo_reconstruction']['total']['packed'][0]+=1e-8
positive('wrong-raw-MO-reconstruction-rejected',any(c['check']=='DENSITY-MO-SAME-INPUT-total' and c['status']=='fail' for c in run('corrupt-reconstruction-negative',bad)))
bad=copy.deepcopy(density);bad['actual']['spin']['status']='missing-input'
positive('missing-state-cannot-be-known-zero',any(c['check']=='DENSITY-ACTUAL-AVAILABILITY-spin' and c['status']=='fail' for c in run('wrong-availability-negative',bad)))
atomic_json(out/'summary.json',dict(checks_passed=len(tests),tests=tests,measured_native_density_checks=len(actual),
    measured_case='OLD-001',measured_round_identity=read(pilot/'manifest.json')['round_identity'],
    measured_events_sha256=hashlib.sha256(events_path.read_bytes()).hexdigest(),
    source_sha256=hashlib.sha256((repo/'tests/review_density_matrices.py').read_bytes()).hexdigest(),
    scope='Checker preflight only; retain NUM007 pilot collection gaps. The new full-round checklist is not retroactively applied to that frozen round.',formal_case_passes=0))
print(json.dumps(dict(output=str(out),checker_checks=len(tests),measured_native_checks=len(actual))),flush=True)
