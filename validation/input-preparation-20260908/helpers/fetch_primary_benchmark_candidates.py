"""Fetch selected geometry records from the authors' pinned public repository."""
from pathlib import Path, PurePosixPath
import ctypes
import hashlib
import json
import sys
import time
import urllib.request

sys.path.insert(0,r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import physical_core_masks, atomic_json
kernel=ctypes.WinDLL('kernel32',use_last_error=True);kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
workspace=Path(__file__).resolve().parent.parent
base=workspace/'outputs/cov-complete-validation-20260906'
root=base/'benchmark-geometry-preparation-v1'
branch=json.loads((base/'external-primary-sources-v1/GMTKN55-v1-branch-api.json').read_text())
commit=branch['commit']['sha']
assert branch['name']=='v1' and len(commit)==40
screen=json.loads((root/'local-formula-screen.json').read_text(encoding='utf-8'))
selected=[r for r in screen['records'] if r['pubchem_status']!='computed-starting-conformer-retrieved' and len(r['formula_charge_matches'])==1]
results=[]
def fetch(url,path):
    if path.exists():return path.read_bytes()
    req=urllib.request.Request(url,headers={'User-Agent':'COV-primary-geometry-preparation/1.0','Accept':'application/vnd.github+json'})
    with urllib.request.urlopen(req,timeout=20) as response:
        raw=response.read(2*1024*1024+1)
    if len(raw)>2*1024*1024:raise ValueError('Response exceeds geometry record bound')
    path.write_bytes(raw)
    return raw
for candidate in selected:
    out=root/'primary-records'/candidate['candidate_id'];out.mkdir(parents=True,exist_ok=True)
    receipt=out/'receipt.json'
    if receipt.exists():
        result=json.loads(receipt.read_text())
        if result['commit']!=commit:raise ValueError('Do not mix primary-source commits')
        results.append(result);continue
    local=candidate['formula_charge_matches'][0]
    relative=PurePosixPath(local['relative_path'])
    directory=relative.with_suffix('').as_posix()
    result=dict(candidate_id=candidate['candidate_id'],query=candidate['query'],commit=commit,
                local_candidate=local,source_directory=directory,started_epoch=time.time(),files=[],
                formal_case=False,gaussian_started=False,accepted_reference_geometry=False)
    try:
        url=f'https://api.github.com/repos/grimme-lab/GMTKN55/contents/{directory}?ref={commit}'
        raw=fetch(url,out/'directory.json')
        result.update(directory_url=url,directory_sha256=hashlib.sha256(raw).hexdigest())
        entries=json.loads(raw)
        result['entry_names']=[entry['name'] for entry in entries]
        for entry in entries:
            if entry['type']!='file' or entry['name'] not in ['coord','coord.xyz','.CHRG','.UHF']:
                continue
            if entry['size']>1024*1024:raise ValueError('Unusually large geometry record')
            payload=fetch(entry['download_url'],out/entry['name'])
            git_blob=hashlib.sha1(b'blob '+str(len(payload)).encode()+b'\0'+payload).hexdigest()
            if git_blob!=entry['sha']:raise ValueError('Downloaded blob differs from pinned primary tree')
            result['files'].append(dict(name=entry['name'],url=entry['download_url'],
                git_blob_sha1=git_blob,sha256=hashlib.sha256(payload).hexdigest(),bytes=len(payload)))
        result['status']='primary-records-collected'
    except Exception as error:
        result.update(status='source-retrieval-failed-retained',error=f'{type(error).__name__}: {error}')
    result['finished_epoch']=time.time();atomic_json(receipt,result);results.append(result)
    print(json.dumps({k:result.get(k) for k in ['candidate_id','status','entry_names','error']},ensure_ascii=False),flush=True)
atomic_json(root/'primary-collection.json',dict(commit=commit,selected_count=len(selected),terminal_count=len(results),records=results,
    formal_external_cases=0,gaussian_started=False,accepted_reference_geometries=0,
    primary_repository='https://github.com/grimme-lab/GMTKN55',primary_doi='10.1039/C7CP04913G'))
