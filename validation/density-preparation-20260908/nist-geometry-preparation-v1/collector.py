"""Archive five observed NIST source URLs as preparation; no quantum jobs."""
from pathlib import Path
import ctypes
import hashlib
import json
import sys
import time
import urllib.request

workspace=Path(__file__).resolve().parent.parent
repo=Path(r'F:\Dev\cov-native-validation-20260905')
sys.path.insert(0,str(repo/'tests'))
from validation_process import atomic_json,physical_core_masks
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
root=workspace/'outputs/cov-complete-validation-20260906/nist-geometry-preparation-v1'
sources=[
 dict(candidate_id='PREP-002',name='germane',elements={'Ge':1,'H':4},charge=0,
      url='https://cccbdb.nist.gov/energy3x.asp?basis=21&casno=7782652&method=57',
      scope='computed HSEh1PBE/6-311G** geometry; singlet minimum in source; starting input only'),
 dict(candidate_id='PREP-007',name='chlorine trifluoride',elements={'Cl':1,'F':3},charge=0,
      url='https://cccbdb.nist.gov/exp2x.asp?casno=7790912&charge=0',
      scope='experimental-data index; initial state row D3h is nonminimum; explicit C2v minimum is another conformer; no Cartesian values observed in web extraction'),
 dict(candidate_id='PREP-008',name='bromine pentafluoride',elements={'Br':1,'F':5},charge=0,
      url='https://cccbdb.nist.gov/exp2x.asp?casno=7789302&charge=0',
      scope='experimental-data index; coordinate payload and original references to be reviewed'),
 dict(candidate_id='PREP-009',name='iodine heptafluoride',elements={'I':1,'F':7},charge=0,
      url='https://cccbdb.nist.gov/alldata2x.asp?casno=16921963&charge=0',
      scope='species index; follow actual geometry links after source retrieval'),
 dict(candidate_id='PREP-042',name='phenoxy radical',elements={'C':6,'H':5,'O':1},charge=0,
      url='https://cccbdb.nist.gov/energy3x.asp?basis=0&casno=2122465&method=64',
      scope='G3B3 result page whose geometry is explicitly B3LYP/6-31G*; doublet minimum; starting input only')]
if root.exists():raise FileExistsError('Keep the existing finite source collection')
root.mkdir()
shutil_source=root/'collector.py';shutil_source.write_bytes(Path(__file__).read_bytes())
atomic_json(root/'manifest.json',dict(expected_sources=len(sources),sources=sources,timeout_seconds=20,
    maximum_bytes=2*1024*1024,concurrent_requests=1,automatic_retries=0,formal_external_cases=0,
    accepted_quality_references=0,gaussian_started=False,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
records=[]
for source in sources:
    record=dict(source,started_epoch=time.time());out=root/source['candidate_id'];out.mkdir()
    try:
        request=urllib.request.Request(source['url'],headers={'User-Agent':'COV-scientific-source-preparation/1.0','Accept':'text/html'})
        with urllib.request.urlopen(request,timeout=20) as response:
            raw=response.read(2*1024*1024+1)
            if len(raw)>2*1024*1024:raise ValueError('Source response exceeds declared bound')
            record.update(http_status=response.status,content_type=response.headers.get('Content-Type'),
                          charset=response.headers.get_content_charset(),resolved_url=response.geturl())
        (out/'source.html').write_bytes(raw)
        record.update(status='html-retrieved-content-review-pending',bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
    except Exception as error:
        record.update(status='retrieval-failure-retained',error=f'{type(error).__name__}: {error}')
    record['finished_epoch']=time.time();records.append(record);atomic_json(out/'receipt.json',record)
    atomic_json(root/'progress.json',dict(expected_sources=5,terminal_sources=len(records),records=records))
    print(json.dumps({k:record.get(k) for k in ['candidate_id','status','bytes','error']}),flush=True)
    time.sleep(.7)
atomic_json(root/'collection-complete.json',dict(expected_sources=5,terminal_sources=len(records),all_terminal=True,
    records=records,formal_external_cases=0,accepted_quality_references=0,gaussian_started=False))
