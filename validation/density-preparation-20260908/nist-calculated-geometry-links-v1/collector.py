"""Follow public NIST geometry links using each molecule's normal anonymous session."""
from pathlib import Path
from html.parser import HTMLParser
import ctypes
import hashlib
import http.cookiejar
import json
import sys
import time
import urllib.parse
import urllib.request

workspace=Path(r'F:\Codex\2026-09-05\branch-15')
sys.path.insert(0,r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json,physical_core_masks
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
root=workspace/'outputs/cov-complete-validation-20260906/nist-calculated-geometry-links-v1'
if root.exists():raise FileExistsError('Keep previous source navigation evidence')
root.mkdir();(root/'collector.py').write_bytes(Path(__file__).read_bytes())
sources=[('PREP-007','7790912'),('PREP-009','16921963'),('PREP-045','2143580')]
atomic_json(root/'manifest.json',dict(sources=[dict(candidate_id=i,cas_digits=c) for i,c in sources],
    http_timeout_seconds=20,max_bytes_per_response=2*1024*1024,retries=0,parallel_requests=1,
    session='Anonymous per-molecule CookieJar; cookies are not persisted or exported',
    formal_external_cases=0,gaussian_started=False))

class Links(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True);self.links=[];self.anchor=None;self.text=[];self.title=[];self.in_title=False
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if tag=='a' and attrs.get('href'):self.anchor=dict(href=attrs['href'],text=[],preceding_text=' '.join(self.text)[-220:])
        if tag=='title':self.in_title=True
    def handle_data(self,data):
        clean=' '.join(data.split())
        if clean:self.text.append(clean)
        if self.anchor is not None:self.anchor['text'].append(data)
        if self.in_title:self.title.append(data)
    def handle_endtag(self,tag):
        if tag=='a' and self.anchor is not None:
            self.anchor['text']=' '.join(''.join(self.anchor['text']).split());self.links.append(self.anchor);self.anchor=None
        if tag=='title':self.in_title=False

records=[]
for ident,cas in sources:
    directory=root/ident;directory.mkdir();record=dict(candidate_id=ident,cas_digits=cas,requests=[])
    opener=urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    def fetch(url,stage):
        request=urllib.request.Request(url,headers={'User-Agent':'COV-scientific-source-preparation/1.0','Accept':'text/html'})
        before=time.time()
        with opener.open(request,timeout=20) as response:
            raw=response.read(2*1024*1024+1)
            if len(raw)>2*1024*1024:raise ValueError('Response exceeds source-navigation bound')
            charset=response.headers.get_content_charset() or 'utf-8';resolved=response.url
            receipt=dict(url=url,resolved_url=resolved,status=response.status,content_type=response.headers.get('Content-Type'),
                         bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),started_epoch=before,finished_epoch=time.time())
        (directory/(stage+'.html')).write_bytes(raw)
        document=raw.decode(charset,errors='strict')
        if any(term in document.lower() for term in ('recaptcha','checking your browser','akamai','access denied')):
            raise ValueError('Access-verification page retained; no further navigation')
        parsed=Links();parsed.feed(document)
        for link in parsed.links:link['url']=urllib.parse.urljoin(resolved,link['href'])
        receipt.update(title=' '.join(parsed.title),stage=stage);record['requests'].append(receipt)
        return parsed
    try:
        initial=fetch(f'https://cccbdb.nist.gov/exp2x.asp?casno={cas}&charge=0','experimental-index')
        matches=[l for l in initial.links if urllib.parse.urlparse(l['url']).path.lower().endswith('/geom2x.asp')]
        assert len({l['url'] for l in matches})==1,'No unique source-provided calculated-geometry link'
        record['followed_link']=matches[0]
        geometry=fetch(matches[0]['url'],'calculated-geometry-index')
        useful=[l for l in geometry.links if any(term in (l['url']+' '+l['text']).lower() for term in ('geom3','energy3','basis=','method=','conform'))]
        record.update(status='public-geometry-index-retrieved',title=' '.join(geometry.title),links=useful,
                      content_tail=' '.join(geometry.text)[-13000:])
    except Exception as error:record.update(status='navigation-failure-retained',error=f'{type(error).__name__}: {error}')
    atomic_json(directory/'receipt.json',record);records.append(record)
    print(json.dumps(dict(candidate_id=ident,status=record['status'],links=len(record.get('links',[])),error=record.get('error'))),flush=True)
atomic_json(root/'collection-complete.json',dict(all_terminal=True,terminal_sources=3,records=records,
    formal_external_cases=0,accepted_quality_references=0,gaussian_started=False))
