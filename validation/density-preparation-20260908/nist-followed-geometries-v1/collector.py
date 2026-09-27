"""Retrieve two already observed model-chemistry links with ordinary molecule sessions."""
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
base=workspace/'outputs/cov-complete-validation-20260906'
previous=base/'nist-calculated-geometry-links-v1'
root=base/'nist-followed-geometries-v1'
if root.exists():raise FileExistsError('Retain existing source retrievals')
root.mkdir();(root/'collector.py').write_bytes(Path(__file__).read_bytes())
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources=[]
for ident in ('PREP-007','PREP-045'):
    receipt=read(previous/ident/'receipt.json')
    observed=[l for l in receipt['links'] if urllib.parse.urlparse(l['url']).path=='/geom3x.asp'
        and urllib.parse.parse_qs(urllib.parse.urlparse(l['url']).query)=={'method':['64'],'basis':['0']}]
    assert len({l['url'] for l in observed})==1
    sources.append(dict(candidate_id=ident,initial_url=receipt['requests'][0]['url'],followed_link=observed[0],
                        source_link_receipt_sha256=sha(previous/ident/'receipt.json')))
atomic_json(root/'manifest.json',dict(sources=sources,request_timeout_seconds=20,retries=0,
    selection='Observed G3B3 geometry link for both molecules; actual geometry method and conformation must be reviewed in payload.',
    formal_external_cases=0,gaussian_started=False))

class Extract(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True);self.tables=[];self.stack=[];self.text=[];self.links=[];self.anchor=None
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if tag=='table':
            table=dict(index=len(self.tables),rows=[],row=None,cell=None,preceding_text=' '.join(self.text)[-450:])
            self.tables.append(table);self.stack.append(table)
        if self.stack:
            if tag=='tr':self.stack[-1]['row']=[]
            if tag in ('td','th'):self.stack[-1]['cell']=[]
        if tag=='a' and attrs.get('href'):self.anchor=dict(href=attrs['href'],text=[])
    def handle_data(self,value):
        clean=' '.join(value.split())
        if clean:self.text.append(clean)
        if self.stack and self.stack[-1]['cell'] is not None:self.stack[-1]['cell'].append(value)
        if self.anchor is not None:self.anchor['text'].append(value)
    def handle_endtag(self,tag):
        if self.stack:
            t=self.stack[-1]
            if tag in ('td','th') and t['cell'] is not None:
                if t['row'] is not None:t['row'].append(' '.join(''.join(t['cell']).split()))
                t['cell']=None
            if tag=='tr' and t['row'] is not None:t['rows'].append(t['row']);t['row']=None
            if tag=='table':self.stack.pop()
        if tag=='a' and self.anchor is not None:
            self.anchor['text']=' '.join(''.join(self.anchor['text']).split());self.links.append(self.anchor);self.anchor=None
records=[]
for source in sources:
    ident=source['candidate_id'];out=root/ident;out.mkdir();record=dict(source,requests=[])
    opener=urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    try:
        for stage,url in [('experimental-index',source['initial_url']),('geometry',source['followed_link']['url'])]:
            before=time.time()
            with opener.open(urllib.request.Request(url,headers={'User-Agent':'COV-scientific-source-preparation/1.0','Accept':'text/html'}),timeout=20) as response:
                raw=response.read(2*1024*1024+1)
                if len(raw)>2*1024*1024:raise ValueError('Source exceeds retrieval bound')
                charset=response.headers.get_content_charset() or 'utf-8';resolved=response.url
                entry=dict(url=url,resolved_url=resolved,http_status=response.status,content_type=response.headers.get('Content-Type'),
                           stage=stage,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),started_epoch=before,finished_epoch=time.time())
            (out/(stage+'.html')).write_bytes(raw);record['requests'].append(entry)
            document=raw.decode(charset,errors='strict')
            if any(t in document.lower() for t in ('recaptcha','checking your browser','access denied')):
                raise ValueError('Access verification retained; no further navigation')
            if stage=='geometry':
                parsed=Extract();parsed.feed(document)
                for link in parsed.links:link['url']=urllib.parse.urljoin(resolved,link['href'])
                extracted=dict(candidate_id=ident,tables=[{k:v for k,v in t.items() if k not in ('row','cell')} for t in parsed.tables],
                    related_links=[l for l in parsed.links if any(t in l['url'].lower() for t in ('geom3','energy3','conform','casno='))],
                    content_tail=' '.join(parsed.text)[-17000:],source_sha256=entry['sha256'])
                atomic_json(out/'tables.json',extracted)
                record.update(status='geometry-page-retrieved-awaiting-review',tables=len(parsed.tables),normalized_sha256=sha(out/'tables.json'))
    except Exception as error:record.update(status='retrieval-failure-retained',error=f'{type(error).__name__}: {error}')
    atomic_json(out/'receipt.json',record);records.append(record)
    print(json.dumps(dict(candidate_id=ident,status=record['status'],error=record.get('error'))),flush=True)
atomic_json(root/'collection-complete.json',dict(all_terminal=True,terminal_sources=2,records=records,
    formal_external_cases=0,accepted_quality_references=0,gaussian_started=False))
