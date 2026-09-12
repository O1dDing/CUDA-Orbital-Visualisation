"""Extract literal source tables without fitting or inventing geometry data."""
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
import hashlib
import json
import math
import re
import urllib.parse

root=Path(__file__).resolve().parent.parent/'outputs/cov-complete-validation-20260906/nist-geometry-preparation-v1'
collection=json.loads((root/'collection-complete.json').read_text(encoding='utf-8'))
assert collection['all_terminal'] and collection['terminal_sources']==5
out=root/'table-review-v1'
if out.exists():raise FileExistsError('Retain the existing table extraction')
out.mkdir()
(out/'extractor.py').write_bytes(Path(__file__).read_bytes())

class Tables(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables=[];self.stack=[];self.text=[];self.links=[];self.anchor=None;self.headings=[];self.heading=None
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if tag=='table':
            self.stack.append(dict(index=len(self.tables),rows=[],row=None,cell=None,
                                   preceding_text=' '.join(self.text)[-500:]))
            self.tables.append(self.stack[-1])
        if self.stack:
            if tag=='tr':self.stack[-1]['row']=[]
            if tag in ('td','th'):self.stack[-1]['cell']=[]
        if tag=='a' and 'href' in attrs:self.anchor=dict(href=attrs['href'],text=[])
        if tag in ('h1','h2','h3'):self.heading=dict(tag=tag,text=[])
    def handle_data(self,value):
        clean=' '.join(value.split())
        if clean:self.text.append(clean)
        if self.stack and self.stack[-1]['cell'] is not None:self.stack[-1]['cell'].append(value)
        if self.anchor is not None:self.anchor['text'].append(value)
        if self.heading is not None:self.heading['text'].append(value)
    def handle_endtag(self,tag):
        if self.stack:
            table=self.stack[-1]
            if tag in ('td','th') and table['cell'] is not None:
                if table['row'] is not None:table['row'].append(' '.join(''.join(table['cell']).split()))
                table['cell']=None
            if tag=='tr' and table['row'] is not None:
                table['rows'].append(table['row']);table['row']=None
            if tag=='table':self.stack.pop()
        if tag=='a' and self.anchor is not None:
            self.anchor['text']=' '.join(''.join(self.anchor['text']).split());self.links.append(self.anchor);self.anchor=None
        if tag in ('h1','h2','h3') and self.heading is not None:
            self.heading['text']=' '.join(''.join(self.heading['text']).split());self.headings.append(self.heading);self.heading=None

records=[]
for item in collection['records']:
    record=dict(candidate_id=item['candidate_id'],status=item['status'],source_url=item['url'],
                source_scope=item['scope'],formal_case=False,accepted_reference_geometry=False)
    if item['status']=='html-retrieved-content-review-pending':
        source=root/item['candidate_id']/'source.html';raw=source.read_bytes()
        assert hashlib.sha256(raw).hexdigest()==item['sha256']
        document=raw.decode(item.get('charset') or 'utf-8',errors='strict')
        parsed=Tables();parsed.feed(document)
        geometries=[]
        for table in parsed.tables:
            rows=[]
            for cells in table['rows']:
                if len(cells)!=4:continue
                match=re.fullmatch(r'([A-Z][a-z]?)([1-9]\d*)',cells[0])
                if not match:continue
                try:xyz=[float(x) for x in cells[1:]]
                except ValueError:continue
                if not all(math.isfinite(x) for x in xyz):raise ValueError('Nonfinite source coordinates')
                rows.append(dict(label=cells[0],element=match[1],atom_number=int(match[2]),
                                 literal_coordinates=cells[1:],coordinates_angstrom=xyz))
            if not rows:continue
            count=dict(Counter(row['element'] for row in rows))
            ids=[row['atom_number'] for row in rows]
            if ids!=list(range(1,len(rows)+1)):raise ValueError('Review nonsequential source atom identifiers')
            min_distance=min(math.dist(a['coordinates_angstrom'],b['coordinates_angstrom']) for i,a in enumerate(rows) for b in rows[:i])
            geometries.append(dict(source_table_index=table['index'],preceding_text=table['preceding_text'],
                rows=rows,element_counts=count,expected_element_counts_match=count==item['elements'],
                minimum_pair_distance_angstrom=min_distance,
                units='angstrom as explicitly headed in source table',
                precision='Literal printed decimals retained; no extra accuracy implied or coordinate fitting performed.'))
        useful_links=[dict(text=l['text'],url=urllib.parse.urljoin(item['resolved_url'],l['href']))
                      for l in parsed.links if l['text'] and any(key in (l['href']+' '+l['text']).lower()
                      for key in ['casno=','cartesian','more geometry','input/output','geom3','doi.org'])]
        detailed=dict(candidate_id=item['candidate_id'],source_url=item['url'],raw_sha256=item['sha256'],
                      headings=parsed.headings,tables=[{k:v for k,v in t.items() if k not in ['row','cell']} for t in parsed.tables],
                      related_links=useful_links,coordinate_tables=geometries,
                      formal_case=False,gaussian_started=False,accepted_reference_geometry=False)
        target=out/(item['candidate_id']+'-tables.json');target.write_text(json.dumps(detailed,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        record.update(status='literal-tables-extracted-awaiting-state-and-source-review',
                      coordinate_table_count=len(geometries),headings=parsed.headings,
                      coordinate_element_counts=[g['element_counts'] for g in geometries],
                      normalized_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    records.append(record)
summary=dict(all_5_sources_terminal_and_examined=True,records=records,formal_external_cases=0,
             accepted_quality_references=0,gaussian_started=False)
(out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False))
