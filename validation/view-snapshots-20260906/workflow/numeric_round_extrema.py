from pathlib import Path
import hashlib
import json
import argparse

BASE=Path(r'F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906')


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


parser=argparse.ArgumentParser()
parser.add_argument('--round',type=Path,required=True)
parser.add_argument('--review-id',default='review-v1')
args=parser.parse_args()
for round_root in (args.round.resolve(),):
    name=round_root.name
    root=round_root/args.review_id
    barrier=read(root/'review-complete.json')
    assert barrier['all_terminal']
    values={}
    unavailable=[]
    def maximum(key,value,case_id,mo=None):
        if key not in values or value>values[key]['value']:
            values[key]={'value':value,'case_id':case_id,'mo_zero_based':mo}
    def minimum(key,value,case_id,mo=None):
        if key not in values or value<values[key]['value']:
            values[key]={'value':value,'case_id':case_id,'mo_zero_based':mo}
    for terminal in barrier['cases']:
        if 'review' not in terminal:
            unavailable.append(terminal['case_id'])
            continue
        case=terminal['case_id']; report_path=root/terminal['review']; report=read(report_path)
        for check in report['checks']:
            data=check['observed'];key=check['check']
            if key=='AO-INTEGRAL':maximum('AO_integral_max_absolute_error',data['max_absolute_error'],case)
            if key.startswith('MO-GRAM-'):
                for field in ('reference_gram_max_error','production_gram_max_error','production_reference_gram_max_error'):
                    maximum(field,data[field],case)
            if key=='MO-RECONSTRUCTED-ELECTRON-TRACE':
                maximum('electron_trace_error',abs(data['trace']-data['expected_electrons']),case)
            if key=='DECOMPOSITION-CONSERVATION':
                values.setdefault('available_chemistry_mos',{'count':0})['count']+=data['available_mos']
        for file,kind in (('sampled-textures.json','sampled'),('complete-textures.json','complete')):
            if not (report_path.parent/file).exists():
                unavailable.append(case+':'+file)
                continue
            for record in read(report_path.parent/file):
                if kind=='sampled' and record['full_grid']:continue
                for key in ('nrms','relative_peak_error','max_absolute_error'):
                    if key in record: maximum(kind+'_'+key,record[key],case,record['mo'])
                if 'abs_cosine' in record: minimum(kind+'_abs_cosine',record['abs_cosine'],case,record['mo'])
                if kind=='sampled' and 'points' in record:minimum('sampled_unique_points',record['points'],case,record['mo'])
    output=root/'unified-review'/'extrema.json'
    if output.exists():raise FileExistsError(output)
    result={'round_identity':barrier['round_identity'],'extrema':values,
            'review_barrier_sha256':hashlib.sha256((root/'review-complete.json').read_bytes()).hexdigest(),
            'formal_case_passes':0,'unavailable_reports':unavailable,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'round':name,**values}),flush=True)
