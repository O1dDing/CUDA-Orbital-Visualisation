"""Review every collection gap only after the complete frozen case barrier."""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json
import shlex
import time

parser=argparse.ArgumentParser();parser.add_argument('--round',type=Path,required=True);args=parser.parse_args()
root=args.round.resolve()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
barrier=load(root/'collection-complete.json');manifest=load(root/'manifest.json')
assert barrier['all_terminal'] and barrier['terminal_cases']==manifest['case_count']
assert barrier['round_identity']==manifest['round_identity']
destination=root/'collection-gap-review-v1';destination.mkdir(exist_ok=False)
cases=[];checks=[]
for record in barrier['cases']:
    if record['collection_status']=='complete':continue
    evidence=root/record['evidence_directory']
    inventory_path=root/record['evidence_manifest'];assert sha(inventory_path)==record['evidence_manifest_sha256']
    inventory={Path(k).as_posix():v for k,v in load(inventory_path).items()}
    consumed={}
    def consume(name):
        path=evidence/name
        actual=sha(path);assert actual==inventory[Path(name).as_posix()]['sha256']
        consumed[name]=actual
        return path
    plan=[shlex.split(line) for line in consume('native/plan.txt').read_text(encoding='utf-8').splitlines()[1:]
          if line and not line.startswith('#')]
    contexts={}
    for action in record['failed_actions']:
        command=plan[action['command']]
        assert command[:2]==[action['op'],action['id']]
        capture=next((c[1] for c in plan[action['command']+1:] if c[0]=='capture'),None)
        key=(action['id'],capture)
        if key in contexts:continue
        context={'target':action['id'],'next_capture':capture,'failed_action':action}
        if capture:
            shot=load(consume('native/'+capture+'.ui.json'));consume('native/'+capture+'.png')
            target=next((t for t in shot['targets'] if t['id']==action['id']),None)
            context.update(state=shot['state'],target_capture=target)
            if action['id'].startswith('diagram.mo.'):
                mo=int(action['id'].split('.')[-1])
                drawn=next((x['data'] for x in shot['draw_trace'] if x['kind']=='draw.level' and x['data'].get('used_internal_mo')==mo),None)
                context['drawn_target']=drawn
                context['target_y_inside_drawn_clip']=bool(drawn and drawn['clip_y'][0]<=drawn['y']<drawn['clip_y'][1])
                context['actual_selection_matches_target']=shot['state']['applied_mo']==mo
        contexts[key]=context
        checks.append(context)
    cases.append({'case_id':record['case_id'],'collection_status':record['collection_status'],
                  'failed_actions':record['failed_actions'],'missing':record['missing'],'integrity_errors':record['integrity_errors'],
                  'contexts':list(contexts.values()),'consumed':consumed,'evidence_directory':str(evidence)})
result={'round_identity':manifest['round_identity'],'production_commit':manifest['git_commit'],'created_epoch':time.time(),
        'collection_barrier_sha256':sha(root/'collection-complete.json'),'review_script_sha256':sha(Path(__file__)),
        'all_case_collection_terminal':True,'case_count':manifest['case_count'],'gap_case_count':len(cases),
        'failed_actions':sum(len(c['failed_actions']) for c in cases),'unique_target_capture_contexts':len(checks),
        'missing_files':sum(len(c['missing']) for c in cases),'integrity_errors':sum(len(c['integrity_errors']) for c in cases),
        'action_details':dict(Counter(a['detail'] for c in cases for a in c['failed_actions'])),
        'contexts_with_target_y_inside_clip':sum(c.get('target_y_inside_drawn_clip',False) for c in checks),
        'contexts_target_reported_invisible':sum(c.get('target_capture',{}).get('visible') is False for c in checks),
        'contexts_selection_did_not_change_to_target':sum(c.get('actual_selection_matches_target') is False for c in checks),
        'cases':cases,'formal_case_passes':0,'algorithm_changed':False,
        'scope':'Collection navigation gap analysis; numerical and complete export review still run independently.'}
(destination/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='cases'},ensure_ascii=False,indent=2))
