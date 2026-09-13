from pathlib import Path
import argparse
import hashlib
import json
import shutil
import time
parser=argparse.ArgumentParser();parser.add_argument('--round',type=Path,required=True)
args=parser.parse_args();root=args.round.resolve();work=Path(__file__).resolve().parent
out=root/'orchestration';out.mkdir(exist_ok=False)
names=['run_numeric_collection.py','run_numeric_review_batch.py','run_view_snapshot_review.py',
       'continue_frozen_round.py','summarize_numeric_review.py','numeric_round_extrema.py',
       'summarize_view_snapshot_round.py','summarize_collection_navigation_gaps.py',
       'freeze_numeric_round.py','freeze_round_orchestration.py']
for name in names:shutil.copy2(work/name,out/name)
record=dict(created_epoch=time.time(),round_identity=json.loads((root/'manifest.json').read_text(encoding='utf-8'))['round_identity'],
            helpers={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in names})
(root/'orchestration-identity.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(directory=str(out),files=len(names),round_identity=record['round_identity'])),flush=True)
