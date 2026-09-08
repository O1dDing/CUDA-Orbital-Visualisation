from pathlib import Path
import hashlib
import json
import time
import xml.etree.ElementTree as ET
base=Path(r'F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906');fixes=base.parent/'general-fixes-20260906'
read=lambda p:json.loads(p.read_text(encoding='utf-8'));sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
program='aaabc3edb13c2768f6bc9063ee933a8881e7f57f';round_root=fixes/'NUM-FIX-007C-original'
manifest=read(round_root/'manifest.json');assert manifest['case_count']==273 and manifest['git_commit']==program
preflight=fixes/'native-evidence-preflight-1788856542626153100';complete=read(preflight/'complete.json')
assert complete['commit']==program and complete['all_stages_passed']
components={}
for mode in ('on','off'):
    records=[json.loads(line) for line in (preflight/f'tests-{mode}.log').read_text(encoding='utf-8').splitlines() if line.startswith('{')]
    receipt=next(r for r in records if 'output' in r)
    xml=Path(receipt['output'])/'ctest.xml';cases=list(ET.parse(xml).getroot().iter('testcase'))
    assert len(cases)==37 and not any(c.find('failure') is not None or c.find('skipped') is not None or c.find('error') is not None for c in cases)
    components[mode]=dict(tests=37,failures=0,skipped=0,junit=str(xml),junit_sha256=sha(xml))
pilots=[]
for suffix in ('pilot','diagnostic'):
    root=fixes/('NUM-FIX-007C-'+suffix);m=read(root/'manifest.json');barrier=read(root/'collection-complete.json')
    numeric=read(root/'review-v1/unified-review/summary.json');view=read(root/'view-review-v1/summary.json')
    assert m['git_commit']==program and barrier['collection_counts']==dict(complete=1,completed_with_gaps=0,error=0)
    assert numeric['numeric_subset_counts']=={'numeric_subset_pass':1} and view['status_counts']=={'view_subset_pass':1}
    cid=m['cases'][0]['case_id'];terminal=read(root/'cases'/cid/'terminal.json');native=root/terminal['evidence_directory']/'native'
    assert terminal['native_session']['failed_commands']==0
    pics=['orbital-details-bottom.png','orbital-details-reopened.png'] if suffix=='pilot' else ['member-0058-tooltip.png']
    pilots.append(dict(kind=suffix,case_id=cid,round_identity=m['round_identity'],
        collection_seconds=barrier['wall_seconds_this_session'],numeric_case_seconds=numeric['timing']['independent_review_per_case']['sum_seconds'],
        view_review_seconds=view['wall_seconds'],sampled_mos=numeric['sampled_mos'],full_frontier_grids=numeric['complete_frontier_textures'],
        export_bundles=view['bundle_count'],details_frames=view['details_control_frame_count'],
        actual_images_visually_inspected={p:sha(native/p) for p in pics},formal_case_pass=False))
out=base/'num007c-launch-checkpoint';out.mkdir(exist_ok=False)
record=dict(created_epoch=time.time(),program_commit=program,runner_commit=manifest['runner_source_commit'],
    current_round=str(round_root),current_round_identity=manifest['round_identity'],manifest_sha256=sha(round_root/'manifest.json'),
    expected_original_cases=273,round_status='collection-running; unified comparison waits for all terminal states',
    preflight=str(preflight),preflight_complete_sha256=sha(preflight/'complete.json'),component_results=components,
    pilots=pilots,resource_budget=manifest['resources'],formal_case_passes=0,gaussian_control_actions_performed=0,
    previous_pilot_failures_retained=['NUM-FIX-007-pilot','NUM-FIX-007B-pilot'],
    next='Collect all 273, run frozen numeric/density/view reviewers, then adjudicate every failure together.')
def atomic(p,d):
    temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n');temp.replace(p)
atomic(out/'summary.json',record)
title='## NUM-FIX-007C：试运行与长提示诊断通过，已冻结启动全 273'
note='\n\n源码 `'+program+'`，取证/复核源码 `'+manifest['runner_source_commit']+'`。ON/OFF 各 37 项组件测试通过，0 失败、0 跳过。NUM-FIX-007 和 007B 的试运行缺口继续保留；修复关闭动作完成判定并让悬停通过真实滚轮完整显示可容纳的目标后，007C 的首个案例及 OLD-013 长提示诊断全部已冻结子检查通过。\n\n'
for row in pilots:
    note+=f"- {row['case_id']}：取证 {row['collection_seconds']:.6f} 秒，独立数值检查本体 {row['numeric_case_seconds']:.6f} 秒，视图复核 {row['view_review_seconds']:.6f} 秒；{row['sampled_mos']} 个 MO、{row['full_frontier_grids']} 个完整前线网格、{row['export_bundles']} 组导出、{row['details_frames']} 张详情控件实帧通过。\n"
note+='\n全库轮身份 `'+manifest['round_identity']+'`。273 个输入及其旧编号全部冻结，2 个 COV 进程并行，每个 2 线程；整体 4 个物理核槽、64 GiB 进程树上限。只有全体采集终止且监督进程退出后才开始冻结的独立数值、密度和视图复核。不会重跑 NUM-FIX-005，也不会在这轮中改变算法、判据或输入。\n\n实际查看了首个案例的详情底部/重新打开帧以及 OLD-013 的 MO 58 短提示帧；指定内容完整可见。此范围不代表全部字形、科学解释或所有像素已通过。DEN-001 至 DEN-005、UI-004 及本次取证问题均等待全库回归裁决；全部正式案例完整通过数仍为 0。Gaussian 由用户独立队列管理，未执行控制操作。\n'
for path in (Path(r'F:\CalChem\COV\TO-DO\validation_20260905\ISSUES.md'),base/'执行检查点.md'):
    text=path.read_text(encoding='utf-8');assert title not in text;path.write_text(text+'\n\n'+title+note,encoding='utf-8',newline='\n')
state=read(base/'阶段状态.json');state['num007c_launch']=record;state['snapshot_epoch']=record['created_epoch'];atomic(base/'阶段状态.json',state)
print(json.dumps(dict(output=str(out),round_identity=manifest['round_identity'],case_count=273,pilots=pilots),ensure_ascii=False))
