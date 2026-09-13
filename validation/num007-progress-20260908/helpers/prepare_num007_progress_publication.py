"""Archive concrete completed evidence and the new frozen full-round identity."""
from pathlib import Path
import ctypes
import hashlib
import json
import shutil
import sys
import time
workspace=Path(r'F:\Codex\2026-09-05\branch-15');repo=Path(r'F:\Dev\cov-native-validation-20260905')
base=workspace/'outputs/cov-complete-validation-20260906';fixes=workspace/'outputs/general-fixes-20260906'
sys.path.insert(0,str(repo/'tests'))
from validation_process import physical_core_masks
k=ctypes.WinDLL('kernel32',use_last_error=True);k.GetCurrentProcess.restype=ctypes.c_void_p;k.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert k.SetProcessAffinityMask(k.GetCurrentProcess(),physical_core_masks(12)[8])
read=lambda p:json.loads(p.read_text(encoding='utf-8'));sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
target=base/'publication-prepared-num007';target.mkdir(exist_ok=False)
references=[]
def copy(source,relative):
    source=Path(source);destination=target/relative;destination.parent.mkdir(parents=True,exist_ok=True)
    if destination.exists():assert sha(destination)==sha(source);return
    shutil.copy2(source,destination)
    references.append(dict(relative=str(relative).replace('\\','/'),local_source=str(source),source_sha256=sha(source)))
def selected_tree(source,prefix,suffixes):
    for path in sorted(source.rglob('*')):
        if path.is_file() and path.suffix.lower() in suffixes and '__pycache__' not in path.parts:
            copy(path,Path(prefix)/path.relative_to(source))
for name in ('阶段状态.json','执行检查点.md','external-fifty-starting-geometries-v1.json','external-fifty-identity-screen-v1.json'):
    copy(base/name,'campaign/'+name)
copy(Path(r'F:\CalChem\COV\TO-DO\validation_20260905\ISSUES.md'),'campaign/ISSUES.md')
for name in ('density-details-fifty-checkpoint-v1','num007c-launch-checkpoint','external-fifty-formula-collision-review-v1'):
    selected_tree(base/name,'campaign/'+name,{'.json','.py','.md'})
for name in ('rcsb-and-cisplatin-coordinate-preparation-v1','mor41-geometry-preparation-v2/selection-v1'):
    selected_tree(base/name,'external/'+name,{'.json','.xyz'})
for name in ('regression-fixtures/density-source-v1/baseline-production-21036b2-v3',
             'regression-fixtures/density-source-v1/fixed-production-e2e835a7-1788854219897649100'):
    selected_tree(base/name,'density/'+Path(name).name,{'.json','.log'})
for name in ('NUM-FIX-006-original','NUM-FIX-007-pilot','NUM-FIX-007B-pilot','NUM-FIX-007C-pilot','NUM-FIX-007C-diagnostic'):
    root=fixes/name;assert read(root/'collection-complete.json')['all_terminal']
    for item in ('manifest.json','orchestration-identity.json','collection-complete.json','review-v1/review-complete.json',
                 'review-v1/unified-review/summary.json','review-v1/unified-review/extrema.json','view-review-v1/summary.json'):
        copy(root/item,Path('rounds')/name/item)
    if name=='NUM-FIX-006-original':
        copy(root/'unified-view-numeric-review/summary.json',Path('rounds')/name/'unified-view-numeric-review/summary.json')
    else:
        cid=read(root/'manifest.json')['cases'][0]['case_id'];terminal=read(root/'cases'/cid/'terminal.json')
        copy(root/'cases'/cid/'terminal.json',Path('rounds')/name/cid/'terminal.json')
        copy(root/'view-review-v1'/(cid+'.json'),Path('rounds')/name/'view-review-v1'/(cid+'.json'))
        evidence=root/terminal['evidence_directory'];native=evidence/'native'
        for item in ('case.plan','expected-evidence.json'):
            copy(evidence/item,Path('rounds')/name/cid/item)
        for item in ('actions.jsonl','session.json'):
            copy(native/item,Path('rounds')/name/cid/'native'/item)
        for path in native.glob('orbital-details-*.ui.json'):copy(path,Path('rounds')/name/cid/'native'/path.name)
        pics=['orbital-details-bottom.png','orbital-details-reopened.png']
        if cid=='OLD-013':pics+=['member-0058-tooltip.png']
        for item in pics:copy(native/item,Path('rounds')/name/cid/'native'/item)
current=fixes/'NUM-FIX-007C-original'
for item in ('manifest.json','orchestration-identity.json'):copy(current/item,Path('rounds/NUM-FIX-007C-original')/item)
for name in ('density-preflight-1788853484469266400','tooltip-preflight-1788854415993140600',
             'native-evidence-preflight-1788855763189363000','native-evidence-preflight-1788856174260660400',
             'native-evidence-preflight-1788856542626153100'):
    for item in ('identity.json','complete.json','process.json'):copy(fixes/name/item,Path('preflight')/name/item)
for name in ('details-plan-preflight-1788855447276760300','details-plan-preflight-1788855735271638400',
             'density-matrix-checker-preflight-1788856141656420600'):
    copy(fixes/name/'summary.json',Path('preflight')/name/'summary.json')
for name in ('checkpoint_num007c_launch.py','checkpoint_density_details_and_fifty.py','resolve_fifty_formula_collisions.py',
             'check_details_plan_coverage.py','check_density_comparison_v1.py','prepare_num007_progress_publication.py'):
    copy(workspace/'work'/name,Path('helpers')/name)
(target/'README.md').write_text('''# COV 进度证据：密度通解、轨道详情与 NUM-FIX-007C

这是正在执行的验证任务的证据检查点，正式案例完整通过数仍为 **0**，不是最终可交付版本。

- NUM-FIX-006 原 273 案例全部终止，31,871 个 MO 数值子集、1,943 组导出及 2,771 个线性轴成员目标通过；上一轮 24 案例的 72 个失败目标均恢复。原问题台账及限定范围保留。
- 密度旧库 20 个控制中 12 个失败；共同根因修复后，普通版和验证版分别 20/20 通过，同一检查器及 2e-12 控制矩阵界限保持不变。数学控制不计入分子数。
- NUM-FIX-007 和 007B 的试运行失败证据保留。007C 的首个案例 OLD-001 及 OLD-013 诊断全部冻结子检查通过，包含实际 P/Q、来源/可用状态、完整显示的详情控件及真实关闭动作。
- 007C 程序源码 aaabc3edb13c2768f6bc9063ee933a8881e7f57f，取证/复核源码 93e294b0bca5d920494cfbb1ab76f4c218aaf898；同一版本的全 273 轮已冻结启动。该轮须全部终止后统一分析，此包不声称它已经通过。
- 50 个外部起始候选已按结构去重，彼此及相对原库均不同。结构身份确认不等于几何、电子态或高质量参考验收；正式外部案例仍为 0。

本目录保留选定报告、动作记录、实帧、来源事实及各自 SHA-256。完整原生网格、二进制程序和全库原始文件仍在本地不可变证据目录，未宣称包含在此包内。Gaussian 由用户的独立 Resume 队列管理，本轮 COV 调度不控制该队列；COV 使用保留的 4 核、64 GiB 配额。

外部坐标来源与适用范围见各 source/state 报告。Cisplatin 采用 Johnston 等 2012 小分子晶体结构的坐标数据，DOI 10.1107/S1600536812024014；RCSB CCD 的实验模型和计算 ideal 坐标严格分开。新增 MOR41 候选来自 Dohm 等 2018 作者数据，DOI 10.1021/acs.jctc.7b01183，固定提交 089ca8ea4b350e2024a2ed013aed4927cda72b9f；仅包含选定坐标事实与来源记录，不分发完整作者代码库。先前 44 个候选的来源、原始失败响应和许可证留在相邻 input/density preparation 检查点。
''',encoding='utf-8',newline='\n')
(target/'.gitattributes').write_text('* -text\n',encoding='ascii',newline='\n')
(target/'source-records.json').write_text(json.dumps(dict(created_epoch=time.time(),sources=references),ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
files=sorted(p for p in target.rglob('*') if p.is_file())
assert not any(p.suffix.lower() in {'.exe','.dll','.pdf','.html','.npz','.pyc'} for p in files)
size=sum(p.stat().st_size for p in files);assert size<64*1024*1024
manifest={str(p.relative_to(target)).replace('\\','/'):sha(p) for p in files}
(target/'files-sha256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(prepared=str(target),advertised_files=len(files),bytes=size,contains_final_acceptance=False),ensure_ascii=False),flush=True)
