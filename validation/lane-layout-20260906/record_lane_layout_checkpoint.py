from pathlib import Path
import hashlib
import json
import shutil
import time

workspace=Path(r'F:\Codex\2026-09-05\branch-15')
base=workspace/'outputs/general-fixes-20260906'
pilot=base/'NUM-FIX-005-pilot'
root=base/'NUM-FIX-005-original'
target=Path(r'F:\Dev\cov-native-validation-20260905\validation\lane-layout-20260906')
def read(p): return json.loads(p.read_text(encoding='utf-8'))
def atomic(p,d):
    temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');temp.replace(p)
numeric=read(pilot/'review-v1/unified-review/summary.json')
record={'snapshot_epoch':time.time(),'production_commit':'e1e118d6dc34839932c2f28a0fd8379a40e4f3b4',
    'previous_round':'NUM-FIX-004-original','previous_view_failures':4,'issue_id':'EXP-002',
    'implementation':'Shared screen/SVG/PNG interval lane packing; physical energies/y/member identities preserved; real horizontal scrolling; disjoint member hit regions.',
    'component_tests':{'ON':{'passed':35,'failed':0,'skipped':0,'evidence':str(base/'component-tests-on-1788707689780709800/ctest.xml')},
                       'OFF':{'passed':35,'failed':0,'skipped':0,'evidence':str(base/'component-tests-off-1788707777172218400/ctest.xml')}},
    'pilot':{'numerical':numeric,'view':read(pilot/'view-review-v1/summary.json')},
    'round_identity':read(root/'manifest.json')['round_identity'],
    'full_round_progress':read(root/'progress.json'),
    'full_round_directory':str(root),'pipeline_progress':str(root/'pipeline-progress.json'),
    'formal_case_passes':0,'external_molecule_count':0,'REF-001':'paused-by-user',
    'scientific_rule_changes':'None in this iteration; symmetry/local projection/topology/pi-field issues remain open.',
    'closure':'Await complete 273 collection, numerical review, actual-frame and PNG review.',
    'remaining':['Global/local symmetry scope and local metric projection','Topology and scoped interaction repairs',
        'Pi partners vs crystal-field gaps and ligand priors','Full selected-MO details/tooltip/clipping review',
        'Required invariants, stored density and ordinary-build interaction checks','Paused physical references and external molecules']}
atomic(base/'lane-layout-progress.json',record)
state_path=workspace/'outputs/cov-complete-validation-20260906/阶段状态.json'
state=read(state_path);state['lane_layout_implementation']=record;state['snapshot_epoch']=time.time();atomic(state_path,state)
target.mkdir(exist_ok=True)
target.joinpath('.gitattributes').write_text('* -text\n',encoding='utf-8')
copies={base/'lane-layout-progress.json':'status.json',root/'manifest.json':'original-frozen-manifest.json',
    pilot/'manifest.json':'pilot-frozen-manifest.json',pilot/'collection-complete.json':'pilot-collection-complete.json',
    pilot/'review-v1/unified-review/summary.json':'pilot-numerical-summary.json',
    pilot/'view-review-v1/summary.json':'pilot-view-summary.json',pilot/'view-review-v1/OLD-001.json':'pilot-view-review.json',
    base/'component-tests-on-1788707689780709800/ctest.xml':'component-tests-on.xml',
    base/'component-tests-off-1788707777172218400/ctest.xml':'component-tests-off.xml'}
for source,name in copies.items(): shutil.copy2(source,target/name)
shutil.copy2(Path(__file__),target/Path(__file__).name)
body='''## NUM-FIX-005 近邻能级布局通解与原 273 新一轮

构建源码 `e1e118d6dc34839932c2f28a0fd8379a40e4f3b4`。上一轮完整分析保留的 EXP-002 已实施共同根因修复：按实际绘制与点击范围进行横向通道分配，屏幕与 PNG/SVG 使用同一布局算法。能量纵坐标、原始成员和真实能量均未改变；不会把像素接近当作严格简并。需要额外宽度时屏幕画布支持真实横向滚动，导出画布扩大，PNG/SVG 尺寸一致。

普通版与验证版各 35 项组件检查通过，无跳过。新增检查从实际 ImGui 绘制顶点定位相邻 MO，分别通过鼠标点击验证选择，并用真实横向滚轮检查溢出通道。保留完整文字/像素验收边界。

首个 OLD-001 试运行已完成：单案例取证 12.569970 秒，独立数值复核 1.992701 秒；7 组多格式导出通过，全部线性轴成员点击通过。原库 NUM-FIX-005 已冻结并启动，全部 273 终止后自动进行整库数值与视图复核。本轮增加逐成员线性点击、实际命中范围重叠、导出线条重叠与画布尺寸检查，不降低原数值或像素阈值。

整库结果未出，EXP-002 尚未关闭。SYM-004/005、局部 S 度量投影、TOPO、POL、其余 UI 和独立物理参考仍待处理。REF-001 保持暂停，正式完整通过 0，外部分子 0。

活动目录：`F:\\Codex\\2026-09-05\\branch-15\\outputs\\general-fixes-20260906\\NUM-FIX-005-original`。`progress.json` 保存取证进度；`pipeline-progress.json` 保存顺序复核状态。不得仅因对话服务中断而重复启动已有进程。此目录只包含紧凑恢复记录；不声称本轮全部纹理/截图已上传。
'''
(target/'README.md').write_text(body,encoding='utf-8')
for path in (Path(r'F:\CalChem\COV\TO-DO\validation_20260905\ISSUES.md'),workspace/'outputs/cov-complete-validation-20260906/执行检查点.md'):
    previous=path.read_text(encoding='utf-8')
    if body.splitlines()[0] not in previous:path.write_text(previous+'\n\n'+body,encoding='utf-8')
hashes={str(p.relative_to(target)):hashlib.sha256(p.read_bytes()).hexdigest() for p in target.rglob('*')
        if p.is_file() and p.name!='files-sha256.json'}
atomic(target/'files-sha256.json',hashes)
print(json.dumps({'checkpoint':str(base/'lane-layout-progress.json'),
    'terminal_cases':record['full_round_progress'].get('terminal_cases'),'publication_directory':str(target)}))
