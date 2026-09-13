from pathlib import Path
import json
import time
root=Path(r'F:\Codex\2026-09-05\branch-15')
base=root/'outputs/general-fixes-20260906'
round_root=base/'NUM-FIX-004-original'
pilot=base/'NUM-FIX-004-pilot'
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def atomic(path,data):
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    temp.replace(path)
record={
    'snapshot_epoch':time.time(),'production_commit':'70e72ef9d4c956b374a36efb4a478489127bca90',
    'production_changes_commit':'5bbc23ee5bd92a957a5c65936f48e18fdcfcdf58',
    'symmetry_changes_commit':'94447b057be12f3126f6cbef8743bb894c622144',
    'original_round_identity':read(round_root/'manifest.json')['round_identity'],
    'original_progress':read(round_root/'progress.json'),
    'pilot_collection':read(pilot/'collection-complete.json'),
    'pilot_numerical_review':{k:read(pilot/'review-v1/cases/OLD-001/attempt-001/result/review.json')[k]
        for k in ('status','wall_seconds','formal_case_pass')},
    'pilot_view_review':read(pilot/'view-review-v1/summary.json'),
    'component_tests':{'on':{'passed':34,'failed':0,'skipped':0,
        'evidence':str(base/'component-tests-on-1788702776911683800/ctest.xml')},
        'off':{'passed':34,'failed':0,'skipped':0,
        'evidence':str(base/'component-tests-off-1788702847481874000/ctest.xml')}},
    'issue_dispositions':{
        'EXP-001':'Common snapshot/member/export-frame implementation complete; 273 numerical and 1943 snapshot identity/member bundles reviewed. Four PNG overlap failures retained as EXP-002.',
        'EXP-002':'New: close but unequal energies occupy overlapping strokes and hit regions. Four selected PNG stroke failures in NUM-FIX-004; generic lane packing pending.',
        'SYM-001/003/006':'Component-verified fixes included in a numerically passing full round; independent interpretation closure remains separate.',
        'SYM-004/005':'Global/local label-scope separation still open.',
        'UI-003/006/007':'Export now retains actual counterpart/member identities and energy ranges; full interactive semantics still open.',
        'UI-004/005/008':'Retained. A scrollable canvas is not itself a clipping defect; unscrollable tooltip and complete pixel checks remain separate.'},
    'formal_case_passes':0,'external_molecule_count':0,'REF-001':'paused-by-user',
    'original_numerical_review':read(round_root/'review-v1/unified-review/summary.json'),
    'original_view_review':read(round_root/'view-review-v1/summary.json'),
    'unified_review':str(round_root/'unified-view-numeric-review/summary.json'),
    'remaining':['Repair unequal-energy stroke/hit-region overlap found by the complete view review.',
        'Global/local symmetry namespaces and missing finite-group derivation.',
        'Topology, interaction scope and ligand-prior common-root fixes.',
        'Full text clipping/interaction/ordinary-build behavior and invariants.',
        'Paused physical reference calculations and external molecule campaign.'],
    'github_verified_prior_head':'57faaa10b09f4d92e834ed317d5da7bcb72384f2',
    'github_current_source_push':read(base/'github-source-receipt-70e72ef.json')
        if (base/'github-source-receipt-70e72ef.json').exists()
        else 'pending receipt; do not infer upload from local commits'}
atomic(base/'view-snapshot-progress.json',record)
stage=root/'outputs/cov-complete-validation-20260906/阶段状态.json'
data=read(stage);data['view_snapshot_implementation']=record;data['snapshot_epoch']=time.time();atomic(stage,data)
heading='## NUM-FIX-004 当前视图快照与全量回归'
body=f'''{heading}

源码构建身份 `70e72ef9d4c956b374a36efb4a478489127bca90`，包含上一阶段对称性修复及当前视图快照通解。验证版与普通版均已编译，34 项组件检查各自全部通过。

- EXP-001 已实施：真实导出按钮传递已经绘制的不可变快照；原始 MO、决定行选择的 anchor 和当前实际检查的 MO 分开记录。后续控件变化不会回写旧快照。
- PNG/SVG 使用完整实际成员及电子箭头，β 对应关系和选中成员保留；JSON/CSV 记录实际成员原始值、行布局值、成员能量范围及可见行映射。机读浮点输出保留 double 往返精度。
- COV_VALIDATION 1 保持兼容；新增只控制证据文件名的 export-name 指令，导出仍经过真实按钮。当前轮分别保留至少 7 个视图状态，含 β 对应成员的案例额外取证。
- 首个 OLD-001 完整取证 11.7034204 秒，独立数值检查 1.9584581 秒；7 组导出视图检查全部通过。该项仅通过相应子集，不等于完整化学或界面验收。
- 已查看真实 PNG 和 compact-top 截图。能级图仍不能在当前卡片内一次完整显示，摘要也有裁切；保留 UI 的相关裁决。绘制调用、选中 PNG 线条探针和快照身份检查不能代替全部文字/像素验收。
- 已启动原 273 全量取证，轮身份 `{record['original_round_identity']}`。所有案例终止后再做本轮统一复核；批次中不修改生产代码、输入或阈值。
- REF-001 继续暂停；正式完整通过 0，新增外部分子 0。既有 NUM-FIX-003 的 273 数值通过证据仍属于 a9bf30a，不能归到本版。

机读进度：`{base/'view-snapshot-progress.json'}`。活动轮进度：`{round_root/'progress.json'}`。
'''
for path in (Path(r'F:\CalChem\COV\TO-DO\validation_20260905\ISSUES.md'),
             root/'outputs/cov-complete-validation-20260906/执行检查点.md'):
    current=path.read_text(encoding='utf-8')
    if heading not in current: path.write_text(current+'\n\n'+body,encoding='utf-8')
print(json.dumps({'checkpoint':str(base/'view-snapshot-progress.json'),
    'original_progress':record['original_progress'].get('terminal_cases'),
    'formal_case_passes':0},ensure_ascii=False))
