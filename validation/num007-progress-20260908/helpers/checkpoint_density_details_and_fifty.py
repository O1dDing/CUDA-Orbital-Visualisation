from pathlib import Path
import hashlib
import json
import time
base=Path(r'F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906')
fixes=base.parent/'general-fixes-20260906'
read=lambda p:json.loads(p.read_text(encoding='utf-8'));sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
fixed=base/'regression-fixtures/density-source-v1/fixed-production-e2e835a7-1788854219897649100/review-complete.json'
density=read(fixed);assert density['total_collected_controls']==40 and density['all_fixed_requirements_satisfied'] and density['ordinary_validation_identical']
external=base/'external-fifty-formula-collision-review-v1/summary.json';identities=read(external)
assert identities['nonidentity_to_originals_established']==50 and identities['all_fifty_pairwise_constitutionally_distinct']
pilot=fixes/'NUM-FIX-007-pilot';barrier=read(pilot/'collection-complete.json')
assert barrier['all_terminal'] and barrier['collection_counts']['completed_with_gaps']==1
out=base/'density-details-fifty-checkpoint-v1';out.mkdir(exist_ok=False)
record=dict(created_epoch=time.time(),density_fixed_control_count=40,ordinary_validation_density_controls_identical=True,
    fixed_density_review=str(fixed),fixed_density_review_sha256=sha(fixed),
    density_issue_decisions={f'DEN-{i:03d}':'已实施通解且固定对照通过；全 273 回归前保留待关闭' for i in range(1,6)},
    tooltip_implementation_commit='6b8b99b9f76c8d1f9d8c30d2e27c21796c2ee4b8',
    native_action_and_clip_fix_commit='20a326365350a356356726e1e9e262dec7f79917',
    actual_density_crosscheck_commit='8b2522716ec8c391cfc744699048b31c2b6d2538',
    first_pilot_round_identity=barrier['round_identity'],first_pilot_status='保留已完成且有缺口的试运行；不得按子检查通过改写为全通过',
    first_pilot_failed_commands=2,first_pilot_collection_seconds=barrier['wall_seconds_this_session'],
    external_constitutionally_distinct_candidates=50,external_identity_evidence=str(external),external_identity_sha256=sha(external),
    formal_external_cases=0,accepted_external_quality_references=0,formal_original_complete_passes=0,
    gaussian_control_actions_performed=0,next='Complete final evidence build, repeat the designated pilot and long-tooltip diagnostic, then freeze a full 273 round with density and actual details evidence.')
atomic=lambda p,d:p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
atomic(out/'summary.json',record)
title='## 密度通解固定对照通过、详情试运行缺口与 50 个外部结构去重'
note='''

密度通解 e2e835a7 已用同一冻结检查器复验普通版 20 项、验证版 20 项。全部 40 项先进入终态，再统一比较；40/40 满足原固定要求，20 组 ON/OFF 实际输出完全一致，原矩阵界限 2e-12 未改变。DEN-001 至 DEN-005 更新为“已实施通解且固定对照通过，待完整 273 回归关闭”，没有把数学控制算作新增分子。

UI-004 已实施短提示与可滚动详情窗口。NUM-FIX-007-pilot 的 OLD-001 全部采集结束，取证 24.196287600003416 秒；数值及导出子集通过，但原生操作有 2 条失败，整次试运行保留为有缺口。真实关闭已生效，自动操作器却继续要求已消失的关闭按钮存在，随后误报目标未绘制。新增 EVID-ACT-001：真实关闭动作的完成判定不应依赖目标继续存在；20a3263 已修复，待后续实际试运行及全库检查。

新增 EVID-CLIP-001：可点击的部分可见目标不能证明说明文字完整可读。已查看实际详情顶部和底部 PNG；第一版末尾说明仍部分裁切。后续计划先通过真实滚轮到达底部关闭按钮再截图，并保存真实 ImGui 裁切矩形；独立检查器核对完整矩形及全部失败动作。该问题保留原反例，不用自报 visible=true 关闭裁切问题。UI-004 此时仍不关闭。

8b25227 增加普通/验证构建共同的实际 P/Q 序列化及独立矩阵检查：从同一 FCHK 的 IOData 系数/占据、原始 producer SCF 密度和独立 AO 约定建立参考；输入矩阵变换要求精确一致，MO 重建只允许由运算次数决定的双精度误差，不借用十进制打印精度。新检查器的 10 项正反例通过，旧试运行实测密度的 17 项诊断通过；这是冻结下一轮前的检查器准备，不改写旧轮判据。

50 个外部起始结构的全部同分子式冲突已解决：噁唑/异噁唑的 O–N 连接数不同；苄基自由基的六元环加支链与 OLD-097 七元环的碳连接度序列不同。各自已冻结的四档距离检查与作者/PubChem 连接证据一致。全部 50 个彼此不同且与原库不同；尚未完成参考质量和电子态验收，正式外部案例仍为 0。

Gaussian 继续由用户的 Resume 流程独立管理。本检查点没有启动、重启、终止该队列，也没有修改其暂停文件。COV 工作沿用物理核槽 8–11，最多 4 核和 64 GiB 的进程树限制。
'''
for path in (Path(r'F:\CalChem\COV\TO-DO\validation_20260905\ISSUES.md'),base/'执行检查点.md'):
    prior=path.read_text(encoding='utf-8');assert title not in prior
    path.write_text(prior+'\n\n'+title+note,encoding='utf-8',newline='\n')
state=read(base/'阶段状态.json');state['density_details_fifty_checkpoint_v1']=record;state['snapshot_epoch']=record['created_epoch']
temp=base/'阶段状态.json.tmp';atomic(temp,state);temp.replace(base/'阶段状态.json')
print(json.dumps({'output':str(out),'fixed_density_controls':40,'distinct_external_starting_structures':50,'formal_case_passes':0},ensure_ascii=False))
