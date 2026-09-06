from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

workspace=Path(r'F:\Codex\2026-09-05\branch-15')
base=workspace/'outputs/general-fixes-20260906'
root=base/'NUM-FIX-004-original'
result=json.loads((root/'unified-view-numeric-review/summary.json').read_text(encoding='utf-8'))
assert result['case_count']==273
target=Path(r'F:\Dev\cov-native-validation-20260905\validation\view-snapshots-20260906')
subprocess.run([sys.executable,workspace/'work/publish_view_snapshot_checkpoint.py'],check=True)
for source,name in ((root/'collection-complete.json','original-collection-complete.json'),
                    (root/'review-v1/review-complete.json','original-numeric-review-complete.json'),
                    (root/'review-v1/unified-review','original-numeric-summary'),
                    (root/'view-review-v1','original-view-review'),
                    (root/'unified-view-numeric-review','original-unified-review')):
    destination=target/name
    if source.is_dir(): shutil.copytree(source,destination,dirs_exist_ok=True)
    else: shutil.copy2(source,destination)
for name in ('summarize_view_snapshot_round.py','record_view_snapshot_checkpoint.py',
             'publish_view_round_completion.py'):
    shutil.copy2(workspace/'work'/name,target/'workflow'/name)
reports=target/'original-numeric-case-reports';reports.mkdir(exist_ok=True)
terminal=json.loads((root/'review-v1/review-complete.json').read_text(encoding='utf-8'))
for case in terminal['cases']:
    if case.get('review'): shutil.copy2(root/'review-v1'/case['review'],reports/(case['case_id']+'.json'))

body='''## NUM-FIX-004 全 273 复核完成与重叠布局问题

程序构建 `70e72ef9d4c956b374a36efb4a478489127bca90`，轮身份 `79a352a6ef9de3413505f2c8c25c06561d6e5da2086240a8f836cd4da3b30693`。全部取证及复核进程已结束；无 Gaussian 新作业。

- 273/273 取证完整，无崩溃、超时或证据缺项。273/273 通过本轮独立数值子集，31,871 个 MO 的实际纹理采样及 1,214 份完整前线纹理均通过。普通版/验证版的完整科学输出与实际 EXE 数值证据在全库相符。
- 多状态视图复核完成 273/273。1,943 组四格式导出均通过快照身份、实际成员、原始能量、电子箭头和导出帧关联检查；其中 4 组另有选中 PNG 内部线条像素失败。因此视图案例状态为 269 通过、4 未通过，不能写成全库导出已全部通过。
- 新增 **EXP-002：近邻能级的线条及点击区域重叠**。OLD-071、OLD-175、OLD-176、OLD-249 的 export-linear 状态中，后绘制能级覆盖了先绘制的选中线条内部。OLD-071 的蓝色外缘仍可见，不将它描述成整个选中标记消失。其四行图仅有两处可分辨的线条位置，近邻行的命中区域也相互重叠。
- 根因：屏幕和导出只把能量四舍五入到 1e-12 Ha 后完全相同的行横向错开。实际相差 1e-5 Ha 左右的行在线性轴上可能远小于一个像素；它们既不会归入同桶，也不能在原横向位置分辨。新修复将根据实际绘制/命中范围分配横向通道，保持原始能量、行分组和纵坐标，不把近邻强制当作严格简并。该布局修复尚未实施/验证。
- EXP-001 的“导出重新构建而丢失当前视图”根因在本版本的全库多状态检查中未再出现；保留 EXP-002 和完整图像/文字验收，按子项裁决。SYM-004/005、POL、TOPO、其余 UI 与物理参考问题仍保留。
- 总取证墙钟 2182.597 秒；独立数值复核 1090.602 秒；视图复核 24.977 秒。单案例取证中位数 25.455 秒、P95 60.295 秒；数值复核中位数 12.481 秒、P95 34.684 秒。
- OLD-015 萘本轮取证 41.415758 秒，独立数值复核 24.582416 秒，合计 65.998174 秒。这不包含 Gaussian 参考计算及尚未完成的化学/完整界面验收，也没有把批量视图复核的 24.977 秒均分成未经测量的单案例耗时。
- REF-001 按用户指示暂停；正式完整通过 0，新增外部分子 0。所有失败与冻结证据保留。对话服务曾出现的 Bad Request 不属于本轮 COV/Gaussian 科学失败。

证据：`original-numeric-summary`、`original-numeric-case-reports`、`original-view-review`、`original-unified-review`。其中每个视图报告均记录消耗证据的哈希。仍需完成后续不依赖 REF 的通解修复及全库循环。
'''
(root/'unified-view-numeric-review/统一复核与裁决.md').write_text(body,encoding='utf-8')
(target/'original-unified-review/统一复核与裁决.md').write_text(body,encoding='utf-8')
for path in (Path(r'F:\CalChem\COV\TO-DO\validation_20260905\ISSUES.md'),
             workspace/'outputs/cov-complete-validation-20260906/执行检查点.md'):
    old=path.read_text(encoding='utf-8')
    if body.splitlines()[0] not in old: path.write_text(old+'\n\n'+body,encoding='utf-8')
(target/'README.md').write_text('# NUM-FIX-004 当前视图快照修复与整库复核\n\n'+body+
    '\n源码、流程、进度、全部紧凑复核报告已在本目录归档。大体积原始纹理、截图和网格仍在本地冻结目录，未声称它们已随本次 Git 提交上传。此前 REF 暂停交接 release 保持完整。\n',encoding='utf-8')
receipt={str(p.relative_to(target)):hashlib.sha256(p.read_bytes()).hexdigest()
         for p in target.rglob('*') if p.is_file() and p.name!='publication-files-sha256.json'}
(target/'publication-files-sha256.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'directory':str(target),'files':len(receipt),'formal_case_passes':0,
                  'view_subset_failures':4,'github':'pending successful push verification'}))
