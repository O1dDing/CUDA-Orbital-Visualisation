"""Archive the next density control and public starting-coordinate preparation."""
from pathlib import Path
import ctypes
import hashlib
import json
import shutil
import sys
import time

workspace=Path(r'F:\Codex\2026-09-05\branch-15')
repo=Path(r'F:\Dev\cov-native-validation-20260905')
base=workspace/'outputs/cov-complete-validation-20260906'
numeric=workspace/'outputs/general-fixes-20260906'
sys.path.insert(0,str(repo/'tests'))
from validation_process import atomic_json,physical_core_masks
kernel=ctypes.WinDLL('kernel32',use_last_error=True)
kernel.GetCurrentProcess.restype=ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes=[ctypes.c_void_p,ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(),physical_core_masks(12)[8])
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
controls=base/'regression-fixtures/density-source-v1'
baseline=controls/'baseline-production-21036b2-v3'
math_test=base/'regression-fixtures/density-print-precision-math-1788846604407094400'
reviewer_test=base/'regression-fixtures/density-reviewer-preflight-1788847493973793400'
for test in (math_test,reviewer_test):assert read(test/'process.json')['exit_code']==0
source_sets=[
    {r['candidate_id'] for r in read(base/'external-coordinate-preparation-v2/summary.json')['records'] if r['status']=='computed-starting-conformer-retrieved'},
    {r['candidate_id'] for r in read(base/'benchmark-geometry-preparation-v1/review-v2/summary.json')['records']},
    {r['candidate_id'] for r in read(base/'nist-geometry-preparation-v1/starting-geometries-v1/summary.json')['records'] if r['status']=='literal-starting-geometry-prepared'},
    {r['candidate_id'] for r in read(base/'nist-followed-geometries-v1/starting-geometries-v1/summary.json')['records']}]
assert [len(s) for s in source_sets]==[35,4,3,2]
assert len(set.union(*source_sets))==44
record=dict(snapshot_epoch=time.time(),phase='density-baseline-queued-and-input-preparation',
    candidates_with_archived_starting_coordinates=44,starting_coordinate_candidate_ids=sorted(set.union(*source_sets)),
    formal_external_cases=0,accepted_quality_references=0,formal_case_passes=0,
    independent_print_precision_tests=dict(tests=7,passed=7,evidence=str(math_test)),
    independent_control_reviewer_tests=dict(tests=7,passed=7,evidence=str(reviewer_test)),
    production_density_controls=20,production_density_baseline=str(baseline),
    production_density_queue_at_snapshot=read(baseline/'post-num006-pipeline/progress.json'),
    running_round_progress_at_snapshot=read(numeric/'NUM-FIX-006-original/progress.json'),
    previous_publication=read(base/'github-input-preparation-git-bytes.json'),
    independent_gaussian='User-managed F:/CalChem/COV/Resume; queue/pause request left untouched',
    gaussian_started_by_this_work=False,production_algorithm_changed=False)
title='## 2026-09-08 密度基线排队及 44 份起始坐标检查点'
note=title+'''

现有 67 个预备候选中，44 个已具可追溯起始坐标：35 份 PubChem 计算构象、4 份作者 GMTKN55 基准几何、5 份 NIST 表格/可视化载荷。新取得的 ClF3 几何明确为 C2v；甲基过氧自由基为 Cs 双重态。它们的 G3B3 页面明确使用 B3LYP/6-31G* 几何。只解码页面中的字面坐标字符串，没有执行外部 JavaScript，也没有对坐标进行拟合。IF7 的已取页面没有可用计算几何行，继续留为待补来源。

独立密度打印误差工具的 7 项测试和密度检查器的 7 项防误通过测试均通过。20 个真实生产库对照仍在等待 NUM-FIX-006 全部终止和自动复核，尚未取得生产密度基线结论。当前有效准备目录为 density-source-v1/baseline-production-21036b2-v3；44 份文件的哈希固定，包含其进程监督依赖。v1/v2 均未执行并保留历史。

前一准备档案已推送 GitHub 并核对远端提交 060c02e2ce64bfd364d0b2fa795eb8284e57fd0d，285 个 Git blob 的 SHA256 全部一致。本检查点没有修改 NUM-FIX-006 的程序、输入或判据。独立 Gaussian 队列仍由用户控制；完整正式通过数和正式外部分子数均为 0。
'''
target=repo/'validation/density-preparation-20260908'
if target.exists():raise FileExistsError('Retain the existing publication checkpoint')
target.mkdir(parents=True)
mapping={}
def copy(source,relative):
    destination=target/relative;destination.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(source,destination);assert sha(source)==sha(destination)
    mapping[str(source)]=relative.replace('\\','/')
for folder in ['nist-geometry-preparation-v1','nist-calculated-geometry-links-v1','nist-followed-geometries-v1']:
    for source in sorted((base/folder).rglob('*')):
        if source.is_file() and source.suffix not in ('.html','.pyc') and '__pycache__' not in source.parts:
            copy(source,str(source.relative_to(base)))
for test in (math_test,reviewer_test):
    for source in sorted(test.rglob('*')):
        if source.is_file() and source.suffix not in ('.pyc',) and '__pycache__' not in source.parts:
            copy(source,str(source.relative_to(base)))
copy(baseline/'identity.json','density-baseline/identity.json')
for source in sorted((baseline/'source').rglob('*')):
    if source.is_file() and source.suffix not in ('.pyc',) and '__pycache__' not in source.parts and 'include' not in source.relative_to(baseline/'source').parts:
        copy(source,'density-baseline/'+str(source.relative_to(baseline)))
queue=baseline/'post-num006-pipeline'
copy(queue/'identity.json','density-baseline/post-num006-pipeline/identity.json')
for source in (queue/'source').iterdir():
    if source.is_file():copy(source,'density-baseline/post-num006-pipeline/source/'+source.name)
copy(numeric/'density-shared-interface-v4.md','density-shared-interface-v4.md')
copy(base/'github-input-preparation-git-bytes.json','previous-publication-verification.json')
copy(workspace/'work/checkpoint_density_preparation_v2.py','helpers/checkpoint_density_preparation_v2.py')
(target/'.gitattributes').write_text('* -text\n',encoding='utf-8',newline='\n')
(target/'README.md').write_text(note+'''

此目录仅为准备证据。完整源网页、生产库二进制及头文件快照仍在本机保留，没有再次发布；生产源码由基线提交定位。`local-to-repository-paths.json` 提供原始证据绝对路径到此处副本的映射。排队器内的本机路径是已执行环境的身份记录，不是通用安装路径。
''',encoding='utf-8',newline='\n')
atomic_json(target/'status.json',record)
atomic_json(target/'local-to-repository-paths.json',mapping)
atomic_json(target/'files-sha256.json',{str(p.relative_to(target)).replace('\\','/'):sha(p) for p in sorted(target.rglob('*')) if p.is_file() and p.name!='files-sha256.json'})
state=read(base/'阶段状态.json');state['density_preparation_20260908_v2']=record;state['snapshot_epoch']=record['snapshot_epoch']
atomic_json(base/'阶段状态.json',state)
checkpoint=base/'执行检查点.md';text=checkpoint.read_text(encoding='utf-8')
assert title not in text
checkpoint.write_text(text+'\n\n'+note,encoding='utf-8',newline='\n')
print(json.dumps(dict(status='preparation-checkpoint-written',directory=str(target),advertised_files=len(read(target/'files-sha256.json')),candidate_coordinates=44,formal_case_passes=0)),flush=True)
