"""Archive bounded input preparation without changing the running COV round."""
from pathlib import Path
import ctypes
import hashlib
import json
import shutil
import sys
import time

workspace = Path(__file__).resolve().parent.parent
repo = Path(r'F:\Dev\cov-native-validation-20260905')
base = workspace / 'outputs/cov-complete-validation-20260906'
numeric = workspace / 'outputs/general-fixes-20260906'
sys.path.insert(0, str(repo / 'tests'))
from validation_process import atomic_json, physical_core_masks

kernel = ctypes.WinDLL('kernel32', use_last_error=True)
kernel.GetCurrentProcess.restype = ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(), physical_core_masks(12)[8])
read = lambda path: json.loads(path.read_text(encoding='utf-8'))
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
external = read(base / 'external-coordinate-preparation-v2/summary.json')
external_review = read(base / 'external-coordinate-preparation-v2/review-v1/summary.json')
benchmark = read(base / 'benchmark-geometry-preparation-v1/review-v2/summary.json')
density = read(base / 'regression-fixtures/density-source-v1/manifest.json')
assert external['all_terminal'] and external['terminal_candidates'] == 67
assert external['counts']['computed-starting-conformer-retrieved'] == 35
assert external_review['raw_identity_integrity_failures'] == 0
assert len(benchmark['records']) == 4
assert all(r['status'] == 'author-coordinate-identity-established' for r in benchmark['records'])
assert density['control_count'] == 20 and not density['cov_executed']
pubchem_ids = {r['candidate_id'] for r in external['records'] if r['status'] == 'computed-starting-conformer-retrieved'}
benchmark_ids = {r['candidate_id'] for r in benchmark['records']}
assert not pubchem_ids.intersection(benchmark_ids)

benchmark_readme = '''# 四个作者基准坐标的补充准备

从作者的 GMTKN55 仓库固定提交 `8d485b37a1ca8837e395042671ca5ba4e0714691` 取得硼嗪、苄基自由基、乙炔基自由基和氢过氧自由基的坐标。它们补充了本次 PubChem 未提供三维构象的四个候选。

每份下载均与固定仓库目录中的 Git blob SHA1 核对，并保存 SHA256。原始 `coord` 的单位为 Bohr，转换系数为 0.529177210903 Å/Bohr。三个自由基有作者明确的电荷 0、未配对电子 1 记录；硼嗪的源目录缺少显式状态字段，局部副本提出的中性单重态仍单独保留为待核实状态。

`review-v1` 保留原来的坐标分量比较。该比较没有建立共同坐标系；`review-v2` 纠正了这项解释，增加同一原子顺序下的全部原子对距离核对，允许坐标系的旋转和平移。四个最大距离差分别约为 3.59e-7、4.19e-7、1.62e-7 和 1.30e-7 Å，没有拟合缩放或按观察结果调整阈值。

这些是作者分发的计算基准几何，不是实验坐标，也不是本项目已验收的高质量参考。未启动 Gaussian，未开始这些候选的 COV 验证，正式外部分子数仍为 0。连接关系、电子态、适用的环境及后续优化、频率、稳定性和全库验收仍需完成。

来源和署名：Goerigk, Hansen, Bauer, Ehrlich, Najibi, Grimme 等，PCCP 2017, 19, 32184–32215，[DOI 10.1039/C7CP04913G](https://doi.org/10.1039/C7CP04913G)，[作者仓库](https://github.com/grimme-lab/GMTKN55)，数据集 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/legalcode)。本次只摘取几何数据和来源记录，未运行作者仓库的软件。上述坐标换算及核对是本项目新增内容。
'''
(base / 'benchmark-geometry-preparation-v1/README.md').write_text(benchmark_readme, encoding='utf-8', newline='\n')

sources_root = base / 'external-primary-sources-v1'
source_review_path = sources_root / 'source-review.json'
if source_review_path.exists():
    raise FileExistsError('Retain the existing source review instead of silently replacing it')
source_review = dict(
    schema_version=1, reviewed_epoch=time.time(), formal_external_cases=0,
    accepted_reference_geometries=0, gaussian_started=False,
    transport_and_content_status_are_distinct=True,
    records=[
        dict(candidate='[Co@Ge10]3-', status='primary-literature-reviewed; complete-coordinates-pending',
             sources=['https://onlinelibrary.wiley.com/doi/10.1002/anie.200805511',
                      'https://mediatum.ub.tum.de/doc/1099289/document.pdf',
                      'https://mediatum.ub.tum.de/doc/1534198/1534198.pdf'],
             findings=['2012 thesis article section contains citation followed by a visually confirmed blank page.',
                       '2020 thesis Table 74 title names CoGe10, but its printed sum formula omits Ge; confirmed in page image.',
                       'Different ammonia environments do not create different core-cluster molecule identities.',
                       'Shape ratios are not proof of all ideal point-group operations.'],
             pending=['Complete coordinates, occupancies, disorder model and justified finite environment.']),
        dict(candidate='[Pd@Bi10]4+', status='primary-index-reviewed; complete-coordinates-pending',
             sources=['https://www.mdpi.com/1422-8599/2025/2/M2020',
                      'https://www.degruyterbrill.com/document/doi/10.1515/znb-2021-0159/html'],
             deposition_lead='CCDC 2448552',
             findings=['HTTP 200 for the CCDC request returned a captcha/terms page, not CIF data.',
                       'HTTP 200 for the local MDPI request returned an automatic-access interstitial, not the article.',
                       'Finite-cluster molecular orbitals and periodic-crystal bands are different comparison scopes.',
                       'The 2021 primary article supplies another supplementary-coordinate lead, not a downloaded coordinate file.'],
             pending=['Actual CIF or coordinate-table payload and independent interpretation of disorder/environment.']),
        dict(candidate='W(CH3)6', status='primary-study-scope-reviewed; coordinates-pending',
             sources=['https://pubs.acs.org/doi/10.1021/ja952231p'],
             findings=['The cited 1996 study is computational; whole-molecule and WC6 skeleton symmetry differ.',
                       'Previously retrieved PubChem entry with charge -6 is retained as a rejected identity source.'],
             pending=['Complete coordinates and independent state/source checks.']),
        dict(candidate='[La(NO3)2(phen)(H2O)4]+', status='primary-abstract-reviewed; coordinates-pending',
             sources=['https://link.springer.com/article/10.1007/s11243-007-0224-4'],
             findings=['The reported surrounding components include hmt, nitrate and crystal water.',
                       'An environment component is not automatically a first-shell La ligand.',
                       'A reported polyhedron name is a hypothesis to compare with coordinates, not a mandated COV output.'],
             pending=['Complete experimental coordinates, shell identification and finite environment choice.'])])
source_review['local_source_artifacts'] = [
    dict(file=p.name, bytes=p.stat().st_size, sha256=sha(p), published_full_content=False)
    for p in sorted(sources_root.iterdir())
    if p.is_file() and (p.suffix in ['.pdf', '.html'] or p.name.endswith('.retrieval.json'))
]
atomic_json(source_review_path, source_review)

record = dict(
    snapshot_epoch=time.time(), phase='input-and-control-preparation',
    external_candidate_count=67, pubchem_retrieval_counts=external['counts'],
    benchmark_geometry_candidates=sorted(benchmark_ids),
    candidates_with_archived_starting_coordinates=len(pubchem_ids | benchmark_ids),
    formal_external_cases=0, accepted_quality_references=0,
    density_controls_prepared=20, density_controls_cov_executed=False,
    gaussian_started_by_this_preparation=False,
    independently_managed_gaussian='F:/CalChem/COV/Resume; queue and pause request untouched',
    running_round='NUM-FIX-006',
    running_round_progress_at_snapshot=read(numeric / 'NUM-FIX-006-original/progress.json'),
    production_or_frozen_round_changes=False,
    evidence={
        'pubchem':'external-coordinate-preparation-v2/review-v1/summary.json',
        'benchmark':'benchmark-geometry-preparation-v1/review-v2/summary.json',
        'rare_geometries':'external-primary-sources-v1/source-review.json',
        'density_controls':'regression-fixtures/density-source-v1/manifest.json'})
checkpoint = '''## 2026-09-08 输入坐标与密度控制补充准备

67 个外部预备候选已全部完成本轮公开三维构象检索：取得 35 个计算起始构象，19 个无可用构象，13 个身份待解决。另从作者固定版本基准库取得四个不重叠候选的坐标，现有 39 个预备候选带可追溯起始坐标；这不是 39 个已验收新分子。原始响应、失败收据、元素/电荷/连接图记录和坐标校验值均保留。

补充了 20 份有限、解析归一化 AO 密度控制输入，覆盖缺失占据、自旋来源、非顺序占据、分数占据、相同电子数但非零自旋密度、重排/变号及虚轨道或占据空间截断。独立解析矩阵和预定阈值已固定，尚未运行 COV。现有 IOData 对共享轨道的自旋拆分含模型启发式，不能仅凭程序间一致就把不唯一的自旋密度当作已知。

NUM-FIX-006 的程序、输入和判据保持冻结，全部 273 例终止后统一分析。上述准备没有启动 Gaussian，也没有控制用户在 F:/CalChem/COV/Resume 独立运行的队列。正式完整通过数和正式外部分子数均仍为 0。
'''
state_path = base / '阶段状态.json'
state = read(state_path)
state['input_preparation_20260908'] = record
state['snapshot_epoch'] = record['snapshot_epoch']
atomic_json(state_path, state)
checkpoint_path = base / '执行检查点.md'
previous = checkpoint_path.read_text(encoding='utf-8')
if checkpoint.splitlines()[0] not in previous:
    checkpoint_path.write_text(previous + '\n\n' + checkpoint, encoding='utf-8', newline='\n')

target = repo / 'validation/input-preparation-20260908'
if target.exists():
    raise FileExistsError('Use a new publication checkpoint instead of replacing this one')
target.mkdir(parents=True)
path_map = {}
def copy(source, relative):
    destination = target / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    assert sha(source) == sha(destination)
    path_map[str(source)] = relative.replace('\\', '/')

for folder in ['external-coordinate-preparation-v1', 'external-coordinate-preparation-v2',
               'benchmark-geometry-preparation-v1', 'regression-fixtures/density-source-v1']:
    for source in sorted((base / folder).rglob('*')):
        if source.is_file():
            copy(source, str(source.relative_to(base)))
for name in ['README.md', 'source-review.json', 'GMTKN55-v1-branch-api.json']:
    copy(sources_root / name, 'external-primary-sources-v1/' + name)
for name in ['density-source-contract-v2.md', 'density-model-addendum-v3.md']:
    copy(numeric / name, 'density-design/' + name)
for name in ['checkpoint_input_preparation.py', 'prepare_external_coordinates.py',
             'review_external_coordinate_preparation.py', 'prepare_density_controls.py',
             'index_local_benchmark_candidates.py', 'fetch_primary_benchmark_candidates.py',
             'review_primary_benchmark_candidates.py']:
    copy(workspace / 'work' / name, 'helpers/' + name)
(target / '.gitattributes').write_text('* -text\n', encoding='utf-8', newline='\n')
(target / 'README.md').write_text(checkpoint + '''

`status.json` 是准备阶段快照；正式验收不从这里推导。`local-to-repository-paths.json` 将原始证据内的本机绝对路径映射到本目录的副本，原始字节保持不变。PubChem 构象来源见逐条 URL；GMTKN55 坐标署名、单位转换及修改说明见相应 README。完整版权论文、逐页全文和访问验证页面没有随此档案发布；只公开来源说明和哈希。

密度控制与通解设计尚待 NUM-FIX-006 全批分析结束后执行。运行中的生产版本及冻结运行器不依赖本目录。
''', encoding='utf-8', newline='\n')
atomic_json(target / 'status.json', record)
atomic_json(target / 'local-to-repository-paths.json', path_map)
atomic_json(target / 'files-sha256.json', {
    str(p.relative_to(target)).replace('\\', '/'):sha(p)
    for p in sorted(target.rglob('*')) if p.is_file() and p.name != 'files-sha256.json'})
print(json.dumps(dict(status='preparation-checkpoint-written', directory=str(target),
                     files=len(list(target.rglob('*'))), record=record), ensure_ascii=False))
