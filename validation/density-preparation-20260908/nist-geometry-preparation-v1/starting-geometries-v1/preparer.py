"""Keep three literal public coordinates as starting geometries, never references."""
from pathlib import Path
import ctypes
import hashlib
import json
import math
import sys
import time

workspace = Path(r'F:\Codex\2026-09-05\branch-15')
sys.path.insert(0, r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json, physical_core_masks
kernel = ctypes.WinDLL('kernel32', use_last_error=True)
kernel.GetCurrentProcess.restype = ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(), physical_core_masks(12)[8])
root = workspace / 'outputs/cov-complete-validation-20260906/nist-geometry-preparation-v1'
read = lambda p: json.loads(p.read_text(encoding='utf-8'))
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
collection = read(root / 'collection-complete.json')
assert collection['all_terminal'] and collection['terminal_sources'] == 5
out = root / 'starting-geometries-v1'
if out.exists():
    raise FileExistsError('Keep the prepared source interpretation and geometry bytes')
out.mkdir()
(out / 'preparer.py').write_bytes(Path(__file__).read_bytes())
interpretations = {
    'PREP-002': dict(name='germane', multiplicity=1, state='1A1', point_group='Td',
        coordinate_method='HSEh1PBE/6-311G**', source_type='NIST archived computed geometry',
        state_table_index=1, expected_state_row=['1','1','yes','Td','1A1'],
        limitations=['Coordinates have three decimal places in angstrom.',
                     'The more precise printed bond length and rounded Cartesians need not coincide exactly.']),
    'PREP-008': dict(name='bromine pentafluoride', multiplicity=1, state='1A1', point_group='C4v',
        coordinate_method=None, source_type='NIST curated experimental-geometry compilation',
        state_table_index=3, expected_state_row=['1A1','C4V'],
        geometry_reference='Gurvich, Veyts, Alcock, Thermodynamic Properties of Individual Substances, fourth edition, Hemisphere, New York, 1989',
        limitations=['The original experimental coordinate derivation has not been obtained.',
                     'The Cartesian table has four decimal places; this is not the experimental uncertainty.',
                     'The 1965 vibrational-spectrum paper listed on the same page is not automatically the geometry source.']),
    'PREP-042': dict(name='phenoxy radical', multiplicity=2, state='2B1', point_group='C2v',
        coordinate_method='B3LYP/6-31G*', source_type='NIST archived computed geometry',
        energy_page_method='G3B3', state_table_index=1,
        expected_state_row=['1','1','yes','C2V','2B1'],
        limitations=['Coordinates have three decimal places in angstrom.',
                     'The G3B3 energy-page label does not identify the geometry optimization method.'])}
records = []
for source in collection['records']:
    ident = source['candidate_id']
    table_file = root / 'table-review-v1' / (ident + '-tables.json')
    data = read(table_file)
    assert data['raw_sha256'] == sha(root / ident / 'source.html') == source['sha256']
    record = dict(candidate_id=ident, source_url=source['url'], raw_source_sha256=source['sha256'],
                  literal_table_sha256=sha(table_file), formal_case=False,
                  accepted_quality_reference=False, gaussian_started=False, cov_executed=False)
    if ident not in interpretations:
        assert not data['coordinate_tables']
        record.update(status='no-cartesian-coordinate-payload-in-this-source',
                      notes=source['scope'])
        records.append(record)
        continue
    interpretation = interpretations[ident]
    state_table = next(t for t in data['tables'] if t['index'] == interpretation['state_table_index'])
    assert interpretation['expected_state_row'] in state_table['rows']
    assert len(data['coordinate_tables']) == 1
    geometry = data['coordinate_tables'][0]
    assert geometry['expected_element_counts_match']
    if interpretation['coordinate_method']:
        assert interpretation['coordinate_method'] in geometry['preceding_text']
    if ident == 'PREP-008':
        assert 'Gurvich' in geometry['preceding_text']
    rows = geometry['rows']
    assert all(math.isfinite(v) for row in rows for v in row['coordinates_angstrom'])
    directory = out / ident
    directory.mkdir()
    lines = [str(len(rows)), f'{ident}; literal NIST Cartesians in angstrom; starting input only']
    lines += [row['element'] + ' ' + ' '.join(row['literal_coordinates']) for row in rows]
    path = directory / 'starting.xyz'
    path.write_text('\n'.join(lines) + '\n', encoding='ascii', newline='\n')
    record.update(status='literal-starting-geometry-prepared', charge=source['charge'],
        interpretation=interpretation, source_table_index=geometry['source_table_index'],
        atom_count=len(rows), element_counts=geometry['element_counts'], units='angstrom',
        xyz_file=str(path.relative_to(out)), xyz_sha256=sha(path),
        transformations='None; source element order and literal coordinate tokens retained.',
        geometry_quality='Starting coordinates only; future optimization, frequency, stability, environment and identity checks remain required.')
    atomic_json(directory / 'source-and-state.json', record)
    records.append(record)
summary = dict(created_epoch=time.time(), sources_examined=5, starting_geometries_prepared=3,
               no_cartesian_payload=2, records=records, formal_external_cases=0,
               accepted_quality_references=0, gaussian_started=False, cov_executed=False)
atomic_json(out / 'summary.json', summary)
(out / 'README.md').write_text('''# NIST 三份起始结构

保留 GeH4、BrF5、苯氧自由基三份公开表格的原子顺序和原始坐标小数，没有对称化、拟合或补造坐标。GeH4 是 HSEh1PBE/6-311G** 计算几何；苯氧自由基页面的能量标签是 G3B3，但几何明确来自 B3LYP/6-31G*；BrF5 的几何由 NIST 汇编，表内引用 Gurvich 等 1989 年专著，原始实验推导仍待补查。

ClF3 本页有不同构象的信息：首个状态行的 D3h 不能替代几何部分的 C2v 最小构象。ClF3 和 IF7 本轮保存的页面均没有实际笛卡尔坐标，继续记录为待补来源。

逐项来源 URL、原始响应哈希、表格索引、状态和局限见 `summary.json` 及各候选的 `source-and-state.json`。XYZ 是后续计算的起点，未启动 Gaussian 或外部 COV 验证，不计作已验收参考。页面的小数位数不等同于物理误差。
''', encoding='utf-8', newline='\n')
print(json.dumps(dict(output=str(out), starting_geometries=3, formal_external_cases=0)))
