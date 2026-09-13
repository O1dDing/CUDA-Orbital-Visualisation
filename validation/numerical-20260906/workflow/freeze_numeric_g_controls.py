"""Freeze two existing Gaussian g-shell fixtures against a frozen COV build.

No calculation is started and no molecular identity is added to the formal
external set. The original all-MO/frame/grid review rules are kept unchanged.
"""
from pathlib import Path
import hashlib
import json
import shutil
import time

BASE = Path(r'F:\Codex\2026-09-05\branch-15\outputs\general-fixes-20260906')
SOURCE = BASE / 'NUM-FIX-003-original'
TARGET = BASE / 'NUM-FIX-003-g-controls'
FIXTURES = Path(r'F:\Codex\2026-09-05\branch-15\outputs\cov-complete-validation-20260906\regression-fixtures\g-shell-v1')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


if TARGET.exists():
    raise FileExistsError('Existing g-control round is not overwritten')
manifest = read(SOURCE / 'manifest.json')
original_manifest_sha = sha(SOURCE / 'manifest.json')
fixture_manifest = read(FIXTURES / 'manifest.json')
TARGET.mkdir()
for directory, key in (('programs', 'programs'), ('runner-source', 'runner_sources'),
                       ('image-deps', 'image_dependency_artifacts')):
    for name, wanted in manifest[key].items():
        source = SOURCE / directory / name
        if sha(source) != wanted:
            raise RuntimeError('Frozen source differs: ' + str(source))
        target = TARGET / directory / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
cases = []
artifacts = {}
for fixture in fixture_manifest['fixtures']:
    old = FIXTURES / fixture['name']
    calculation = read(old / 'calculation.json')
    if calculation['status'] != 'collected' or calculation['fixture_identity'] != fixture['identity']:
        raise RuntimeError('Existing Gaussian fixture is incomplete or has a different identity')
    if sha(old / 'wavefunction.fch') != calculation['fchk_sha256']:
        raise RuntimeError('Existing Gaussian wavefunction differs')
    case_id = 'AUX-G-' + fixture['name'].upper()
    directory = TARGET / 'inputs' / case_id
    directory.mkdir(parents=True)
    names = {'wavefunction.fch': 'wavefunction.fch', 'gaussian.log': 'wavefunction.log',
             'input.gjf': 'original-input.gjf', 'calculation.json': 'original-calculation.json'}
    for source, destination in names.items():
        shutil.copy2(old / source, directory / destination)
        artifacts[str((directory / destination).relative_to(TARGET))] = sha(directory / destination)
    fchk = directory / 'wavefunction.fch'
    info = fchk.stat()
    cases.append({'case_id': case_id, 'input': str(fchk), 'relative_path': case_id + '/wavefunction.fch',
                  'source_size': info.st_size, 'source_mtime_ns': info.st_mtime_ns,
                  'log_candidates': [str(directory / 'wavefunction.log')], 'sha256': calculation['fchk_sha256'],
                  'original_snapshot': str(old / 'wavefunction.fch'),
                  'source_case_id': fixture['source_case_id'], 'fixture_identity': fixture['identity'],
                  'counts_as_new_external_molecule': False})
manifest.update(name=TARGET.name, kind='g-controls', case_count=len(cases), cases=cases,
                input_artifacts=artifacts, created_epoch=time.time(), formal_case_passes=0,
                source_original_manifest_sha256=original_manifest_sha,
                source_fixture_manifest_sha256=sha(FIXTURES / 'manifest.json'),
                auxiliary_representation_controls_only=True, new_external_molecules=0)
manifest['resources']['workers'] = 2
manifest.pop('round_identity')
manifest['round_identity'] = hashlib.sha256(json.dumps(manifest, sort_keys=True, separators=(',', ':'),
                                                     ensure_ascii=False).encode()).hexdigest()
(TARGET / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(json.dumps({'round': str(TARGET), 'case_count': len(cases), 'round_identity': manifest['round_identity'],
                  'new_external_molecules': 0, 'gaussian_calculations_started': 0}))
