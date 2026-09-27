"""Package actual campaign data in bounded, checksum-addressed GitHub assets."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import os
import sys
import time
import zipfile

BASE = Path(r'F:\Codex\2026-09-05\branch-15')
ROOT = BASE / 'outputs/cov-complete-validation-20260906'
DEST = BASE / 'outputs/github-handoff-20260906'
sys.path.insert(0, r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import atomic_json, physical_core_masks, pin_current
pin_current(sum(physical_core_masks(4)))
DEST.mkdir(exist_ok=True)
ASSETS = DEST / 'assets'
ASSETS.mkdir(exist_ok=True)
SCOPES = [
    ('campaign', ROOT),
    ('jobs/REF-001', Path(r'F:\Dev\cov-cycle-20260906\jobs\REF-001')),
    ('jobs/R001', Path(r'F:\Dev\cov-cycle-20260906\jobs\R001')),
    ('jobs/S-001', Path(r'F:\Dev\cov-cycle-20260906\jobs\S-001')),
    ('original-collection', BASE / 'outputs/fchk-collection-20260905'),
]
EXCLUDE_SUFFIXES = {'.exe', '.dll', '.obj', '.pyc', '.pyd', '.lock'}
files = []
excluded = []
for scope, source in SCOPES:
    for directory, dirs, names in os.walk(source):
        dirs.sort()
        names.sort()
        for name in names:
            path = Path(directory) / name
            relative = path.relative_to(source)
            record = {'source_path': str(path), 'path': scope + '/' + relative.as_posix(),
                      'bytes': path.stat().st_size}
            if path.suffix.lower() in EXCLUDE_SUFFIXES or '__pycache__' in relative.parts:
                excluded.append(dict(record, reason='generated tool/runtime/build/cache file; scientific data and source retained'))
            else:
                files.append(record)
assert len([r for r in files if r['path'].startswith('original-collection/cases/') and r['path'].endswith('/source.fch')]) == 273
assert len([r for r in files if r['path'].startswith('campaign/references/REF-001/reference-inputs/') and r['path'].endswith('/initial.gjf')]) == 273
groups = []
current = []
size = 0
for record in files:
    assert record['bytes'] < 2 * 1024**3 - 1024**2, record['path']
    if current and size + record['bytes'] > 768 * 1024**2:
        groups.append(current)
        current = []
        size = 0
    current.append(record)
    size += record['bytes']
if current:
    groups.append(current)
atomic_json(DEST / 'packaging-plan.json', {'created_epoch': time.time(), 'archive_count': len(groups),
    'file_count': len(files), 'uncompressed_bytes': sum(r['bytes'] for r in files),
    'licensed_gaussian_binaries_included': False, 'exclusions': excluded})

def sha_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024**2), b''):
            digest.update(chunk)
    return digest.hexdigest()

def pack(entry):
    index, group = entry
    name = f'cov-paused-data-{index+1:03d}-of-{len(groups):03d}.zip'
    target = ASSETS / name
    sidecar = ASSETS / (name + '.json')
    if target.exists() and sidecar.exists():
        old = json.loads(sidecar.read_text())
        assert old['sha256'] == sha_file(target)
        return old
    part = target.with_suffix('.zip.partial')
    records = []
    started = time.perf_counter()
    with zipfile.ZipFile(part, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=1, allowZip64=True) as archive:
        for record in group:
            source = Path(record['source_path'])
            before = source.stat()
            digest = hashlib.sha256()
            with source.open('rb') as stream, archive.open(record['path'], 'w', force_zip64=True) as out:
                for chunk in iter(lambda: stream.read(1024**2), b''):
                    digest.update(chunk)
                    out.write(chunk)
            after = source.stat()
            assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), record['path']
            records.append({'path': record['path'], 'bytes': before.st_size, 'sha256': digest.hexdigest()})
    os.replace(part, target)
    value = {'name': name, 'bytes': target.stat().st_size, 'sha256': sha_file(target),
             'file_count': len(records), 'files': records, 'wall_seconds': time.perf_counter()-started}
    assert value['bytes'] < 2 * 1024**3
    atomic_json(sidecar, value)
    print(json.dumps({k: value[k] for k in ['name', 'bytes', 'file_count', 'wall_seconds']}), flush=True)
    return value

with ThreadPoolExecutor(max_workers=4) as pool:
    archives = list(pool.map(pack, enumerate(groups)))
manifest = {'schema_version': 1, 'created_epoch': time.time(),
    'purpose': 'Paused REF-001 recovery and established COV fixes; not an accepted software release',
    'counts': {'original_cases': 273, 'reference_candidates_collected': 13,
               'reference_cases_interrupted': 2, 'reference_cases_not_started': 258,
               'formal_scientific_passes': 0, 'formal_external_molecules': 0},
    'archives': archives, 'excluded_generated_tools': excluded,
    'includes': ['all original source FCHK and original native evidence', 'all 273 independent reproduction results',
                 'all 273 actual prepared reference inputs and basis/ECP data',
                 'all REF-001 complete and interrupted attempts, checkpoints, scratch and logs',
                 'independent numerical evidence, failed diagnostics, source snapshots and progress'],
    'gaussian_software_included': False}
atomic_json(DEST / 'data-manifest.json', manifest)
print(json.dumps({'complete': True, 'archives': len(archives), 'archive_bytes': sum(a['bytes'] for a in archives),
                  'files': sum(a['file_count'] for a in archives)}), flush=True)
