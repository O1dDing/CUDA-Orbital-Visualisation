"""Strict Gaussian MO cube reader and independent full-grid comparisons.

Only same-grid, same-input evidence is comparable. Global orbital sign is
reported and may be reversed; no amplitude fitting, rotation, or interpolation
is performed. A passing grid is not a formal molecule acceptance verdict.
"""
from __future__ import annotations

import math
from pathlib import Path

from independent_grid_reference import np, sha256


def number(token):
    value = float(token.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(value):
        raise ValueError('Non-finite cube field')
    return value


def read_mo_cube(path):
    with Path(path).open(encoding='ascii') as stream:
        comments = [stream.readline().rstrip('\n'), stream.readline().rstrip('\n')]
        row = stream.readline().split()
        if len(row) not in (4, 5):
            raise ValueError('Malformed cube origin header')
        natom = int(row[0])
        if natom >= 0:
            raise ValueError('An orbital cube with explicit dataset identities is required')
        if len(row) == 5 and int(row[4]) != 1:
            raise ValueError('Orbital cube header NVAL must be one or omitted')
        origin = np.array([number(value) for value in row[1:4]])
        shape, axes = [], []
        for _ in range(3):
            row = stream.readline().split()
            if len(row) != 4 or int(row[0]) < 2:
                raise ValueError('Positive Bohr axes with at least two points are required')
            shape.append(int(row[0]))
            axes.append([number(value) for value in row[1:]])
        axes = np.asarray(axes)
        if abs(np.linalg.det(axes)) < 1e-20:
            raise ValueError('Cube axes do not span three dimensions')
        atoms = []
        for _ in range(abs(natom)):
            row = stream.readline().split()
            if len(row) != 5:
                raise ValueError('Malformed cube atom record')
            atoms.append({'atomic_number': int(row[0]), 'nuclear_charge': number(row[1]),
                          'coordinates_bohr': [number(value) for value in row[2:]]})
        tokens = stream.read().split()
    if not tokens:
        raise ValueError('Missing orbital cube dataset header')
    datasets = int(tokens[0])
    if datasets < 1:
        raise ValueError('At least one orbital dataset is required')
    ids = [int(token) for token in tokens[1:1+datasets]]
    if len(ids) != datasets or len(set(ids)) != datasets:
        raise ValueError('Missing or duplicate cube orbital identity')
    count = math.prod(shape)
    body = tokens[1+datasets:]
    if len(body) != count*datasets:
        raise ValueError('Truncated or extra cube grid values')
    # Gaussian interleaves datasets at each point. Preserve C-order (z fastest).
    values = np.fromiter((number(token) for token in body), dtype=np.float64).reshape(count, datasets).T
    return {'path': str(Path(path).resolve()), 'sha256': sha256(path), 'comments': comments,
            'origin_bohr': origin, 'shape': shape, 'step_vectors_bohr': axes,
            'atoms': atoms, 'orbital_ids': ids, 'values': values}


def grid_metrics(actual, reference, thresholds):
    actual, reference = np.asarray(actual, float), np.asarray(reference, float)
    if actual.shape != reference.shape or not actual.size:
        raise ValueError('Nonempty, equal-sized grid arrays are required')
    if not np.isfinite(actual).all() or not np.isfinite(reference).all():
        raise ValueError('Non-finite orbital grid values')
    actual, reference = actual.ravel(), reference.ravel()
    scale = max(float(np.max(np.abs(actual))), float(np.max(np.abs(reference))))
    if scale == 0:
        return {'status': 'insufficient', 'pass': False, 'reason': 'Both sampled fields are zero',
                'points': actual.size}
    a, r = actual/scale, reference/scale
    aa, rr, ar = float(a@a), float(r@r), float(a@r)
    if aa == 0 or rr == 0:
        return {'status': 'insufficient', 'pass': False, 'reason': 'One sampled field has zero norm',
                'points': actual.size}
    phase = 1 if ar >= 0 else -1
    error = phase*a-r
    nrms = math.sqrt(float(error@error)/rr)
    cosine = min(1., abs(ar)/math.sqrt(aa*rr))
    relative_peak = float(np.max(np.abs(error)))/float(np.max(np.abs(r)))
    passed = (nrms <= thresholds['nrms_max'] and cosine >= thresholds['abs_cosine_min']
              and relative_peak <= thresholds['relative_peak_error_max'])
    return {'status': 'pass' if passed else 'fail', 'pass': passed, 'points': actual.size,
            'nrms': nrms, 'abs_cosine': cosine, 'relative_peak_error': relative_peak, 'phase': phase,
            'actual_peak': float(np.max(np.abs(actual))), 'reference_peak': float(np.max(np.abs(reference))),
            'max_absolute_error': float(np.max(np.abs(error)))*scale,
            'thresholds': dict(thresholds), 'amplitude_fitting': False}


def compare_cube_to_reference(cube, record, directory, column, expected_orbital, thresholds):
    """Compare one explicitly identified spin/MO column, after its batch barrier."""
    from iodata import load_one
    directory = Path(directory)
    input_path = Path(record['source_input'])
    if sha256(input_path) != record['input_sha256']:
        raise ValueError('Independent reference input identity changed')
    grid = record['grid']
    if cube['shape'] != grid['shape'] or not np.allclose(cube['origin_bohr'], grid['origin_bohr'], atol=1e-12, rtol=0):
        raise ValueError('Different cube/reference point grids')
    if not np.allclose(cube['step_vectors_bohr'], grid['step_vectors_bohr'], atol=1e-12, rtol=0):
        raise ValueError('Different cube/reference step vectors')
    mol = load_one(str(input_path))
    if len(cube['atoms']) != len(mol.atnums) or any(row['atomic_number'] != int(z) for row, z in zip(cube['atoms'], mol.atnums)):
        raise ValueError('Cube/reference atomic identities differ')
    coordinates = np.array([row['coordinates_bohr'] for row in cube['atoms']])
    if not np.allclose(coordinates, mol.atcoords, atol=6e-7, rtol=0):
        raise ValueError('Cube/reference atom coordinate frames differ beyond cube rounding precision')
    if cube['orbital_ids'] != [expected_orbital]:
        raise ValueError('Cube dataset does not contain exactly the requested orbital')
    if column < 0 or column >= record['nmo']:
        raise ValueError('MO column outside independent reference')
    mo_path = directory / 'mo_values.npy'
    if sha256(mo_path) != record['artifacts'][mo_path.name]['sha256']:
        raise ValueError('Independent grid artifact identity changed')
    values = np.load(mo_path, mmap_mode='r', allow_pickle=False)
    result = grid_metrics(cube['values'][0], values[column], thresholds)
    result.update(cube_sha256=cube['sha256'], reference_input_sha256=record['input_sha256'],
                  orbital_identity=record['orbitals'][column], scientific_case_pass=False)
    return result
