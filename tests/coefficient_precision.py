"""Independent coefficient-rounding error in the full-space AO metric."""
import os
from pathlib import Path
import sys

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
sys.path.insert(0, str(Path(r'F:\Dev\cov-validation-20260905\reference-deps')))
import numpy as np


def measure_coefficient_rounding(coefficients, overlap):
    coefficients = np.asarray(coefficients, dtype=np.float64)
    overlap = np.asarray(overlap, dtype=np.float64)
    if coefficients.ndim != 2 or not all(coefficients.shape):
        raise ValueError('A nonempty AO-by-MO matrix is required')
    if overlap.shape != (coefficients.shape[0], coefficients.shape[0]):
        raise ValueError('AO metric dimension mismatch')
    if not np.isfinite(coefficients).all() or not np.isfinite(overlap).all():
        raise ValueError('Finite coefficients and AO metric are required')
    symmetry_tolerance = 8*np.finfo(np.float64).eps*max(1.0, float(np.abs(overlap).max()))
    if not np.allclose(overlap, overlap.T, rtol=0, atol=symmetry_tolerance):
        raise ValueError('AO metric is not symmetric at arithmetic precision')
    with np.errstate(over='raise', invalid='raise'):
        rounded = coefficients.astype(np.float32).astype(np.float64)
    # Compute the difference itself in the AO metric. Subtracting two large
    # self-norms would lose the small error in cancellation-sensitive cases.
    chol = np.linalg.cholesky(overlap)
    original_metric_vectors = chol.T @ coefficients
    difference_metric_vectors = chol.T @ (rounded-coefficients)
    norms = np.linalg.norm(original_metric_vectors, axis=0)
    errors = np.linalg.norm(difference_metric_vectors, axis=0)
    if not np.isfinite(norms).all() or not np.isfinite(errors).all() or not (norms > 0).all():
        raise ValueError('Finite nonzero original orbital norms are required')
    return {'original_norms': norms, 'error_norms': errors, 'relative_errors': errors/norms,
            'maximum_absolute_coefficient_rounding_error': float(np.abs(rounded-coefficients).max()),
            'scope': 'Coefficient-only IEEE float32 roundtrip; not a complete CPU/CUDA or acceptance result'}
