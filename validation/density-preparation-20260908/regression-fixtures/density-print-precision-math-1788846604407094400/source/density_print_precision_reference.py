"""Independent matrix envelopes for explicitly declared decimal-print rounding.

This is preparation for a separate source-consistency check. It does not infer
occupations or spin models, modify COV, or set a rendering tolerance. A matrix
inside an entrywise envelope is merely not contradicted by that envelope; this
does not prove that all entries share a single underlying exact wavefunction.
"""
from decimal import Decimal, localcontext
import math
import os
from pathlib import Path
import re
import sys

for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
sys.path.insert(0,str(Path(r'F:\Dev\cov-validation-20260905\reference-deps')))
import numpy as np

TOKEN=re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$')
UNIT_ROUNDOFF=np.finfo(np.float64).eps/2
MIN_SUBNORMAL=np.nextafter(np.float64(0),np.float64(1))

def gamma(operation_count):
    if not isinstance(operation_count,int) or isinstance(operation_count,bool) or operation_count<0:
        raise ValueError('An integer nonnegative operation count is required')
    x=operation_count*UNIT_ROUNDOFF
    if x>=.1:raise ValueError('Arithmetic bound is not useful at this operation count')
    return x/(1-x)

def printed_value_interval(token, *, rounding_mode):
    """Retain the actual final decimal place; never guess the print convention.

    nearest: half a last-place unit for nearest decimal formatting.
    one-last-place: a whole unit, if the declared writer may truncate.
    Both include binary64 parsing uncertainty and retain nonzero subnormal bounds.
    """
    if rounding_mode not in ('nearest','one-last-place'):
        raise ValueError('Explicitly declare the assumed decimal writer rounding mode')
    if not isinstance(token,str) or not TOKEN.fullmatch(token):
        raise ValueError('A finite plain Fortran real token is required')
    with localcontext() as ctx:
        ctx.prec=max(80,len(token)*2)
        value=Decimal(token.replace('D','E').replace('d','e'))
        if not value.is_finite():raise ValueError('A finite decimal value is required')
        quantum=Decimal(1).scaleb(value.as_tuple().exponent)
        decimal_radius=quantum*(Decimal('.5') if rounding_mode=='nearest' else Decimal(1))
        binary=float(value)
        if not math.isfinite(binary):raise ValueError('Decimal value is outside finite binary64 range')
        parse_error=abs(Decimal.from_float(binary)-value)
        radius=float(decimal_radius+parse_error)
        if not math.isfinite(radius):raise ValueError('Print uncertainty is outside finite binary64 range')
    # Outward rounding of the bound, including underflow to zero.
    radius=float(np.nextafter(radius,np.inf))
    return dict(value=binary,radius=radius,decimal_quantum=str(quantum),
                exact_decimal_radius=str(decimal_radius),binary_parse_error=str(parse_error),
                rounding_mode=rounding_mode,literal=token)

def outer_sum_envelope(coefficients, weights, coefficient_radii=None, weight_radii=None):
    """Enclose sum_k w_k c_k c_k.T for caller-supplied independent intervals.

    Use zero radii for a strict implementation comparison from the *same parsed
    values*. Use printed-token radii only for a separate producer-density versus
    producer-MO consistency diagnosis. Signed weights support alpha-minus-beta.
    The arithmetic allowance is fixed from operation counts, never fitted to data.
    """
    c=np.asarray(coefficients,dtype=np.float64);w=np.asarray(weights,dtype=np.float64)
    if c.ndim!=2 or not c.shape[0] or w.shape!=(c.shape[1],):
        raise ValueError('A compatible AO-by-MO matrix and weights are required')
    r=np.zeros_like(c) if coefficient_radii is None else np.asarray(coefficient_radii,dtype=np.float64)
    rw=np.zeros_like(w) if weight_radii is None else np.asarray(weight_radii,dtype=np.float64)
    if (r.shape!=c.shape or rw.shape!=w.shape or any(not np.isfinite(v).all() for v in (c,w,r,rw))
            or np.any(r<0) or np.any(rw<0)):
        raise ValueError('Finite values and compatible nonnegative radii are required')
    k=c.shape[1]
    if not k:
        zeros=np.zeros((c.shape[0],c.shape[0]))
        return dict(nominal=zeros,source_radius=zeros.copy(),arithmetic_radius=zeros.copy(),
                    total_radius=zeros.copy(),terms=0,arithmetic_operation_bound=0)
    a=np.abs(c);aw=np.abs(w)
    with np.errstate(over='raise',invalid='raise'):
        nominal=(c*w)@c.T
        magnitude=(a*aw)@a.T
        source=(a*aw)@r.T+(r*aw)@a.T+(r*aw)@r.T
        if np.any(rw):source+=((a+r)*rw)@(a+r).T
        # A deliberately conservative forward bound for the weighted products,
        # sums, positive uncertainty products, and their final additions.
        # Fused operations require no larger allowance than unfused operations.
        operation_bound=12*k+64
        arithmetic=gamma(operation_bound)*(magnitude+source)+MIN_SUBNORMAL*operation_bound
        source_upper=np.nextafter(source+gamma(operation_bound)*source,np.inf)
        total=np.nextafter(source_upper+arithmetic,np.inf)
    if any(not np.isfinite(v).all() for v in (nominal,source_upper,arithmetic,total)):
        raise ValueError('Overflowed density envelope cannot be accepted')
    return dict(nominal=nominal,source_radius=source_upper,arithmetic_radius=arithmetic,
                total_radius=total,terms=k,arithmetic_operation_bound=operation_bound)

def compare_matrix_to_envelope(observed,envelope,observed_radii=None):
    observed=np.asarray(observed,dtype=np.float64)
    nominal=envelope['nominal'];radius=envelope['total_radius']
    ro=np.zeros_like(observed) if observed_radii is None else np.asarray(observed_radii,dtype=np.float64)
    if (observed.shape!=nominal.shape or ro.shape!=nominal.shape
            or not np.isfinite(observed).all() or not np.isfinite(ro).all() or np.any(ro<0)):
        raise ValueError('Finite dimensionally compatible observed values/radii are required')
    with np.errstate(over='raise',invalid='raise'):
        difference=np.abs(observed-nominal)
        budget=radius+ro
        budget=np.nextafter(budget+gamma(8)*(np.abs(observed)+np.abs(nominal)+budget)+8*MIN_SUBNORMAL,np.inf)
        excess=np.maximum(difference-budget,0)
    if not np.isfinite(budget).all():raise ValueError('Overflowed comparison budget is invalid')
    return dict(entrywise_compatible=bool(np.all(excess==0)),
                violating_entries=int(np.count_nonzero(excess)),
                maximum_absolute_difference=float(np.max(difference)),
                maximum_absolute_budget=float(np.max(budget)),
                maximum_excess=float(np.max(excess)),
                formal_case_pass=False,
                scope='Entrywise necessary source-rounding compatibility only; no spin/model/physical/rendering verdict')
