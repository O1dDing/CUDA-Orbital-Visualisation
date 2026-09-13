"""Prepare independent finite-AO density controls; do not run or alter COV."""
from pathlib import Path
from decimal import Decimal, getcontext
import ctypes
import hashlib
import json
import sys

sys.path.insert(0, r'F:\Dev\cov-native-validation-20260905\tests')
from validation_process import physical_core_masks
kernel = ctypes.WinDLL('kernel32', use_last_error=True)
kernel.GetCurrentProcess.restype = ctypes.c_void_p
kernel.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
assert kernel.SetProcessAffinityMask(kernel.GetCurrentProcess(), physical_core_masks(12)[8])
root = Path(__file__).resolve().parent.parent / 'outputs/cov-complete-validation-20260906/regression-fixtures/density-source-v1'
if root.exists():
    raise FileExistsError('Retain prepared fixtures and expected values; create a new version for revisions')
root.mkdir(parents=True)
getcontext().prec = 70
D = Decimal
xyz = [[D(v) for v in p] for p in [('0','0','0'),('1.25','0','0'),('0.3','1.7','0')]]
s = [[(-sum((a-b)**2 for a,b in zip(p,q))/2).exp() for q in xyz] for p in xyz]
def inner(a,b):
    return sum(a[i]*s[i][j]*b[j] for i in range(3) for j in range(3))
c = []
for i in range(3):
    v = [D(int(j==i)) for j in range(3)]
    for u in c:
        overlap = inner(u,v)
        v = [a-overlap*b for a,b in zip(v,u)]
    norm = inner(v,v).sqrt()
    c.append([a/norm for a in v])
gram_error = max(abs(inner(a,b)-int(i==j)) for i,a in enumerate(c) for j,b in enumerate(c))
assert gram_error < D('1e-65')
def matrix(vectors, weights):
    return [[sum(D(w)*v[i]*v[j] for v,w in zip(vectors,weights)) for j in range(3)] for i in range(3)]
def pack(mat):
    return [str(mat[i][j]) for i in range(3) for j in range(i+1)]
def trace_s(mat):
    return str(sum(mat[i][j]*s[j][i] for i in range(3) for j in range(3)))
def block(v, occ, spin='Alpha', energy='-0.5'):
    return dict(coefficients=v,occupation=occ,spin=spin,energy=energy)
def molden(blocks, omit_last_coefficient=False):
    lines = ['[Molden Format]', '[Atoms] AU']
    lines += [f'He {i+1} 2 '+ ' '.join(str(a) for a in p) for i,p in enumerate(xyz)]
    lines += ['[GTO]']
    for i in range(3):
        lines += [f'{i+1} 0','s 1 1.0','1.0 1.0','']
    lines += ['[MO]']
    for k,b in enumerate(blocks):
        lines += ['Ene= '+b['energy']]
        if b['spin'] is not None:
            lines += ['Spin= '+b['spin']]
        if b['occupation'] is not None:
            lines += ['Occup= '+str(b['occupation'])]
        lines += [f'{i+1} {value:.17E}' for i,value in enumerate(b['coefficients']) if not(omit_last_coefficient and k==len(blocks)-1 and i==2)]
    return '\n'.join(lines)+'\n'
def fchk(vectors, na, nb, stored_p=None, stored_q=None):
    lines = ['Finite-AO mathematical density input control; not a Gaussian calculation','SP        ROHF     explicit-s-control']
    def scalar(name,value):
        lines.append(f'{name:<40}I{value:16d}')
    def array(name,typ,values):
        lines.append(f'{name:<40}{typ}   N={len(values):12d}')
        for i in range(0,len(values),4):
            lines.append(' '.join(str(v) if typ=='I' else f'{D(v):.17E}' for v in values[i:i+4]))
    for name,value in [('Number of atoms',3),('Charge',6-na-nb),('Multiplicity',abs(na-nb)+1),
                       ('Number of electrons',na+nb),('Number of alpha electrons',na),('Number of beta electrons',nb),
                       ('Number of basis functions',3),('Number of independent functions',len(vectors))]:
        scalar(name,value)
    array('Atomic numbers','I',[2,2,2])
    array('Current cartesian coordinates','R',[v for p in xyz for v in p])
    array('Shell types','I',[0,0,0])
    array('Number of primitives per shell','I',[1,1,1])
    array('Shell to atom map','I',[1,2,3])
    array('Primitive exponents','R',[1,1,1])
    array('Contraction coefficients','R',[1,1,1])
    array('Alpha Orbital Energies','R',[D('-.8')+D('.5')*i for i in range(len(vectors))])
    array('Alpha MO coefficients','R',[v for col in vectors for v in col])
    if stored_p is not None:
        array('Total SCF Density','R',pack(stored_p))
    if stored_q is not None:
        array('Spin SCF Density','R',pack(stored_q))
    return '\n'.join(lines)+'\n'
records = []
def add(ident, name, text, p=None, q=None, p_status='available', q_status='available', notes='', relation=None):
    path = root / (ident+'-'+name)
    path.write_text(text,encoding='ascii',newline='\n')
    record = dict(control_id=ident, file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  expected_total_status=p_status, expected_spin_status=q_status, notes=notes,
                  cov_execution_status='not-run; wait for NUM-FIX-006 full collection and unified analysis',
                  formal_case=False, new_external_molecule=False)
    if p is not None:
        record.update(total_density_lower_packed_decimal=pack(p), total_electrons_trace_ps=trace_s(p))
    if q is not None:
        record.update(spin_density_lower_packed_decimal=pack(q), spin_electrons_trace_qs=trace_s(q))
    if relation:
        record['metamorphic_relation'] = relation
    records.append(record)

zero = matrix(c,[0,0,0])
shared = [block(c[i],occ,energy=str(D('-.7')+D('.5')*i)) for i,occ in enumerate([1,2,0])]
p = matrix(c,[1,2,0]); q = matrix(c,[1,0,0])
add('DEN-C01','missing-occupation.molden',molden([block(c[0],2),block(c[1],None),block(c[2],0)]),
    p_status='missing-input',q_status='missing-input',notes='Unknown occupation cannot be replaced by zero. Complete density and total Mayer evidence must remain unavailable.')
add('DEN-C02','shared-open-unsorted.molden',molden(shared),p,q,
    q_status='available-under-explicitly-recorded-integer-shared-determinant-interpretation',
    notes='With conventional integer shared-orbital interpretation, the actual singly occupied direction is c0. A stricter unknown-model result must expose its reason; c1 spin density is incorrect.')
add('DEN-C03','shared-open-reordered.molden',molden([shared[1],shared[2],shared[0]]),p,q,
    q_status='same-interpretation-as-DEN-C02',relation={'control':'DEN-C02','old_mo_for_new':[1,2,0],'total':'equal','spin':'equal','occupations':'follow-coefficient-vector'})
for ident,occs in [('DEN-C04',['.5','.5',0]),('DEN-C05',['1.5','.5',0])]:
    add(ident,'fractional-shared.molden',molden([block(v,n) for v,n in zip(c,occs)]),matrix(c,occs),
        q_status='missing-input-unless-external-model-supplied',notes='Total Occup does not encode a unique alpha-minus-beta partition. Aggregated integer counts do not remove this ambiguity.')
beta = [[D('.6')*a+D('.8')*b for a,b in zip(c[0],c[1])],[-D('.8')*a+D('.6')*b for a,b in zip(c[0],c[1])],c[2]]
explicit = [block(v,1 if i==0 else 0,'Alpha') for i,v in enumerate(c)] + [block(v,1 if i==0 else 0,'Beta') for i,v in enumerate(beta)]
up = matrix(c+beta,[1,0,0,1,0,0]); uq = matrix(c+beta,[1,0,0,-1,0,0])
assert max(abs(v) for row in uq for v in row)>D('.1') and abs(D(trace_s(uq)))<D('1e-65')
add('DEN-C06','explicit-spins-equal-counts.molden',molden(explicit),up,uq,
    notes='Explicit alpha/beta spaces each contain one electron, but Q is nonzero. Equal electron counts do not prove a closed shell.')
flipped = [{**b,'spin':'Beta' if b['spin']=='Alpha' else 'Alpha'} for b in explicit]
add('DEN-C07','explicit-spins-reversed.molden',molden(flipped),up,[[-v for v in row] for row in uq],
    relation={'control':'DEN-C06','total':'equal','spin':'negated'})
phase_order = [5,0,3,2,1,4]
reordered = [{**explicit[i],'coefficients':[(-v if k%2 else v) for v in explicit[i]['coefficients']]} for k,i in enumerate(phase_order)]
add('DEN-C08','explicit-spins-reordered-phases.molden',molden(reordered),up,uq,
    relation={'control':'DEN-C06','old_mo_for_new':phase_order,'total':'equal','spin':'equal'})
add('DEN-C09','closed-zero-spin.molden',molden([block(v,n) for v,n in zip(c,[2,0,0])]),matrix(c,[2,0,0]),zero,
    notes='An explicit complete shared closed-shell occupation has an available all-zero spin matrix.')
add('DEN-C10','missing-spin-field.molden',molden([block(v,n,None) for v,n in zip(c,[1,0,0])]),matrix(c,[1,0,0]),
    q_status='missing-input-unless-explicit-model-supplied',notes='Default Alpha enum is not a producer Spin field. Preserve independent total-density usability.')
for ident,occ in [('DEN-C11','NaN'),('DEN-C12','-.5'),('DEN-C13','2.5')]:
    add(ident,'invalid-ordinary-occupation.molden',molden([block(c[0],occ),block(c[1],0),block(c[2],0)]),
        p_status='invalid-input-or-unsupported-input-purpose',q_status='invalid-input-or-unsupported-input-purpose',
        notes='Do not output fabricated zero/full density. Signed difference-density files may have another legitimate purpose, but are outside this ordinary electronic-state control.')
bad_explicit = [{**b,'occupation':2 if i==0 else b['occupation']} for i,b in enumerate(explicit)]
add('DEN-C14','invalid-explicit-spin-occupation.molden',molden(bad_explicit),
    p_status='invalid-input',q_status='invalid-input',notes='Explicit-spin occupation exceeds Pauli range [0,1]. Electron-count unavailability alone is insufficient if complete density is still published.')
add('DEN-C15','missing-virtual-coefficient.molden',molden([block(v,n) for v,n in zip(c,[2,0,0])],omit_last_coefficient=True),
    p_status='parser-rejects-incomplete-coefficient-block',q_status='parser-rejects-incomplete-coefficient-block',
    notes='Existing parser protection control: malformed coefficient dimensions must not be hidden by a zero occupation.')
cp = matrix(c,[2,1,0]); cq=matrix(c,[0,1,0])
add('DEN-C16','canonical-full.fchk',fchk(c,2,1),cp,cq,
    notes='Canonical Gaussian-style count-to-MO rule is appropriate for this explicitly defined complete shared determinant.')
add('DEN-C17','canonical-virtual-truncated.fchk',fchk(c[:2],2,1),cp,cq,
    relation={'control':'DEN-C16','total':'equal','spin':'equal','removed_mo':'unoccupied-only'})
add('DEN-C18','canonical-occupied-truncated.fchk',fchk(c[:1],2,1),
    p_status='missing-input-occupied-subspace',q_status='missing-input-occupied-subspace',notes='Alpha count 2 with only one supplied canonical direction: cannot reconstruct the absent occupied orbital.')
add('DEN-C19','stored-density-with-truncated-mos.fchk',fchk(c[:1],2,1,cp,cq),cp,cq,
    notes='Actual stored producer-format P/Q are complete, but MO-reconstruction diagnostics must report missing occupied subspace; do not overwrite the stored matrices.')
add('DEN-C20','explicit-fractional-spins.molden',molden([block(v,n,'Alpha') for v,n in zip(c,['.5','.5',0])]+[block(v,n,'Beta') for v,n in zip(beta,['.25','.25',0])]),
    matrix(c+beta,['.5','.5',0,'.25','.25',0]),matrix(c+beta,['.5','.5',0,'-.25','-.25',0]),
    notes='Explicit alpha/beta fractional occupations define P/Q even when integer electron counters cannot represent their trace; density availability must not depend on integer count storage.')

manifest = dict(description='Finite normalized s-AO algebra controls, not molecular reference calculations',
                generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scientific_source_commit_to_probe='21036b2990ef86dcdc419c65cdebd2a03bf4a12b',
                source_ao_definition='Normalized exponent-1 s Gaussians; Sij=exp(-|Ri-Rj|^2/2)',
                atomic_numbers=[2,2,2], coordinates_bohr=[[str(v) for v in p] for p in xyz],
                overlap_decimal=[[str(v) for v in row] for row in s],
                coefficients_by_mo_decimal=[[str(v) for v in col] for col in c],
                decimal_precision=70, analytic_gram_max_error=str(gram_error),
                fixed_matrix_absolute_tolerance=2e-12, coefficient_print_significant_digits=18,
                threshold_status='set before any execution against COV; not adjusted to observed failures',
                control_count=len(records), formal_case_count=0, new_external_molecule_count=0,
                cov_executed=False, gaussian_executed=False, controls=records)
(root/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
(root/'generator.py').write_bytes(Path(__file__).read_bytes())
print(json.dumps(dict(directory=str(root),control_count=len(records),analytic_gram_max_error=str(gram_error),cov_executed=False)),flush=True)
