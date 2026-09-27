# Conservative molecular-counting evidence

`structure_identity.py` and `structure_identity_batch.py` are independent validation
infrastructure. They do not replace COV bond perception or certify a scientific
validation pass. Every original case remains in the formal validation set.

The comparison graph retains all explicit nuclei and their elements. Charge,
spin, bond order, stereochemistry, coordinates and basis choice cannot by
themselves produce novelty credit. This deliberately conservative convention can
exclude distinctions that matter elsewhere in chemistry. It must not be described
as a complete chemical-identity model.

Exact element-preserving graph isomorphism supplies a checked atom bijection.
Weisfeiler-Lehman hashes are diagnostic prefilters only; equal hashes do not prove
equal graphs. A radius-derived proximity edge is not an assigned chemical bond.
Record cutoff sensitivity and borderline distances before reviewing connectivity.

Before requesting novelty credit, review the primary molecule and its environment
separately. Pass the complete original exclusion catalog, including possible
component scopes, to `assess_primary_molecule_novelty`. Then add every accepted
external primary molecule to that catalog before judging the next candidate.
Any unresolved baseline with the same nuclear composition blocks automatic
credit. Never set `chemical_scope_reviewed` merely because a proximity graph was
generated or an input filename named a molecule.

## Frozen original collection

The batch reader independently checks the original FCHK hash, nucleus counts,
elements and coordinates against the frozen campaign. It records all 273 terminal
states before whole-set comparison. Collection failure, radius sensitivity,
unreviewed chemical scope and final scientific acceptance are separate states.

GI-001 used radius scales 1.15, 1.25 and 1.35, a 1e-6 Angstrom numerical boundary
tolerance, and one fixed rigid transform plus one deterministic atom permutation
per case. Its 819 graph checks are sampled controls, not exhaustive molecular,
orbital or rendering invariance proofs. The manifest and per-case evidence must
remain immutable after collection; follow-up interpretation belongs in a new
review artifact.

The OLD-098/099 control illustrates why one distance graph is insufficient:
at scale 1.35, the square input's two C-C diagonals are included, while the
rectangle input's diagonals are excluded. Keep this recorded control disagreement;
do not tune the cutoff to turn the review into a perfect-score test.

## Reproducing the checks

The collector requires Python, NumPy, NetworkX and RDKit. GI-001 recorded Python
3.12, NumPy 2.5.2, NetworkX 3.6.1 and RDKit 2026.03.5. Its isolated package lock,
wheel hashes and imported module paths are in the round evidence. These packages
are separate from the Gaussian/IOData/GBasis numerical-reference environment.
Reproduction must use the recorded runtime and source hashes, not a blind upgrade
of another running round's dependencies.

Run `python -m unittest test_structure_identity -v` from this directory in that
isolated environment. A prepared, version-frozen round runs
`structure_identity_batch.py preflight --root <round>` followed by
`structure_identity_batch.py collect --root <round>` through the common resource
supervisor. Whole-set review follows only after every original case terminates.
The nine tests cover element matching, a hash collision, explicit hydrogen
connectivity, ring-block junctions, coordinate transformations, cutoff ambiguity,
invalid geometry, charge-state exclusion and unresolved-scope exclusion.

References: [NetworkX exact isomorphism](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.isomorphism.is_isomorphic.html),
[RDKit molecular input and explicit hydrogens](https://www.rdkit.org/docs/GettingStartedInPython.html).
