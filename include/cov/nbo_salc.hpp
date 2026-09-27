#pragma once
#include "cov/nbo_integration.hpp"
#include "cov/symmetry.hpp"

namespace cov {
// A fixed, orthonormal basis derived from geometry and verified NAOs, never
// from the currently selected canonical MO. No display cutoff is applied here.
struct NboSalcOptions {
    double metric_tolerance=5e-5;
    double symmetry_tolerance=2e-4;
    double energy_tolerance_hartree=2e-5;
    double eigenvalue_cluster_tolerance=1e-7;
    std::size_t maximum_group_order=128;
};
struct NboSalcEnergyEvidence {
    NboSpin spin=NboSpin::Total;
    bool available=false;
    std::string status="unavailable", detail, input_units="hartree";
    double hermiticity_error=0, canonical_residual=0, eigenvalue_error_hartree=0;
    double projected_residual=0, nullspace_residual=0, fock_symmetry_error=0;
    double density_symmetry_error=0;
    bool density_symmetry_checked=false;
    std::size_t effective_rank=0;
    std::size_t overlap_numerical_rank=0,canonical_effective_rank=0,canonical_null_directions=0;
    double overlap_rank_threshold=0,outside_canonical_nao_norm=0,outside_canonical_fock_coupling=0;
    bool electronic_symmetry_verified=false;
    std::size_t canonical_columns_checked=0;
    NboSource source;
};
struct NboSalcFragment {
    std::string id, label;
    std::vector<std::size_t> atoms;
    // Layout suggestion only. Zero is left; one is right. Not a bond claim.
    unsigned side=0;
};
struct NboSalcSubspace {
    std::string id, label, fragment_id, detail;
    NboSpin spin=NboSpin::Total;
    std::vector<std::size_t> orbital_indices;
    std::vector<double> characters; // same order as model.operations
    std::size_t dimension=0, irrep_dimension=0, multiplicity=0;
    bool symmetry_verified=false;
    double closure_error=0, orthogonality_error=0, character_norm=0;
};
struct NboSalcOrbital {
    std::string id, label, fragment_id, subspace_id, type, angular, detail;
    NboSpin spin=NboSpin::Total;
    std::vector<std::size_t> atoms;
    std::vector<NboOrbitalTerm> terms; // existing NAO references; render Combination
    bool symmetry_adapted=false;
    std::size_t partner_index=0, partner_dimension=1;
    std::optional<double> energy_hartree, occupation;
    std::string energy_semantics="unavailable";
};
struct NboSalcLink {
    std::size_t side_index=0, canonical_index=0;
    double coefficient=0, weight=0;
};
struct NboSalcCoverage {
    std::size_t canonical_index=0;
    double weight_sum=0, residual_norm=0;
    bool available=false;
};
struct NboSalcModel {
    std::string dataset_id, canonical_fingerprint, cache_key;
    std::string point_group, used_group, status="unavailable", detail;
    bool available=false, group_verified=false;
    double group_closure_error=0, representation_error=0, orthogonality_error=0;
    std::vector<SymmetryOperation> operations;
    std::vector<NboSalcFragment> fragments;
    std::vector<NboSalcSubspace> subspaces;
    std::vector<NboSalcOrbital> orbitals;
    std::vector<NboSalcLink> links;
    std::vector<NboSalcCoverage> coverage;
    std::vector<NboSalcEnergyEvidence> energies;
    std::vector<std::string> diagnostics;
};
// Call once when attaching/changing the immutable dataset and cache the result.
// model.cache_key records dataset+canonical identity and options. No UI state,
// selected MO, energy window or reading-budget argument enters this function.
NboSalcModel build_nbo_salc_model(const Wavefunction& canonical,
    const NboIntegration& integration,const NboSalcOptions& options={});
NboOrbitalSelection nbo_salc_selection(const NboSalcModel& model,std::size_t side_index);
NboOrbitalSelection nbo_salc_component_selection(const NboSalcModel& model,const NboSalcLink& link);
std::string serialize_nbo_salc_json(const NboSalcModel& model);

// Reuses the existing AO angular transform; coefficients are row-major n x k.
// No dense n x n operation is retained, and empty means unsupported/invalid.
std::vector<double> apply_orbital_symmetry_operation(const Wavefunction& canonical,
    const SymmetryOperation& operation,const std::vector<double>& coefficients,
    std::size_t columns);
} // namespace cov
