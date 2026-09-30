#pragma once
#include "cov/nbo_integration.hpp"
#include <utility>

namespace cov {
// All indices are zero-based except literal producer NAO IDs. The local
// projectors, not a chosen SVD vector gauge, define each coupled channel.
struct NboPiCanonicalGroup {
    std::vector<std::size_t> members;
    NboSpin spin=NboSpin::Total;
    double centre_weight=0, ligand_weight=0;
    std::optional<double> occupation_per_mo;
    double cross_fock_min_hartree=0, cross_fock_max_hartree=0,
           cross_fock_mean_hartree=0;
    std::string character="unresolved";
};
struct NboPiAngularEvidence {
    std::size_t atom=0;
    std::string method="AO S-metric local angular projector in bond-axis frame";
    std::size_t rank=2;
    double p_metric_min_eigenvalue=0, max_sigma_leakage=0,
           centre_projection=0, partition_residual=0;
    std::array<double,3> pi_generalized_eigenvalues{};
};
struct NboPiDirectionProjection {
    std::string e2_id;
    std::size_t donor_nbo_id=0, acceptor_nbo_id=0;
    double donor_ligand_weight=0, donor_centre_weight=0,
           acceptor_ligand_weight=0, acceptor_centre_weight=0;
    NboSource source;
};
struct NboPiCoupling {
    std::string id, channel="pi", direction="unresolved", direction_evidence;
    NboSpin spin=NboSpin::Total;
    std::vector<std::size_t> centre_atoms, ligand_atoms;
    std::vector<std::size_t> centre_nao_ids, ligand_nao_ids;
    std::vector<double> singular_values_hartree;
    std::size_t coupled_rank=0;
    double centre_onsite_hartree=0, ligand_onsite_hartree=0;
    std::optional<double> centre_occupation, ligand_occupation;
    std::array<double,2> centre_onsite_range_hartree{}, ligand_onsite_range_hartree{};
    std::optional<std::array<double,2>> centre_occupation_range, ligand_occupation_range;
    std::string occupation_status="insufficient_evidence", occupation_reason;
    double operator_max_error_hartree=0, nao_orthogonality_error=0;
    double angular_leakage=0, minimum_centre_projection=0;
    std::optional<double> direction_donor_weight, direction_acceptor_weight;
    std::vector<NboPiAngularEvidence> angular_projector_evidence;
    std::vector<NboPiDirectionProjection> direction_projection_evidence;
    std::vector<NboPiCanonicalGroup> groups;
    // Empty for the broad atom-pi projector. Nonempty families retain full
    // ligand support and independently assigned internal pi/antipi character.
    std::string ligand_family;
    std::vector<std::size_t> ligand_family_nbo_ids;
    bool localized_family_verified=false;
    std::string operator_kind="canonical-same-operator";
    double canonical_operator_residual_hartree=0;
    double operator_validation_tolerance_hartree=2e-5;
    bool canonical_members_are_verified_shared_spatial=false;
    NboSource source;
};
struct NboPiCouplingAnalysis {
    std::string status="insufficient_evidence", reason;
    std::vector<NboPiCoupling> couplings;
    std::vector<NboSource> evidence;
};
// Same-operator, rotation-covariant local p-pi projections. No canonical or
// producer orbital is modified. A failed operator gate yields no pseudozero.
NboPiCouplingAnalysis analyse_nbo_pi_couplings(const Wavefunction& canonical,
    const NboIntegration& integration,
    const std::vector<std::pair<std::size_t,std::size_t>>& strong_connectivity);
} // namespace cov
