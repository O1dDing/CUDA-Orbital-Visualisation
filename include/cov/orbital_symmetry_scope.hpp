#pragma once

#include "cov/model.hpp"
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace cov {
struct LocalIrrepAssignment;
struct PiPartnerAssessment;
struct LocalAngularProjection;
class LocalAngularProjectionWorkspace;
struct LigandFieldEnvironment;

enum class OrbitalSymmetryOrigin {
    Unavailable,
    Producer,
    MolecularOperations,
    LocalMetricProjection,
    LocalDimensionCandidate,
    PiPartnerCandidate,
    SpinCounterpartCandidate,
    MixedMembers,
    LocalMetricDecomposition,
};

// A label and its scope form one value. A local/template/candidate label is
// never a replacement for the full-orbital label stored on MolecularOrbital.
// Candidate sources are immutable snapshots; target MO identities remain here.
struct OrbitalSymmetryExplanation {
    std::string label;
    OrbitalSymmetryOrigin origin = OrbitalSymmetryOrigin::Unavailable;
    std::string point_group;
    std::string point_group_basis = "unresolved";
    std::vector<std::size_t> orbital_indices;
    std::vector<std::size_t> atom_indices;
    std::array<double,9> rotation_reference_to_input{};
    bool axes_available = false;
    std::string geometry_id;
    double geometry_angular_rms = std::numeric_limits<double>::quiet_NaN();
    double geometry_shape_measure = std::numeric_limits<double>::quiet_NaN();
    std::shared_ptr<const LocalIrrepAssignment> local_assignment;
    std::shared_ptr<const LocalAngularProjection> local_decomposition;
    std::string assessment_detail;
    std::shared_ptr<const PiPartnerAssessment> pi_partner_evidence;
    std::shared_ptr<const DerivedOrbitalSymmetryAssignment> molecular_assignment;
    std::shared_ptr<const OrbitalSymmetryExplanation> candidate_source;
    double candidate_score = std::numeric_limits<double>::quiet_NaN();
    std::string source_path;
    std::size_t source_line_begin = 0;
    std::size_t source_line_end = 0;
    std::size_t source_job_segment = 0;
    // Context reported by the producer is kept without assuming that a printed
    // MO label belongs to either the detected group or its Abelian subgroup.
    std::string producer_detected_group;
    std::string producer_abelian_group;
};

[[nodiscard]] const char* orbital_symmetry_origin_name(OrbitalSymmetryOrigin) noexcept;
[[nodiscard]] bool orbital_symmetry_is_local(const OrbitalSymmetryExplanation&) noexcept;
[[nodiscard]] bool orbital_symmetry_is_candidate(const OrbitalSymmetryExplanation&) noexcept;
[[nodiscard]] OrbitalSymmetryExplanation molecular_orbital_symmetry(
    const Wavefunction&, std::size_t orbital_index);
[[nodiscard]] OrbitalSymmetryExplanation local_orbital_symmetry(
    const LocalIrrepAssignment&, const LigandFieldEnvironment&,
    std::span<const std::size_t> orbital_indices);
[[nodiscard]] OrbitalSymmetryExplanation evaluate_local_orbital_symmetry(
    const LocalAngularProjectionWorkspace&,const LigandFieldEnvironment&,
    std::span<const std::size_t> orbital_indices,int admitted_shell_family);
[[nodiscard]] OrbitalSymmetryExplanation candidate_orbital_symmetry(
    const OrbitalSymmetryExplanation& source, OrbitalSymmetryOrigin,
    std::span<const std::size_t> target_indices, double score);
[[nodiscard]] bool compatible_symmetry_scopes(
    const OrbitalSymmetryExplanation&, const OrbitalSymmetryExplanation&) noexcept;
[[nodiscard]] OrbitalSymmetryExplanation aggregate_molecular_symmetry(
    const Wavefunction&, std::span<const std::size_t> orbital_indices);
// Human-facing compact text always marks local and candidate explanations.
// More detailed values use the same explanation object in UI and exports.
[[nodiscard]] std::string orbital_symmetry_compact_text(const OrbitalSymmetryExplanation&);
[[nodiscard]] std::string orbital_symmetry_json(const OrbitalSymmetryExplanation&);
}
