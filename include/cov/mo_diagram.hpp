#pragma once

#include "cov/model.hpp"
#include "cov/orbital_view.hpp"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace cov {

enum class MODiagramMode {
    ValenceCentral = 0,
    DelocalisedPiFamilyOnly,
    MulticentreActiveSpaceOnly,
};

// Stable machine and human labels shared by the live view and every export
// format.  The built diagram is authoritative because a requested compact
// mode may legitimately fall back to the valence view when no active space is
// supported by the wavefunction.
[[nodiscard]] const char* mo_diagram_mode_name(MODiagramMode mode) noexcept;
[[nodiscard]] const char* mo_diagram_mode_title(MODiagramMode mode) noexcept;

// Compact main-group conjugated systems are clearer as their complete
// delocalised-pi active space.  Transition-metal ligand-field cases retain the
// valence-central diagram even if a ligand happens to contain a pi family.
[[nodiscard]] MODiagramMode preferred_compact_mo_diagram_mode(
    const Wavefunction& wavefunction,
    bool compact) noexcept;

enum class EnergyAxisMode {
    Linear = 0,
    NonlinearFocus,
};

struct EnergyAxisKnot {
    double energy_hartree = 0.0;
    double coordinate = 0.0;
};

struct EnergyTransform {
    EnergyAxisMode mode = EnergyAxisMode::NonlinearFocus;
    double focus_hartree = 0.0;
    double scale_hartree = 1.0e-4;
    double minimum_gap_weight = 0.070;
    std::vector<EnergyAxisKnot> knots;
};

[[nodiscard]] const char* energy_axis_mode_name(EnergyAxisMode mode) noexcept;
[[nodiscard]] const char* energy_transform_name(EnergyAxisMode mode) noexcept;
[[nodiscard]] EnergyTransform build_energy_transform(
    const std::vector<double>& energies_hartree,
    EnergyAxisMode mode,
    double minimum_gap_weight = 0.070);
[[nodiscard]] double energy_display_coordinate(
    double energy_hartree,
    const EnergyTransform& transform) noexcept;
[[nodiscard]] double energy_from_display_coordinate(
    double coordinate,
    const EnergyTransform& transform) noexcept;

enum class AnnotationSource {
    Direct,
    ParsedLabel,
    Derived,
    Heuristic,
    Unavailable,
};

enum class BondingClass {
    Unclassified,
    Bonding,
    Nonbonding,
    Antibonding,
};

struct MulticentreDescriptor {
    bool available = false;
    std::size_t centres = 0;
    double electrons = 0.0;
    std::vector<std::size_t> atom_indices;
    std::string label;
    std::size_t channel_count = 0;
    std::string source_subspace_id;
    double source_subspace_electron_count = 0.0;
    AnnotationSource source = AnnotationSource::Unavailable;
    double confidence = 0.0;
    bool heuristic = false;
};

struct DelocalisedPiDescriptor {
    bool available = false;
    std::size_t participating_atoms = 0;
    double participating_electrons = 0.0;
    std::vector<std::size_t> atom_indices;
    std::vector<std::size_t> orbital_indices;
    std::string family_id;
    // `available` says that a delocalised-pi family is known.  Topology is a
    // separate claim: parsed producer labels can establish the family without
    // establishing channel directions or cyclicity.
    bool topology_available = false;
    DelocalisedPiTopology topology = DelocalisedPiTopology::Unknown;
    std::size_t orientation_channels = 0;
    bool cyclic_topology = false;
    std::vector<PiOrientationChannel> orientation_channel_details;
    std::string label;
    AnnotationSource source = AnnotationSource::Unavailable;
    double confidence = 0.0;
    bool heuristic = false;
};

// The same compact topology suffix is used by the live diagram and bitmap/
// vector exports. An empty suffix means that the family is known but a
// concrete topology is not.
[[nodiscard]] inline const char* compact_pi_topology_code(
    const DelocalisedPiTopology topology) noexcept {
    switch (topology) {
        case DelocalisedPiTopology::Path: return "P";
        case DelocalisedPiTopology::Cycle: return "C";
        case DelocalisedPiTopology::BranchedResonance: return "B";
        case DelocalisedPiTopology::Spiro: return "S";
        case DelocalisedPiTopology::HapticMetal: return "H";
        case DelocalisedPiTopology::SymmetryDirectSum: return "D";
        default: return "?";
    }
}
[[nodiscard]] inline std::string compact_pi_topology_suffix(
    const DelocalisedPiDescriptor& descriptor) {
    if (!descriptor.topology_available) return {};
    return std::to_string(descriptor.orientation_channels)+"ch "+
           compact_pi_topology_code(descriptor.topology);
}

struct OrbitalAnnotation {
    // Machine-friendly canonical family names: sigma / pi / delta / phi / unavailable.
    std::string family = "unavailable";
    BondingClass bonding_class = BondingClass::Unclassified;
    AnnotationSource family_source = AnnotationSource::Unavailable;
    AnnotationSource bonding_source = AnnotationSource::Unavailable;
    double family_confidence = 0.0;
    double bonding_confidence = 0.0;
    MulticentreDescriptor multicentre;
    DelocalisedPiDescriptor delocalised_pi;
    bool heuristic = false;
};

struct SymmetryNotation {
    std::string base;
    std::string subscript;
    std::string superscript;
    std::string raw;
};

[[nodiscard]] SymmetryNotation parse_symmetry_notation(std::string_view raw);
[[nodiscard]] std::string format_symmetry_unicode(std::string_view raw);

struct DiagramSelectionPlan {
    std::vector<std::size_t> included_indices;
    std::size_t hidden_count = 0;
    std::size_t valence_occupied_count = 0;
    std::size_t frontier_virtual_count = 0;
    // Complete protected manifolds may legitimately exceed the visual row
    // target. This records only that unavoidable excess; it is never a
    // licence to keep arbitrary unprotected rows.
    std::size_t protected_overflow_count = 0;
    std::string summary;
};

struct MODiagramOptions {
    MODiagramMode mode = MODiagramMode::ValenceCentral;
    EnergyUnit energy_unit = EnergyUnit::Hartree;
    EnergyAxisMode energy_axis_mode = EnergyAxisMode::NonlinearFocus;
    DegeneracySettings degeneracy{};
    OrbitalFilterSettings filter{};
    std::size_t selected_index = 0;

    std::size_t neighbourhood = 12;
    std::size_t max_levels = 0;
    std::size_t max_virtual_levels = 10;
    bool hide_ligand_centred_intermediates = false;

    double nonlinear_minimum_gap_weight = 0.070;

    // A small same-symmetry pi splitting is shown as an approximately
    // nonbonding weak-coupling level instead of forcing a donor/acceptor
    // assignment.  The retained member is the one with the larger metal
    // valence contribution, so the reduced ligand-field diagram keeps the
    // chemically relevant d-level rather than an arbitrary canonical MO.
    double weak_pi_split_hartree = 0.020;
    double weak_metal_ligand_overlap = 0.025;

    int width = 1200;
    int height = 900;
    bool include_hidden_in_metadata = true;
};

enum class PiInteractionKind {
    Donor,
    Acceptor,
    Coupled,
    WeakNearNonbonding,
};

struct PiInteractionDescriptor {
    std::size_t lower_level = 0;
    std::size_t upper_level = 0;
    std::vector<std::size_t> lower_orbitals;
    std::vector<std::size_t> upper_orbitals;
    std::string symmetry;
    PiInteractionKind kind = PiInteractionKind::Coupled;
    double splitting_hartree = 0.0;
    double confidence = 0.0;
    bool lower_visible = true;
    bool upper_visible = true;
    std::size_t retained_level = 0;
};

[[nodiscard]] const char* pi_interaction_kind_name(
    PiInteractionKind kind) noexcept;

struct MODiagramLevel {
    OrbitalMetadata metadata;
    OrbitalAnnotation annotation;
    OrbitalChemistry chemistry;
    ElectronGlyphs electrons;
    bool homo = false;
    bool lumo = false;
    double layout_energy_hartree = 0.0;

    // One row represents a complete degenerate subspace.  Group-level
    // quantities are averaged per member and therefore remain invariant to a
    // rotation of canonical MOs inside an exactly degenerate subspace.
    std::vector<std::size_t> member_indices;
    std::vector<ElectronGlyphs> member_electrons;
    // For a UDFT spatial-row view, each majority-spin member can retain the
    // matched minority-spin canonical MO.  This keeps selection highlighting
    // and exact member metadata available without drawing a duplicate row.
    std::vector<std::size_t> member_spin_counterparts;
    double energy_spread_hartree = 0.0;
    double total_occupation = 0.0;
    double metal_s_weight = 0.0;
    double metal_p_weight = 0.0;
    double metal_d_weight = 0.0;
    // p population on the atoms directly coordinated to the selected metal.
    // ligand_p_weight may additionally include the rest of those ligand
    // fragments (for example O in CO or N in CN).
    double direct_ligand_p_weight = 0.0;
    double ligand_p_weight = 0.0;
    double sigma_fraction = 0.0;
    double pi_fraction = 0.0;
    double metal_ligand_overlap = 0.0;
    bool raw_data_fallback = false;
    bool approximate_nonbonding = false;
};

// One visible row may represent an exactly-degenerate canonical-MO set and,
// after UDFT spatial collapse, its matched opposite-spin counterparts.
[[nodiscard]] inline bool mo_diagram_level_covers_orbital(
    const MODiagramLevel& level,
    const std::size_t orbital_index) noexcept {
    if (std::find(level.member_indices.begin(),level.member_indices.end(),
                  orbital_index)!=level.member_indices.end()) {
        return true;
    }
    if (std::find(level.member_spin_counterparts.begin(),
                  level.member_spin_counterparts.end(),orbital_index)!=
            level.member_spin_counterparts.end()) {
        return true;
    }
    return level.member_indices.empty() &&
           level.metadata.orbital_index==orbital_index;
}

struct LocalGeometryDiagramDescriptor {
    std::size_t centre_atom = 0;
    std::string geometry_id;
    std::string geometry_name;
    std::string point_group;
    std::vector<std::size_t> neighbour_atoms;
    double confidence = 0.0;
    double angular_rms = 0.0;
    double shape_measure = 0.0;
    double radial_cv = 0.0;
};

struct ElectronicStateDiagramMetadata {
    WavefunctionSource source = WavefunctionSource::Unknown;
    std::int32_t charge = 0;
    std::uint32_t multiplicity = 0;
    std::uint32_t alpha_electrons = 0;
    std::uint32_t beta_electrons = 0;
    DataProvenance charge_provenance = DataProvenance::Unavailable;
    DataProvenance multiplicity_provenance = DataProvenance::Unavailable;
    DataProvenance electron_counts_provenance = DataProvenance::Unavailable;
    ScfConvergenceStatus scf_convergence = ScfConvergenceStatus::Unavailable;
    DataProvenance scf_convergence_provenance = DataProvenance::Unavailable;
    WavefunctionStabilityStatus stability =
        WavefunctionStabilityStatus::Unavailable;
    DataProvenance stability_provenance = DataProvenance::Unavailable;
    std::string stability_detail;
    double spin_squared_before = 0.0;
    double spin_squared_after = 0.0;
    DataProvenance spin_squared_provenance = DataProvenance::Unavailable;
    std::vector<double> atomic_partial_charges;
    std::string atomic_partial_charge_scheme;
    DataProvenance atomic_partial_charge_provenance =
        DataProvenance::Unavailable;
    std::string point_group_detected;
    std::string point_group_used;
    DataProvenance point_group_provenance = DataProvenance::Unavailable;
    std::string source_title;
    std::string source_route;
    std::string enrichment_source;
};

struct MODiagramData {
    DiagramPlan plan;
    FrontierOrbitals frontier;
    DiagramSelectionPlan selection;
    std::vector<MODiagramLevel> levels;
    std::vector<OrbitalMetadata> metadata;
    std::vector<OrbitalAnnotation> annotations;
    std::vector<PiInteractionDescriptor> pi_interactions;
    ElectronicStateDiagramMetadata electronic_state;
    // Every unambiguous Mayer-supported CN2--CN10 centre, including main-group
    // and non-metal centres.  This is structural metadata and never creates a
    // transition-metal ligand-field diagram by itself.
    std::vector<LocalGeometryDiagramDescriptor> local_geometries;
    MODiagramMode mode = MODiagramMode::ValenceCentral;
    EnergyTransform energy_transform;

    // Local first-shell ligand-field information is kept separate from the
    // full molecular point group.  For an unrestricted calculation the
    // diagram can use majority-spin spatial representatives while retaining
    // the paired alpha/beta occupation and member metadata.
    std::string ligand_field_point_group;
    std::string ligand_field_geometry_id;
    std::string ligand_field_geometry_name;
    std::size_t ligand_field_coordination_number = 0;
    std::size_t ligand_field_metal_atom = 0;
    std::vector<std::size_t> ligand_field_ligand_atoms;
    double ligand_field_confidence = 0.0;
    double ligand_field_angular_rms = 0.0;
    double ligand_field_shape_measure = 0.0;
    double ligand_field_radial_cv = 0.0;
    bool spin_counterparts_collapsed = false;
    bool spin_counterparts_partial = false;
    std::size_t spin_counterpart_pair_count = 0;
    std::size_t spin_counterpart_unmatched_visible = 0;
};

[[nodiscard]] const char* annotation_source_name(AnnotationSource source) noexcept;
[[nodiscard]] const char* bonding_class_name(BondingClass value) noexcept;
[[nodiscard]] OrbitalAnnotation annotate_orbital(const MolecularOrbital& orbital);
[[nodiscard]] DiagramSelectionPlan build_valence_selection_plan(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::vector<OrbitalMetadata>& metadata);
[[nodiscard]] MODiagramData build_mo_diagram_data(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options);

struct MODiagramExportResult {
    bool svg = false;
    bool png = false;
    bool json = false;
    bool csv = false;
    std::filesystem::path svg_path;
    std::filesystem::path png_path;
    std::filesystem::path json_path;
    std::filesystem::path csv_path;
    std::string error;
};

[[nodiscard]] MODiagramExportResult export_mo_diagram_bundle(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::filesystem::path& base_path);
[[nodiscard]] bool write_mo_diagram_svg(
    const MODiagramData& data,
    const MODiagramOptions& options,
    const std::filesystem::path& path,
    std::string* error = nullptr);
[[nodiscard]] bool write_mo_diagram_png(
    const MODiagramData& data,
    const MODiagramOptions& options,
    const std::filesystem::path& path,
    std::string* error = nullptr);
[[nodiscard]] bool write_mo_diagram_json(
    const MODiagramData& data,
    const MODiagramOptions& options,
    const std::filesystem::path& path,
    std::string* error = nullptr);
[[nodiscard]] bool write_mo_diagram_csv(
    const MODiagramData& data,
    const MODiagramOptions& options,
    const std::filesystem::path& path,
    std::string* error = nullptr);

} // namespace cov
