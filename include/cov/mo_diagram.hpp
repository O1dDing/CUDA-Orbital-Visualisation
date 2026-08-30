#pragma once

#include "cov/model.hpp"
#include "cov/orbital_view.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace cov {

enum class MODiagramMode {
    ValenceCentral = 0,
};

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
    std::string label;
    AnnotationSource source = AnnotationSource::Unavailable;
    double confidence = 0.0;
    bool heuristic = false;
};

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

struct MODiagramData {
    DiagramPlan plan;
    FrontierOrbitals frontier;
    DiagramSelectionPlan selection;
    std::vector<MODiagramLevel> levels;
    std::vector<OrbitalMetadata> metadata;
    std::vector<OrbitalAnnotation> annotations;
    std::vector<PiInteractionDescriptor> pi_interactions;
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
