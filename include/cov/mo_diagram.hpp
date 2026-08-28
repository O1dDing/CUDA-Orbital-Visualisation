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

    double nonlinear_minimum_gap_weight = 0.070;

    int width = 1200;
    int height = 900;
    bool include_hidden_in_metadata = true;
};

struct MODiagramLevel {
    OrbitalMetadata metadata;
    OrbitalAnnotation annotation;
    OrbitalChemistry chemistry;
    ElectronGlyphs electrons;
    bool homo = false;
    bool lumo = false;
    double layout_energy_hartree = 0.0;
};

struct MODiagramData {
    DiagramPlan plan;
    FrontierOrbitals frontier;
    DiagramSelectionPlan selection;
    std::vector<MODiagramLevel> levels;
    std::vector<OrbitalMetadata> metadata;
    std::vector<OrbitalAnnotation> annotations;
    MODiagramMode mode = MODiagramMode::ValenceCentral;
    EnergyTransform energy_transform;
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
