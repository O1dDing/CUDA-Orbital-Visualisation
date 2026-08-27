#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace cov {

enum class EnergyUnit {
    Hartree = 0,
    ElectronVolt,
    JoulePerMol,
    KilojoulePerMol,
    CaloriePerMol,
    KilocaloriePerMol,
};

[[nodiscard]] double convert_hartree(double energy_hartree, EnergyUnit unit) noexcept;
[[nodiscard]] const char* energy_unit_symbol(EnergyUnit unit) noexcept;
[[nodiscard]] std::string format_energy(double energy_hartree,
                                        EnergyUnit unit,
                                        int precision = 6);

enum class OrbitalFilterMode {
    AutoReasonable = 0,
    All,
    Occupied,
    Virtual,
    Core,
    Valence,
};

struct OrbitalFilterSettings {
    OrbitalFilterMode mode = OrbitalFilterMode::AutoReasonable;
    double occupation_threshold = 1.0e-4;
    double virtual_window_hartree = 1.5;
    double core_energy_cutoff_hartree = -1.5;
};

struct FrontierOrbitals {
    std::optional<std::size_t> homo;
    std::optional<std::size_t> lumo;
    std::optional<std::size_t> alpha_homo;
    std::optional<std::size_t> alpha_lumo;
    std::optional<std::size_t> beta_homo;
    std::optional<std::size_t> beta_lumo;
    bool separate_spin_sets = false;
};

[[nodiscard]] FrontierOrbitals find_frontier_orbitals(
    const std::vector<MolecularOrbital>& orbitals,
    double occupation_threshold = 1.0e-4);

enum class OrbitalRegion {
    Core,
    Valence,
    Virtual,
};

[[nodiscard]] OrbitalRegion classify_orbital_region(
    const MolecularOrbital& orbital,
    const OrbitalFilterSettings& settings) noexcept;

[[nodiscard]] bool orbital_visible(
    const std::vector<MolecularOrbital>& orbitals,
    std::size_t index,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings) noexcept;

[[nodiscard]] std::vector<std::size_t> visible_orbital_indices(
    const std::vector<MolecularOrbital>& orbitals,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings);

struct DegeneracySettings {
    double tolerance_hartree = 1.0e-5;
    bool require_same_spin = true;
};

struct OrbitalLabel {
    std::size_t orbital_index = 0;      // internal zero-based index
    std::size_t raw_mo_number = 1;      // Molden/user-facing one-based number
    std::size_t group_base_number = 1;  // first raw MO number in a degenerate set
    std::size_t group_member_index = 0; // 0=a, 1=b, ...
    std::size_t group_size = 1;
    std::string display_label;          // e.g. "17-a" or "18"
};

[[nodiscard]] std::vector<OrbitalLabel> build_orbital_labels(
    const std::vector<MolecularOrbital>& orbitals,
    const DegeneracySettings& settings = {});

[[nodiscard]] std::string degeneracy_suffix(std::size_t member_index);

struct SymmetrySummary {
    bool has_meaningful_labels = false;
    double labelled_fraction = 0.0;
    std::vector<std::string> distinct_labels;
};

[[nodiscard]] SymmetrySummary analyse_symmetry_labels(
    const std::vector<MolecularOrbital>& orbitals);

enum class DiagramClassification {
    Simple,
    SymmetryGrouped,
    SalcUnavailable,
};

struct DiagramPlan {
    DiagramClassification classification = DiagramClassification::Simple;
    bool complex_system = false;
    bool strict_salc_available = false;
    std::string machine_reason;
};

[[nodiscard]] DiagramPlan choose_diagram_plan(
    const Wavefunction& wavefunction,
    std::size_t complex_atom_threshold = 12,
    double symmetry_coverage_threshold = 0.60);

struct ElectronGlyphs {
    int alpha = 0;
    int beta = 0;
};

[[nodiscard]] ElectronGlyphs electron_glyphs_for_orbital(
    const MolecularOrbital& orbital,
    bool separate_spin_sets) noexcept;

struct OrbitalMetadata {
    std::size_t orbital_index = 0;
    std::size_t raw_mo_number = 1;
    std::string display_label;
    double energy_hartree = 0.0;
    float occupation = 0.0f;
    Spin spin = Spin::Alpha;
    std::string symmetry;
    std::size_t degeneracy_size = 1;
    OrbitalRegion region = OrbitalRegion::Virtual;
    bool visible = true;
    bool selected = false;
};

[[nodiscard]] std::vector<OrbitalMetadata> build_orbital_metadata(
    const Wavefunction& wavefunction,
    std::size_t selected_index,
    const DegeneracySettings& degeneracy,
    const OrbitalFilterSettings& filter);

} // namespace cov
