#pragma once

#include "cov/model.hpp"
#include "cov/orbital_view.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace cov {

struct MODiagramOptions {
    EnergyUnit energy_unit = EnergyUnit::Hartree;
    DegeneracySettings degeneracy{};
    OrbitalFilterSettings filter{};
    std::size_t selected_index = 0;
    std::size_t neighbourhood = 16;
    int width = 1200;
    int height = 900;
    bool include_hidden_in_metadata = true;
};

struct MODiagramLevel {
    OrbitalMetadata metadata;
    ElectronGlyphs electrons;
    bool homo = false;
    bool lumo = false;
};

struct MODiagramData {
    DiagramPlan plan;
    FrontierOrbitals frontier;
    std::vector<MODiagramLevel> levels;
    std::vector<OrbitalMetadata> metadata;

    // Presentation metadata. `textbook_diatomic` permits a qualitative side-AO
    // reference while central MO energies remain quantitative wavefunction data.
    bool textbook_diatomic = false;
    bool textbook_h2 = false;
    std::string species_label;
};

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
