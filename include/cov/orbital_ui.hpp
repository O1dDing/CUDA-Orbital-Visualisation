#pragma once

#include "cov/model.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/orbital_view.hpp"
#include "cov/ui.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace cov::ui {

struct OrbitalUIDiagramCache {
    const Wavefunction* wavefunction = nullptr;
    const MolecularOrbital* orbital_data = nullptr;
    std::size_t atom_count = 0;
    std::size_t orbital_count = 0;
    std::optional<MODiagramOptions> options;
    std::optional<MODiagramData> data;
};

struct OrbitalUIBrowserCache {
    const Wavefunction* wavefunction = nullptr;
    const MolecularOrbital* orbital_data = nullptr;
    std::size_t atom_count = 0;
    std::size_t orbital_count = 0;
    std::optional<DegeneracySettings> degeneracy;
    std::optional<OrbitalFilterSettings> filter;
    std::optional<FrontierOrbitals> frontier;
    std::vector<OrbitalMetadata> metadata;
};

struct OrbitalUIState {
    EnergyUnit energy_unit = EnergyUnit::Hartree;
    EnergyAxisMode energy_axis_mode = EnergyAxisMode::NonlinearFocus;
    DegeneracySettings degeneracy{};
    OrbitalFilterSettings filter{};
    bool grouped_labels = true;
    bool hide_ligand_centred_intermediates = true;
    int diagram_neighbourhood = 12;
    std::array<char, 96> search{};
    OrbitalUIBrowserCache browser_cache{};
    OrbitalUIDiagramCache diagram_cache{};
};

struct OrbitalUIActions {
    std::optional<std::size_t> select_orbital;
    bool export_diagram = false;
};

void draw_orbital_browser(const Wavefunction& wavefunction,
                          std::size_t selected_index,
                          OrbitalUIState& state,
                          Language language,
                          float ui_scale,
                          OrbitalUIActions& actions);

void draw_energy_diagram(const Wavefunction& wavefunction,
                         std::size_t selected_index,
                         OrbitalUIState& state,
                         Language language,
                         float ui_scale,
                         OrbitalUIActions& actions);

} // namespace cov::ui
