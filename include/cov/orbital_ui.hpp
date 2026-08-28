#pragma once

#include "cov/model.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/orbital_view.hpp"
#include "cov/ui.hpp"

#include <array>
#include <cstddef>
#include <optional>

namespace cov::ui {

struct OrbitalUIState {
    EnergyUnit energy_unit = EnergyUnit::Hartree;
    EnergyAxisMode energy_axis_mode = EnergyAxisMode::NonlinearFocus;
    DegeneracySettings degeneracy{};
    OrbitalFilterSettings filter{};
    bool grouped_labels = true;
    int diagram_neighbourhood = 12;
    std::array<char, 96> search{};
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
