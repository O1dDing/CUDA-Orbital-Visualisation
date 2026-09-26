#pragma once

#include "cov/nbo.hpp"
#include "cov/ui.hpp"
#include <array>
#include <filesystem>
#include <optional>
#include <string>

namespace cov::ui {

struct NboUIState {
    std::array<char, 2048> path{};
    std::array<char, 2048> archive47{};
    std::array<char, 2048> aonbo{};
    std::array<char, 2048> nbomo{};
    std::array<char, 2048> export_path{};
    std::optional<NboDataset> dataset;
    std::string error;
    std::string export_status;
    std::size_t selected_orbital = 0;
    int analysis_segment = -1;
    double contribution_threshold = 0.01;
};

struct NboUIActions {
    bool attach = false;
    bool canonical_set = false;
    bool nbo_set = false;
    bool export_bundle = false;
    std::optional<std::size_t> selected_orbital;
};

NboUIActions draw_nbo_panel(NboUIState& state, Language language,
                            bool canonical_loaded, bool nbo_renderable,
                            bool nbo_active, float scale,
                            const Wavefunction* canonical, std::size_t canonical_index);
std::string nbo_glyph_seed(Language language);

// Exports one immutable dataset snapshot and a matched SVG/PNG report view.
// The report uses occupation only; diagonal Fock values retain their own label.
void export_nbo_bundle(const NboDataset& dataset, std::size_t selected_orbital,
                       const std::filesystem::path& base,
                       const Wavefunction* canonical, std::size_t canonical_index,
                       double contribution_threshold, bool nbo_active,
                       std::size_t rendered_orbital_index);

} // namespace cov::ui
