#pragma once

#include "cov/nbo.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/ui.hpp"
#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <set>
#include <vector>

namespace cov::ui {

struct NboFocusAtomGroup {
    std::size_t id = 0;
    std::set<std::size_t> atoms; // explicit producer one-based atom IDs
};

// Inspection controls only. These never change the central canonical MO diagram.
struct NboFocusUIState {
    std::optional<std::size_t> canonical_index;
    std::set<std::size_t> visible_atoms; // NBO producer one-based atom IDs
    std::set<std::size_t> ligand_atoms;  // unsaved group editor only
    std::vector<NboFocusAtomGroup> ligand_groups;
    std::size_t next_group_id = 1;
    std::set<std::string> hidden_shells;
    std::set<std::string> explicitly_included_shells;
    bool group_by_atom = false;
    bool group_ligands_by_l = false;
    bool show_core = false;
    bool show_rydberg = false;
    bool show_full_details = false;
    std::string diagram_id;
    std::optional<std::size_t> last_inspected;
    std::optional<std::size_t> pending_canonical_selection;
};

struct NboUIState {
    std::array<char, 2048> path{};
    std::array<char, 2048> archive47{};
    std::array<char, 2048> aonbo{};
    std::array<char, 2048> nbomo{};
    std::array<char, 2048> naomo{};
    std::array<char, 2048> aonao{};
    std::array<char, 2048> naonbo{};
    std::array<char, 2048> export_path{};
    std::optional<NboDataset> dataset;
    std::string error;
    std::string export_status;
    std::size_t selected_orbital = 0;
    int analysis_segment = -1;
    double contribution_threshold = 0.01;
    NboFocusUIState focus;
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
                            const Wavefunction* canonical, std::size_t canonical_index,
                            const MODiagramViewSnapshot* diagram);
std::string nbo_glyph_seed(Language language);

void draw_nbo_focus_view(NboFocusUIState& focus, const NboDataset& dataset,
                         const Wavefunction& canonical,
                         const MODiagramViewSnapshot& diagram,
                         Language language, float scale);

// Exports one immutable dataset snapshot and a matched SVG/PNG report view.
// The report uses occupation only; diagonal Fock values retain their own label.
void export_nbo_bundle(const NboDataset& dataset, std::size_t selected_orbital,
                       const std::filesystem::path& base,
                       const Wavefunction* canonical, std::size_t canonical_index,
                       double contribution_threshold, bool nbo_active,
                       std::size_t rendered_orbital_index,
                       const MODiagramViewSnapshot* diagram,
                       const NboFocusUIState* focus);

} // namespace cov::ui
