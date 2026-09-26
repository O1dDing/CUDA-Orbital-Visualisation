#pragma once

#include "cov/nbo_integration.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/ui.hpp"
#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace cov::ui {

struct NboAomoFragmentGroup {
    std::size_t id=0;
    std::set<std::size_t> atoms; // zero-based canonical atom identities
    bool suggested=false;
};

struct NboAomoNode {
    std::string id,label,detail,energy_semantics;
    std::optional<NboOrbitalRef> orbital;
    std::optional<std::size_t> canonical_index;
    std::optional<std::size_t> fragment_group_id;
    std::vector<std::size_t> atoms;
    std::optional<double> energy_hartree,occupation;
    std::optional<double> metric_norm2;
    float x=0,y=0,width=160,height=24;
    bool group_header=false,available=true,composition_available=true;
};
struct NboAomoEdge {
    std::string id;
    std::size_t source_node=0,target_node=0;
    double coefficient=0;
    std::optional<double> weight;
    NboSource source;
};
struct NboAomoViewSnapshot {
    std::string id,integration_id,mo_snapshot_id,capability_status,capability_detail;
    std::string mo_energy_axis_mode,mo_energy_axis_detail;
    NboOrbitalKind basis_kind=NboOrbitalKind::NAO;
    std::size_t focused_canonical_index=0;
    std::optional<double> focused_projection_weight,focused_projection_residual_norm;
    std::vector<std::size_t> central_mo_indices;
    std::vector<NboAomoNode> nodes;
    std::vector<NboAomoEdge> edges;
    std::optional<NboOrbitalSelection> selection;
    std::vector<NboAomoFragmentGroup> fragment_groups;
    std::vector<std::string> sum_component_ids;
    std::size_t hidden_basis_count=0,hidden_mo_count=0;
    std::size_t hidden_class_count=0;
    bool show_core=false,show_rydberg=false;
    float zoom=1,pan_x=0,pan_y=0;
};

struct NboAomoUIState {
    NboOrbitalKind basis_kind=NboOrbitalKind::NAO;
    float zoom=1,pan_x=0,pan_y=0;
    std::set<std::size_t> collapsed_atoms; // zero-based canonical
    std::set<std::size_t> collapsed_levels; // index into the frozen central levels
    bool show_core=false,show_rydberg=false;
    std::set<std::size_t> draft_fragment_atoms;
    std::optional<std::size_t> editing_fragment_id;
    std::set<std::string> sum_component_ids;
    std::optional<std::size_t> sum_canonical_index;
    NboOrbitalKind sum_basis_kind=NboOrbitalKind::NAO;
    std::vector<std::string> ambiguous_edge_ids;
    std::vector<NboAomoFragmentGroup> fragment_groups;
    bool suggested_fragments_initialized=false;
    std::size_t next_fragment_id=1;
    std::optional<std::size_t> focused_canonical_index;
    std::optional<std::size_t> last_inspected;
    std::optional<NboOrbitalSelection> pending_selection;
    std::optional<NboOrbitalSelection> selection; // last applied 3D identity; root updates
    bool export_requested=false;
    bool show_full_numeric=false;
    std::string source_id,status,export_status;
    std::array<char,2048> export_path{};
    std::shared_ptr<const NboAomoViewSnapshot> drawn_snapshot;
    std::uint64_t revision=0;
};

// Returns true only when the whole AO/NAO--MO view is scientifically available.
// Otherwise the caller should draw its existing Gaussian-only MO diagram.
bool draw_nbo_aomo_diagram(NboAomoUIState& state,const NboIntegration& integration,
                           const Wavefunction& canonical,const MODiagramViewSnapshot& diagram,
                           Language language,float ui_scale);

struct NboAomoExportResult {
    bool json=false,csv=false,svg=false,png=false;
    std::filesystem::path json_path,csv_path,svg_path,png_path;
    std::string error;
};
NboAomoExportResult export_nbo_aomo_bundle(const NboAomoViewSnapshot& snapshot,
    const NboIntegration& integration,const std::filesystem::path& base);

// Picks an occupied, non-core NBO associated with the current central MO when
// matrix evidence exists. No selection is made when that association is absent.
std::optional<std::size_t> suggest_initial_nbo_index(const NboIntegration& integration,
    const MODiagramViewSnapshot* diagram);

} // namespace cov::ui
