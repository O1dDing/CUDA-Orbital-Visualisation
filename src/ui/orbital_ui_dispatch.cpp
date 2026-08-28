#include "cov/orbital_ui.hpp"

#include <imgui.h>

#include <cstddef>

namespace cov::ui {

void draw_orbital_browser_legacy(const Wavefunction& wavefunction,
                                 std::size_t selected_index,
                                 OrbitalUIState& state,
                                 Language language,
                                 float ui_scale,
                                 OrbitalUIActions& actions);

void draw_energy_diagram_legacy(const Wavefunction& wavefunction,
                                std::size_t selected_index,
                                OrbitalUIState& state,
                                Language language,
                                float ui_scale,
                                OrbitalUIActions& actions);

namespace {

struct ProvenanceCounts {
    std::size_t producer=0;
    std::size_t derived=0;
    std::size_t unavailable=0;
};

ProvenanceCounts count_symmetry(const Wavefunction& wf) {
    ProvenanceCounts out;
    for (const auto& mo:wf.orbitals) {
        switch (mo.symmetry_provenance) {
            case DataProvenance::Producer: ++out.producer; break;
            case DataProvenance::Derived: ++out.derived; break;
            default: ++out.unavailable; break;
        }
    }
    return out;
}

ProvenanceCounts count_occupation(const Wavefunction& wf) {
    ProvenanceCounts out;
    for (const auto& mo:wf.orbitals) {
        switch (mo.occupation_provenance) {
            case DataProvenance::Producer: ++out.producer; break;
            case DataProvenance::Derived: ++out.derived; break;
            default: ++out.unavailable; break;
        }
    }
    return out;
}

void provenance_strip(const Wavefunction& wf) {
    const auto sym=count_symmetry(wf);
    const auto occ=count_occupation(wf);
    ImGui::TextDisabled(
        "%s | symmetry P/D/U %zu/%zu/%zu | occupation P/D/U %zu/%zu/%zu",
        wavefunction_source_name(wf.source),
        sym.producer,sym.derived,sym.unavailable,
        occ.producer,occ.derived,occ.unavailable);

    ImGui::TextDisabled(
        "density %s | overlap %s | bond order %s%s%s",
        data_provenance_name(wf.total_density_provenance),
        data_provenance_name(wf.ao_overlap_provenance),
        data_provenance_name(wf.bond_order_provenance),
        wf.point_group_detected.empty()?"":" | point group ",
        wf.point_group_detected.empty()?"":wf.point_group_detected.c_str());

    if (!wf.enrichment_source.empty()) {
        ImGui::TextDisabled("Gaussian LOG/OUT enrichment attached");
    }
    ImGui::Spacing();
}

} // namespace

void draw_orbital_browser(const Wavefunction& wavefunction,
                          std::size_t selected_index,
                          OrbitalUIState& state,
                          Language language,
                          float ui_scale,
                          OrbitalUIActions& actions) {
    provenance_strip(wavefunction);
    draw_orbital_browser_legacy(wavefunction,selected_index,state,language,ui_scale,actions);
}

void draw_energy_diagram(const Wavefunction& wavefunction,
                         std::size_t selected_index,
                         OrbitalUIState& state,
                         Language language,
                         float ui_scale,
                         OrbitalUIActions& actions) {
    draw_energy_diagram_legacy(wavefunction,selected_index,state,language,ui_scale,actions);
}

} // namespace cov::ui
