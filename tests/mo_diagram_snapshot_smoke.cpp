#include "cov/mo_diagram.hpp"
#include "cov/orbital_ui.hpp"

#include <imgui.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

namespace {
void require(bool ok, const char* message) {
    if (!ok) { std::cerr << message << '\n'; std::exit(1); }
}
std::string read(const std::filesystem::path& path) {
    std::ifstream in(path,std::ios::binary);
    return {std::istreambuf_iterator<char>(in),{}};
}
std::size_t count(const std::string& value,const std::string& needle) {
    std::size_t result=0,pos=0;
    while ((pos=value.find(needle,pos))!=std::string::npos) { ++result; pos+=needle.size(); }
    return result;
}
}

int main() {
    // A display contract fixture, not a physical electronic-state reference.
    cov::MODiagramData data;
    const double energy[]={-0.20,-0.19999,-0.18,-0.17999,-3.0};
    for (std::size_t i=0;i<5;++i) {
        cov::OrbitalMetadata raw;
        raw.orbital_index=i; raw.raw_mo_number=i+1;
        raw.energy_hartree=energy[i]; raw.occupation=1;
        raw.spin=i>=2?cov::Spin::Beta:cov::Spin::Alpha;
        data.metadata.push_back(raw);
    }
    data.annotations.resize(5);
    data.selection.included_indices={0,1};
    cov::MODiagramLevel level;
    level.metadata=data.metadata[0]; level.metadata.degeneracy_size=2;
    level.member_indices={0,1}; level.member_spin_counterparts={2,3};
    level.member_electrons={{1,1},{1,0}};
    level.layout_energy_hartree=-0.199995;
    level.energy_spread_hartree=0.00001;
    data.levels={level};
    data.energy_transform=cov::build_energy_transform({-0.4,0.1},cov::EnergyAxisMode::Linear);
    cov::MODiagramOptions options;
    options.selected_index=5; options.energy_axis_mode=cov::EnergyAxisMode::Linear;
    options.energy_unit=cov::EnergyUnit::ElectronVolt;
    options.include_hidden_in_metadata=false;
    auto snapshot=cov::make_mo_diagram_view_snapshot(data,options,3,"interactive-canvas");
    const auto members=cov::mo_diagram_member_views(snapshot.data,snapshot.data.levels[0]);
    require(members.size()==2 && !members[0].selected && members[1].selected &&
        members[1].inspected_orbital_index==3 && members[1].orbital_index==1,
        "beta inspection must select exactly its spatial member while retaining both identities");
    require(!snapshot.data.metadata[1].selected && snapshot.data.metadata[3].selected &&
        snapshot.data.metadata[3].energy_hartree==energy[3] &&
        snapshot.options.selected_index==5,"inspection must not alter raw metadata or the row-selection anchor");
    data.metadata[3].energy_hartree=99;
    data.levels.clear(); options.energy_unit=cov::EnergyUnit::Hartree;
    require(snapshot.data.metadata[3].energy_hartree==energy[3] &&
        snapshot.data.levels.size()==1 && snapshot.options.energy_unit==cov::EnergyUnit::ElectronVolt,
        "an export snapshot must own its frame values despite subsequent UI/cache mutations");

    const auto root=std::filesystem::temp_directory_path()/
        ("cov_snapshot_"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(root);
    const auto result=cov::export_mo_diagram_bundle(snapshot,root/"captured");
    require(result.svg&&result.png&&result.json&&result.csv,"snapshot export failed");
    const auto svg=read(result.svg_path),png=read(result.png_path),json=read(result.json_path),csv=read(result.csv_path);
    const auto& id=snapshot.data.view->id;
    require(svg.find(id)!=std::string::npos && png.find(id)!=std::string::npos &&
        json.find(id)!=std::string::npos && csv.find(id)!=std::string::npos,
        "all four files must identify the same snapshot");
    require(count(svg,"class=\"mo-member\"")==2 && count(svg,"data-selected=\"true\"")==1 &&
        svg.find("data-inspected-orbital-index=\"3\"")!=std::string::npos,
        "SVG must preserve real member count and beta inspection");
    require(json.find("\"spin_counterpart\":3")!=std::string::npos &&
        json.find("\"index\": 3")!=std::string::npos &&
        json.find("\"index\": 4")==std::string::npos &&
        json.find("\"all_members_energy_min_hartree\":-0.2")!=std::string::npos &&
        json.find("\"all_members_energy_max_hartree\":-0.17999")!=std::string::npos,
        "compact machine output must retain counterpart metadata and the actual energy range");
    auto hidden=cov::make_mo_diagram_view_snapshot(snapshot.data,snapshot.options,4);
    require(!hidden.data.levels[0].metadata.selected && hidden.data.metadata[4].selected,
        "inspecting a hidden MO must not manufacture a visible row");
    require(cov::write_mo_diagram_json(hidden.data,hidden.options,root/"hidden.json"),"hidden export failed");
    require(read(root/"hidden.json").find("\"index\": 4")!=std::string::npos,
        "the inspected hidden MO must retain raw metadata even with compact metadata export");
    auto invalid=cov::make_mo_diagram_view_snapshot(snapshot.data,snapshot.options,999);
    require(!invalid.data.view->inspected_orbital_index,"invalid inspection must not become a false index");
    auto raw=level; raw.member_indices.clear(); raw.member_spin_counterparts.clear();
    raw.metadata.degeneracy_size=3;
    require(cov::mo_diagram_member_views(snapshot.data,raw).size()==1,
        "a degeneracy-size field alone must not invent canonical members");

    IMGUI_CHECKVERSION(); ImGui::CreateContext();
    auto& io=ImGui::GetIO(); io.IniFilename=nullptr;
    io.DisplaySize={1280,1800}; io.DeltaTime=1.0f/60;
    io.Fonts->AddFontDefault(); io.Fonts->Build();
    cov::Wavefunction wf; wf.atoms.resize(1); wf.atoms[0].atomic_number=1;
    for (int i=0;i<8;++i) {
        cov::MolecularOrbital mo;
        mo.energy_hartree=-0.5+0.11*(i%4); mo.occupation=i%4<2?1.0f:0.0f;
        mo.spin=i>=4?cov::Spin::Beta:cov::Spin::Alpha;
        wf.orbitals.push_back(mo);
    }
    cov::ui::OrbitalUIState state;
    cov::ui::OrbitalUIActions actions;
    ImVec2 linear_target;
    auto frame=[&](std::size_t selected) {
        ImGui::NewFrame();
        ImGui::SetNextWindowPos({10,10},ImGuiCond_Always);
        ImGui::SetNextWindowSize({1200,1700},ImGuiCond_Always);
        ImGui::Begin("snapshot controls",nullptr,ImGuiWindowFlags_NoSavedSettings);
        const auto start=ImGui::GetCursorScreenPos();
        linear_target={start.x+ImGui::CalcTextSize(cov::ui::tr(cov::ui::Text::EnergyScale,
            cov::ui::Language::English)).x+ImGui::GetStyle().ItemSpacing.x+8,start.y+8};
        cov::ui::draw_energy_diagram(wf,selected,state,cov::ui::Language::English,1.0f,actions);
        ImGui::End(); ImGui::Render();
        require(actions.drawn_diagram!=nullptr,"live UI must provide its drawn snapshot");
    };
    frame(0); const auto first=actions.drawn_diagram;
    frame(0); require(first==actions.drawn_diagram,"unchanged view should reuse its immutable snapshot");
    frame(7); const auto beta=actions.drawn_diagram;
    require(beta!=first && beta->data.view->inspected_orbital_index==7 &&
        first->data.view->inspected_orbital_index==0,"frame selection must not mutate earlier snapshots");
    // Exercise the actual ImGui radio button. The same click frame must have
    // consistent widget mode, transform and export options.
    io.AddMousePosEvent(linear_target.x,linear_target.y); frame(7);
    io.AddMouseButtonEvent(0,true); frame(7);
    io.AddMouseButtonEvent(0,false); frame(7);
    require(state.energy_axis_mode==cov::EnergyAxisMode::Linear &&
        actions.drawn_diagram->options.energy_axis_mode==state.energy_axis_mode &&
        actions.drawn_diagram->data.energy_transform.mode==state.energy_axis_mode,
        "actual axis click must freeze the new transform in the same drawn frame");
    state.energy_unit=cov::EnergyUnit::ElectronVolt; frame(7);
    require(actions.drawn_diagram->options.energy_unit==cov::EnergyUnit::ElectronVolt,
        "energy-unit change must produce a matching export snapshot");
    ImGui::DestroyContext();
    // Only this uniquely created test directory is removed.
    std::filesystem::remove_all(root);
    std::cout << "mo_diagram_snapshot_smoke ok\n";
}
