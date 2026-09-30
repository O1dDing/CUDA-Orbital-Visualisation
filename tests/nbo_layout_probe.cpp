#include "cov/chemistry_route.hpp"
#include "cov/nbo_aomo_ui.hpp"
#include "cov/wavefunction_io.hpp"
#include <imgui.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

// Compact observation of the production draw path for external case replays.
// This is not a scientific oracle: raw-file checks and ordinary-window checks
// remain separate from these exported layout measurements.
int main(int argc,char** argv) {try {
    if(argc<4)throw std::runtime_error("Usage: cov_nbo_layout_probe canonical.fchk analysis_directory output.json [export_base]");
    const auto start=std::chrono::steady_clock::now();
    auto wf=cov::parse_wavefunction(argv[1]);
    const auto original=wf;
    const auto found=cov::discover_nbo_inputs({argv[2]});
    if(found.candidates.size()!=1)throw std::runtime_error("Expected one analysis candidate");
    const auto integration=cov::read_nbo_integration(wf,found.candidates.front());
    const auto route=cov::route_chemistry(wf,&integration);
    cov::ui::NboAomoUIState state;
    if(!cov::ui::prepare_nbo_aomo_state(state,integration,wf))
        throw std::runtime_error("Unified graph unavailable: "+state.status);
    cov::MODiagramOptions options;
    options.routed=&route;options.routed_identity=integration.id;
    for(std::size_t i=0;i<wf.orbitals.size();++i)
        if(wf.orbitals[i].occupation>0)options.selected_index=i;
    if(state.names)for(const auto& name:state.names->canonical)
        options.display_names.push_back(name.label);
    std::cerr<<"Building central diagram\n";
    const auto data=cov::build_mo_diagram_data(wf,options);
    const auto graph=cov::make_mo_diagram_view_snapshot(data,options,options.selected_index,"case-layout-probe");
    IMGUI_CHECKVERSION();ImGui::CreateContext();
    auto& io=ImGui::GetIO();io.IniFilename=nullptr;io.DisplaySize={1280,1800};io.DeltaTime=1.0f/60;
    std::cerr<<"Loading fonts\n";
    if(!cov::ui::configure_fonts())throw std::runtime_error("Font setup failed");
    const auto draw=[&](cov::ui::Language language) {
        ImGui::NewFrame();ImGui::SetNextWindowPos({10,10},ImGuiCond_Always);
        ImGui::SetNextWindowSize({1200,1700},ImGuiCond_Always);
        ImGui::Begin("case replay",nullptr,ImGuiWindowFlags_NoSavedSettings);
        const bool ok=cov::ui::draw_nbo_aomo_diagram(state,integration,wf,graph,language,1.0f);
        ImGui::End();ImGui::Render();
        if(!ok||!state.drawn_snapshot)throw std::runtime_error("Draw did not provide a snapshot");
    };
    std::cerr<<"Drawing unified snapshot\n";
    draw(cov::ui::Language::English);
    bool preserved=wf.ao_overlap==original.ao_overlap&&wf.orbitals.size()==original.orbitals.size();
    for(std::size_t i=0;i<wf.orbitals.size();++i) {
        const auto& a=wf.orbitals[i];const auto& b=original.orbitals[i];
        preserved=preserved&&a.coefficients==b.coefficients&&a.gaussian_source_coefficients==b.gaussian_source_coefficients
            &&a.energy_hartree==b.energy_hartree&&a.occupation==b.occupation&&a.spin==b.spin
            &&a.source_orbital_index==b.source_orbital_index;
    }
    const auto& view=*state.drawn_snapshot;
    std::filesystem::path output(argv[3]);
    if(output.has_parent_path())std::filesystem::create_directories(output.parent_path());
    std::ofstream out(output);out<<std::setprecision(17)<<std::boolalpha;
    out<<"{\"canonical_preserved\":"<<preserved<<",\"canonical_count\":"<<wf.orbitals.size()
       <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
       <<",\"width\":"<<view.canvas_width<<",\"height\":"<<view.canvas_height
       <<",\"numeric_span\":"<<view.numeric_span<<",\"hidden_class_count\":"<<view.hidden_class_count
       <<",\"energies\":[";
    bool comma=false;
    for(const auto& e:state.salc_model->energies) {
        if(comma)out<<',';comma=true;
        out<<"{\"spin\":"<<std::quoted(cov::nbo_spin_name(e.spin))<<",\"status\":"<<std::quoted(e.status)
           <<",\"available\":"<<e.available<<",\"canonical_same_operator\":"<<e.canonical_same_operator
           <<",\"printed_operator_verified\":"<<e.printed_operator_verified<<'}';
    }
    out<<"],\"nodes\":[";comma=false;
    for(const auto& n:view.nodes) {
        if(comma)out<<',';comma=true;
        out<<"{\"id\":"<<std::quoted(n.id)<<",\"label\":"<<std::quoted(n.label)
           <<",\"lane\":"<<static_cast<int>(n.lane)<<",\"x\":"<<n.x<<",\"y\":"<<n.y
           <<",\"width\":"<<n.width<<",\"height\":"<<n.height
           <<",\"label_x\":"<<n.label_x<<",\"label_y\":"<<n.label_y
           <<",\"label_width\":"<<n.label_width<<",\"label_height\":"<<n.label_height
           <<",\"quantitative\":"<<n.quantitative_energy<<",\"group_header\":"<<n.group_header
           <<",\"group\":"<<std::quoted(n.display_group_id)<<",\"members\":"<<n.shell_member_count
           <<",\"shell_label\":"<<std::quoted(n.shell_label)<<",\"offset\":"<<n.display_offset_y
           <<",\"energy\":";
        if(n.energy_hartree)out<<*n.energy_hartree;else out<<"null";
        out<<",\"display_energy\":";
        if(n.display_energy_hartree)out<<*n.display_energy_hartree;else out<<"null";
        out<<",\"salc_index\":";
        if(n.salc_index)out<<*n.salc_index;else out<<"null";
        out<<'}';
    }
    out<<"],\"captions\":[";comma=false;
    for(const auto& c:view.captions){if(comma)out<<',';comma=true;out<<std::quoted(c.text);}
    out<<"]}";out.close();if(!out)throw std::runtime_error("Output write failed");
    if(argc>4) {
        const cov::ui::Language languages[]={cov::ui::Language::English,cov::ui::Language::ChineseSimplified,
            cov::ui::Language::Japanese,cov::ui::Language::French};
        for(std::size_t i=0;i<4;++i) {
            draw(languages[i]);
            const auto exported=cov::ui::export_nbo_aomo_bundle(*state.drawn_snapshot,integration,
                std::filesystem::path(std::string(argv[4])+"-"+std::to_string(i)));
            if(!exported.svg||!exported.png||!exported.json||!exported.csv)
                throw std::runtime_error("Bundle export failed: "+exported.error);
        }
    }
    ImGui::DestroyContext();return preserved?0:2;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
