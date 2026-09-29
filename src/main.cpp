#include "cov/cuda_orbital.hpp"
#include "cov/numerical_diagnostics.hpp"
#include "cov/pi_topology_evidence.hpp"
#include <sstream>
#include "cov/file_dialog.hpp"
#include "cov/gl_api.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/molden_parser.hpp"
#include "cov/molecule_style.hpp"
#include "cov/nbo_ui.hpp"
#include "cov/nbo_molecular_overlay.hpp"
#include "cov/nbo_channels.hpp"
#include "cov/chemistry_route.hpp"
#include "cov/orbital_tracking.hpp"
#include "cov/nbo_aomo_labels.hpp"
#include "cov/orbital_ui.hpp"
#include "cov/ui.hpp"
#include "cov/volume_renderer.hpp"
#include "cov/validation.hpp"
#include "cov/viewer_layout.hpp"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl2.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::filesystem::path> g_dropped_paths;

enum class StatusKind {
    Ready,
    Parsing,
    Loaded,
    GridUpdated,
    Exported,
    Error,
};

void drop_callback(GLFWwindow*, const int count, const char** paths) {
    if (!paths) return;
    g_dropped_paths.clear();
    for (int i=0;i<count;++i) if(paths[i])
        g_dropped_paths.push_back(std::filesystem::u8path(paths[i]));
}

cov::GridBox make_grid_box(const cov::Wavefunction& wf, const float padding_bohr = 4.0f) {
    if (wf.atoms.empty()) return {};

    float min_x = static_cast<float>(wf.atoms.front().x);
    float min_y = static_cast<float>(wf.atoms.front().y);
    float min_z = static_cast<float>(wf.atoms.front().z);
    float max_x = min_x;
    float max_y = min_y;
    float max_z = min_z;

    for (const auto& atom : wf.atoms) {
        min_x = std::min(min_x, static_cast<float>(atom.x));
        min_y = std::min(min_y, static_cast<float>(atom.y));
        min_z = std::min(min_z, static_cast<float>(atom.z));
        max_x = std::max(max_x, static_cast<float>(atom.x));
        max_y = std::max(max_y, static_cast<float>(atom.y));
        max_z = std::max(max_z, static_cast<float>(atom.z));
    }

    const float cx = 0.5f * (min_x + max_x);
    const float cy = 0.5f * (min_y + max_y);
    const float cz = 0.5f * (min_z + max_z);
    const float extent = std::max({max_x - min_x, max_y - min_y, max_z - min_z});
    const float half = 0.5f * extent + padding_bohr;

    return {cx-half, cy-half, cz-half, cx+half, cy+half, cz+half};
}

std::size_t initial_orbital(const cov::Wavefunction& wf) {
    const auto frontier=cov::find_frontier_orbitals(wf.orbitals);
    return frontier.homo.value_or(0u);
}

const char* status_label(const StatusKind status, const cov::ui::Language language) {
    using cov::ui::Text;
    switch (status) {
        case StatusKind::Parsing: return cov::ui::tr(Text::Parsing, language);
        case StatusKind::Loaded: return cov::ui::tr(Text::Loaded, language);
        case StatusKind::GridUpdated: return cov::ui::tr(Text::GridUpdated, language);
        case StatusKind::Exported: return cov::ui::tr(Text::Exported, language);
        case StatusKind::Error: return cov::ui::tr(Text::Error, language);
        default: return cov::ui::tr(Text::Ready, language);
    }
}

cov::ui::Tone status_tone(const StatusKind status) {
    switch (status) {
        case StatusKind::Loaded:
        case StatusKind::GridUpdated:
        case StatusKind::Exported: return cov::ui::Tone::Success;
        case StatusKind::Error: return cov::ui::Tone::Danger;
        case StatusKind::Parsing: return cov::ui::Tone::Accent;
        default: return cov::ui::Tone::Neutral;
    }
}

std::filesystem::path path_from_utf8(const std::string& value) {
#ifdef _WIN32
    return std::filesystem::u8path(value);
#else
    return std::filesystem::path(value);
#endif
}

std::string path_to_utf8(const std::filesystem::path& path) {
#if defined(__cpp_lib_char8_t)
    const auto value = path.u8string();
    return std::string(reinterpret_cast<const char*>(value.data()), value.size());
#else
    return path.u8string();
#endif
}

void copy_path_to_buffer(const std::filesystem::path& path,
                         std::array<char, 2048>& buffer) {
    const std::string value = path_to_utf8(path);
    std::snprintf(buffer.data(), buffer.size(), "%s", value.c_str());
}

void disabled_wrapped(const char* text) {
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    ImGui::TextWrapped("%s", text);
    ImGui::PopStyleColor();
}

void metric_row(const char* label, const char* value) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    disabled_wrapped(label);
    ImGui::TableNextColumn();
    ImGui::TextWrapped("%s", value);
}

void push_recent(std::vector<std::filesystem::path>& recent,
                 const std::filesystem::path& path) {
    const auto normalized = path.lexically_normal();
    recent.erase(std::remove_if(recent.begin(), recent.end(), [&](const auto& existing) {
        return existing.lexically_normal() == normalized;
    }), recent.end());
    recent.insert(recent.begin(), normalized);
    if (recent.size() > 8) recent.resize(8);
}

const char* molecule_style_name(const cov::MoleculeStyle style,
                                const cov::ui::Language language) {
    return cov::ui::tr(style == cov::MoleculeStyle::StickDelocalisation
                           ? cov::ui::Text::StickDelocalisation
                           : cov::ui::Text::MediumBallStick,
                       language);
}

struct OrbitalAppearanceText {
    const char* material;
    const char* standard;
    const char* glass;
    const char* surface;
    const char* solid;
    const char* wire;
    const char* solid_wire;
    const char* auto_light;
};

OrbitalAppearanceText orbital_appearance_text(const cov::ui::Language language) {
    switch (language) {
        case cov::ui::Language::ChineseSimplified:
            return {"轨道材质", "标准", "玻璃", "表面模式", "实体", "线框", "实体 + 线框", "柔和自动打光"};
        case cov::ui::Language::Japanese:
            return {"軌道マテリアル", "標準", "ガラス", "表示モード", "ソリッド", "ワイヤー", "ソリッド + ワイヤー", "ソフト自動照明"};
        case cov::ui::Language::French:
            return {"Matériau orbital", "Standard", "Verre", "Mode de surface", "Solide", "Filaire", "Solide + filaire", "Éclairage automatique doux"};
        default:
            return {"Orbital material", "Standard", "Glass", "Surface mode", "Solid", "Wire", "Solid + Wire", "Soft automatic lighting"};
    }
}

const char* orbital_material_name(const cov::OrbitalMaterial material,
                                  const OrbitalAppearanceText& text) {
    return material == cov::OrbitalMaterial::Glass ? text.glass : text.standard;
}

const char* scene_text(cov::ui::Language language, const char* en, const char* zh,
                      const char* ja, const char* fr) {
    switch(language) {
        case cov::ui::Language::ChineseSimplified: return zh;
        case cov::ui::Language::Japanese: return ja;
        case cov::ui::Language::French: return fr;
        default: return en;
    }
}

const char* orbital_surface_name(const cov::OrbitalSurfaceMode mode,
                                 const OrbitalAppearanceText& text) {
    switch (mode) {
        case cov::OrbitalSurfaceMode::Wire: return text.wire;
        case cov::OrbitalSurfaceMode::SolidWire: return text.solid_wire;
        default: return text.solid;
    }
}

} // namespace

int main(int argc, char** argv) {
    try { cov::validation::configure(argc, argv); }
    catch (const std::exception& e) { std::fprintf(stderr,"Validation: %s\n",e.what()); return 2; }
    if (!glfwInit()) {
        std::fprintf(stderr, "GLFW initialisation failed\n");
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
    if (cov::validation::background()) {
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        glfwWindowHint(GLFW_FOCUSED, GLFW_FALSE);
    }

    GLFWwindow* window = glfwCreateWindow(
        1500, 940, "Chemical Orbital Visualiser", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "Unable to create OpenGL window\n");
        glfwTerminate();
        return 1;
    }

    glfwSetWindowSizeLimits(window, 640, 360, GLFW_DONT_CARE, GLFW_DONT_CARE);
    glfwMakeContextCurrent(window);
    if (cov::validation::active()) {
        glfwSetWindowSize(window, cov::validation::window_width(), cov::validation::window_height());
        glfwSetWindowTitle(window, "COV native validation");
    }
    glfwSwapInterval(1);
    glfwSetDropCallback(window, drop_callback);

    float x_scale = 1.0f;
    float y_scale = 1.0f;
    glfwGetWindowContentScale(window, &x_scale, &y_scale);
    const float ui_scale = std::clamp(std::max(x_scale, y_scale), 1.0f, 1.75f);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& startup_io = ImGui::GetIO();
    if (cov::validation::active()) {
        startup_io.IniFilename = nullptr;
        // Each native-plan step supplies an ordered event batch for this
        // frame. Do not defer part of it into the next GLFW polling batch.
        startup_io.ConfigInputTrickleEventQueue = false;
    }
    startup_io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    cov::ui::apply_theme(ui_scale);
    cov::ui::configure_fonts(16.5f * ui_scale);
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    int exit_code = 0;
    try {
        cov::VolumeRenderer renderer;
        cov::OrbitCamera camera;

        std::unique_ptr<cov::Wavefunction> wavefunction;
        std::optional<cov::Wavefunction> nbo_wavefunction;
        std::optional<cov::NboIntegration> integration;
        std::optional<cov::RoutedAnalysis> routed;
        std::unique_ptr<cov::NboSelectionView> inspection;
        cov::InteractionGraph semantic_graph;
        // Each overlay retains a separate signed field and CUDA evaluation.
        struct AdditionalField {
            std::unique_ptr<cov::VolumeRenderer> renderer;
            std::unique_ptr<cov::CudaOrbitalEvaluator> evaluator;
            std::size_t index=0;
        };
        std::vector<AdditionalField> additional_fields;
        cov::ui::NboUIState nbo_ui;
        bool nbo_active = false;
        std::size_t canonical_mo_index = 0;
        std::size_t nbo_mo_index = 0;
        std::optional<cov::OrbitalTrackingResult> frame_tracking;
        std::unique_ptr<cov::CudaOrbitalEvaluator> evaluator;
        cov::GridBox grid_box;

        cov::ui::Language language = cov::ui::Language::English;
        cov::ui::OrbitalUIState orbital_ui;
        orbital_ui.nbo_ui = &nbo_ui;
        cov::MoleculeRenderSettings molecule_render;
        cov::OrbitalMaterial orbital_material = cov::OrbitalMaterial::Standard;
        cov::OrbitalSurfaceMode orbital_surface_mode = cov::OrbitalSurfaceMode::Solid;
        StatusKind status = StatusKind::Ready;
        std::string status_detail;

        std::size_t mo_index = 0;
        std::optional<std::size_t> pending_mo_index;
        float isovalue = 0.03f;
        int resolution = 128;
        std::array<char, 2048> path_buffer{};
        std::filesystem::path current_file;
        std::vector<std::filesystem::path> recent_files;

        bool recompute = false;
        bool resize_and_recompute = false;

        auto active_wavefunction = [&]() -> const cov::Wavefunction* {
            if(inspection) return &inspection->wavefunction;
            return nbo_active && nbo_wavefunction ? &*nbo_wavefunction
                                                   : wavefunction ? &*wavefunction : nullptr;
        };
        auto active_view = [&]() {
            cov::ActiveOrbitalView view;
            if(inspection) {
                view.kind=cov::ActiveOrbitalKind::Inspection;
                view.source_id=inspection->selection.dataset_id;
                view.label=inspection->label;
                view.source_label=inspection->selection.label;
                view.display_name_evidence="verified selected combination and source terms";
                view.semantic_kind=inspection->selection.semantic_kind;
                view.group_id=inspection->selection.group_id;
                view.selection=inspection->selection;
                view.canonical_index=inspection->selection.target_canonical_index;
                view.rendered_index=mo_index;
                if(!inspection->selection.terms.empty())
                    view.spin=inspection->selection.terms.front().orbital.spin;
                view.source_spin=view.spin;
                view.spin_semantics="actual selected source orbital spin";
            }else if(nbo_active && nbo_ui.dataset &&
                     mo_index<nbo_ui.dataset->orbitals.size()) {
                view.kind=cov::ActiveOrbitalKind::NboSet;
                const auto& row=nbo_ui.dataset->orbitals[mo_index];
                view.source_id=integration?integration->id:nbo_ui.dataset->source.path;
                view.label=row.label;
                view.source_label=row.label;
                view.display_name_evidence="associated producer NBO orbital";
                view.semantic_kind="nbo";
                view.spin=row.spin;
                view.source_spin=row.spin;
                view.spin_semantics="producer NBO spin block";
                view.group_id="nbo:"+std::string(cov::nbo_spin_name(row.spin))+":"+
                    std::to_string(row.id);
                if(integration && row.id) {
                    const cov::NboOrbitalRef ref{cov::NboOrbitalKind::NBO,row.spin,row.id-1};
                    if(cov::nbo_orbital(*integration,ref))
                        view.selection=cov::nbo_single_selection(*integration,ref);
                }
                view.rendered_index=mo_index;
            }else if(wavefunction && mo_index<wavefunction->orbitals.size()) {
                view.kind=cov::ActiveOrbitalKind::Canonical;
                view.source_id=routed?routed->canonical_fingerprint:path_to_utf8(current_file);
                view.source_label=cov::ui::canonical_mo_source_label(*wavefunction,mo_index);
                const cov::ui::NboAomoName* display_name=nullptr;
                if(nbo_ui.aomo.names && mo_index<nbo_ui.aomo.names->canonical.size()){
                    const auto& name=nbo_ui.aomo.names->canonical[mo_index];
                    display_name=&name;
                    if(name.verified)view.display_name_evidence=name.detail;
                }
                view.label=cov::ui::canonical_mo_display_label(*wavefunction,mo_index,display_name);
                view.semantic_kind="canonical";
                view.spin=wavefunction->orbital_occupation_model==cov::OrbitalOccupationModel::CanonicalShared?
                    cov::NboSpin::Total:wavefunction->orbitals[mo_index].spin==cov::Spin::Beta?
                    cov::NboSpin::Beta:cov::NboSpin::Alpha;
                view.source_spin=wavefunction->orbitals[mo_index].spin==cov::Spin::Beta?
                    cov::NboSpin::Beta:cov::NboSpin::Alpha;
                view.spin_semantics=view.spin==cov::NboSpin::Total?
                    "shared spatial-orbital display; alpha is the producer channel field":
                    "explicit canonical spin orbital";
                view.canonical_index=mo_index;
                view.rendered_index=mo_index;
            }
            return view;
        };
        auto export_analysis_companions = [&](std::filesystem::path base) {
            if(!routed || !wavefunction)return;
            base.replace_extension();
            const auto write=[&](const char* suffix,const std::string& payload) {
                auto path=base;path+=suffix;
                std::ofstream out(path,std::ios::binary);
                if(!out)throw std::runtime_error("Cannot write routed analysis companion");
                out<<payload;
                if(!out)throw std::runtime_error("Routed analysis companion write failed");
            };
            const auto analysis=cov::serialize_routed_analysis_json(*routed);
            const auto active=cov::serialize_active_orbital_view_json(active_view());
            write(".analysis.json",analysis);
            write(".active-view.json",active);
            if(integration)
                write(".integration.json",cov::serialize_nbo_integration_json(*integration));
            else
                write(".integration.json",
                    "{\"status\":\"not_analysed\",\"reason\":\"No NBO integration attached\"}");
            if(nbo_ui.aomo.salc_model)
                write(".salc.json",cov::serialize_nbo_salc_json(*nbo_ui.aomo.salc_model));
            else write(".salc.json","{\"status\":\"not_analysed\",\"reason\":\"No verified NBO SALC model attached\"}");
            std::ostringstream names;
            names<<"{\"schema\":\"cov.orbital.display-names.v1\",\"canonical\":[";
            if(nbo_ui.aomo.names)
                for(std::size_t i=0;i<nbo_ui.aomo.names->canonical.size();++i){
                    if(i)names<<',';const auto& name=nbo_ui.aomo.names->canonical[i];
                    names<<"{\"index\":"<<i<<",\"label\":"<<cov::validation::quote(name.label)
                         <<",\"irrep\":"<<cov::validation::quote(name.irrep)
                         <<",\"ordinal\":"<<name.ordinal
                         <<",\"verified\":"<<(name.verified?"true":"false")
                         <<",\"detail\":"<<cov::validation::quote(name.detail)<<'}';
                }
            names<<"],\"salc\":[";
            if(nbo_ui.aomo.names)
                for(std::size_t i=0;i<nbo_ui.aomo.names->salc.size();++i){
                    if(i)names<<',';const auto& name=nbo_ui.aomo.names->salc[i];
                    names<<"{\"index\":"<<i<<",\"label\":"<<cov::validation::quote(name.label)
                         <<",\"irrep\":"<<cov::validation::quote(name.irrep)
                         <<",\"ordinal\":"<<name.ordinal
                         <<",\"verified\":"<<(name.verified?"true":"false")
                         <<",\"detail\":"<<cov::validation::quote(name.detail)<<'}';
                }
            names<<"]}";
            write(".display-names.json",names.str());
            cov::validation::record("chemistry.export",
                "{\"analysis\":"+analysis+",\"active_view\":"+active+"}");
        };

        auto identity = [&]() {
            const auto* active = active_wavefunction();
            const auto* orbital = active && mo_index < active->orbitals.size()
                                    ? &active->orbitals[mo_index] : nullptr;
            const auto spin = inspection && !inspection->selection.terms.empty() ? inspection->selection.terms.front().orbital.spin :
                              nbo_active && nbo_ui.dataset &&
                              mo_index<nbo_ui.dataset->orbitals.size()
                ? nbo_ui.dataset->orbitals[mo_index].spin
                : orbital && orbital->spin == cov::Spin::Beta
                    ? cov::NboSpin::Beta : cov::NboSpin::Alpha;
            const cov::NboCanonicalEvidence* evidence=nullptr;
            if(nbo_ui.dataset)for(const auto& record:nbo_ui.dataset->association.canonical_evidence)
                if(record.spin==spin || (spin==cov::NboSpin::Alpha && record.spin==cov::NboSpin::Total)){
                    evidence=&record;break;
                }
            const bool single_inspection=inspection && inspection->selection.terms.size()==1;
            const auto* descriptor=single_inspection && integration?
                cov::nbo_orbital(*integration,inspection->selection.terms[0].orbital):nullptr;
            const bool canonical_inspection=descriptor && descriptor->ref.kind==cov::NboOrbitalKind::Canonical;
            const bool direct_canonical=wavefunction && wavefunction->source==cov::WavefunctionSource::Fchk &&
                (!inspection || canonical_inspection) && !nbo_active;
            cov::validation::orbital_identity(
                inspection ? (inspection->selection.mode==cov::NboSelectionMode::Orbital && descriptor?
                    cov::nbo_orbital_kind_name(descriptor->ref.kind):"signed-combination") : nbo_active ? "nbo" : "canonical",
                inspection && integration ? integration->id : nbo_active && nbo_ui.dataset ? nbo_ui.dataset->source.path
                                             : path_to_utf8(current_file),
                cov::nbo_spin_name(spin),
                descriptor ? descriptor->ref.index : inspection ? std::numeric_limits<std::size_t>::max() : orbital ? orbital->source_orbital_index
                        : std::numeric_limits<std::size_t>::max(),
                nbo_ui.dataset ? nbo_ui.dataset->association.status : "not_attached",
                descriptor ? descriptor->source.path : inspection ? "verified signed orbital terms" :
                    direct_canonical ? path_to_utf8(current_file) : evidence ? evidence->coefficient_source : "not_available",
                direct_canonical,
                evidence && evidence->density_verified);
        };

        auto activate_set = [&](bool use_nbo) {
            if (use_nbo == nbo_active && !inspection) return;
            if (use_nbo && !nbo_wavefunction) throw std::runtime_error("NBO coefficients are not available for rendering");
            const cov::Wavefunction& target = use_nbo ? *nbo_wavefunction : *wavefunction;
            if (target.orbitals.empty()) throw std::runtime_error("Selected orbital set is empty");
            auto next_evaluator = std::make_unique<cov::CudaOrbitalEvaluator>(target);
            const std::size_t next_index = std::min(use_nbo ? nbo_mo_index : canonical_mo_index,
                                                    target.orbitals.size()-1);
            if (evaluator) evaluator->detach_gl_texture();
            const bool resized=resize_and_recompute || renderer.nx()!=resolution;
            try {
                if(resized)renderer.resize_volume(resolution,resolution,resolution);
                next_evaluator->attach_gl_texture(renderer.volume_texture());
                next_evaluator->evaluate(next_index, grid_box,
                                         resolution, resolution, resolution);
            } catch (...) {
                next_evaluator->detach_gl_texture();
                if (evaluator) {
                    evaluator->attach_gl_texture(renderer.volume_texture());
                    if(resized)evaluator->evaluate(mo_index,grid_box,
                                                   resolution,resolution,resolution);
                }
                throw;
            }
            evaluator = std::move(next_evaluator);
            additional_fields.clear();
            inspection.reset();
            nbo_ui.aomo.selection.reset();
            nbo_active = use_nbo;
            mo_index = next_index;
            pending_mo_index.reset();
            orbital_ui.browser_cache={};
            orbital_ui.diagram_cache={};
            renderer.invalidate_geometry_cache();
            resize_and_recompute=false;
            identity();
            cov::validation::evaluated(mo_index,"set-switch",evaluator->last_kernel_ms());
            status=StatusKind::GridUpdated;
            status_detail=(use_nbo?"NBO ":"Canonical MO ")+std::to_string(next_index+1);
            recompute = false;
        };

        auto evaluate_now = [&]() {
            const auto* active = active_wavefunction();
            if (!active || !evaluator || active->orbitals.empty()) return;
            if (resize_and_recompute || renderer.nx() != resolution) {
                evaluator->detach_gl_texture();
                renderer.resize_volume(resolution, resolution, resolution);
                evaluator->attach_gl_texture(renderer.volume_texture());
                resize_and_recompute = false;
            }
            evaluator->evaluate(mo_index, grid_box,
                                resolution, resolution, resolution);
            for(auto& field:additional_fields){
                if(field.renderer->nx()!=resolution){
                    field.evaluator->detach_gl_texture();
                    field.renderer->resize_volume(resolution,resolution,resolution);
                    field.evaluator->attach_gl_texture(field.renderer->volume_texture());
                }
                field.evaluator->evaluate(field.index,grid_box,resolution,resolution,resolution);
            }
            identity();
            cov::validation::evaluated(mo_index,"selection-or-grid",evaluator->last_kernel_ms());
            status = StatusKind::GridUpdated;
            status_detail = evaluator->device_name();
            recompute = false;
        };

        auto load_file = [&](const std::filesystem::path& path) {
            try {
                status = StatusKind::Parsing;
                status_detail = path_to_utf8(path);

                cov::MoldenParseOptions options;
                options.max_atoms = 100;
                options.require_orbitals = true;
                auto next_wavefunction=std::make_unique<cov::Wavefunction>(cov::parse_molden(path, options));
                auto& wf=*next_wavefunction;
                if (cov::validation::active()) {
                    std::ostringstream diagnostics;
                    cov::write_numerical_diagnostics_json(diagnostics,wf);
                    cov::validation::record("input.numerical_diagnostics",diagnostics.str());
                    std::ostringstream density_evidence;
                    cov::write_density_evidence_json(density_evidence,wf);
                    cov::validation::record("input.density_evidence",density_evidence.str());
                    std::ostringstream topology_evidence;
                    cov::write_pi_topology_assignments_json(topology_evidence,wf);
                    cov::validation::record("input.pi_topology_evidence",topology_evidence.str());
                }
                const auto new_mo = initial_orbital(wf);
                const auto new_box = make_grid_box(wf);
                std::optional<cov::OrbitalTrackingResult> new_tracking;
                if (wavefunction) {
                    // Cross-frame identity is descriptive state only. Both
                    // canonical wavefunctions remain immutable, and loading a
                    // new frame still resets selection to that frame's own HOMO.
                    // ✳ TODO: Profile slow frame opens and bound search work per conflict group.
                    new_tracking = cov::track_orbital_subspaces(*wavefunction, wf);
                }

                auto next_route=cov::route_chemistry(wf);
                auto next_graph=*next_route.interaction_graph.value;
                auto next_evaluator=std::make_unique<cov::CudaOrbitalEvaluator>(wf);
                if (evaluator) evaluator->detach_gl_texture();
                try {
                    renderer.resize_volume(resolution,resolution,resolution);
                    next_evaluator->attach_gl_texture(renderer.volume_texture());
                    next_evaluator->evaluate(new_mo,new_box,resolution,resolution,resolution);
                } catch(...) {
                    next_evaluator->detach_gl_texture();
                    if(evaluator){evaluator->attach_gl_texture(renderer.volume_texture());
                        evaluator->evaluate(mo_index,grid_box,resolution,resolution,resolution);}
                    throw;
                }
                evaluator=std::move(next_evaluator);
                additional_fields.clear();
                inspection.reset();
                nbo_ui.integration=nullptr;
                integration.reset();
                routed.reset();
                nbo_ui.aomo={};
                nbo_ui.selected_atoms.clear();
                nbo_ui.selected_structure.reset();
                nbo_ui.atom_colour_mode=0;
                nbo_ui.show_bond_indices=false;nbo_ui.show_e2=false;
                nbo_ui.input_discovery.reset();nbo_ui.pending_candidate.reset();
                nbo_active = false;
                nbo_wavefunction.reset();
                nbo_ui.dataset.reset();
                nbo_ui.focus = {};
                nbo_ui.error.clear();
                nbo_ui.export_status.clear();
                wavefunction = std::move(next_wavefunction);
                cov::ui::invalidate_canonical_mo_names_cache();
                routed=std::move(next_route);
                nbo_ui.routed=&*routed;
                semantic_graph=std::move(next_graph);
                frame_tracking = std::move(new_tracking);
                renderer.invalidate_geometry_cache();
                mo_index = new_mo;
                canonical_mo_index = new_mo;
                nbo_mo_index = 0;
                pending_mo_index.reset();
                orbital_ui.browser_cache={};
                orbital_ui.diagram_cache={};
                if(cov::validation::active())
                    cov::validation::record("chemistry.route",cov::serialize_routed_analysis_json(*routed));
                grid_box = new_box;
                current_file = path;

                identity();
                cov::validation::evaluated(mo_index,"input-load",evaluator->last_kernel_ms());
                copy_path_to_buffer(path, path_buffer);
                push_recent(recent_files, path);
                status = StatusKind::Loaded;
                status_detail = path_to_utf8(path.filename());
            } catch (const std::exception& e) {
                status = StatusKind::Error;
                status_detail = e.what();
            }
        };

        auto apply_selection = [&](const cov::NboOrbitalSelection& selection) {
            if(!integration || !wavefunction) throw std::runtime_error("No verified orbital data is attached");
            auto next=std::make_unique<cov::NboSelectionView>(cov::make_nbo_selection_view(*integration,*wavefunction,selection));
            if(!next->available || next->wavefunction.orbitals.empty())
                throw std::runtime_error(next->detail.empty()?next->status:next->detail);
            auto next_evaluator=std::make_unique<cov::CudaOrbitalEvaluator>(next->wavefunction);
            std::vector<AdditionalField> extra;
            for(std::size_t i=1;i<next->wavefunction.orbitals.size();++i){
                AdditionalField f;
                f.renderer=std::make_unique<cov::VolumeRenderer>();
                f.renderer->resize_volume(resolution,resolution,resolution);
                f.evaluator=std::make_unique<cov::CudaOrbitalEvaluator>(next->wavefunction);
                f.evaluator->attach_gl_texture(f.renderer->volume_texture());
                f.evaluator->evaluate(i,grid_box,resolution,resolution,resolution);
                f.index=i;extra.push_back(std::move(f));
            }
            if(evaluator)evaluator->detach_gl_texture();
            try {
                if(renderer.nx()!=resolution)renderer.resize_volume(resolution,resolution,resolution);
                next_evaluator->attach_gl_texture(renderer.volume_texture());
                next_evaluator->evaluate(0,grid_box,resolution,resolution,resolution);
            } catch(...) {
                next_evaluator->detach_gl_texture();
                if(evaluator){evaluator->attach_gl_texture(renderer.volume_texture());
                    evaluator->evaluate(mo_index,grid_box,resolution,resolution,resolution);}
                throw;
            }
            evaluator=std::move(next_evaluator);
            additional_fields=std::move(extra);
            inspection=std::move(next);
            nbo_active=false;mo_index=0;
            nbo_ui.aomo.selection=selection;
            status_detail=inspection->label;
            // Match the scene's canonical identity rather than an ambiguous
            // spin-local descriptor number (for example beta MO 3 vs MO 31).
            if(selection.terms.size()==1 &&
               selection.terms.front().orbital.kind==cov::NboOrbitalKind::Canonical) {
                const auto& ref=selection.terms.front().orbital;
                status_detail="MO "+std::to_string(ref.index+1)+" ["+
                    cov::nbo_spin_name(ref.spin)+"]";
            }
            nbo_ui.aomo.status=status_detail;
            nbo_ui.error.clear();status=StatusKind::GridUpdated;
            nbo_ui.selected_atoms={inspection->atoms.begin(),inspection->atoms.end()};
            if(selection.target_canonical_index)canonical_mo_index=*selection.target_canonical_index;
            else if(selection.terms.size()==1 && selection.terms.front().orbital.kind==cov::NboOrbitalKind::Canonical)
                canonical_mo_index=selection.terms.front().orbital.index;
            pending_mo_index.reset();recompute=false;resize_and_recompute=false;
            identity();
            if(cov::validation::active())
                cov::validation::record("aomo.selection",cov::serialize_nbo_selection_json(*inspection));
            cov::validation::evaluated(0,"typed-orbital-selection",evaluator->last_kernel_ms());
        };

        auto attach_integration = [&](cov::NboIntegration next,bool keep_discovery=false) {
            if(!wavefunction)throw std::runtime_error("Load the matching Gaussian wavefunction first");
            cov::annotate_nbo_bond_channels(next,*wavefunction);
            if(inspection || nbo_active)activate_set(false);
            nbo_wavefunction.reset();
            integration=std::move(next);
            nbo_ui.integration=&*integration;
            nbo_ui.dataset=integration->dataset;
            routed=cov::route_chemistry(*wavefunction,&*integration);
            nbo_ui.routed=&*routed;
            semantic_graph=*routed->interaction_graph.value;
            renderer.invalidate_geometry_cache();
            if(cov::validation::active())
                cov::validation::record("chemistry.route",cov::serialize_routed_analysis_json(*routed));
            nbo_ui.focus={};nbo_ui.aomo={};nbo_ui.selected_atoms.clear();nbo_ui.selected_structure.reset();
            nbo_ui.show_bond_indices=false;nbo_ui.show_e2=false;
            nbo_ui.pending_candidate.reset();
            if(!keep_discovery)nbo_ui.input_discovery.reset();
            nbo_ui.input_status=integration->dataset.source.path;
            nbo_ui.error.clear();
            if(const auto* c=cov::nbo_capability(*integration,"nbo");c && c->available()) {
                try { nbo_wavefunction=cov::make_nbo_wavefunction(integration->dataset,*wavefunction); }
                catch(const std::exception& e){nbo_ui.error=e.what();}
            }
            const auto suggested=cov::ui::suggest_initial_nbo_index(*integration,orbital_ui.diagram_cache.snapshot.get());
            nbo_ui.selected_orbital=suggested.value_or(std::numeric_limits<std::size_t>::max());
            nbo_mo_index=nbo_ui.selected_orbital;
            status=nbo_ui.error.empty()?StatusKind::Loaded:StatusKind::Error;
            status_detail=nbo_ui.error.empty()?integration->dataset.source.path:nbo_ui.error;
            // Do not build large diagnostic payloads during ordinary viewing.
            if(cov::validation::active())
                cov::validation::record("nbo.integration",cov::serialize_nbo_integration_json(*integration));
            cov::validation::record("nbo.attach","{\"source\":"+cov::validation::quote(integration->dataset.source.path)+
                ",\"association\":"+cov::validation::quote(integration->dataset.association.status)+
                ",\"renderable\":"+(nbo_wavefunction?"true":"false")+"}");
        };
        auto clear_integration = [&] {
            if(inspection || nbo_active)activate_set(false);
            nbo_wavefunction.reset();nbo_ui.integration=nullptr;integration.reset();
            nbo_ui.dataset.reset();nbo_ui.aomo={};nbo_ui.focus={};
            nbo_ui.selected_atoms.clear();nbo_ui.selected_structure.reset();
            nbo_ui.selected_orbital=std::numeric_limits<std::size_t>::max();
            nbo_ui.atom_colour_mode=0;nbo_ui.show_bond_indices=false;nbo_ui.show_e2=false;
            nbo_ui.pending_candidate.reset();
            orbital_ui.browser_cache={};orbital_ui.diagram_cache={};
            routed.reset();nbo_ui.routed=nullptr;semantic_graph={};
            if(wavefunction){
                routed=cov::route_chemistry(*wavefunction);
                nbo_ui.routed=&*routed;semantic_graph=*routed->interaction_graph.value;
                if(cov::validation::active())
                    cov::validation::record("chemistry.route",cov::serialize_routed_analysis_json(*routed));
            }
            renderer.invalidate_geometry_cache();
        };
        auto apply_candidate = [&](const cov::NboInputCandidate& candidate) {
            if(!candidate.canonical.empty() && candidate.canonical.lexically_normal()!=current_file.lexically_normal()) {
                load_file(candidate.canonical);
                if(status==StatusKind::Error)throw std::runtime_error(status_detail);
            }
            if(!wavefunction)throw std::runtime_error("The package has no loaded canonical wavefunction");
            if(candidate.report.empty()){
                clear_integration();
                nbo_ui.input_status="NBO data missing — Gaussian molecular orbitals remain available";
                status=StatusKind::Loaded;
                status_detail=nbo_ui.input_status;
                return;
            }
            attach_integration(cov::read_nbo_integration(*wavefunction,candidate),true);
            nbo_ui.input_status=candidate.label;
            if(const auto* matched=cov::nbo_capability(*integration,"source_association");
               matched && matched->state==cov::NboCapabilityState::Rejected){
                nbo_ui.error="Gaussian / NBO mismatch: "+matched->detail;
                nbo_ui.input_status=nbo_ui.error+". Gaussian molecular orbitals remain available.";
                status=StatusKind::Error;status_detail=nbo_ui.input_status;
            }
        };
        auto load_inputs = [&](const std::vector<std::filesystem::path>& paths) {
            try {
                nbo_ui.error.clear();
                auto found=cov::discover_nbo_inputs(paths);
                nbo_ui.pending_candidate.reset();
                if(found.candidates.size()==1 && !found.selection_required) {
                    apply_candidate(found.candidates.front());
                } else if(found.candidates.empty()) {
                    // Existing Molden and standalone Gaussian use the same loader.
                    if(paths.size()==1 && !std::filesystem::is_directory(paths.front())){
                        load_file(paths.front());
                        if(status==StatusKind::Error)throw std::runtime_error(status_detail);
                    }
                    nbo_ui.input_status="No unique NBO calculation found. Gaussian view is retained.";
                } else {
                    nbo_ui.input_status="Multiple calculations found — choose the matching calculation below";
                }
                nbo_ui.input_discovery=std::move(found);
            } catch(const std::exception& e){
                clear_integration();
                nbo_ui.input_discovery.reset();
                nbo_ui.error=e.what();nbo_ui.input_status="NBO association failed; Gaussian view is retained";
                status=StatusKind::Error;status_detail=nbo_ui.error;
                cov::validation::record("input.package.error","{\"reason\":"+cov::validation::quote(e.what())+"}");
            }
        };

        if (argc >= 2) {
            const std::string p = argv[1];
            std::snprintf(path_buffer.data(), path_buffer.size(), "%s", p.c_str());
            load_inputs({path_from_utf8(p)});
        }
        if (cov::validation::active() && !wavefunction) {
            throw std::runtime_error("Native validation input failed: "+status_detail);
        }

        bool scene_drag_active = false;
        ImVec2 scene_press{};
        bool scene_was_dragged=false;
        bool choose_nbo_input=false;

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            // Open at the next frame boundary, after previous draw references
            // are released. Reuse the package loader and native picker.
            if(choose_nbo_input){
                choose_nbo_input=false;
                const auto dialog=cov::open_molden_file_dialog();
                if(dialog.selected())load_inputs({dialog.path});
                else if(!dialog.cancelled&&!dialog.error.empty()){
                    status=StatusKind::Error;
                    status_detail=dialog.supported?dialog.error:
                        cov::ui::tr(cov::ui::Text::OpenDialogUnsupported,language);
                }
            }
            cov::validation::begin_frame(camera,molecule_render,isovalue,resolution,resize_and_recompute);
            if (resize_and_recompute) recompute = true;
            if(auto paths=cov::validation::take_dropped_paths();!paths.empty())load_inputs(paths);

            if (!g_dropped_paths.empty()) {
                const auto paths=std::move(g_dropped_paths);
                g_dropped_paths.clear();
                load_inputs(paths);
            }

            if (cov::validation::active()) {
                int window_width=0, window_height=0;
                glfwGetWindowSize(window,&window_width,&window_height);
                if (window_width!=cov::validation::window_width() ||
                    window_height!=cov::validation::window_height()) {
                    glfwSetWindowSize(window,cov::validation::window_width(),cov::validation::window_height());
                    glfwPollEvents();
                }
            }
            int fb_w = 0, fb_h = 0;
            glfwGetFramebufferSize(window, &fb_w, &fb_h);
            if (cov::validation::active() && (fb_w != cov::validation::window_width() ||
                                             fb_h != cov::validation::window_height())) {
                throw std::runtime_error("Validation framebuffer differs from the requested size");
            }
            if (fb_w <= 0 || fb_h <= 0) {
                glfwWaitEventsTimeout(0.05);
                continue;
            }
            ImGui_ImplOpenGL2_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            cov::validation::input_frame();
            ImGui::NewFrame();
            ImGuiIO& io = ImGui::GetIO();
            bool canonical_requested=false;
            const auto* aomo_cap=integration?cov::nbo_capability(*integration,"aomo"):nullptr;
            const bool aomo_available=aomo_cap && aomo_cap->available();
            const auto layout = cov::viewer_layout(io.DisplaySize.x, io.DisplaySize.y,
                                                    fb_w, fb_h, ui_scale,aomo_available);
            const auto& viewport = layout.framebuffer;
            const bool over_scene = layout.scene.contains(io.MousePos.x, io.MousePos.y);
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                scene_drag_active = over_scene && !io.WantCaptureMouse;
                scene_press=io.MousePos;scene_was_dragged=false;
            }
            if(scene_drag_active && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                if(!scene_was_dragged && over_scene && integration){
                    const auto picked=renderer.pick_geometry((io.MousePos.x-layout.scene.x)/layout.scene.width,
                        (io.MousePos.y-layout.scene.y)/layout.scene.height);
                    if(picked){
                        if(picked->kind==cov::GeometryTargetKind::Atom){
                            nbo_ui.selected_structure.reset();
                            if(!io.KeyCtrl)nbo_ui.selected_atoms.clear();
                            if(nbo_ui.selected_atoms.count(picked->index))nbo_ui.selected_atoms.erase(picked->index);
                            else nbo_ui.selected_atoms.insert(picked->index);
                        }else if(picked->index<integration->structure.size()){
                            nbo_ui.selected_structure=picked->index;
                            const auto& e=integration->structure[picked->index];
                            nbo_ui.selected_atoms={e.atoms.begin(),e.atoms.end()};
                            if(e.kind=="donor_acceptor" && e.orbitals.size()>=2){
                                cov::NboOrbitalSelection s;s.dataset_id=integration->id;s.label=e.label;
                                s.mode=cov::NboSelectionMode::Overlay;
                                s.terms={{e.orbitals[0],1},{e.orbitals[1],1}};nbo_ui.aomo.pending_selection=s;
                            }
                        }
                        cov::validation::record("scene.pick","{\"kind\":"+
                            std::to_string(static_cast<int>(picked->kind))+",\"index\":"+std::to_string(picked->index)+"}");
                    }
                }
                scene_drag_active=false;
            }
            if(scene_drag_active && std::hypot(io.MousePos.x-scene_press.x,io.MousePos.y-scene_press.y)>4*ui_scale)
                scene_was_dragged=true;
            if (scene_drag_active && scene_was_dragged && ImGui::IsMouseDown(ImGuiMouseButton_Left) && over_scene && !io.WantCaptureMouse) {
                camera.yaw += io.MouseDelta.x * 0.007f;
                camera.pitch = std::clamp(camera.pitch + io.MouseDelta.y * 0.007f, -1.45f, 1.45f);
            }
            if (over_scene && !io.WantCaptureMouse && std::abs(io.MouseWheel) > 0.0f) {
                camera.distance *= std::pow(0.88f, io.MouseWheel);
                camera.distance = std::clamp(camera.distance, 1.1f, 6.0f);
            }
            glViewport(0, 0, fb_w, fb_h);
            glClearColor(0.025f, 0.031f, 0.043f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glViewport(viewport.x, viewport.y, viewport.width, viewport.height);
            identity();
            cov::validation::field("scene.active_view",
                cov::serialize_active_orbital_view_json(active_view()));
            if(routed) {
                cov::validation::field("chemistry.route.charge",
                    std::string(cov::routed_status_name(routed->total_atomic_charge.status))+":"+
                    cov::routed_provider_name(routed->total_atomic_charge.provider)+":"+
                    routed->total_atomic_charge.fallback_reason);
                cov::validation::field("chemistry.route.spin",
                    std::string(cov::routed_status_name(routed->atomic_spin.status))+":"+
                    cov::routed_provider_name(routed->atomic_spin.provider)+":"+
                    routed->atomic_spin.reason);
            }
            std::optional<cov::MoleculeOverlay> overlay;
            if(integration && wavefunction)overlay=cov::make_nbo_molecule_overlay(*integration,semantic_graph,
                wavefunction->atoms.size(),{nbo_ui.selected_atoms.begin(),nbo_ui.selected_atoms.end()},
                nbo_ui.selected_structure,static_cast<cov::AtomScalarMode>(nbo_ui.atom_colour_mode),
                nbo_ui.show_bond_indices,nbo_ui.show_e2,routed?&*routed:nullptr);
            if(overlay)cov::validation::field("overlay.scalars",
                cov::serialize_molecule_overlay_scalars_json(*overlay,routed?&*routed:nullptr));
            if (const auto* active = active_wavefunction();
                active && viewport.width > 0 && viewport.height > 0) {
                // Geometry and the actual orbital texture share one viewport,
                // projection and depth buffer, outside the control panel.
                renderer.render_geometry(*wavefunction, grid_box, viewport.width, viewport.height,
                                         camera, molecule_render,overlay?&*overlay:nullptr,
                                         &semantic_graph);
                renderer.render_volume(viewport.width, viewport.height, isovalue, camera,
                                       molecule_render.orbital_opacity,
                                       orbital_material, orbital_surface_mode);
                for(auto& field:additional_fields)field.renderer->render_volume(viewport.width,viewport.height,
                    isovalue,camera,molecule_render.orbital_opacity,orbital_material,orbital_surface_mode,1);
                cov::validation::after_scene(renderer, grid_box, mo_index);
                for(std::size_t f=0;f<additional_fields.size();++f)
                    cov::validation::after_scene(*additional_fields[f].renderer,grid_box,mo_index,f+1);
            }
            cov::validation::scene_view(layout, camera);
            glViewport(0, 0, fb_w, fb_h);

            // A non-intercepting scene legend also supplies validation hit
            // locations. Both ordinary clicks and test clicks use pick_geometry.
            ImGui::SetNextWindowPos(ImVec2(layout.scene.x,layout.scene.y));
            ImGui::SetNextWindowSize(ImVec2(layout.scene.width,layout.scene.height));
            ImGui::SetNextWindowBgAlpha(0);
            ImGui::Begin("##scene_evidence",nullptr,ImGuiWindowFlags_NoInputs|ImGuiWindowFlags_NoDecoration|
                ImGuiWindowFlags_NoSavedSettings|ImGuiWindowFlags_NoBringToFrontOnFocus);
            if(wavefunction){
                const auto spin_caption=[&](cov::NboSpin spin)->const char* {
                    return spin==cov::NboSpin::Alpha?"α":spin==cov::NboSpin::Beta?"β":
                        scene_text(language,"total","总","全体","total");
                };
                const auto mo_caption=[&](std::size_t index) {
                    const auto* name=nbo_ui.aomo.names && index<nbo_ui.aomo.names->canonical.size()
                        ?&nbo_ui.aomo.names->canonical[index]:nullptr;
                    return cov::ui::canonical_mo_display_label(*wavefunction,index,name);
                };
                std::string label;
                if(inspection){
                    label=inspection->label;
                    if(inspection->selection.terms.size()==1){const auto& ref=inspection->selection.terms[0].orbital;
                        label+=" | "+std::string(cov::nbo_orbital_kind_name(ref.kind))+" "+std::to_string(ref.index+1)+" / "+spin_caption(ref.spin);
                    }else {
                        label+=" | "+std::to_string(inspection->selection.terms.size())+
                            scene_text(language," signed terms"," 个带系数项"," 個の符号付き成分"," termes signés");
                        std::vector<std::string> source_kinds;
                        for(const auto& term:inspection->selection.terms){
                            const auto identity=std::string(cov::nbo_orbital_kind_name(term.orbital.kind))+
                                " / "+spin_caption(term.orbital.spin);
                            if(std::find(source_kinds.begin(),source_kinds.end(),identity)==source_kinds.end())
                                source_kinds.push_back(identity);
                        }
                        for(const auto& identity:source_kinds)label+=" | "+identity;
                    }
                }
                else if(nbo_active && nbo_ui.dataset && mo_index<nbo_ui.dataset->orbitals.size())
                    label="NBO "+std::to_string(mo_index+1)+" / "+spin_caption(nbo_ui.dataset->orbitals[mo_index].spin);
                else label=mo_caption(canonical_mo_index);
                // Translate presentation only; inspection identities and the actual signed field stay intact.
                if(inspection) {
                    const auto& selected=inspection->selection;
                    if(selected.mode==cov::NboSelectionMode::Orbital && selected.terms.size()==1 &&
                       selected.terms[0].orbital.kind==cov::NboOrbitalKind::Canonical)
                        label=mo_caption(selected.terms[0].orbital.index);
                    else if(selected.target_canonical_index)
                        label=std::string(scene_text(language,"Orbital components of ","轨道组成：","軌道成分：","Composantes orbitales de "))+
                            mo_caption(*selected.target_canonical_index);
                    else if(selected.mode==cov::NboSelectionMode::Combination && nbo_ui.aomo.salc_model &&
                            selected.dataset_id==nbo_ui.aomo.salc_model->dataset_id) {
                        const auto& model=*nbo_ui.aomo.salc_model;
                        for(std::size_t i=0;i<model.orbitals.size();++i) {
                            const auto& side=model.orbitals[i];
                            if(side.terms.size()!=selected.terms.size())continue;
                            bool same=true;
                            for(std::size_t j=0;j<side.terms.size();++j)
                                same=same && side.terms[j].orbital==selected.terms[j].orbital &&
                                    side.terms[j].coefficient==selected.terms[j].coefficient;
                            if(!same)continue;
                            label=side.atoms.size()==1?"NAO":
                                ((nbo_ui.aomo.names && i<nbo_ui.aomo.names->salc.size())?
                                    nbo_ui.aomo.names->salc[i].label+" · SALC":"SALC");
                            label+=" (";
                            for(std::size_t j=0;j<side.atoms.size();++j) {
                                if(j)label+=", ";const auto atom=side.atoms[j];
                                if(atom<wavefunction->atoms.size())label+=wavefunction->atoms[atom].symbol+std::to_string(atom+1);
                            }
                            label+=") / ";label+=spin_caption(side.spin);break;
                        }
                    }
                }
                ImGui::TextWrapped("%s",label.c_str());
                if(inspection && inspection->selection.mode==cov::NboSelectionMode::Overlay){
                    for(std::size_t i=0;i<inspection->selection.terms.size();++i){
                        const auto& term=inspection->selection.terms[i];
                        ImGui::TextWrapped("%zu: %s %zu / %s; c=%+.6g",i+1,
                            cov::nbo_orbital_kind_name(term.orbital.kind),term.orbital.index+1,
                            cov::nbo_spin_name(term.orbital.spin),term.coefficient);
                    }
                }
                ImGui::TextColored(ImVec4(.96f,.5f,.45f,1),"%s",scene_text(language,
                    "Phase: red + / blue -","相位：红 + / 蓝 −","位相：赤 + / 青 −","Phase : rouge + / bleu -"));
                if(!additional_fields.empty())ImGui::TextColored(ImVec4(.3f,.92f,.75f,1),"%s",
                    scene_text(language,"Overlay: green + / gold - (separate field)",
                        "叠加轨道：绿 + / 金 −（独立场）","重ね合わせ：緑 + / 金 −（独立場）",
                        "Superposition : vert + / or - (champ distinct)"));
                if(overlay && std::any_of(overlay->bonds.begin(),overlay->bonds.end(),[](const auto& bond){
                    return bond.style==cov::OverlayBondStyle::Unresolved;}))
                    ImGui::TextWrapped("%s",scene_text(language,
                        "Grey dotted links: connectivity is evidenced; integer bond order is unresolved. Click for orbital details and indices.",
                        "灰色点划线：连接证据存在，整数键型未确定；点选查看轨道详情和键级指数。",
                        "灰色の点線：結合の証拠はありますが、整数結合次数は未確定です。クリックして軌道と指数を確認できます。",
                        "Traits pointillés gris : connexion attestée, ordre entier non résolu. Cliquez pour les orbitales et les indices."));
                if(overlay && overlay->colour_mode!=cov::AtomScalarMode::Element){
                    ImGui::Text("%s: -%.3f ... 0 ... +%.3f",
                        overlay->colour_mode==cov::AtomScalarMode::NaturalCharge?
                            scene_text(language,"NPA charge / e","NPA 电荷 / e","NPA 電荷 / e","Charge NPA / e"):
                            scene_text(language,"Spin population / e","自旋布居 / e","スピン分布 / e","Population de spin / e"),overlay->scalar_range,overlay->scalar_range);
                    ImGui::TextDisabled("%s",scene_text(language,
                        "Atoms: purple (-), grey (0), gold (+); missing stays uncoloured",
                        "原子色标：紫（负）→灰（零）→金（正）；缺失不着色",
                        "原子の色：紫（負）→灰（零）→金（正）；欠損値は着色しません",
                        "Atomes : violet (-), gris (0), or (+) ; sans couleur si absent"));
                }
                for(const auto& t:renderer.geometry_targets()){
                    const float px=layout.scene.x+(t.segment?(t.x+t.end_x)/2:t.x)*layout.scene.width;
                    const float py=layout.scene.y+(t.segment?(t.y+t.end_y)/2:t.y)*layout.scene.height;
                    const std::string kind=t.kind==cov::GeometryTargetKind::Atom?"atom":
                        t.kind==cov::GeometryTargetKind::Bond?"bond":t.kind==cov::GeometryTargetKind::Relation?"relation":"multicentre";
                    cov::validation::hit("scene."+kind+"."+std::to_string(t.index),ImVec2(px-3,py-3),ImVec2(px+3,py+3));
                    if(overlay && nbo_ui.show_bond_indices && t.kind==cov::GeometryTargetKind::Bond &&
                       integration && t.index<integration->structure.size()){
                        const auto& e=integration->structure[t.index];
                        if(e.wiberg){char value[64];std::snprintf(value,sizeof(value),"WBI %.3f",*e.wiberg);
                            ImGui::GetWindowDrawList()->AddText(ImVec2(px+8,py+8),IM_COL32(255,225,145,255),value);}
                    }
                }
                if(integration && nbo_ui.selected_structure && *nbo_ui.selected_structure<integration->structure.size()){
                    const auto& evidence=integration->structure[*nbo_ui.selected_structure];
                    ImGui::TextWrapped("%s",evidence.label.c_str());
                    if(evidence.value)ImGui::Text("%.5g %s",*evidence.value,evidence.units.c_str());
                }
                if(overlay)for(auto atom:nbo_ui.selected_atoms)if(atom<overlay->atom_values.size() && overlay->atom_values[atom])
                    ImGui::Text("Atom %zu: %+.6f",atom+1,*overlay->atom_values[atom]);
                cov::validation::field("scene.orbital_label",label);
                cov::validation::field("scene.overlay_fields",std::to_string(additional_fields.size()+1));
            }
            ImGui::End();

            if(wavefunction && layout.scene.width>220 && layout.scene.height>160){
                const float toolbar_width=std::min(layout.scene.width-24.0f,520.0f*ui_scale);
                const char* opacity_help=scene_text(language,"Orbital opacity (lower to see atoms and bonds)",
                    "轨道不透明度（降低可看清原子与键）","軌道の不透明度（下げると原子と結合が見えます）",
                    "Opacité orbitale (réduire pour voir atomes et liaisons)");
                const char* threshold_help=scene_text(language,
                    "Display threshold only; orbital coefficients and amplitudes stay unchanged.",
                    "仅调整显示阈值；轨道系数与幅度不变。",
                    "表示しきい値のみ調整します。軌道係数と振幅は変わりません。",
                    "Seuil d’affichage seul ; coefficients et amplitudes restent inchangés.");
                const auto& toolbar_style=ImGui::GetStyle();
                const float text_width=std::max(1.0f,toolbar_width-2*toolbar_style.WindowPadding.x);
                const char* reveal_label=scene_text(language,"Reveal bonds","看清骨架","骨格を表示","Voir les liaisons");
                const char* iso_label=scene_text(language,"Isovalue","等值面","等値面","Isovaleur");
                const char* fit_label=scene_text(language,"Fit component","适合当前成分","成分に合わせる","Adapter");
                const float fit_width=ImGui::CalcTextSize(fit_label).x+2*toolbar_style.FramePadding.x;
                const float reveal_width=ImGui::CalcTextSize(reveal_label).x+2*toolbar_style.FramePadding.x;
                const float value_width=ImGui::CalcTextSize("1.2345e-07").x+2*toolbar_style.FramePadding.x;
                const bool stacked=text_width<ImGui::CalcTextSize(iso_label).x+value_width+
                    fit_width+2*toolbar_style.ItemSpacing.x;
                const float toolbar_height=ImGui::CalcTextSize(opacity_help,nullptr,false,text_width).y+
                    ImGui::CalcTextSize(threshold_help,nullptr,false,text_width).y+
                    (stacked?4:2)*ImGui::GetFrameHeight()+(stacked?6:4)*toolbar_style.ItemSpacing.y+
                    2*toolbar_style.WindowPadding.y;
                ImGui::SetNextWindowPos(ImVec2(layout.scene.x+12.0f,
                    std::max(layout.scene.y+12.0f,layout.scene.y+layout.scene.height-toolbar_height-12.0f)));
                ImGui::SetNextWindowSize(ImVec2(toolbar_width,toolbar_height));
                ImGui::SetNextWindowBgAlpha(.82f);
                ImGui::Begin("##scene_display_controls",nullptr,ImGuiWindowFlags_NoDecoration|
                    ImGuiWindowFlags_NoMove|ImGuiWindowFlags_NoSavedSettings);
                ImGui::TextWrapped("%s",opacity_help);
                ImGui::SetNextItemWidth(stacked?text_width:std::max(70.0f,text_width-reveal_width-toolbar_style.ItemSpacing.x));
                ImGui::SliderFloat("##scene_orbital_opacity",&molecule_render.orbital_opacity,.02f,1.0f,"%.2f");
                cov::validation::item("scene.opacity");
                if(!stacked)ImGui::SameLine();
                if(ImGui::Button(reveal_label))molecule_render.orbital_opacity=.24f;
                cov::validation::item("scene.reveal_bonds");
                ImGui::TextUnformatted(iso_label);ImGui::SameLine();
                ImGui::SetNextItemWidth(std::max(65.0f,ImGui::GetContentRegionAvail().x-
                    (stacked?0:fit_width+toolbar_style.ItemSpacing.x)));
                ImGui::SliderFloat("##scene_isovalue",&isovalue,1e-7f,.2f,"%.5g",ImGuiSliderFlags_Logarithmic);
                cov::validation::item("scene.isovalue");if(!stacked)ImGui::SameLine();
                if(ImGui::Button(fit_label)){
                    double norm2=1;
                    if(inspection && !inspection->selection.normalize && !inspection->metric_norm2.empty())
                        norm2=*std::max_element(inspection->metric_norm2.begin(),inspection->metric_norm2.end());
                    isovalue=std::clamp(static_cast<float>(.03*std::sqrt(std::max(0.0,norm2))),1e-7f,.2f);
                }
                cov::validation::item("scene.fit_component");
                ImGui::PushStyleColor(ImGuiCol_Text,ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
                ImGui::TextWrapped("%s",threshold_help);
                ImGui::PopStyleColor();
                ImGui::End();
            }

            ImGui::SetNextWindowPos(ImVec2(layout.controls.x, layout.controls.y), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(layout.controls.width, layout.controls.height), ImGuiCond_Always);
            ImGui::SetNextWindowBgAlpha(0.965f);
            constexpr ImGuiWindowFlags panel_flags =
                ImGuiWindowFlags_NoTitleBar |
                ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoSavedSettings;

            ImGui::Begin("##cov_control_panel", nullptr, panel_flags);
            const auto panel_position = ImGui::GetWindowPos();
            const auto panel_size = ImGui::GetWindowSize();
            cov::validation::hit("layout.control-panel", panel_position,
                ImVec2(panel_position.x + panel_size.x, panel_position.y + panel_size.y));

            // Keep the complete panel reachable when the window is short.
            ImGui::BeginChild("##cov_panel_scroll", ImVec2(0, 0), false,
                              ImGuiWindowFlags_None);
            const auto brand = [&] {
                ImGui::TextWrapped("%s", cov::ui::tr(cov::ui::Text::AppTitle, language));
                disabled_wrapped(cov::ui::tr(cov::ui::Text::Tagline, language));
            };
            const auto language_control = [&] {
                disabled_wrapped(cov::ui::tr(cov::ui::Text::LanguageLabel, language));
                int language_index = static_cast<int>(language);
                ImGui::SetNextItemWidth(-1.0f);
                if (ImGui::Combo("##language_combo", &language_index,
                                 "English\0简体中文\0日本語\0Français\0")) {
                    language = static_cast<cov::ui::Language>(language_index);
                    glfwSetWindowTitle(window,
                        cov::ui::tr(cov::ui::Text::AppTitle, language));
                }
                cov::validation::item("language");
            };
            const float header_width = ImGui::CalcTextSize(
                cov::ui::tr(cov::ui::Text::AppTitle, language)).x +
                142.0f * ui_scale + 4.0f * ImGui::GetStyle().ItemSpacing.x;
            if (ImGui::GetContentRegionAvail().x < header_width) {
                brand();
                ImGui::Spacing();
                language_control();
            } else if (ImGui::BeginTable("##cov_header", 2,
                                  ImGuiTableFlags_SizingStretchProp |
                                  ImGuiTableFlags_NoSavedSettings)) {
                ImGui::TableSetupColumn("##brand", ImGuiTableColumnFlags_WidthStretch);
                ImGui::TableSetupColumn("##language", ImGuiTableColumnFlags_WidthFixed,
                                        142.0f * ui_scale);
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                brand();
                ImGui::TableNextColumn();
                language_control();
                ImGui::EndTable();
            }

            ImGui::Spacing();
            cov::ui::status_badge(status_label(status, language), status_tone(status));
            if (!status_detail.empty()) {
                disabled_wrapped(status_detail.c_str());
                cov::validation::field("status.detail",status_detail);
            } else {
                disabled_wrapped(cov::ui::tr(cov::ui::Text::IdleHint, language));
            }
            ImGui::Separator();
            ImGui::Spacing();

            const bool show_input_panels=!aomo_available || ImGui::CollapsingHeader(
                scene_text(language,"Input files and calculation details","输入文件与计算详情",
                    "入力ファイルと計算の詳細","Fichiers et détails du calcul"));
            cov::validation::item("input.details");
            if(show_input_panels){
            cov::ui::begin_card("##file_card", 192.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::FileSection, language));
            std::optional<std::filesystem::path> recent_to_load;
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::OpenFile, language),
                              ImVec2(150.0f * ui_scale, 0.0f))) {
                const cov::FileDialogResult dialog = cov::open_molden_file_dialog();
                if (dialog.selected()) {
                    load_inputs({dialog.path});
                } else if (!dialog.cancelled && !dialog.error.empty()) {
                    status = StatusKind::Error;
                    status_detail = dialog.supported
                                        ? dialog.error
                                        : cov::ui::tr(cov::ui::Text::OpenDialogUnsupported, language);
                }
            }
            if (!current_file.empty()) {
                const std::string file_label = std::string(
                    cov::ui::tr(cov::ui::Text::CurrentFile, language)) + ": " +
                    path_to_utf8(current_file.filename());
                disabled_wrapped(file_label.c_str());
            }

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoldenPath, language));
            const float load_width = 76.0f * ui_scale;
            ImGui::SetNextItemWidth(std::max(100.0f,
                ImGui::GetContentRegionAvail().x - load_width - 8.0f));
            ImGui::InputText("##molden_path", path_buffer.data(), path_buffer.size());
            ImGui::SameLine();
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::Load, language),
                              ImVec2(load_width, 0.0f))) {
                load_inputs({path_from_utf8(path_buffer.data())});
            }

            if (!recent_files.empty()) {
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::RecentFiles, language));
                const std::string preview = path_to_utf8(recent_files.front().filename());
                ImGui::SetNextItemWidth(-1.0f);
                if (ImGui::BeginCombo("##recent_files", preview.c_str())) {
                    for (std::size_t i = 0; i < recent_files.size(); ++i) {
                        ImGui::PushID(static_cast<int>(i));
                        const std::string label = path_to_utf8(recent_files[i].filename());
                        if (ImGui::Selectable(label.c_str(), i == 0)) {
                            recent_to_load = recent_files[i];
                        }
                        ImGui::PopID();
                    }
                    ImGui::EndCombo();
                }
            }
            cov::ui::end_card();
            if (recent_to_load) load_inputs({*recent_to_load});
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##wavefunction_card", 270.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::WavefunctionSection, language));
            if (wavefunction) {
                if (ImGui::BeginTable("##wavefunction_metrics", 2,
                                      ImGuiTableFlags_SizingStretchProp |
                                      ImGuiTableFlags_NoSavedSettings)) {
                    const std::string atoms = std::to_string(wavefunction->atoms.size()) + " / 100";
                    const std::string shells = std::to_string(wavefunction->shells.size());
                    const std::string basis = std::to_string(wavefunction->basis_count);
                    const std::string orbitals = std::to_string(wavefunction->orbitals.size());
                    const std::string convention =
                        std::string("D=") + (wavefunction->pure_d ? "5D" : "6D") +
                        "  F=" + (wavefunction->pure_f ? "7F" : "10F") +
                        "  G=" + (wavefunction->pure_g ? "9G" : "15G");
                    std::string state = "—";
                    if (wavefunction->charge_provenance != cov::DataProvenance::Unavailable ||
                        wavefunction->multiplicity_provenance != cov::DataProvenance::Unavailable) {
                        state.clear();
                        if (wavefunction->charge_provenance != cov::DataProvenance::Unavailable) {
                            if (wavefunction->charge > 0) state += "+";
                            state += std::to_string(wavefunction->charge);
                        } else {
                            state += "?";
                        }
                        state += " / ";
                        state += wavefunction->multiplicity_provenance !=
                                     cov::DataProvenance::Unavailable
                                     ? std::to_string(wavefunction->multiplicity) : "?";
                    }
                    const std::string electron_split =
                        wavefunction->electron_counts_provenance ==
                                cov::DataProvenance::Unavailable
                            ? "— / —"
                            : std::to_string(wavefunction->alpha_electrons) + " / " +
                                  std::to_string(wavefunction->beta_electrons);
                    std::string diagnostics = "—";
                    if (wavefunction->scf_convergence !=
                            cov::ScfConvergenceStatus::Unavailable ||
                        wavefunction->stability !=
                            cov::WavefunctionStabilityStatus::Unavailable) {
                        const char* scf = wavefunction->scf_convergence ==
                                                  cov::ScfConvergenceStatus::Converged
                                              ? cov::ui::tr(cov::ui::Text::Converged,language)
                                              : wavefunction->scf_convergence ==
                                                        cov::ScfConvergenceStatus::Failed
                                                    ? cov::ui::tr(cov::ui::Text::Failed,language)
                                                    : "—";
                        const char* stability = wavefunction->stability ==
                                                        cov::WavefunctionStabilityStatus::Stable
                                                    ? cov::ui::tr(cov::ui::Text::Stable,language)
                                                    : wavefunction->stability ==
                                                              cov::WavefunctionStabilityStatus::Unstable
                                                          ? cov::ui::tr(cov::ui::Text::Unstable,language)
                                                          : "—";
                        diagnostics = std::string(scf) + " / " + stability;
                    }
                    std::string spin_squared = "—";
                    if (wavefunction->spin_squared_provenance !=
                        cov::DataProvenance::Unavailable) {
                        char value[64]{};
                        std::snprintf(value, sizeof(value), "%.4f / %.4f",
                                      wavefunction->spin_squared_before_annihilation,
                                      wavefunction->spin_squared_after_annihilation);
                        spin_squared = value;
                    }
                    metric_row(cov::ui::tr(cov::ui::Text::Atoms, language), atoms.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::Shells, language), shells.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::BasisFunctions, language), basis.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::Orbitals, language), orbitals.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::ShellConvention, language), convention.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::ChargeMultiplicity, language),
                               state.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::AlphaBetaElectrons, language),
                               electron_split.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::SCFStability, language),
                               diagnostics.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::SpinSquared, language),
                               spin_squared.c_str());
                    ImGui::EndTable();
                }
            } else {
                ImGui::TextDisabled("—");
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##frame_tracking_card", 156.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::FrameTracking,
                                               language));
            if (frame_tracking) {
                if (ImGui::BeginTable("##frame_tracking_metrics", 2,
                                      ImGuiTableFlags_SizingStretchProp |
                                      ImGuiTableFlags_NoSavedSettings)) {
                    std::size_t matched_members = 0u;
                    for (const auto& match : frame_tracking->matches) {
                        matched_members += match.from_members.size();
                    }
                    const std::string matched =
                        std::to_string(frame_tracking->matches.size()) +
                        " (" + std::to_string(matched_members) + " MO)";
                    const std::string unmatched =
                        std::to_string(frame_tracking->unmatched_from.size()) +
                        " / " +
                        std::to_string(frame_tracking->unmatched_to.size());
                    metric_row(cov::ui::tr(
                                   cov::ui::Text::AtomMappingCompatibility,
                                   language),
                               cov::ui::tr(
                                   frame_tracking->atom_mapping_compatible
                                       ? cov::ui::Text::Compatible
                                       : cov::ui::Text::Incompatible,
                                   language));
                    metric_row(cov::ui::tr(cov::ui::Text::MatchedSubspaces,
                                           language),
                               matched.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::UnmatchedSubspaces,
                                           language),
                               unmatched.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::TrackingOptimisation,
                                           language),
                               cov::ui::tr(
                                   frame_tracking->composite_optimisation_truncated
                                       ? cov::ui::Text::ConservativeFallback
                                       : cov::ui::Text::ExactOrNotNeeded,
                                   language));
                    ImGui::EndTable();
                }
            } else {
                ImGui::TextDisabled("%s", cov::ui::tr(
                    cov::ui::Text::NoPreviousFrame, language));
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            }
            {
                cov::ui::OrbitalUIActions orbital_actions;
                const bool show_browser=!aomo_available || ImGui::CollapsingHeader(
                    scene_text(language,"Full MO browser","完整 MO 浏览器","全 MO 一覧","Liste complète des OM"));
                cov::validation::item("browser.expand");
                if(show_browser){
                cov::validation::anchor("panel.browser");
                cov::ui::begin_card("##orbital_browser_card", 0);
                cov::ui::section_title(cov::ui::tr(cov::ui::Text::OrbitalBrowser, language));
                if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                    cov::ui::draw_orbital_browser(*wavefunction, canonical_mo_index, orbital_ui,
                                                  language, ui_scale, orbital_actions);
                } else {
                    ImGui::TextDisabled("—");
                }
                cov::ui::end_card();
                ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

                }

                cov::validation::anchor("panel.diagram");
                cov::ui::begin_card("##energy_diagram_card", 0);
                cov::ui::section_title(cov::ui::tr(cov::ui::Text::EnergyDiagram, language));
                cov::ui::OrbitalUIActions diagram_actions;
                if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                    cov::ui::draw_energy_diagram(*wavefunction, canonical_mo_index, orbital_ui,
                                                 language, ui_scale, diagram_actions);
                } else {
                    ImGui::TextDisabled("—");
                }
                cov::ui::end_card();
                ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

                if (orbital_actions.select_orbital) pending_mo_index = orbital_actions.select_orbital;
                if (diagram_actions.select_orbital) pending_mo_index = diagram_actions.select_orbital;
                canonical_requested=orbital_actions.select_orbital.has_value() || diagram_actions.select_orbital.has_value();
                if (nbo_ui.focus.pending_canonical_selection) {
                    pending_mo_index = nbo_ui.focus.pending_canonical_selection;
                    canonical_requested=true;
                    nbo_ui.focus.pending_canonical_selection.reset();
                }
                const bool export_requested = orbital_actions.export_diagram || diagram_actions.export_diagram;
                if (export_requested && wavefunction) {
                    std::filesystem::path base = current_file.empty()
                                                     ? std::filesystem::current_path() / "mo_diagram"
                                                     : current_file;
                    base = cov::validation::export_base(base);
                    const auto snapshot=diagram_actions.drawn_diagram;
                    cov::MODiagramExportResult result;
                    if (snapshot) result=cov::export_mo_diagram_bundle(*snapshot,base);
                    else result.error="No current diagram view is available for export";
    #ifdef COV_ENABLE_VALIDATION
                    cov::validation::record("export.actual","{\"base\":"+cov::validation::quote(path_to_utf8(base))+
                        ",\"snapshot_id\":"+(snapshot?cov::validation::quote(snapshot->data.view->id):"null")+
                        ",\"mode\":"+(snapshot?std::to_string(static_cast<int>(snapshot->data.mode)):"null")+
                        ",\"selected_index\":"+(snapshot && snapshot->data.view->inspected_orbital_index
                            ?std::to_string(*snapshot->data.view->inspected_orbital_index):"null")+
                        ",\"success\":"+((result.svg&&result.png&&result.json&&result.csv)?"true":"false")+"}");
    #endif
                    if (result.svg && result.png && result.json && result.csv) {
                        export_analysis_companions(base);
                        status = StatusKind::Exported;
                        status_detail = path_to_utf8(result.svg_path.parent_path() /
                            result.svg_path.stem()) + ".{png,svg,json,csv}";
                    } else {
                        status = StatusKind::Error;
                        status_detail = result.error.empty()
                                            ? cov::ui::tr(cov::ui::Text::ExportFailed, language)
                                            : result.error;
                    }
                }
            }

            cov::validation::anchor("panel.nbo");
            cov::ui::begin_card("##nbo_card", 0);
            const auto nbo_actions = cov::ui::draw_nbo_panel(
                nbo_ui, language, static_cast<bool>(wavefunction),
                nbo_wavefunction.has_value(), active_view().kind, ui_scale,
                wavefunction ? &*wavefunction : nullptr, canonical_mo_index,
                orbital_ui.diagram_cache.snapshot.get());
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));
            const bool attach_requested = nbo_actions.attach;
            if(nbo_actions.choose_input)choose_nbo_input=true;
            std::optional<bool> requested_set;
            if (nbo_actions.canonical_set) requested_set=false;
            if (nbo_actions.nbo_set) requested_set=true;
            if(requested_set && *requested_set && nbo_mo_index==std::numeric_limits<std::size_t>::max()){
                requested_set.reset();nbo_ui.error="Select a specific available NBO orbital first";
            }
            if (nbo_actions.selected_orbital && nbo_ui.dataset) {
                nbo_ui.selected_orbital=*nbo_actions.selected_orbital;
                nbo_mo_index=nbo_ui.selected_orbital;
                if (nbo_active) pending_mo_index=nbo_mo_index;
            }
            if (nbo_actions.export_bundle && nbo_ui.dataset) {
                try {
                    auto base=nbo_ui.export_path[0]
                        ? path_from_utf8(nbo_ui.export_path.data())
                        : path_from_utf8(nbo_ui.dataset->source.path);
                    base=cov::validation::export_base(base);
                    cov::ui::export_nbo_bundle(*nbo_ui.dataset,nbo_ui.selected_orbital,base,
                                               wavefunction ? &*wavefunction : nullptr,
                                               canonical_mo_index,nbo_ui.contribution_threshold,
                                               active_view(),integration?&*integration:nullptr,
                                               orbital_ui.diagram_cache.snapshot.get(),&nbo_ui.focus);
                    export_analysis_companions(base);
                    nbo_ui.export_status=path_to_utf8(base)+".{nbo.json,npa.csv,nao.csv,nbo.csv,wiberg.csv,e2.csv,e2-sections.csv,view.json,view.svg,view.png,focus.json,focus.csv,focus.groups.csv,focus.svg,focus.png}";
                    cov::validation::record("nbo.export","{\"base\":"+
                        cov::validation::quote(path_to_utf8(base))+
                        ",\"selected_index\":"+std::to_string(nbo_ui.selected_orbital)+
                        ",\"source\":"+cov::validation::quote(nbo_ui.dataset->source.path)+"}");
                } catch (const std::exception& e) { nbo_ui.export_status=e.what(); }
            }

            cov::ui::begin_card("##render_card", 545.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::RenderingSection, language));
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoleculeStyle, language));
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::BeginCombo("##molecule_style",
                                  molecule_style_name(molecule_render.style, language))) {
                for (const cov::MoleculeStyle style : {
                         cov::MoleculeStyle::MediumBallAndStick,
                         cov::MoleculeStyle::StickDelocalisation}) {
                    const bool selected = molecule_render.style == style;
                    if (ImGui::Selectable(molecule_style_name(style, language), selected)) {
                        molecule_render.style = style;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
            if (molecule_render.style == cov::MoleculeStyle::StickDelocalisation) {
                ImGui::TextDisabled("%s",
                    cov::ui::tr(cov::ui::Text::DelocalisationHeuristic, language));
            }

            const OrbitalAppearanceText appearance = orbital_appearance_text(language);
            ImGui::TextDisabled("%s", appearance.material);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::BeginCombo("##orbital_material",
                                  orbital_material_name(orbital_material, appearance))) {
                for (const cov::OrbitalMaterial material : {
                         cov::OrbitalMaterial::Standard,
                         cov::OrbitalMaterial::Glass}) {
                    const bool selected = orbital_material == material;
                    if (ImGui::Selectable(orbital_material_name(material, appearance), selected)) {
                        orbital_material = material;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }

            ImGui::TextDisabled("%s", appearance.surface);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::BeginCombo("##orbital_surface",
                                  orbital_surface_name(orbital_surface_mode, appearance))) {
                for (const cov::OrbitalSurfaceMode mode : {
                         cov::OrbitalSurfaceMode::Solid,
                         cov::OrbitalSurfaceMode::Wire,
                         cov::OrbitalSurfaceMode::SolidWire}) {
                    const bool selected = orbital_surface_mode == mode;
                    if (ImGui::Selectable(orbital_surface_name(mode, appearance), selected)) {
                        orbital_surface_mode = mode;
                    }
                    if (selected) ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }
            ImGui::TextDisabled("%s", appearance.auto_light);

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::AtomSize, language));
            ImGui::SliderFloat("##atom_size", &molecule_render.atom_scale, 0.55f, 1.8f, "%.2f");
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::BondSize, language));
            ImGui::SliderFloat("##bond_size", &molecule_render.bond_scale, 0.5f, 2.0f, "%.2f");
            ImGui::Checkbox(cov::ui::tr(cov::ui::Text::ShowHydrogens, language),
                            &molecule_render.show_hydrogens);
            ImGui::Checkbox(cov::ui::tr(cov::ui::Text::ShowCoordinationContacts, language),
                            &molecule_render.show_coordination_contacts);
            ImGui::Checkbox(cov::ui::tr(cov::ui::Text::ShowMulticentreSupport, language),
                            &molecule_render.show_multicentre_support);
            ImGui::Checkbox(cov::ui::tr(
                                cov::ui::Text::ShowPolyhedralCageSupport,language),
                            &molecule_render.show_polyhedral_cage_support);
            ImGui::Checkbox(cov::ui::tr(cov::ui::Text::ShowWeakInteractions, language),
                            &molecule_render.show_weak_interactions);
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("%s", cov::ui::tr(cov::ui::Text::WeakInteractionsHint,
                                                     language));
            }

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoleculeOpacity, language));
            ImGui::SliderFloat("##molecule_opacity", &molecule_render.molecule_opacity,
                               0.15f, 1.0f, "%.2f");
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::OrbitalOpacity, language));
            ImGui::SliderFloat("##orbital_opacity", &molecule_render.orbital_opacity,
                               0.02f, 1.0f, "%.2f");

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::Isovalue, language));
            ImGui::SliderFloat("##isovalue", &isovalue, 0.002f, 0.12f, "%.4f",
                               ImGuiSliderFlags_Logarithmic);

            constexpr int resolutions[] = {64, 128, 256, 512};
            int resolution_index = 1;
            for (int i = 0; i < 4; ++i) {
                if (resolutions[i] == resolution) resolution_index = i;
            }
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::Grid, language));
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::Combo("##grid", &resolution_index,
                             "64³\0" "128³\0" "256³\0" "512³\0")) {
                resolution = resolutions[resolution_index];
                resize_and_recompute = true;
                recompute = true;
            }

            const float button_gap = ImGui::GetStyle().ItemSpacing.x;
            const float half_button = (ImGui::GetContentRegionAvail().x - button_gap) * 0.5f;
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::RecomputeGrid, language),
                              ImVec2(half_button, 0.0f))) {
                recompute = true;
            }
            ImGui::SameLine();
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::ResetCamera, language),
                              ImVec2(half_button, 0.0f))) {
                camera = {};
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##performance_card", 190.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::PerformanceSection, language));
            if (evaluator) {
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::CUDADevice, language));
                ImGui::TextUnformatted(evaluator->device_name());
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::LastKernel, language));
                ImGui::Text("%.3f ms", evaluator->last_kernel_ms());
                cov::ui::status_badge(cov::ui::tr(cov::ui::Text::GPUResident, language),
                                      cov::ui::Tone::Success);
            } else {
                ImGui::TextDisabled("CUDA —");
            }
            ImGui::TextDisabled("%s: %s",
                                cov::ui::tr(cov::ui::Text::FontStatus, language),
                                cov::ui::font_status());
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::InteractionHint, language));
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::IsovalueHint, language));
            cov::ui::end_card();

            ImGui::EndChild();
            ImGui::End();
            cov::validation::field("language",std::to_string(static_cast<int>(language)));
            identity();
            cov::validation::ui_frame(mo_index,pending_mo_index.value_or(mo_index));

            ImGui::Render();
            ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
            cov::validation::end_frame(fb_w,fb_h,mo_index,orbital_ui,active_wavefunction());

            glfwSwapBuffers(window);
            if (attach_requested) {
                try {
                    cov::NboReadOptions options;
                    if (nbo_ui.analysis_segment >= 0)
                        options.analysis_segment = static_cast<std::size_t>(nbo_ui.analysis_segment);
                    if (nbo_ui.archive47[0]) options.archive47=path_from_utf8(nbo_ui.archive47.data());
                    if (nbo_ui.aonbo[0]) options.aonbo=path_from_utf8(nbo_ui.aonbo.data());
                    if (nbo_ui.nbomo[0]) options.nbomo=path_from_utf8(nbo_ui.nbomo.data());
                    if (nbo_ui.naomo[0]) options.naomo=path_from_utf8(nbo_ui.naomo.data());
                    if (nbo_ui.aonao[0]) options.aonao=path_from_utf8(nbo_ui.aonao.data());
                    if (nbo_ui.naonbo[0]) options.naonbo=path_from_utf8(nbo_ui.naonbo.data());
                    auto dataset=cov::read_nbo(path_from_utf8(nbo_ui.path.data()),options);
                    attach_integration(cov::integrate_nbo(*wavefunction,dataset));
                } catch (const std::exception& e) {
                    nbo_ui.error=e.what();
                    cov::validation::record("nbo.attach.error","{\"reason\":"+
                        cov::validation::quote(nbo_ui.error)+"}");
                }
            }

            if(nbo_ui.pending_candidate && nbo_ui.input_discovery){
                const auto selected=*nbo_ui.pending_candidate;nbo_ui.pending_candidate.reset();
                if(selected<nbo_ui.input_discovery->candidates.size()){
                    const auto candidate=nbo_ui.input_discovery->candidates[selected];
                    try {apply_candidate(candidate);}catch(const std::exception& e){
                        clear_integration();nbo_ui.input_discovery.reset();nbo_ui.error=e.what();
                        nbo_ui.input_status="NBO association failed; Gaussian view is retained";
                        status=StatusKind::Error;status_detail=nbo_ui.error;
                    }
                }
            }
            if(nbo_ui.aomo.pending_selection){
                const auto selection=*nbo_ui.aomo.pending_selection;nbo_ui.aomo.pending_selection.reset();
                try{apply_selection(selection);}catch(const std::exception& e){nbo_ui.aomo.status=e.what();
                    nbo_ui.error=e.what();status=StatusKind::Error;status_detail=e.what();
                    cov::validation::record("aomo.selection.error","{\"reason\":"+cov::validation::quote(e.what())+"}");}
            }
            if(nbo_ui.aomo.export_requested){
                nbo_ui.aomo.export_requested=false;
                if(integration && nbo_ui.aomo.drawn_snapshot){
                    const auto base=cov::validation::export_base(nbo_ui.aomo.export_path[0]?
                        path_from_utf8(nbo_ui.aomo.export_path.data()):current_file);
                    const auto result=cov::ui::export_nbo_aomo_bundle(*nbo_ui.aomo.drawn_snapshot,*integration,base);
                    if(result.json && result.svg && result.png && result.csv)
                        export_analysis_companions(base);
                    nbo_ui.aomo.export_status=result.error.empty()?path_to_utf8(result.svg_path):result.error;
                }
            }
            if(canonical_requested){
                const auto selected=pending_mo_index;
                try{activate_set(false);pending_mo_index=selected;}
                catch(const std::exception& e){nbo_ui.error=e.what();}
            }
            if (requested_set) {
                try { activate_set(*requested_set); }
                catch (const std::exception& e) { nbo_ui.error=e.what(); }
            }
            // Selection debounce: at most the latest requested orbital is evaluated
            // once at the end of this frame. Browser hover/filtering never launches CUDA.
            if (pending_mo_index && active_wavefunction() &&
                *pending_mo_index < active_wavefunction()->orbitals.size()) {
                if (*pending_mo_index != mo_index) {
                    mo_index = *pending_mo_index;
                    if (nbo_active) nbo_mo_index=mo_index;
                    else canonical_mo_index=mo_index;
                    recompute = true;
                }
                pending_mo_index.reset();
            }

            if (recompute) {
                try {
                    evaluate_now();
                } catch (const std::exception& e) {
                    status = StatusKind::Error;
                    status_detail = e.what();
                    recompute = false;
                }
            }

            if (cov::validation::done()) {exit_code=cov::validation::result();break;}
        }

        if (evaluator) evaluator->detach_gl_texture();
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Fatal error: %s\n", e.what());
        exit_code = 1;
    }

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return exit_code;
}
