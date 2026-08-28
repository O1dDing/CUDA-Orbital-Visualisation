#include "cov/cuda_orbital.hpp"
#include "cov/file_dialog.hpp"
#include "cov/gl_api.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/molden_parser.hpp"
#include "cov/molecule_style.hpp"
#include "cov/orbital_ui.hpp"
#include "cov/ui.hpp"
#include "cov/volume_renderer.hpp"

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
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string g_dropped_path;

enum class StatusKind {
    Ready,
    Parsing,
    Loaded,
    GridUpdated,
    Exported,
    Error,
};

void drop_callback(GLFWwindow*, const int count, const char** paths) {
    if (count > 0 && paths && paths[0]) {
        g_dropped_path = paths[0];
    }
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
    std::size_t selected = 0;
    for (std::size_t i = 0; i < wf.orbitals.size(); ++i) {
        if (wf.orbitals[i].occupation > 1.0e-4f) selected = i;
    }
    return selected;
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

void metric_row(const char* label, const char* value) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextDisabled("%s", label);
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(value);
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
    if (!glfwInit()) {
        std::fprintf(stderr, "GLFW initialisation failed\n");
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);

    GLFWwindow* window = glfwCreateWindow(
        1500, 940, "CUDA Orbital Visualisation", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "Unable to create OpenGL window\n");
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glfwSetDropCallback(window, drop_callback);

    float x_scale = 1.0f;
    float y_scale = 1.0f;
    glfwGetWindowContentScale(window, &x_scale, &y_scale);
    const float ui_scale = std::clamp(std::max(x_scale, y_scale), 1.0f, 1.75f);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& startup_io = ImGui::GetIO();
    startup_io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    cov::ui::apply_theme(ui_scale);
    cov::ui::configure_fonts(16.5f * ui_scale);
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    int exit_code = 0;
    try {
        cov::VolumeRenderer renderer;
        cov::OrbitCamera camera;

        std::optional<cov::Wavefunction> wavefunction;
        std::unique_ptr<cov::CudaOrbitalEvaluator> evaluator;
        cov::GridBox grid_box;

        cov::ui::Language language = cov::ui::Language::English;
        cov::ui::OrbitalUIState orbital_ui;
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

        auto evaluate_now = [&]() {
            if (!wavefunction || !evaluator || wavefunction->orbitals.empty()) return;
            if (resize_and_recompute || renderer.nx() != resolution) {
                evaluator->detach_gl_texture();
                renderer.resize_volume(resolution, resolution, resolution);
                evaluator->attach_gl_texture(renderer.volume_texture());
                resize_and_recompute = false;
            }
            evaluator->evaluate(mo_index, grid_box,
                                resolution, resolution, resolution);
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
                auto wf = cov::parse_molden(path, options);
                const auto new_mo = initial_orbital(wf);
                const auto new_box = make_grid_box(wf);

                if (evaluator) evaluator->detach_gl_texture();
                evaluator.reset();
                wavefunction = std::move(wf);
                evaluator = std::make_unique<cov::CudaOrbitalEvaluator>(*wavefunction);
                mo_index = new_mo;
                pending_mo_index.reset();
                grid_box = new_box;

                renderer.resize_volume(resolution, resolution, resolution);
                evaluator->attach_gl_texture(renderer.volume_texture());
                evaluator->evaluate(mo_index, grid_box,
                                    resolution, resolution, resolution);
                current_file = path;
                copy_path_to_buffer(path, path_buffer);
                push_recent(recent_files, path);
                status = StatusKind::Loaded;
                status_detail = path_to_utf8(path.filename());
            } catch (const std::exception& e) {
                status = StatusKind::Error;
                status_detail = e.what();
            }
        };

        if (argc >= 2) {
            const std::string p = argv[1];
            std::snprintf(path_buffer.data(), path_buffer.size(), "%s", p.c_str());
            load_file(path_from_utf8(p));
        }

        double last_x = 0.0;
        double last_y = 0.0;
        glfwGetCursorPos(window, &last_x, &last_y);

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();

            if (!g_dropped_path.empty()) {
                load_file(path_from_utf8(g_dropped_path));
                g_dropped_path.clear();
            }

            int fb_w = 1, fb_h = 1;
            glfwGetFramebufferSize(window, &fb_w, &fb_h);
            glViewport(0, 0, fb_w, fb_h);
            glClearColor(0.025f, 0.031f, 0.043f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            if (wavefunction) {
                // Geometry is rendered first into colour + depth. The implicit
                // orbital surface then depth-tests against it and alpha-blends
                // over geometry behind the lobe. This preserves correct front/
                // back occlusion while allowing a genuinely transparent glass skin.
                renderer.render_geometry(*wavefunction, grid_box, fb_w, fb_h, camera,
                                         molecule_render);
                renderer.render_volume(fb_w, fb_h, isovalue, camera,
                                       molecule_render.orbital_opacity,
                                       orbital_material, orbital_surface_mode);
            }

            ImGui_ImplOpenGL2_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            ImGuiIO& io = ImGui::GetIO();
            double mx = 0.0, my = 0.0;
            glfwGetCursorPos(window, &mx, &my);
            if (!io.WantCaptureMouse &&
                glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
                camera.yaw += static_cast<float>((mx - last_x) * 0.007);
                camera.pitch += static_cast<float>((my - last_y) * 0.007);
                camera.pitch = std::clamp(camera.pitch, -1.45f, 1.45f);
            }
            if (!io.WantCaptureMouse && std::abs(io.MouseWheel) > 0.0f) {
                camera.distance *= std::pow(0.88f, io.MouseWheel);
                camera.distance = std::clamp(camera.distance, 1.1f, 6.0f);
            }
            last_x = mx;
            last_y = my;

            const float margin = 14.0f * ui_scale;
            const float panel_width = std::min(540.0f * ui_scale,
                                               std::max(370.0f, io.DisplaySize.x * 0.46f));
            const float panel_height = std::max(320.0f, io.DisplaySize.y - margin * 2.0f);

            ImGui::SetNextWindowPos(ImVec2(margin, margin), ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(panel_width, panel_height), ImGuiCond_Always);
            ImGui::SetNextWindowBgAlpha(0.965f);
            constexpr ImGuiWindowFlags panel_flags =
                ImGuiWindowFlags_NoTitleBar |
                ImGuiWindowFlags_NoMove |
                ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoSavedSettings;

            ImGui::Begin("##cov_control_panel", nullptr, panel_flags);

            if (ImGui::BeginTable("##cov_header", 2,
                                  ImGuiTableFlags_SizingStretchProp |
                                  ImGuiTableFlags_NoSavedSettings)) {
                ImGui::TableSetupColumn("##brand", ImGuiTableColumnFlags_WidthStretch);
                ImGui::TableSetupColumn("##language", ImGuiTableColumnFlags_WidthFixed,
                                        142.0f * ui_scale);
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(cov::ui::tr(cov::ui::Text::AppTitle, language));
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::Tagline, language));
                ImGui::TableNextColumn();
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::LanguageLabel, language));
                int language_index = static_cast<int>(language);
                ImGui::SetNextItemWidth(-1.0f);
                if (ImGui::Combo("##language_combo", &language_index,
                                 "English\0简体中文\0日本語\0Français\0")) {
                    language = static_cast<cov::ui::Language>(language_index);
                    glfwSetWindowTitle(window,
                        cov::ui::tr(cov::ui::Text::AppTitle, language));
                }
                ImGui::EndTable();
            }

            ImGui::Spacing();
            cov::ui::status_badge(status_label(status, language), status_tone(status));
            if (!status_detail.empty()) {
                ImGui::SameLine();
                ImGui::TextDisabled("%s", status_detail.c_str());
            } else {
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::IdleHint, language));
            }
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::ExperimentalNote, language));
            ImGui::Separator();
            ImGui::Spacing();

            ImGui::BeginChild("##cov_panel_scroll", ImVec2(0, 0), false,
                              ImGuiWindowFlags_None);

            cov::ui::begin_card("##file_card", 192.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::FileSection, language));
            std::optional<std::filesystem::path> recent_to_load;
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::OpenFile, language),
                              ImVec2(150.0f * ui_scale, 0.0f))) {
                const cov::FileDialogResult dialog = cov::open_molden_file_dialog();
                if (dialog.selected()) {
                    load_file(dialog.path);
                } else if (!dialog.cancelled && !dialog.error.empty()) {
                    status = StatusKind::Error;
                    status_detail = dialog.supported
                                        ? dialog.error
                                        : cov::ui::tr(cov::ui::Text::OpenDialogUnsupported, language);
                }
            }
            if (!current_file.empty()) {
                ImGui::SameLine();
                ImGui::TextDisabled("%s: %s",
                    cov::ui::tr(cov::ui::Text::CurrentFile, language),
                    path_to_utf8(current_file.filename()).c_str());
            }

            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoldenPath, language));
            const float load_width = 76.0f * ui_scale;
            ImGui::SetNextItemWidth(std::max(100.0f,
                ImGui::GetContentRegionAvail().x - load_width - 8.0f));
            ImGui::InputText("##molden_path", path_buffer.data(), path_buffer.size());
            ImGui::SameLine();
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::Load, language),
                              ImVec2(load_width, 0.0f))) {
                load_file(path_from_utf8(path_buffer.data()));
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
            if (recent_to_load) load_file(*recent_to_load);
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##wavefunction_card", 174.0f * ui_scale);
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
                    metric_row(cov::ui::tr(cov::ui::Text::Atoms, language), atoms.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::Shells, language), shells.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::BasisFunctions, language), basis.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::Orbitals, language), orbitals.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::ShellConvention, language), convention.c_str());
                    ImGui::EndTable();
                }
            } else {
                ImGui::TextDisabled("—");
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##orbital_browser_card", 620.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::OrbitalBrowser, language));
            cov::ui::OrbitalUIActions orbital_actions;
            if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                cov::ui::draw_orbital_browser(*wavefunction, mo_index, orbital_ui,
                                              language, ui_scale, orbital_actions);
            } else {
                ImGui::TextDisabled("—");
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##energy_diagram_card", 430.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::EnergyDiagram, language));
            cov::ui::OrbitalUIActions diagram_actions;
            if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                cov::ui::draw_energy_diagram(*wavefunction, mo_index, orbital_ui,
                                             language, ui_scale, diagram_actions);
            } else {
                ImGui::TextDisabled("—");
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            if (orbital_actions.select_orbital) pending_mo_index = orbital_actions.select_orbital;
            if (diagram_actions.select_orbital) pending_mo_index = diagram_actions.select_orbital;
            const bool export_requested = orbital_actions.export_diagram || diagram_actions.export_diagram;
            if (export_requested && wavefunction) {
                cov::MODiagramOptions options;
                options.energy_unit = orbital_ui.energy_unit;
                options.energy_axis_mode = orbital_ui.energy_axis_mode;
                options.degeneracy = orbital_ui.degeneracy;
                options.filter = orbital_ui.filter;
                options.selected_index = mo_index;
                options.neighbourhood = static_cast<std::size_t>(std::max(2, orbital_ui.diagram_neighbourhood));
                std::filesystem::path base = current_file.empty()
                                                 ? std::filesystem::current_path() / "mo_diagram"
                                                 : current_file;
                const auto result = cov::export_mo_diagram_bundle(*wavefunction, options, base);
                if (result.svg && result.png && result.json && result.csv) {
                    status = StatusKind::Exported;
                    if (base.has_extension()) base.replace_extension();
                    status_detail = path_to_utf8(base) + ".mo.{png,svg,json,csv}";
                } else {
                    status = StatusKind::Error;
                    status_detail = result.error.empty()
                                        ? cov::ui::tr(cov::ui::Text::ExportFailed, language)
                                        : result.error;
                }
            }

            cov::ui::begin_card("##render_card", 470.0f * ui_scale);
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

            // Selection debounce: at most the latest requested orbital is evaluated
            // once at the end of this frame. Browser hover/filtering never launches CUDA.
            if (pending_mo_index && wavefunction &&
                *pending_mo_index < wavefunction->orbitals.size()) {
                if (*pending_mo_index != mo_index) {
                    mo_index = *pending_mo_index;
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

            ImGui::Render();
            ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());

            glfwSwapBuffers(window);
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
