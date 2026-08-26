#include "cov/cuda_orbital.hpp"
#include "cov/gl_api.hpp"
#include "cov/molden_parser.hpp"
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

namespace {

std::string g_dropped_path;

enum class StatusKind {
    Ready,
    Parsing,
    Loaded,
    GridUpdated,
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
        case StatusKind::Error: return cov::ui::tr(Text::Error, language);
        default: return cov::ui::tr(Text::Ready, language);
    }
}

cov::ui::Tone status_tone(const StatusKind status) {
    switch (status) {
        case StatusKind::Loaded:
        case StatusKind::GridUpdated: return cov::ui::Tone::Success;
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

void metric_row(const char* label, const char* value) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextDisabled("%s", label);
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(value);
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
        1440, 900, "CUDA Orbital Visualisation", nullptr, nullptr);
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
        StatusKind status = StatusKind::Ready;
        std::string status_detail;

        std::size_t mo_index = 0;
        float isovalue = 0.03f;
        int resolution = 128;
        std::array<char, 1024> path_buffer{};

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
                status_detail = path.string();

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
                grid_box = new_box;

                renderer.resize_volume(resolution, resolution, resolution);
                evaluator->attach_gl_texture(renderer.volume_texture());
                evaluator->evaluate(mo_index, grid_box,
                                    resolution, resolution, resolution);
                status = StatusKind::Loaded;
                status_detail = path.filename().string();
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
                std::snprintf(path_buffer.data(), path_buffer.size(), "%s",
                              g_dropped_path.c_str());
                load_file(path_from_utf8(g_dropped_path));
                g_dropped_path.clear();
            }

            int fb_w = 1, fb_h = 1;
            glfwGetFramebufferSize(window, &fb_w, &fb_h);
            glViewport(0, 0, fb_w, fb_h);
            glClearColor(0.025f, 0.031f, 0.043f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            if (wavefunction) {
                renderer.render_volume(fb_w, fb_h, isovalue, camera);
                renderer.render_geometry(*wavefunction, grid_box, fb_w, fb_h, camera);
            }

            ImGui_ImplOpenGL2_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            ImGuiIO& io = ImGui::GetIO();
            double mx = 0.0, my = 0.0;
            glfwGetCursorPos(window, &mx, &my);
            if (!io.WantCaptureMouse && glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
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
            const float panel_width = std::min(430.0f * ui_scale,
                                               std::max(320.0f, io.DisplaySize.x * 0.42f));
            const float panel_height = std::max(300.0f, io.DisplaySize.y - margin * 2.0f);

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

            cov::ui::begin_card("##file_card", 126.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::FileSection, language));
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoldenPath, language));
            const float load_width = 76.0f * ui_scale;
            ImGui::SetNextItemWidth(std::max(80.0f, ImGui::GetContentRegionAvail().x - load_width - 8.0f));
            ImGui::InputText("##molden_path", path_buffer.data(), path_buffer.size());
            ImGui::SameLine();
            if (ImGui::Button(cov::ui::tr(cov::ui::Text::Load, language),
                              ImVec2(load_width, 0.0f))) {
                load_file(path_from_utf8(path_buffer.data()));
            }
            ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::IdleHint, language));
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##wavefunction_card", 176.0f * ui_scale);
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

            cov::ui::begin_card("##orbital_card", 214.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::OrbitalSection, language));
            if (wavefunction && evaluator && !wavefunction->orbitals.empty()) {
                ImGui::TextDisabled("%s", cov::ui::tr(cov::ui::Text::MoldenMO, language));
                int mo_ui = static_cast<int>(mo_index) + 1;
                if (ImGui::SliderInt("##molden_mo", &mo_ui, 1,
                                     static_cast<int>(wavefunction->orbitals.size()))) {
                    mo_index = static_cast<std::size_t>(mo_ui - 1);
                    recompute = true;
                }

                const auto& mo = wavefunction->orbitals[mo_index];
                if (ImGui::BeginTable("##orbital_metrics", 2,
                                      ImGuiTableFlags_SizingStretchProp |
                                      ImGuiTableFlags_NoSavedSettings)) {
                    const std::string index = std::to_string(mo_index);
                    char energy[64]{};
                    char occupation[64]{};
                    std::snprintf(energy, sizeof(energy), "%.8f Ha", mo.energy_hartree);
                    std::snprintf(occupation, sizeof(occupation), "%.3f", mo.occupation);
                    metric_row(cov::ui::tr(cov::ui::Text::InternalIndex, language), index.c_str());
                    metric_row(cov::ui::tr(cov::ui::Text::Energy, language), energy);
                    metric_row(cov::ui::tr(cov::ui::Text::Occupation, language), occupation);
                    metric_row(cov::ui::tr(cov::ui::Text::Spin, language),
                               cov::ui::tr(mo.spin == cov::Spin::Beta ? cov::ui::Text::Beta
                                                                      : cov::ui::Text::Alpha,
                                           language));
                    metric_row(cov::ui::tr(cov::ui::Text::Symmetry, language),
                               mo.symmetry.empty() ? cov::ui::tr(cov::ui::Text::NoneValue, language)
                                                   : mo.symmetry.c_str());
                    ImGui::EndTable();
                }
            } else {
                ImGui::TextDisabled("—");
            }
            cov::ui::end_card();
            ImGui::Dummy(ImVec2(0, 7.0f * ui_scale));

            cov::ui::begin_card("##render_card", 192.0f * ui_scale);
            cov::ui::section_title(cov::ui::tr(cov::ui::Text::RenderingSection, language));
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
