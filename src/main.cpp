#include "cov/cuda_orbital.hpp"
#include "cov/molden_parser.hpp"
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

const char* spin_name(const cov::Spin spin) {
    return spin == cov::Spin::Beta ? "Beta" : "Alpha";
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

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    int exit_code = 0;
    try {
        cov::VolumeRenderer renderer;
        cov::OrbitCamera camera;

        std::optional<cov::Wavefunction> wavefunction;
        std::unique_ptr<cov::CudaOrbitalEvaluator> evaluator;
        cov::GridBox grid_box;

        std::size_t mo_index = 0;
        float isovalue = 0.03f;
        int resolution = 128;
        std::string status = "Drop a .molden file here, pass it on the command line, or enter a path.";
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
            status = "MO grid updated on " + std::string(evaluator->device_name()) +
                     " in " + std::to_string(evaluator->last_kernel_ms()) + " ms";
            recompute = false;
        };

        auto load_file = [&](const std::filesystem::path& path) {
            try {
                status = "Parsing " + path.string() + "...";
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
                status = "Loaded " + path.filename().string() +
                         " | " + std::to_string(wavefunction->atoms.size()) + " atoms | " +
                         std::to_string(wavefunction->basis_count) + " basis functions | " +
                         std::to_string(wavefunction->orbitals.size()) + " MOs | " +
                         std::to_string(evaluator->last_kernel_ms()) + " ms";
            } catch (const std::exception& e) {
                status = std::string("ERROR: ") + e.what();
            }
        };

        if (argc >= 2) {
            const std::string p = argv[1];
            std::snprintf(path_buffer.data(), path_buffer.size(), "%s", p.c_str());
            load_file(p);
        }

        double last_x = 0.0;
        double last_y = 0.0;
        glfwGetCursorPos(window, &last_x, &last_y);

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();

            if (!g_dropped_path.empty()) {
                std::snprintf(path_buffer.data(), path_buffer.size(), "%s",
                              g_dropped_path.c_str());
                load_file(g_dropped_path);
                g_dropped_path.clear();
            }

            int fb_w = 1, fb_h = 1;
            glfwGetFramebufferSize(window, &fb_w, &fb_h);
            glViewport(0, 0, fb_w, fb_h);
            glClearColor(0.025f, 0.028f, 0.035f, 1.0f);
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

            ImGui::SetNextWindowSize(ImVec2(420, 560), ImGuiCond_FirstUseEver);
            ImGui::Begin("Orbital controls");

            ImGui::TextWrapped("%s", status.c_str());
            ImGui::Separator();

            ImGui::InputText("Molden path", path_buffer.data(), path_buffer.size());
            ImGui::SameLine();
            if (ImGui::Button("Load")) {
                load_file(path_buffer.data());
            }

            if (wavefunction && evaluator) {
                ImGui::Separator();
                ImGui::Text("CUDA: %s", evaluator->device_name());
                ImGui::Text("Atoms: %zu / 100", wavefunction->atoms.size());
                ImGui::Text("Shells: %zu", wavefunction->shells.size());
                ImGui::Text("Basis functions: %u", wavefunction->basis_count);
                ImGui::Text("Orbitals: %zu", wavefunction->orbitals.size());
                ImGui::Text("Spherical: D=%s F=%s G=%s",
                            wavefunction->pure_d ? "5D" : "6D",
                            wavefunction->pure_f ? "7F" : "10F",
                            wavefunction->pure_g ? "9G" : "15G");

                int mo_ui = static_cast<int>(mo_index) + 1;
                if (ImGui::SliderInt("Molden MO (1-based)", &mo_ui, 1,
                                     static_cast<int>(wavefunction->orbitals.size()))) {
                    mo_index = static_cast<std::size_t>(mo_ui - 1);
                    recompute = true;
                }

                const auto& mo = wavefunction->orbitals[mo_index];
                ImGui::Text("Internal index: %zu", mo_index);
                ImGui::Text("Energy: %.8f Ha", mo.energy_hartree);
                ImGui::Text("Occupation: %.3f", mo.occupation);
                ImGui::Text("Spin: %s", spin_name(mo.spin));
                ImGui::Text("Symmetry: %s", mo.symmetry.empty() ? "(none)" : mo.symmetry.c_str());

                if (ImGui::SliderFloat("Isovalue", &isovalue, 0.002f, 0.12f, "%.4f",
                                       ImGuiSliderFlags_Logarithmic)) {
                }

                constexpr int resolutions[] = {64, 128, 256, 512};
                int resolution_index = 1;
                for (int i = 0; i < 4; ++i) {
                    if (resolutions[i] == resolution) resolution_index = i;
                }
                if (ImGui::Combo("Grid", &resolution_index,
                                 "64^3\0" "128^3\0" "256^3\0" "512^3\0")) {
                    resolution = resolutions[resolution_index];
                    resize_and_recompute = true;
                    recompute = true;
                }

                if (ImGui::Button("Recompute grid")) recompute = true;
                ImGui::SameLine();
                if (ImGui::Button("Reset camera")) camera = {};

                ImGui::Text("Last CUDA kernel: %.3f ms", evaluator->last_kernel_ms());
                ImGui::TextWrapped(
                    "Mouse: left-drag to orbit, wheel to zoom. "
                    "Changing isovalue does not recompute the CUDA grid.");
            }

            ImGui::End();

            if (recompute) {
                try {
                    evaluate_now();
                } catch (const std::exception& e) {
                    status = std::string("CUDA ERROR: ") + e.what();
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
