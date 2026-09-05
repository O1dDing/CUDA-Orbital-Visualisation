#pragma once
#include "cov/orbital_ui.hpp"
#include "cov/volume_renderer.hpp"
#include <imgui.h>
#include <filesystem>
#include <string>

namespace cov::validation {
#ifdef COV_ENABLE_VALIDATION
bool configure(int argc, char** argv);
bool active();
bool done();
int result();
void begin_frame(OrbitCamera&, MoleculeRenderSettings&, float&, int&, bool&);
void input_frame();
void evaluated(std::size_t mo, const char* reason, float milliseconds);
void ui_frame(std::size_t drawn, std::size_t requested);
void after_scene(const VolumeRenderer&, const GridBox&, std::size_t mo);
void end_frame(int width, int height, std::size_t applied,
               const ui::OrbitalUIState&, const Wavefunction*);
void item(const std::string& id);
void hit(const std::string& id, ImVec2 lo, ImVec2 hi);
void anchor(const std::string& id);
void record(const std::string& kind, const std::string& json);
void field(const std::string& label, const std::string& value);
std::string quote(const std::string& value);
std::filesystem::path export_base(const std::filesystem::path& original);
#else
inline bool configure(int, char**) { return false; }
inline bool active() { return false; }
inline bool done() { return false; }
inline int result() { return 0; }
inline void begin_frame(OrbitCamera&, MoleculeRenderSettings&, float&, int&, bool&) {}
inline void input_frame() {}
inline void evaluated(std::size_t, const char*, float) {}
inline void ui_frame(std::size_t, std::size_t) {}
inline void after_scene(const VolumeRenderer&, const GridBox&, std::size_t) {}
inline void end_frame(int, int, std::size_t, const ui::OrbitalUIState&, const Wavefunction*) {}
inline void item(const std::string&) {}
inline void hit(const std::string&, ImVec2, ImVec2) {}
inline void anchor(const std::string&) {}
inline void record(const std::string&, const std::string&) {}
inline void field(const std::string&, const std::string&) {}
inline std::string quote(const std::string&) { return {}; }
inline std::filesystem::path export_base(const std::filesystem::path& p) { return p; }
#endif
}
