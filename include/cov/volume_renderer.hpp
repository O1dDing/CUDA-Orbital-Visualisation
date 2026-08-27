#pragma once

#include "cov/model.hpp"
#include "cov/molecule_style.hpp"

#include <cstdint>

namespace cov {

struct OrbitCamera {
    float yaw = 0.65f;
    float pitch = 0.35f;
    float distance = 2.2f;
    float fov_degrees = 42.0f;
};

class VolumeRenderer {
public:
    VolumeRenderer();
    ~VolumeRenderer();

    VolumeRenderer(const VolumeRenderer&) = delete;
    VolumeRenderer& operator=(const VolumeRenderer&) = delete;

    void resize_volume(int nx, int ny, int nz);

    [[nodiscard]] unsigned int volume_texture() const noexcept { return texture_; }
    [[nodiscard]] int nx() const noexcept { return nx_; }
    [[nodiscard]] int ny() const noexcept { return ny_; }
    [[nodiscard]] int nz() const noexcept { return nz_; }

    void render_volume(int framebuffer_width,
                       int framebuffer_height,
                       float isovalue,
                       const OrbitCamera& camera,
                       float opacity = 1.0f);

    void render_geometry(const Wavefunction& wavefunction,
                         const GridBox& box,
                         int framebuffer_width,
                         int framebuffer_height,
                         const OrbitCamera& camera,
                         const MoleculeRenderSettings& settings = {});

private:
    unsigned int texture_ = 0;
    unsigned int program_ = 0;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
};

} // namespace cov
