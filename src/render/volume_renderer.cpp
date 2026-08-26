#include "cov/volume_renderer.hpp"

#include "cov/gl_api.hpp"

#ifdef _WIN32
#include <Windows.h>
#endif
#include <GL/gl.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef GL_TEXTURE_3D
#define GL_TEXTURE_3D 0x806F
#endif
#ifndef GL_TEXTURE_WRAP_R
#define GL_TEXTURE_WRAP_R 0x8072
#endif
#ifndef GL_CLAMP_TO_EDGE
#define GL_CLAMP_TO_EDGE 0x812F
#endif
#ifndef GL_R32F
#define GL_R32F 0x822E
#endif
#ifndef GL_RED
#define GL_RED 0x1903
#endif
#ifndef GL_TEXTURE0
#define GL_TEXTURE0 0x84C0
#endif
#ifndef GL_VERTEX_SHADER
#define GL_VERTEX_SHADER 0x8B31
#endif
#ifndef GL_FRAGMENT_SHADER
#define GL_FRAGMENT_SHADER 0x8B30
#endif
#ifndef GL_COMPILE_STATUS
#define GL_COMPILE_STATUS 0x8B81
#endif
#ifndef GL_LINK_STATUS
#define GL_LINK_STATUS 0x8B82
#endif
#ifndef GL_INFO_LOG_LENGTH
#define GL_INFO_LOG_LENGTH 0x8B84
#endif

namespace cov {
namespace {

constexpr float kPi = 3.14159265358979323846f;

struct Vec3 {
    float x, y, z;
};

Vec3 operator+(const Vec3 a, const Vec3 b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
Vec3 operator-(const Vec3 a, const Vec3 b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
Vec3 operator*(const Vec3 a, const float s) { return {a.x*s,a.y*s,a.z*s}; }

float dot(const Vec3 a, const Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

Vec3 cross(const Vec3 a, const Vec3 b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

Vec3 normalise(const Vec3 v) {
    const float n = std::sqrt(std::max(1.0e-20f, dot(v,v)));
    return v * (1.0f/n);
}

struct CameraBasis {
    Vec3 position;
    Vec3 forward;
    Vec3 right;
    Vec3 up;
};

CameraBasis camera_basis(const OrbitCamera& camera) {
    const Vec3 centre{0.5f, 0.5f, 0.5f};
    const float cp = std::cos(camera.pitch);
    const Vec3 offset{
        camera.distance * cp * std::cos(camera.yaw),
        camera.distance * std::sin(camera.pitch),
        camera.distance * cp * std::sin(camera.yaw)
    };
    CameraBasis b{};
    b.position = centre + offset;
    b.forward = normalise(centre - b.position);

    Vec3 world_up{0.0f, 1.0f, 0.0f};
    if (std::abs(dot(b.forward, world_up)) > 0.98f) {
        world_up = {0.0f, 0.0f, 1.0f};
    }
    b.right = normalise(cross(b.forward, world_up));
    b.up = normalise(cross(b.right, b.forward));
    return b;
}

GLuint compile_shader(const GLenum type, const char* source) {
    const GLuint shader = gl::CreateShader(type);
    gl::ShaderSource(shader, 1, &source, nullptr);
    gl::CompileShader(shader);

    GLint ok = 0;
    gl::GetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        GLint len = 0;
        gl::GetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
        std::vector<char> log(static_cast<std::size_t>(std::max(1, len)), '\0');
        gl::GetShaderInfoLog(shader, len, nullptr, log.data());
        const std::string message(log.data());
        gl::DeleteShader(shader);
        throw std::runtime_error("OpenGL shader compilation failed: " + message);
    }
    return shader;
}

GLuint make_program() {
    static constexpr const char* kVertex = R"GLSL(
#version 120
varying vec2 vUV;
void main() {
    gl_Position = gl_Vertex;
    vUV = gl_MultiTexCoord0.xy;
}
)GLSL";

    static constexpr const char* kFragment = R"GLSL(
#version 120
uniform sampler3D uVolume;
uniform float uIso;
uniform vec3 uCameraPos;
uniform vec3 uCameraForward;
uniform vec3 uCameraRight;
uniform vec3 uCameraUp;
uniform float uAspect;
uniform float uTanHalfFov;
uniform vec3 uTexel;
varying vec2 vUV;

bool rayBox(vec3 ro, vec3 rd, out float tNear, out float tFar) {
    vec3 inv = 1.0 / rd;
    vec3 t0 = (vec3(0.0) - ro) * inv;
    vec3 t1 = (vec3(1.0) - ro) * inv;
    vec3 lo = min(t0, t1);
    vec3 hi = max(t0, t1);
    tNear = max(max(lo.x, lo.y), lo.z);
    tFar = min(min(hi.x, hi.y), hi.z);
    return tFar >= max(tNear, 0.0);
}

float field(vec3 p) {
    return texture3D(uVolume, clamp(p, vec3(0.0), vec3(1.0))).r;
}

vec3 gradient(vec3 p) {
    float gx = field(p + vec3(uTexel.x,0.0,0.0)) - field(p - vec3(uTexel.x,0.0,0.0));
    float gy = field(p + vec3(0.0,uTexel.y,0.0)) - field(p - vec3(0.0,uTexel.y,0.0));
    float gz = field(p + vec3(0.0,0.0,uTexel.z)) - field(p - vec3(0.0,0.0,uTexel.z));
    return normalize(vec3(gx,gy,gz));
}

void main() {
    vec2 ndc = vUV * 2.0 - 1.0;
    vec3 ro = uCameraPos;
    vec3 rd = normalize(
        uCameraForward +
        ndc.x * uAspect * uTanHalfFov * uCameraRight +
        ndc.y * uTanHalfFov * uCameraUp);

    float t0, t1;
    if (!rayBox(ro, rd, t0, t1)) {
        gl_FragColor = vec4(0.025, 0.028, 0.035, 1.0);
        return;
    }

    float startT = max(t0, 0.0);
    float travel = max(0.0001, t1 - startT);
    float dt = travel / 512.0;
    float t = startT;
    vec3 p = ro + rd * t;
    float prev = field(p);
    float prevF = abs(prev) - uIso;

    for (int i = 0; i < 512; ++i) {
        t += dt;
        if (t > t1) break;
        p = ro + rd * t;
        float cur = field(p);
        float curF = abs(cur) - uIso;

        if (prevF < 0.0 && curF >= 0.0) {
            vec3 hit = p;
            vec3 n = gradient(hit);
            vec3 light = normalize(vec3(0.55, 0.8, 0.35));
            float diffuse = 0.25 + 0.75 * abs(dot(n, light));
            vec3 positive = vec3(0.92, 0.20, 0.17);
            vec3 negative = vec3(0.16, 0.38, 0.95);
            vec3 base = cur >= 0.0 ? positive : negative;
            gl_FragColor = vec4(base * diffuse, 1.0);
            return;
        }
        prev = cur;
        prevF = curF;
    }

    gl_FragColor = vec4(0.025, 0.028, 0.035, 1.0);
}
)GLSL";

    const GLuint vs = compile_shader(GL_VERTEX_SHADER, kVertex);
    const GLuint fs = compile_shader(GL_FRAGMENT_SHADER, kFragment);

    const GLuint program = gl::CreateProgram();
    gl::AttachShader(program, vs);
    gl::AttachShader(program, fs);
    gl::LinkProgram(program);

    GLint ok = 0;
    gl::GetProgramiv(program, GL_LINK_STATUS, &ok);
    if (!ok) {
        GLint len = 0;
        gl::GetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
        std::vector<char> log(static_cast<std::size_t>(std::max(1, len)), '\0');
        gl::GetProgramInfoLog(program, len, nullptr, log.data());
        const std::string message(log.data());
        gl::DeleteShader(vs);
        gl::DeleteShader(fs);
        gl::DeleteProgram(program);
        throw std::runtime_error("OpenGL program link failed: " + message);
    }

    gl::DeleteShader(vs);
    gl::DeleteShader(fs);
    return program;
}

float covalent_radius_angstrom(const int z) {
    switch (z) {
        case 1: return 0.31f;
        case 5: return 0.84f;
        case 6: return 0.76f;
        case 7: return 0.71f;
        case 8: return 0.66f;
        case 9: return 0.57f;
        case 14: return 1.11f;
        case 15: return 1.07f;
        case 16: return 1.05f;
        case 17: return 1.02f;
        case 35: return 1.20f;
        case 53: return 1.39f;
        default: return 0.85f;
    }
}

Vec3 atom_colour(const int z) {
    switch (z) {
        case 1: return {0.95f,0.95f,0.95f};
        case 6: return {0.30f,0.30f,0.32f};
        case 7: return {0.20f,0.35f,0.95f};
        case 8: return {0.92f,0.12f,0.12f};
        case 9: return {0.20f,0.85f,0.30f};
        case 15: return {0.95f,0.55f,0.10f};
        case 16: return {0.95f,0.85f,0.12f};
        case 17: return {0.15f,0.80f,0.22f};
        default: return {0.65f,0.65f,0.70f};
    }
}

Vec3 to_texture(const Atom& atom, const GridBox& box) {
    const float dx = std::max(1.0e-8f, box.max_x - box.min_x);
    const float dy = std::max(1.0e-8f, box.max_y - box.min_y);
    const float dz = std::max(1.0e-8f, box.max_z - box.min_z);
    return {
        (static_cast<float>(atom.x) - box.min_x) / dx,
        (static_cast<float>(atom.y) - box.min_y) / dy,
        (static_cast<float>(atom.z) - box.min_z) / dz
    };
}

bool project(const Vec3 p,
             const CameraBasis& b,
             const float aspect,
             const float tan_half_fov,
             float& x_ndc,
             float& y_ndc) {
    const Vec3 q = p - b.position;
    const float z = dot(q, b.forward);
    if (z <= 0.02f) return false;
    x_ndc = dot(q, b.right) / (z * tan_half_fov * aspect);
    y_ndc = dot(q, b.up) / (z * tan_half_fov);
    return true;
}

} // namespace

VolumeRenderer::VolumeRenderer() {
    if (!gl::load()) {
        throw std::runtime_error(
            "Required OpenGL 2.1 shader/3D-texture entry points are unavailable");
    }
    program_ = make_program();

    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_3D, texture_);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_3D, 0);
}

VolumeRenderer::~VolumeRenderer() {
    if (program_) gl::DeleteProgram(program_);
    if (texture_) glDeleteTextures(1, &texture_);
}

void VolumeRenderer::resize_volume(const int nx, const int ny, const int nz) {
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Invalid 3D texture dimensions");
    }
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;

    glBindTexture(GL_TEXTURE_3D, texture_);
    gl::TexImage3D(GL_TEXTURE_3D, 0, GL_R32F,
                   nx_, ny_, nz_, 0, GL_RED, GL_FLOAT, nullptr);
    glBindTexture(GL_TEXTURE_3D, 0);
}

void VolumeRenderer::render_volume(const int framebuffer_width,
                                   const int framebuffer_height,
                                   const float isovalue,
                                   const OrbitCamera& camera) {
    if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0) return;

    const CameraBasis b = camera_basis(camera);
    const float aspect = framebuffer_height > 0
                             ? static_cast<float>(framebuffer_width) /
                                   static_cast<float>(framebuffer_height)
                             : 1.0f;
    const float tan_half_fov =
        std::tan(camera.fov_degrees * kPi / 360.0f);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    gl::UseProgram(program_);

    gl::Uniform1i(gl::GetUniformLocation(program_, "uVolume"), 0);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uIso"), isovalue);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uCameraPos"),
                  b.position.x, b.position.y, b.position.z);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uCameraForward"),
                  b.forward.x, b.forward.y, b.forward.z);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uCameraRight"),
                  b.right.x, b.right.y, b.right.z);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uCameraUp"),
                  b.up.x, b.up.y, b.up.z);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uAspect"), aspect);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uTanHalfFov"), tan_half_fov);
    gl::Uniform3f(gl::GetUniformLocation(program_, "uTexel"),
                  1.0f / static_cast<float>(nx_),
                  1.0f / static_cast<float>(ny_),
                  1.0f / static_cast<float>(nz_));

    gl::ActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, texture_);

    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);
    glTexCoord2f(1.0f, 0.0f); glVertex2f( 1.0f, -1.0f);
    glTexCoord2f(1.0f, 1.0f); glVertex2f( 1.0f,  1.0f);
    glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f,  1.0f);
    glEnd();

    glBindTexture(GL_TEXTURE_3D, 0);
    gl::UseProgram(0);
}

void VolumeRenderer::render_geometry(const Wavefunction& wavefunction,
                                     const GridBox& box,
                                     const int framebuffer_width,
                                     const int framebuffer_height,
                                     const OrbitCamera& camera) {
    if (wavefunction.atoms.empty()) return;

    const CameraBasis b = camera_basis(camera);
    const float aspect = framebuffer_height > 0
                             ? static_cast<float>(framebuffer_width) /
                                   static_cast<float>(framebuffer_height)
                             : 1.0f;
    const float tan_half_fov =
        std::tan(camera.fov_degrees * kPi / 360.0f);

    std::vector<Vec3> points;
    points.reserve(wavefunction.atoms.size());
    for (const Atom& atom : wavefunction.atoms) {
        points.push_back(to_texture(atom, box));
    }

    glDisable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glLineWidth(2.0f);
    glColor3f(0.72f, 0.72f, 0.76f);
    glBegin(GL_LINES);
    for (std::size_t i = 0; i < wavefunction.atoms.size(); ++i) {
        for (std::size_t j = i + 1; j < wavefunction.atoms.size(); ++j) {
            const Atom& a = wavefunction.atoms[i];
            const Atom& c = wavefunction.atoms[j];
            const double dx = a.x - c.x;
            const double dy = a.y - c.y;
            const double dz = a.z - c.z;
            const double distance_bohr = std::sqrt(dx*dx + dy*dy + dz*dz);
            const double cutoff_bohr =
                1.25 * (covalent_radius_angstrom(a.atomic_number) +
                        covalent_radius_angstrom(c.atomic_number)) *
                kAngstromToBohr;
            if (distance_bohr > cutoff_bohr) continue;

            float x0,y0,x1,y1;
            if (project(points[i], b, aspect, tan_half_fov, x0, y0) &&
                project(points[j], b, aspect, tan_half_fov, x1, y1)) {
                glVertex2f(x0,y0);
                glVertex2f(x1,y1);
            }
        }
    }
    glEnd();

    glPointSize(9.0f);
    glBegin(GL_POINTS);
    for (std::size_t i = 0; i < wavefunction.atoms.size(); ++i) {
        float x,y;
        if (!project(points[i], b, aspect, tan_half_fov, x, y)) continue;
        const Vec3 colour = atom_colour(wavefunction.atoms[i].atomic_number);
        glColor3f(colour.x, colour.y, colour.z);
        glVertex2f(x,y);
    }
    glEnd();
    glColor3f(1.0f,1.0f,1.0f);

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

} // namespace cov
