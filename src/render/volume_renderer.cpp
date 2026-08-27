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
constexpr float kNearPlane = 0.03f;
constexpr float kFarPlane = 8.0f;

struct Vec3 {
    float x, y, z;
};

Vec3 operator+(const Vec3 a, const Vec3 b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
Vec3 operator-(const Vec3 a, const Vec3 b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
Vec3 operator*(const Vec3 a, const float s) { return {a.x*s,a.y*s,a.z*s}; }
Vec3 operator/(const Vec3 a, const float s) { return {a.x/s,a.y/s,a.z/s}; }

float dot(const Vec3 a, const Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

Vec3 cross(const Vec3 a, const Vec3 b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

float length(const Vec3 v) { return std::sqrt(std::max(0.0f, dot(v, v))); }

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
    if (std::abs(dot(b.forward, world_up)) > 0.98f) world_up = {0.0f, 0.0f, 1.0f};
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

    // The ray marcher now writes the real isosurface depth into the shared
    // OpenGL depth buffer. Geometry rendered afterwards is therefore hidden by
    // an orbital lobe when the lobe is physically closer to the camera, while
    // atoms/bonds that genuinely sit in front remain visible.
    static constexpr const char* kFragment = R"GLSL(
#version 120
uniform sampler3D uVolume;
uniform float uIso;
uniform float uOpacity;
uniform vec3 uCameraPos;
uniform vec3 uCameraForward;
uniform vec3 uCameraRight;
uniform vec3 uCameraUp;
uniform float uAspect;
uniform float uTanHalfFov;
uniform vec3 uTexel;
varying vec2 vUV;

const vec3 kBackground = vec3(0.025, 0.031, 0.043);
const float kNear = 0.03;
const float kFar = 8.0;

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

float depthFromViewDistance(float viewDepth) {
    float z = max(viewDepth, kNear + 0.0001);
    float ndc = (kFar + kNear) / (kFar - kNear)
              - (2.0 * kFar * kNear) / ((kFar - kNear) * z);
    return clamp(0.5 * ndc + 0.5, 0.0, 1.0);
}

void main() {
    vec2 ndc = vUV * 2.0 - 1.0;
    vec3 ro = uCameraPos;
    vec3 rd = normalize(
        uCameraForward +
        ndc.x * uAspect * uTanHalfFov * uCameraRight +
        ndc.y * uTanHalfFov * uCameraUp);

    float t0, t1;
    if (!rayBox(ro, rd, t0, t1)) discard;

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
            vec3 shaded = base * diffuse;
            float viewDepth = dot(hit - uCameraPos, uCameraForward);
            gl_FragDepth = depthFromViewDistance(viewDepth);
            gl_FragColor = vec4(mix(kBackground, shaded, clamp(uOpacity, 0.08, 1.0)), 1.0);
            return;
        }
        prev = cur;
        prevF = curF;
    }

    discard;
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

Vec3 atom_colour(const int z) {
    switch (z) {
        case 1: return {0.95f,0.95f,0.95f};
        case 5: return {0.93f,0.62f,0.62f};
        case 6: return {0.30f,0.31f,0.34f};
        case 7: return {0.20f,0.35f,0.95f};
        case 8: return {0.92f,0.12f,0.12f};
        case 9: return {0.20f,0.85f,0.30f};
        case 14: return {0.62f,0.64f,0.67f};
        case 15: return {0.95f,0.55f,0.10f};
        case 16: return {0.95f,0.85f,0.12f};
        case 17: return {0.15f,0.80f,0.22f};
        case 35: return {0.65f,0.16f,0.16f};
        case 53: return {0.45f,0.18f,0.62f};
        default: return {0.65f,0.65f,0.70f};
    }
}

Vec3 to_texture_point(const double x, const double y, const double z, const GridBox& box) {
    const float dx = std::max(1.0e-8f, box.max_x - box.min_x);
    const float dy = std::max(1.0e-8f, box.max_y - box.min_y);
    const float dz = std::max(1.0e-8f, box.max_z - box.min_z);
    return {
        (static_cast<float>(x) - box.min_x) / dx,
        (static_cast<float>(y) - box.min_y) / dy,
        (static_cast<float>(z) - box.min_z) / dz
    };
}

Vec3 to_texture(const Atom& atom, const GridBox& box) {
    return to_texture_point(atom.x, atom.y, atom.z, box);
}

Vec3 shaded_colour(const Vec3 base,
                   const Vec3 normal,
                   const Vec3 point,
                   const CameraBasis& camera) {
    const Vec3 light = normalise(camera.right * -0.35f + camera.up * 0.58f + camera.forward * -0.74f);
    const Vec3 view = normalise(camera.position - point);
    const Vec3 half_vector = normalise(light + view);
    const float diffuse = std::max(0.0f, dot(normal, light));
    const float specular = std::pow(std::max(0.0f, dot(normal, half_vector)), 30.0f);
    const float illumination = 0.28f + 0.72f * diffuse;
    return {
        std::clamp(base.x * illumination + 0.34f * specular, 0.0f, 1.0f),
        std::clamp(base.y * illumination + 0.34f * specular, 0.0f, 1.0f),
        std::clamp(base.z * illumination + 0.34f * specular, 0.0f, 1.0f)
    };
}

void load_3d_camera(const CameraBasis& b,
                    const float aspect,
                    const float tan_half_fov) {
    const double top = static_cast<double>(kNearPlane) * static_cast<double>(tan_half_fov);
    const double right = top * static_cast<double>(std::max(0.01f, aspect));

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-right, right, -top, top, kNearPlane, kFarPlane);

    const GLfloat view[16] = {
        b.right.x, b.up.x, -b.forward.x, 0.0f,
        b.right.y, b.up.y, -b.forward.y, 0.0f,
        b.right.z, b.up.z, -b.forward.z, 0.0f,
        -dot(b.right, b.position), -dot(b.up, b.position), dot(b.forward, b.position), 1.0f
    };
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(view);
}

void draw_sphere(const Vec3 centre,
                 const float radius,
                 const Vec3 base_colour,
                 const float alpha,
                 const CameraBasis& camera) {
    constexpr int latitude_segments = 18;
    constexpr int longitude_segments = 28;

    for (int lat = 0; lat < latitude_segments; ++lat) {
        const float phi0 = -0.5f * kPi + kPi * static_cast<float>(lat) /
                                             static_cast<float>(latitude_segments);
        const float phi1 = -0.5f * kPi + kPi * static_cast<float>(lat + 1) /
                                             static_cast<float>(latitude_segments);
        const float y0 = std::sin(phi0);
        const float y1 = std::sin(phi1);
        const float r0 = std::cos(phi0);
        const float r1 = std::cos(phi1);

        glBegin(GL_TRIANGLE_STRIP);
        for (int lon = 0; lon <= longitude_segments; ++lon) {
            const float theta = 2.0f * kPi * static_cast<float>(lon) /
                                                static_cast<float>(longitude_segments);
            const float ct = std::cos(theta);
            const float st = std::sin(theta);

            const Vec3 n0{r0 * ct, y0, r0 * st};
            const Vec3 p0 = centre + n0 * radius;
            const Vec3 c0 = shaded_colour(base_colour, n0, p0, camera);
            glColor4f(c0.x, c0.y, c0.z, alpha);
            glVertex3f(p0.x, p0.y, p0.z);

            const Vec3 n1{r1 * ct, y1, r1 * st};
            const Vec3 p1 = centre + n1 * radius;
            const Vec3 c1 = shaded_colour(base_colour, n1, p1, camera);
            glColor4f(c1.x, c1.y, c1.z, alpha);
            glVertex3f(p1.x, p1.y, p1.z);
        }
        glEnd();
    }
}

void draw_cylinder(const Vec3 a,
                   const Vec3 b,
                   const float radius,
                   const Vec3 base_colour,
                   const float alpha,
                   const CameraBasis& camera) {
    const Vec3 axis_vector = b - a;
    const float axis_length = std::sqrt(std::max(1.0e-20f, dot(axis_vector, axis_vector)));
    if (axis_length <= 1.0e-6f) return;
    const Vec3 axis = axis_vector * (1.0f / axis_length);
    const Vec3 helper = std::abs(axis.y) < 0.90f ? Vec3{0.0f, 1.0f, 0.0f}
                                                  : Vec3{1.0f, 0.0f, 0.0f};
    const Vec3 u = normalise(cross(axis, helper));
    const Vec3 v = normalise(cross(axis, u));
    constexpr int slices = 18;

    glBegin(GL_TRIANGLE_STRIP);
    for (int i = 0; i <= slices; ++i) {
        const float theta = 2.0f * kPi * static_cast<float>(i) / static_cast<float>(slices);
        const Vec3 normal = u * std::cos(theta) + v * std::sin(theta);
        const Vec3 pa = a + normal * radius;
        const Vec3 pb = b + normal * radius;
        const Vec3 ca = shaded_colour(base_colour, normal, pa, camera);
        const Vec3 cb = shaded_colour(base_colour, normal, pb, camera);
        glColor4f(ca.x, ca.y, ca.z, alpha);
        glVertex3f(pa.x, pa.y, pa.z);
        glColor4f(cb.x, cb.y, cb.z, alpha);
        glVertex3f(pb.x, pb.y, pb.z);
    }
    glEnd();

    const Vec3 cap_a_colour = shaded_colour(base_colour, axis * -1.0f, a, camera);
    glColor4f(cap_a_colour.x, cap_a_colour.y, cap_a_colour.z, alpha);
    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(a.x, a.y, a.z);
    for (int i = 0; i <= slices; ++i) {
        const float theta = -2.0f * kPi * static_cast<float>(i) / static_cast<float>(slices);
        const Vec3 rim = a + (u * std::cos(theta) + v * std::sin(theta)) * radius;
        glVertex3f(rim.x, rim.y, rim.z);
    }
    glEnd();

    const Vec3 cap_b_colour = shaded_colour(base_colour, axis, b, camera);
    glColor4f(cap_b_colour.x, cap_b_colour.y, cap_b_colour.z, alpha);
    glBegin(GL_TRIANGLE_FAN);
    glVertex3f(b.x, b.y, b.z);
    for (int i = 0; i <= slices; ++i) {
        const float theta = 2.0f * kPi * static_cast<float>(i) / static_cast<float>(slices);
        const Vec3 rim = b + (u * std::cos(theta) + v * std::sin(theta)) * radius;
        glVertex3f(rim.x, rim.y, rim.z);
    }
    glEnd();
}

void draw_dashed_cylinder(const Vec3 a,
                          const Vec3 b,
                          const float radius,
                          const Vec3 colour,
                          const float alpha,
                          const CameraBasis& camera) {
    constexpr int dash_count = 8;
    constexpr float dash_fraction = 0.55f;
    const Vec3 delta = b - a;
    for (int i = 0; i < dash_count; ++i) {
        const float t0 = static_cast<float>(i) / static_cast<float>(dash_count);
        const float t1 = std::min(1.0f, t0 + dash_fraction / static_cast<float>(dash_count));
        draw_cylinder(a + delta * t0, a + delta * t1,
                      radius, colour, alpha, camera);
    }
}

Vec3 heavy_atom_centroid(const Wavefunction& wavefunction,
                         const std::vector<Vec3>& points) {
    Vec3 sum{0.0f, 0.0f, 0.0f};
    std::size_t count = 0;
    for (std::size_t i = 0; i < wavefunction.atoms.size(); ++i) {
        if (wavefunction.atoms[i].atomic_number == 1) continue;
        sum = sum + points[i];
        ++count;
    }
    if (count == 0) {
        for (const Vec3 p : points) sum = sum + p;
        count = points.size();
    }
    return count ? sum / static_cast<float>(count) : Vec3{0.5f, 0.5f, 0.5f};
}

Vec3 inward_delocalisation_offset(const Vec3 a,
                                  const Vec3 b,
                                  const Vec3 centre,
                                  const float distance) {
    const Vec3 axis = normalise(b - a);
    const Vec3 midpoint = (a + b) * 0.5f;
    Vec3 inward = centre - midpoint;
    inward = inward - axis * dot(inward, axis);
    if (length(inward) <= 1.0e-6f) return {0.0f, 0.0f, 0.0f};
    return normalise(inward) * distance;
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
                                   const OrbitCamera& camera,
                                   const float opacity) {
    if (nx_ <= 0 || ny_ <= 0 || nz_ <= 0) return;

    const CameraBasis b = camera_basis(camera);
    const float aspect = framebuffer_height > 0
                             ? static_cast<float>(framebuffer_width) /
                                   static_cast<float>(framebuffer_height)
                             : 1.0f;
    const float tan_half_fov = std::tan(camera.fov_degrees * kPi / 360.0f);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    gl::UseProgram(program_);

    gl::Uniform1i(gl::GetUniformLocation(program_, "uVolume"), 0);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uIso"), isovalue);
    gl::Uniform1f(gl::GetUniformLocation(program_, "uOpacity"),
                  std::clamp(opacity, 0.08f, 1.0f));
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

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0f, -1.0f);
    glTexCoord2f(1.0f, 0.0f); glVertex2f( 1.0f, -1.0f);
    glTexCoord2f(1.0f, 1.0f); glVertex2f( 1.0f,  1.0f);
    glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0f,  1.0f);
    glEnd();

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

    glBindTexture(GL_TEXTURE_3D, 0);
    gl::UseProgram(0);
}

void VolumeRenderer::render_geometry(const Wavefunction& wavefunction,
                                     const GridBox& box,
                                     const int framebuffer_width,
                                     const int framebuffer_height,
                                     const OrbitCamera& camera,
                                     const MoleculeRenderSettings& settings) {
    if (wavefunction.atoms.empty()) return;

    const CameraBasis b = camera_basis(camera);
    const float aspect = framebuffer_height > 0
                             ? static_cast<float>(framebuffer_width) /
                                   static_cast<float>(framebuffer_height)
                             : 1.0f;
    const float tan_half_fov = std::tan(camera.fov_degrees * kPi / 360.0f);

    std::vector<Vec3> points;
    points.reserve(wavefunction.atoms.size());
    for (const Atom& atom : wavefunction.atoms) points.push_back(to_texture(atom, box));
    const auto bonds = analyse_bonds(wavefunction);
    const Vec3 molecular_centre = heavy_atom_centroid(wavefunction, points);

    // Do not clear the depth buffer here: the orbital raymarch has already
    // written its true surface depth. Geometry behind a lobe must therefore be
    // occluded by that lobe instead of being painted over it as an overlay.
    gl::UseProgram(0);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glShadeModel(GL_SMOOTH);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    load_3d_camera(b, aspect, tan_half_fov);

    const float opacity = std::clamp(settings.molecule_opacity, 0.05f, 1.0f);
    // The accepted pre.2 screenshot used Atom=1.25 and Bond=1.15. Those
    // proportions are now calibrated into the base radii so the user-facing
    // reference value 1.0 reproduces the same size.
    const float bond_radius = (settings.style == MoleculeStyle::MediumBallAndStick ? 0.00540f : 0.00415f) *
                              std::clamp(settings.bond_scale, 0.35f, 3.0f);
    const Vec3 bond_colour{0.73f, 0.76f, 0.81f};
    const Vec3 delocalised_colour{0.46f, 0.69f, 0.98f};

    for (const auto& bond : bonds) {
        const Atom& atom_a = wavefunction.atoms[bond.atom_a];
        const Atom& atom_b = wavefunction.atoms[bond.atom_b];
        if (!settings.show_hydrogens &&
            (atom_a.atomic_number == 1 || atom_b.atomic_number == 1)) {
            continue;
        }

        const Vec3 a = points[bond.atom_a];
        const Vec3 c = points[bond.atom_b];

        // Standard Molden does not carry authoritative bond orders. Therefore
        // the base structural graph remains a chemically neutral single bond;
        // no double/triple bond is fabricated from distance alone.
        draw_cylinder(a, c, bond_radius, bond_colour, opacity * 0.94f, b);

        if (settings.style == MoleculeStyle::StickDelocalisation && bond.delocalised) {
            // Delocalisation is a visual heuristic: retain the full 3D single
            // bond and add a thinner 3D dashed tube displaced toward the heavy-
            // atom centre. This reads like an inner aromatic/delocalised line
            // rather than replacing the molecular skeleton with floating dashes.
            const Vec3 inward = inward_delocalisation_offset(
                a, c, molecular_centre, bond_radius * 2.35f);
            draw_dashed_cylinder(a + inward, c + inward,
                                 bond_radius * 0.58f,
                                 delocalised_colour, opacity * 0.96f, b);
        }
    }

    if (settings.style == MoleculeStyle::MediumBallAndStick) {
        for (std::size_t i = 0; i < wavefunction.atoms.size(); ++i) {
            const Atom& atom = wavefunction.atoms[i];
            if (!settings.show_hydrogens && atom.atomic_number == 1) continue;

            const float covalent = static_cast<float>(covalent_radius_angstrom(atom.atomic_number));
            const float sphere_radius = std::clamp(
                (0.0275f + 0.0200f * covalent) * settings.atom_scale,
                0.030f, 0.108f);
            draw_sphere(points[i], sphere_radius, atom_colour(atom.atomic_number), opacity, b);
        }
    }

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);

    glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
}

} // namespace cov
