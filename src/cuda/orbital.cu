#include "cov/cuda_orbital.hpp"

#include <cuda_gl_interop.h>
#include <cuda_runtime.h>

#ifdef _WIN32
#include <Windows.h>
#endif
#include <GL/gl.h>

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef GL_TEXTURE_3D
#define GL_TEXTURE_3D 0x806F
#endif

namespace cov {
namespace {

constexpr float kPi = 3.14159265358979323846f;

struct GpuPrimitive {
    float exponent;
    float coefficient;
};

struct GpuShell {
    float cx;
    float cy;
    float cz;
    std::uint32_t primitive_offset;
    std::uint32_t primitive_count;
    std::uint32_t basis_offset;
    std::uint8_t l;
    std::uint8_t pure;
    std::uint16_t pad;
};

void cuda_check(const cudaError_t status, const char* what) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(status));
    }
}

__device__ __forceinline__ float powi(float x, int n) {
    float out = 1.0f;
    for (int i = 0; i < n; ++i) out *= x;
    return out;
}

__device__ __forceinline__ int double_factorial_odd(int n) {
    if (n <= 0) return 1;
    int v = 1;
    for (int k = n; k > 1; k -= 2) v *= k;
    return v;
}

__device__ __forceinline__ float cartesian_primitive_norm(
    const float alpha, const int ax, const int ay, const int az) {
    const int l = ax + ay + az;
    const float base = powf(2.0f * alpha / kPi, 0.75f);
    const float angular = powf(4.0f * alpha, 0.5f * static_cast<float>(l));
    const float denom = sqrtf(static_cast<float>(
        double_factorial_odd(2 * ax - 1) *
        double_factorial_odd(2 * ay - 1) *
        double_factorial_odd(2 * az - 1)));
    return base * angular / denom;
}

__device__ __forceinline__ float gamma_l_plus_three_halves(const int l) {
    const float sqrt_pi = sqrtf(kPi);
    switch (l) {
        case 0: return 0.5f * sqrt_pi;
        case 1: return 0.75f * sqrt_pi;
        case 2: return 1.875f * sqrt_pi;
        case 3: return 6.5625f * sqrt_pi;
        default: return 29.53125f * sqrt_pi;
    }
}

__device__ __forceinline__ float spherical_primitive_norm(const float alpha, const int l) {
    const float power = powf(2.0f * alpha, static_cast<float>(l) + 1.5f);
    return sqrtf(2.0f * power / gamma_l_plus_three_halves(l));
}

__device__ __forceinline__ int factorial_small(int n) {
    int r = 1;
    for (int i = 2; i <= n; ++i) r *= i;
    return r;
}

__device__ float associated_legendre(const int l, const int m, const float x) {
    float pmm = 1.0f;
    if (m > 0) {
        const float somx2 = sqrtf(fmaxf(0.0f, 1.0f - x * x));
        float fact = 1.0f;
        for (int i = 1; i <= m; ++i) {
            pmm *= -fact * somx2;
            fact += 2.0f;
        }
    }
    if (l == m) return pmm;

    float pmmp1 = x * static_cast<float>(2 * m + 1) * pmm;
    if (l == m + 1) return pmmp1;

    float pll = 0.0f;
    for (int ll = m + 2; ll <= l; ++ll) {
        pll = (static_cast<float>(2 * ll - 1) * x * pmmp1 -
               static_cast<float>(ll + m - 1) * pmm) /
              static_cast<float>(ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

__device__ float real_solid_harmonic(
    const int l, const int index, const float x, const float y, const float z) {
    if (l == 0) {
        return 0.28209479177387814f;
    }

    const float r2 = x * x + y * y + z * z;
    if (r2 < 1.0e-20f) return 0.0f;
    const float r = sqrtf(r2);
    const float costheta = fminf(1.0f, fmaxf(-1.0f, z / r));
    const float phi = atan2f(y, x);

    int m = 0;
    bool sine = false;
    if (index > 0) {
        m = (index + 1) / 2;
        sine = (index % 2 == 0);
    }

    const float p = associated_legendre(l, m, costheta);
    const float n = sqrtf(
        (static_cast<float>(2 * l + 1) / (4.0f * kPi)) *
        (static_cast<float>(factorial_small(l - m)) /
         static_cast<float>(factorial_small(l + m))));

    float angular = n * p;
    if (m > 0) {
        angular *= sqrtf(2.0f);
        angular *= sine ? sinf(static_cast<float>(m) * phi)
                        : cosf(static_cast<float>(m) * phi);
    }
    return powi(r, l) * angular;
}

__device__ __forceinline__ void cartesian_exponents(
    const int l, const int index, int& ax, int& ay, int& az) {
    ax = ay = az = 0;
    if (l == 0) return;

    if (l == 1) {
        const int table[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
        ax = table[index][0]; ay = table[index][1]; az = table[index][2];
        return;
    }
    if (l == 2) {
        const int table[6][3] = {
            {2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}
        };
        ax = table[index][0]; ay = table[index][1]; az = table[index][2];
        return;
    }
    if (l == 3) {
        const int table[10][3] = {
            {3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},
            {2,0,1},{1,0,2},{0,1,2},{0,2,1},{1,1,1}
        };
        ax = table[index][0]; ay = table[index][1]; az = table[index][2];
        return;
    }

    const int table[15][3] = {
        {4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},
        {1,3,0},{0,3,1},{1,0,3},{0,1,3},{2,2,0},
        {2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}
    };
    ax = table[index][0]; ay = table[index][1]; az = table[index][2];
}

__device__ __forceinline__ int shell_basis_count_device(const GpuShell& shell) {
    return shell.pure ? (2 * static_cast<int>(shell.l) + 1)
                      : ((static_cast<int>(shell.l) + 1) *
                         (static_cast<int>(shell.l) + 2) / 2);
}

__device__ float evaluate_shell_component(
    const GpuShell& shell,
    const GpuPrimitive* primitives,
    const int component,
    const float dx,
    const float dy,
    const float dz,
    const float r2) {

    if (shell.pure) {
        float radial = 0.0f;
        for (std::uint32_t p = 0; p < shell.primitive_count; ++p) {
            const GpuPrimitive primitive = primitives[shell.primitive_offset + p];
            const float norm = spherical_primitive_norm(
                primitive.exponent, static_cast<int>(shell.l));
            radial += primitive.coefficient * norm *
                      __expf(-primitive.exponent * r2);
        }
        return radial * real_solid_harmonic(
            static_cast<int>(shell.l), component, dx, dy, dz);
    }

    int ax, ay, az;
    cartesian_exponents(static_cast<int>(shell.l), component, ax, ay, az);
    const float monomial = powi(dx, ax) * powi(dy, ay) * powi(dz, az);
    float contracted = 0.0f;
    for (std::uint32_t p = 0; p < shell.primitive_count; ++p) {
        const GpuPrimitive primitive = primitives[shell.primitive_offset + p];
        const float norm = cartesian_primitive_norm(primitive.exponent, ax, ay, az);
        contracted += primitive.coefficient * norm *
                      __expf(-primitive.exponent * r2);
    }
    return contracted * monomial;
}

__global__ void orbital_kernel(
    cudaSurfaceObject_t surface,
    const GpuShell* shells,
    const std::uint32_t shell_count,
    const GpuPrimitive* primitives,
    const float* coefficients,
    const GridBox box,
    const int nx,
    const int ny,
    const int nz) {

    const int ix = blockIdx.x * blockDim.x + threadIdx.x;
    const int iy = blockIdx.y * blockDim.y + threadIdx.y;
    const int iz = blockIdx.z * blockDim.z + threadIdx.z;
    if (ix >= nx || iy >= ny || iz >= nz) return;

    const float tx = nx > 1 ? static_cast<float>(ix) / static_cast<float>(nx - 1) : 0.0f;
    const float ty = ny > 1 ? static_cast<float>(iy) / static_cast<float>(ny - 1) : 0.0f;
    const float tz = nz > 1 ? static_cast<float>(iz) / static_cast<float>(nz - 1) : 0.0f;

    const float x = box.min_x + tx * (box.max_x - box.min_x);
    const float y = box.min_y + ty * (box.max_y - box.min_y);
    const float z = box.min_z + tz * (box.max_z - box.min_z);

    float psi = 0.0f;
    for (std::uint32_t s = 0; s < shell_count; ++s) {
        const GpuShell shell = shells[s];
        const float dx = x - shell.cx;
        const float dy = y - shell.cy;
        const float dz = z - shell.cz;
        const float r2 = dx * dx + dy * dy + dz * dz;

        const int n = shell_basis_count_device(shell);
        for (int c = 0; c < n; ++c) {
            const float basis = evaluate_shell_component(
                shell, primitives, c, dx, dy, dz, r2);
            psi = fmaf(coefficients[shell.basis_offset + c], basis, psi);
        }
    }

    surf3Dwrite(psi, surface,
                static_cast<std::size_t>(ix) * sizeof(float),
                iy, iz);
}

} // namespace

struct CudaOrbitalEvaluator::Impl {
    const Wavefunction* wf = nullptr;
    GpuShell* d_shells = nullptr;
    GpuPrimitive* d_primitives = nullptr;
    float* d_coefficients = nullptr;
    cudaGraphicsResource* texture_resource = nullptr;
    unsigned int texture = 0;
    std::string device_name_storage = "unknown";
    double last_ms = 0.0;

    explicit Impl(const Wavefunction& wavefunction) : wf(&wavefunction) {
        int device = 0;
        cuda_check(cudaGetDevice(&device), "cudaGetDevice");

        cudaDeviceProp prop{};
        cuda_check(cudaGetDeviceProperties(&prop, device), "cudaGetDeviceProperties");
        device_name_storage = prop.name;

        std::vector<GpuShell> shells;
        shells.reserve(wavefunction.shells.size());
        for (const Shell& shell : wavefunction.shells) {
            const Atom& atom = wavefunction.atoms.at(shell.atom_index);
            GpuShell gpu{};
            gpu.cx = static_cast<float>(atom.x);
            gpu.cy = static_cast<float>(atom.y);
            gpu.cz = static_cast<float>(atom.z);
            gpu.primitive_offset = shell.primitive_offset;
            gpu.primitive_count = shell.primitive_count;
            gpu.basis_offset = shell.basis_offset;
            gpu.l = shell.angular_momentum;
            gpu.pure = shell.pure;
            shells.push_back(gpu);
        }

        std::vector<GpuPrimitive> primitives;
        primitives.reserve(wavefunction.primitives.size());
        for (const Primitive& p : wavefunction.primitives) {
            primitives.push_back({p.exponent, p.coefficient});
        }

        cuda_check(cudaMalloc(&d_shells, shells.size() * sizeof(GpuShell)),
                   "cudaMalloc shells");
        cuda_check(cudaMalloc(&d_primitives, primitives.size() * sizeof(GpuPrimitive)),
                   "cudaMalloc primitives");
        cuda_check(cudaMalloc(&d_coefficients,
                              wavefunction.basis_count * sizeof(float)),
                   "cudaMalloc MO coefficients");

        cuda_check(cudaMemcpy(d_shells, shells.data(),
                              shells.size() * sizeof(GpuShell),
                              cudaMemcpyHostToDevice),
                   "cudaMemcpy shells");
        cuda_check(cudaMemcpy(d_primitives, primitives.data(),
                              primitives.size() * sizeof(GpuPrimitive),
                              cudaMemcpyHostToDevice),
                   "cudaMemcpy primitives");
    }

    ~Impl() {
        if (texture_resource) {
            cudaGraphicsUnregisterResource(texture_resource);
        }
        cudaFree(d_coefficients);
        cudaFree(d_primitives);
        cudaFree(d_shells);
    }
};

CudaOrbitalEvaluator::CudaOrbitalEvaluator(const Wavefunction& wavefunction)
    : impl_(std::make_unique<Impl>(wavefunction)) {}

CudaOrbitalEvaluator::~CudaOrbitalEvaluator() = default;
CudaOrbitalEvaluator::CudaOrbitalEvaluator(CudaOrbitalEvaluator&&) noexcept = default;
CudaOrbitalEvaluator& CudaOrbitalEvaluator::operator=(CudaOrbitalEvaluator&&) noexcept = default;

void CudaOrbitalEvaluator::attach_gl_texture(const unsigned int texture) {
    detach_gl_texture();
    impl_->texture = texture;
    cuda_check(cudaGraphicsGLRegisterImage(
                   &impl_->texture_resource,
                   texture,
                   GL_TEXTURE_3D,
                   cudaGraphicsRegisterFlagsSurfaceLoadStore |
                       cudaGraphicsRegisterFlagsWriteDiscard),
               "cudaGraphicsGLRegisterImage");
}

void CudaOrbitalEvaluator::detach_gl_texture() {
    if (impl_ && impl_->texture_resource) {
        cuda_check(cudaGraphicsUnregisterResource(impl_->texture_resource),
                   "cudaGraphicsUnregisterResource");
        impl_->texture_resource = nullptr;
        impl_->texture = 0;
    }
}

void CudaOrbitalEvaluator::evaluate(const std::size_t mo_index,
                                    const GridBox& box,
                                    const int nx,
                                    const int ny,
                                    const int nz) {
    if (!impl_->texture_resource) {
        throw std::runtime_error("No OpenGL 3D texture attached to CUDA evaluator");
    }
    if (mo_index >= impl_->wf->orbitals.size()) {
        throw std::out_of_range("MO index out of range");
    }
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Invalid grid dimensions");
    }

    const auto& mo = impl_->wf->orbitals[mo_index];
    cuda_check(cudaMemcpy(impl_->d_coefficients, mo.coefficients.data(),
                          impl_->wf->basis_count * sizeof(float),
                          cudaMemcpyHostToDevice),
               "cudaMemcpy MO coefficients");

    cuda_check(cudaGraphicsMapResources(1, &impl_->texture_resource, 0),
               "cudaGraphicsMapResources");

    cudaArray_t array = nullptr;
    cuda_check(cudaGraphicsSubResourceGetMappedArray(
                   &array, impl_->texture_resource, 0, 0),
               "cudaGraphicsSubResourceGetMappedArray");

    cudaResourceDesc resource_desc{};
    resource_desc.resType = cudaResourceTypeArray;
    resource_desc.res.array.array = array;

    cudaSurfaceObject_t surface = 0;
    cuda_check(cudaCreateSurfaceObject(&surface, &resource_desc),
               "cudaCreateSurfaceObject");

    cudaEvent_t start{}, stop{};
    cuda_check(cudaEventCreate(&start), "cudaEventCreate start");
    cuda_check(cudaEventCreate(&stop), "cudaEventCreate stop");

    const dim3 block(8, 8, 4);
    const dim3 grid(
        static_cast<unsigned int>((nx + block.x - 1) / block.x),
        static_cast<unsigned int>((ny + block.y - 1) / block.y),
        static_cast<unsigned int>((nz + block.z - 1) / block.z));

    cuda_check(cudaEventRecord(start), "cudaEventRecord start");
    orbital_kernel<<<grid, block>>>(
        surface,
        impl_->d_shells,
        static_cast<std::uint32_t>(impl_->wf->shells.size()),
        impl_->d_primitives,
        impl_->d_coefficients,
        box,
        nx, ny, nz);
    cuda_check(cudaGetLastError(), "orbital_kernel launch");
    cuda_check(cudaEventRecord(stop), "cudaEventRecord stop");
    cuda_check(cudaEventSynchronize(stop), "cudaEventSynchronize");

    float elapsed_ms = 0.0f;
    cuda_check(cudaEventElapsedTime(&elapsed_ms, start, stop),
               "cudaEventElapsedTime");
    impl_->last_ms = elapsed_ms;

    cudaEventDestroy(stop);
    cudaEventDestroy(start);
    cudaDestroySurfaceObject(surface);
    cuda_check(cudaGraphicsUnmapResources(1, &impl_->texture_resource, 0),
               "cudaGraphicsUnmapResources");
}

const char* CudaOrbitalEvaluator::device_name() const noexcept {
    return impl_ ? impl_->device_name_storage.c_str() : "unknown";
}

double CudaOrbitalEvaluator::last_kernel_ms() const noexcept {
    return impl_ ? impl_->last_ms : 0.0;
}

} // namespace cov
