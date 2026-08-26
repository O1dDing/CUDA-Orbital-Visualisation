#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace cov {

constexpr double kAngstromToBohr = 1.8897261254578281;

enum class Spin : std::uint8_t {
    Alpha = 0,
    Beta = 1,
};

struct Atom {
    std::string symbol;
    int atomic_number = 0;
    double x = 0.0; // bohr
    double y = 0.0; // bohr
    double z = 0.0; // bohr
};

struct Primitive {
    float exponent = 0.0f;
    float coefficient = 0.0f;
};

struct Shell {
    std::uint32_t atom_index = 0;
    std::uint32_t primitive_offset = 0;
    std::uint32_t primitive_count = 0;
    std::uint32_t basis_offset = 0;
    std::uint8_t angular_momentum = 0; // s=0, p=1, d=2, ...
    std::uint8_t pure = 0;             // 0 Cartesian, 1 real spherical
};

struct MolecularOrbital {
    double energy_hartree = 0.0;
    float occupation = 0.0f;
    Spin spin = Spin::Alpha;
    std::string symmetry;
    std::vector<float> coefficients;
};

struct Wavefunction {
    std::vector<Atom> atoms;
    std::vector<Primitive> primitives;
    std::vector<Shell> shells;
    std::vector<MolecularOrbital> orbitals;
    std::uint32_t basis_count = 0;
    bool pure_d = false;
    bool pure_f = false;
    bool pure_g = false;
};

struct GridBox {
    float min_x = -5.0f;
    float min_y = -5.0f;
    float min_z = -5.0f;
    float max_x = 5.0f;
    float max_y = 5.0f;
    float max_z = 5.0f;
};

inline std::uint32_t cartesian_basis_count(const std::uint8_t l) {
    return static_cast<std::uint32_t>((l + 1u) * (l + 2u) / 2u);
}

inline std::uint32_t spherical_basis_count(const std::uint8_t l) {
    return static_cast<std::uint32_t>(2u * l + 1u);
}

inline std::uint32_t shell_basis_count(const Shell& shell) {
    return shell.pure ? spherical_basis_count(shell.angular_momentum)
                      : cartesian_basis_count(shell.angular_momentum);
}

} // namespace cov
