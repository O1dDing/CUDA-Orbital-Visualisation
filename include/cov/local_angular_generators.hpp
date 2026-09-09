#pragma once
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace cov {
struct LocalAngularGenerators {
    bool available = false;
    int shell_degree = -1;
    int angular_degree = -1;
    bool pure_shell = false;
    std::size_t basis_count = 0;
    std::size_t columns = 0;
    // Row-major AO coefficients, one generator for each real angular function.
    // The order is m=0, real m=1, imaginary m=1, real m=2, ... with Condon--Shortley phase.
    std::vector<double> coefficients;
    double rotation_orthogonality_error = 0;
    double polynomial_reconstruction_error = 0;
    std::string detail;
};

// Functions are defined in the supplied local frame. Cartesian degree L also
// contains r^(L-l) Y_lm for l=L-2,L-4,...; these directions remain explicit.
// Radial contractions and overlap normalization are handled by the caller's
// common AO metric, not by independent coefficient-squared weights.
LocalAngularGenerators local_angular_generators(
    int shell_degree, bool pure_shell, int angular_degree,
    const std::array<double,9>& rotation_reference_to_input);
}
