#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cov {

// phi_internal[i] = basis_scale * phi_Gaussian[source_index].  Gaussian and
// COV use the same primitive normalization; the scale here is a phase (+/-1).
// Keep coefficient_scale explicit so a future normalization change cannot
// accidentally apply a covariant matrix transform to contravariant data.
struct GaussianAoTransformEntry {
    std::size_t source_index = 0;
    double basis_scale = 1.0;
    double coefficient_scale = 1.0;
};

inline std::vector<GaussianAoTransformEntry> gaussian_ao_transform(
    const int shell_type, const std::size_t source_offset = 0) {
    const int l = shell_type < 0 ? -shell_type : shell_type;
    if (l > 4) throw std::runtime_error("Gaussian AO transform above g is unavailable");
    const bool pure = shell_type <= -2;
    const std::size_t count = shell_type == -1 ? 4u : pure
        ? static_cast<std::size_t>(2*l+1)
        : static_cast<std::size_t>((l+1)*(l+2)/2);
    static constexpr std::size_t cartesian_g[15] = {
        14,4,0,13,12,8,3,5,1,11,9,2,10,7,6};
    std::vector<GaussianAoTransformEntry> result;
    result.reserve(count);
    for (std::size_t i=0; i<count; ++i) {
        // Real solid harmonics use m=0,+1,-1,+2,-2,... . COV's
        // polynomials include the Condon--Shortley (-1)^|m| phase;
        // Gaussian FCHK pure shells omit that phase.
        const std::size_t m = (i+1)/2;
        const double phase = pure && m%2 != 0 ? -1.0 : 1.0;
        const std::size_t source = !pure && l==4 ? cartesian_g[i] : i;
        result.push_back({source_offset+source, phase, 1.0/phase});
    }
    return result;
}

} // namespace cov
