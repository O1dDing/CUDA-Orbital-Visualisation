#include "cov/density.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cov {

std::vector<double> reconstruct_total_density_packed(const Wavefunction& wavefunction) {
    const std::size_t n = wavefunction.basis_count;
    if (n == 0) return {};

    const std::size_t packed_size = n * (n + 1u) / 2u;
    std::vector<double> density(packed_size, 0.0);

    for (const MolecularOrbital& mo : wavefunction.orbitals) {
        if (mo.occupation == 0.0f) continue;
        if (mo.coefficients.size() != n) {
            throw std::runtime_error(
                "Cannot reconstruct density: MO coefficient count does not match basis count");
        }

        const double occupation = static_cast<double>(mo.occupation);
        for (std::size_t i = 0; i < n; ++i) {
            const double ci = static_cast<double>(mo.coefficients[i]);
            const std::size_t row = i * (i + 1u) / 2u;
            for (std::size_t j = 0; j <= i; ++j) {
                density[row + j] += occupation * ci *
                                    static_cast<double>(mo.coefficients[j]);
            }
        }
    }

    return density;
}

} // namespace cov
