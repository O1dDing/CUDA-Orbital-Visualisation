#include "cov/density.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cov {
namespace {

void accumulate_outer(std::vector<double>& density,
                      const MolecularOrbital& mo,
                      const std::size_t n,
                      const double weight) {
    if (weight == 0.0) return;
    if (mo.coefficients.size() != n) {
        throw std::runtime_error(
            "Cannot reconstruct density: MO coefficient count does not match basis count");
    }

    for (std::size_t i = 0; i < n; ++i) {
        const double ci = static_cast<double>(mo.coefficients[i]);
        const std::size_t row = i * (i + 1u) / 2u;
        for (std::size_t j = 0; j <= i; ++j) {
            density[row + j] += weight * ci *
                                static_cast<double>(mo.coefficients[j]);
        }
    }
}

} // namespace

std::vector<double> reconstruct_total_density_packed(const Wavefunction& wavefunction) {
    const std::size_t n = wavefunction.basis_count;
    if (n == 0) return {};

    std::vector<double> density(n * (n + 1u) / 2u, 0.0);
    for (const MolecularOrbital& mo : wavefunction.orbitals) {
        accumulate_outer(density, mo, n, static_cast<double>(mo.occupation));
    }
    return density;
}

std::vector<double> reconstruct_spin_density_packed(const Wavefunction& wavefunction) {
    const std::size_t n = wavefunction.basis_count;
    if (n == 0 || wavefunction.orbitals.empty()) return {};

    std::vector<double> density(n * (n + 1u) / 2u, 0.0);
    bool has_beta_block = false;
    for (const auto& mo : wavefunction.orbitals) {
        if (mo.spin == Spin::Beta) {
            has_beta_block = true;
            break;
        }
    }

    if (has_beta_block) {
        for (const auto& mo : wavefunction.orbitals) {
            const double sign = mo.spin == Spin::Alpha ? 1.0 : -1.0;
            accumulate_outer(density, mo, n,
                             sign * static_cast<double>(mo.occupation));
        }
        return density;
    }

    // Restricted/ROHF-style FCHK stores one canonical orbital block. Recover
    // alpha-beta occupation from the electron counts rather than equating the
    // spin-summed occupation with spin density.
    std::size_t orbital = 0;
    for (const auto& mo : wavefunction.orbitals) {
        const double alpha = orbital < wavefunction.alpha_electrons ? 1.0 : 0.0;
        const double beta = orbital < wavefunction.beta_electrons ? 1.0 : 0.0;
        accumulate_outer(density, mo, n, alpha - beta);
        ++orbital;
    }
    return density;
}

} // namespace cov
