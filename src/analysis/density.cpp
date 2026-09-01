#include "cov/density.hpp"

#include <cmath>
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

    if (wavefunction.alpha_electrons + wavefunction.beta_electrons > 0u) {
        // Restricted/ROHF-style FCHK stores one canonical orbital block.
        // Recover alpha-beta occupation from explicit electron counts.
        std::size_t orbital = 0;
        for (const auto& mo : wavefunction.orbitals) {
            const double alpha = orbital < wavefunction.alpha_electrons ? 1.0 : 0.0;
            const double beta = orbital < wavefunction.beta_electrons ? 1.0 : 0.0;
            accumulate_outer(density, mo, n, alpha - beta);
            ++orbital;
        }
        return density;
    }

    // Molden often has no global alpha/beta electron-count fields but does carry
    // Spin= and Occup= per orbital. For a single Alpha block, occupation 2 is a
    // doubly occupied restricted orbital (zero spin density), occupation 1 is a
    // singly occupied alpha orbital, and 0 is virtual. Fractional occupations
    // between these integer cases are insufficient to reconstruct alpha-beta
    // partitioning uniquely, so leave spin density unavailable instead of
    // inventing a decomposition.
    for (const auto& mo : wavefunction.orbitals) {
        const double occ = static_cast<double>(mo.occupation);
        double spin_occupation = 0.0;
        if (std::abs(occ) <= 1.0e-5 || std::abs(occ - 2.0) <= 1.0e-5) {
            spin_occupation = 0.0;
        } else if (std::abs(occ - 1.0) <= 1.0e-5 && mo.spin == Spin::Alpha) {
            spin_occupation = 1.0;
        } else {
            return {};
        }
        accumulate_outer(density, mo, n, spin_occupation);
    }
    return density;
}

} // namespace cov
