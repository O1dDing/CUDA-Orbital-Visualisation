#include "cov/overlap.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace cov {
namespace {

std::vector<const MolecularOrbital*> complete_spin_block(const Wavefunction& wf,
                                                         const Spin spin) {
    std::vector<const MolecularOrbital*> block;
    block.reserve(wf.basis_count);
    for (const auto& mo : wf.orbitals) {
        if (mo.spin == spin) block.push_back(&mo);
    }
    if (block.size() != wf.basis_count) return {};
    for (const auto* mo : block) {
        if (!mo || mo->coefficients.size() != wf.basis_count) return {};
    }
    return block;
}

bool invert_square(std::vector<double>& a,
                   std::vector<double>& inverse,
                   const std::size_t n,
                   const double relative_tolerance,
                   double& pivot_ratio) {
    inverse.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) inverse[i * n + i] = 1.0;

    double scale = 0.0;
    for (const double value : a) scale = std::max(scale, std::abs(value));
    if (!(scale > 0.0) || !std::isfinite(scale)) return false;

    double min_pivot = std::numeric_limits<double>::infinity();
    double max_pivot = 0.0;

    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot_row = col;
        double pivot_abs = std::abs(a[col * n + col]);
        for (std::size_t row = col + 1; row < n; ++row) {
            const double candidate = std::abs(a[row * n + col]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }
        if (!std::isfinite(pivot_abs) ||
            pivot_abs <= relative_tolerance * scale) {
            return false;
        }

        min_pivot = std::min(min_pivot, pivot_abs);
        max_pivot = std::max(max_pivot, pivot_abs);
        if (pivot_row != col) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(a[col * n + j], a[pivot_row * n + j]);
                std::swap(inverse[col * n + j], inverse[pivot_row * n + j]);
            }
        }

        const double pivot = a[col * n + col];
        for (std::size_t j = 0; j < n; ++j) {
            a[col * n + j] /= pivot;
            inverse[col * n + j] /= pivot;
        }

        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            const double factor = a[row * n + col];
            if (factor == 0.0) continue;
            a[row * n + col] = 0.0;
            for (std::size_t j = col + 1; j < n; ++j) {
                a[row * n + j] -= factor * a[col * n + j];
            }
            for (std::size_t j = 0; j < n; ++j) {
                inverse[row * n + j] -= factor * inverse[col * n + j];
            }
        }
    }

    pivot_ratio = max_pivot > 0.0 ? min_pivot / max_pivot : 0.0;
    return true;
}

} // namespace

OverlapDerivationResult derive_ao_overlap_from_mos(const Wavefunction& wavefunction,
                                                   const double relative_pivot_tolerance,
                                                   const std::size_t maximum_basis) {
    OverlapDerivationResult result;
    const std::size_t n = wavefunction.basis_count;
    if (n == 0 || n > maximum_basis) return result;

    auto block = complete_spin_block(wavefunction, Spin::Alpha);
    if (block.empty()) block = complete_spin_block(wavefunction, Spin::Beta);
    if (block.empty()) return result;

    // C is stored AO-major here: C(mu,i).
    std::vector<double> c(n * n, 0.0);
    for (std::size_t orbital = 0; orbital < n; ++orbital) {
        for (std::size_t basis = 0; basis < n; ++basis) {
            c[basis * n + orbital] =
                static_cast<double>(block[orbital]->coefficients[basis]);
        }
    }

    std::vector<double> work = c;
    std::vector<double> c_inverse;
    if (!invert_square(work, c_inverse, n, relative_pivot_tolerance,
                       result.pivot_ratio)) {
        result.pivot_ratio = 0.0;
        return result;
    }

    // If X=C^{-1}, then S=X^T X. Accumulate only once and mirror explicitly.
    result.matrix.assign(n * n, 0.0);
    for (std::size_t mu = 0; mu < n; ++mu) {
        for (std::size_t nu = 0; nu <= mu; ++nu) {
            double value = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                value += c_inverse[i * n + mu] * c_inverse[i * n + nu];
            }
            result.matrix[mu * n + nu] = value;
            result.matrix[nu * n + mu] = value;
        }
    }

    // Validate C^T S C=I without trusting the algebra alone. This catches
    // non-finite coefficients and future basis-order mistakes at the boundary.
    std::vector<double> sc(n * n, 0.0);
    for (std::size_t mu = 0; mu < n; ++mu) {
        for (std::size_t i = 0; i < n; ++i) {
            double value = 0.0;
            for (std::size_t nu = 0; nu < n; ++nu) {
                value += result.matrix[mu * n + nu] * c[nu * n + i];
            }
            sc[mu * n + i] = value;
        }
    }

    double max_error = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double value = 0.0;
            for (std::size_t mu = 0; mu < n; ++mu) {
                value += c[mu * n + i] * sc[mu * n + j];
            }
            const double target = i == j ? 1.0 : 0.0;
            max_error = std::max(max_error, std::abs(value - target));
        }
    }
    result.max_orthonormality_error = max_error;

    if (!std::isfinite(max_error) || max_error > 5.0e-6) {
        result.matrix.clear();
    }
    return result;
}

} // namespace cov
