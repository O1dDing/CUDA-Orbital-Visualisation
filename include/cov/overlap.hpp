#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <vector>

namespace cov {

struct OverlapDerivationResult {
    std::vector<double> matrix; // row-major basis_count x basis_count
    double max_orthonormality_error = 0.0;
    double pivot_ratio = 0.0;

    [[nodiscard]] bool available() const noexcept { return !matrix.empty(); }
};

// Recover the AO overlap matrix from a complete orthonormal MO coefficient
// matrix C using C^T S C = I, hence S = C^{-T} C^{-1}. This is format-agnostic
// and avoids duplicating basis-integral convention logic in the first FCHK
// analysis path. A complete, nonsingular spin block with exactly basis_count
// orbitals is required; otherwise an unavailable result is returned.
[[nodiscard]] OverlapDerivationResult derive_ao_overlap_from_mos(
    const Wavefunction& wavefunction,
    double relative_pivot_tolerance = 1.0e-11,
    // 1024 keeps the dense recovery workspace comfortably below 128 MiB
    // while covering large, chemically relevant multi-ring/ECP validation
    // systems. Callers analysing still larger matrices may opt in explicitly.
    std::size_t maximum_basis = 1024);

} // namespace cov
