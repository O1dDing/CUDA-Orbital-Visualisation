#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <string>

namespace cov {

struct OrbitalSymmetryOptions {
    double degeneracy_tolerance_hartree = 1.0e-5;
    double character_tolerance = 0.18;
    double minimum_subspace_retention = 0.985;
};

struct OrbitalSymmetryResult {
    std::string point_group;
    std::size_t groups_examined = 0;
    std::size_t groups_labelled = 0;
    std::size_t orbitals_labelled = 0;
    double worst_subspace_retention = 1.0;

    [[nodiscard]] bool applied() const noexcept {
        return orbitals_labelled != 0;
    }
};

// Derive conventional MO irreps from COV's own validated molecular symmetry
// operations. The implementation transforms the actual AO basis, evaluates
// each degenerate MO subspace with the AO-overlap metric, and only assigns a
// label when the resulting characters match a supported point-group irrep.
// Producer-reported MO labels are never overwritten.
OrbitalSymmetryResult derive_orbital_symmetry(
    Wavefunction& wavefunction,
    const OrbitalSymmetryOptions& options = {});

} // namespace cov
