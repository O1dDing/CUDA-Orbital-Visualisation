#pragma once

#include "cov/model.hpp"
#include "cov/point_group_catalog.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>

namespace cov {

// Symmetry assignment made from the central atom's local s/p/d AO content.
// copy_index is retained because lower-symmetry groups can contain multiple
// independent copies carrying the same textual irrep label.  Zero means that
// the irrep label is determined but equivalent copies are internally mixed.
struct LocalIrrepAssignment {
    std::string label;
    MetalAOShell shell = MetalAOShell::S;
    std::uint8_t copy_index = 1;
    double confidence = 0.0;
    std::string basis_functions;
};

// Coefficient-free fallback for a complete local subspace.  It succeeds only
// when all irrep copies with the supplied dimension carry one textual label;
// genuinely different same-dimensional irreps remain unassigned.
[[nodiscard]] std::optional<LocalIrrepAssignment>
classify_local_irrep_by_dimension(
    std::string_view point_group,
    MetalAOShell shell,
    std::size_t dimension);

// rotation_reference_to_input is a row-major orthogonal matrix mapping a
// vector in the point-group reference frame into the input molecular frame.
// The classifier rotates central-metal AO amplitudes back to the reference
// frame, then scores the catalogued s/p/d irreps.  Equivalent repeated copies
// are first aggregated by label; failure to select a copy therefore never
// discards an otherwise rigorous irrep assignment.
[[nodiscard]] std::optional<LocalIrrepAssignment> classify_local_metal_irrep(
    const Wavefunction& wavefunction,
    std::span<const std::size_t> orbital_indices,
    std::size_t metal_atom,
    std::string_view point_group,
    const std::array<double, 9>& rotation_reference_to_input);

} // namespace cov
