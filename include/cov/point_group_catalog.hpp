#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace cov {

enum class PointGroupFamily : std::uint8_t {
    Trivial = 0,
    Reflection,
    Cyclic,
    Dihedral,
    Linear,
    Polyhedral,
};

enum class MetalAOShell : std::uint8_t {
    S = 0,
    P,
    D,
};

struct IrrepDefinition {
    std::string_view label;
    std::uint8_t dimension = 1;
};

struct PointGroupDefinition {
    std::string_view canonical_name;
    PointGroupFamily family = PointGroupFamily::Trivial;
    // Zero denotes an infinite group.
    std::uint16_t order = 1;
    bool linear = false;
    bool centrosymmetric = false;
    std::span<const std::string_view> aliases;
    // Complete finite-group irrep list.  For Dinfh this is the subset needed
    // by central-atom s, p and d orbitals.
    std::span<const IrrepDefinition> irreps;
};

// One occurrence of an irrep in a shell decomposition.  Repeated labels are
// deliberately represented by repeated objects rather than a set or a map:
// C2v d orbitals, for example, contain two distinct A1 copies.
struct IrrepCopy {
    std::string_view label;
    std::uint8_t dimension = 1;
    MetalAOShell shell = MetalAOShell::S;
    // One-based ordinal among equal labels in this shell decomposition.
    std::uint8_t copy_index = 1;
    std::string_view basis_functions;
};

using IrrepDecomposition = std::vector<IrrepCopy>;

struct MetalSPDDecomposition {
    IrrepDecomposition s;
    IrrepDecomposition p;
    IrrepDecomposition d;

    [[nodiscard]] std::size_t total_dimension() const noexcept;
};

// Formal parent groups are listed together with common lower-symmetry groups
// so later geometry/MO layers can share one canonical naming boundary.
[[nodiscard]] std::span<const PointGroupDefinition> point_group_catalog() noexcept;

// Lookup is case-insensitive and accepts common ASCII aliases such as O_h,
// D_inf_h and D-infinity-h.  Returned pointers remain valid for program life.
[[nodiscard]] const PointGroupDefinition* find_point_group(
    std::string_view name) noexcept;

[[nodiscard]] std::optional<IrrepDecomposition> decompose_metal_ao_shell(
    std::string_view point_group,
    MetalAOShell shell);

[[nodiscard]] std::optional<MetalSPDDecomposition> decompose_metal_spd(
    std::string_view point_group);

[[nodiscard]] std::size_t irrep_multiplicity(
    std::span<const IrrepCopy> decomposition,
    std::string_view label) noexcept;

[[nodiscard]] std::size_t irrep_total_dimension(
    std::span<const IrrepCopy> decomposition) noexcept;

} // namespace cov
