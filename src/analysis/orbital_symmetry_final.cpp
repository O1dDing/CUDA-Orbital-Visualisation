#include "cov/orbital_symmetry.hpp"

// Keep the already-tested finite/linear machinery in one translation unit so
// the final dispatcher can reuse its AO-operation helpers without exposing them
// as public API. The included implementation is renamed here; CMake compiles
// this file instead of compiling orbital_symmetry_dispatch.cpp directly.
#define derive_orbital_symmetry derive_orbital_symmetry_dispatch_v1
#include "orbital_symmetry_dispatch.cpp"
#undef derive_orbital_symmetry

namespace cov {
namespace {

SymmetryOperation square_operation(const SymmetryOperation& op) {
    SymmetryOperation squared = op;
    squared.kind = SymmetryOperationKind::ProperRotation;
    squared.order = 2;
    squared.power = 2;

    std::array<double, 9> matrix{};
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            for (int k = 0; k < 3; ++k) {
                matrix[3 * row + col] +=
                    op.matrix[3 * row + k] * op.matrix[3 * k + col];
            }
        }
    }
    squared.matrix = matrix;

    squared.atom_permutation.assign(op.atom_permutation.size(), 0u);
    for (std::size_t i = 0; i < op.atom_permutation.size(); ++i) {
        const std::size_t middle = op.atom_permutation[i];
        if (middle >= op.atom_permutation.size()) {
            squared.atom_permutation.clear();
            break;
        }
        squared.atom_permutation[i] = op.atom_permutation[middle];
    }
    return squared;
}

std::optional<std::string> classify_oh_native(
    const Wavefunction& wavefunction,
    const MolecularSymmetry& symmetry,
    const std::vector<std::size_t>& group,
    const OrbitalSymmetryOptions& options,
    double& retention) {
    if (symmetry.point_group != "Oh") return std::nullopt;

    const auto* c4 = find_operation(
        symmetry, SymmetryOperationKind::ProperRotation, 4);
    const auto* inversion = find_operation(
        symmetry, SymmetryOperationKind::Inversion, 2);
    if (!c4 || !inversion) return std::nullopt;

    // COV's geometry engine may identify Oh from the three C4 axes before a
    // body-diagonal C3 candidate is explicitly materialised.  C4, C4^2 and
    // inversion are nevertheless a complete discriminator for the five O
    // irreps at their known dimensions and also reject the common accidental
    // A+E three-dimensional reducible subspace via the C4^2 character.
    const SymmetryOperation c2_axis = square_operation(*c4);
    if (c2_axis.atom_permutation.size() != wavefunction.atoms.size()) {
        return std::nullopt;
    }

    const auto c4_character = character(wavefunction, *c4, group);
    const auto c2_character = character(wavefunction, c2_axis, group);
    const auto parity_character = character(wavefunction, *inversion, group);
    if (!c4_character.valid || !c2_character.valid ||
        !parity_character.valid) {
        return std::nullopt;
    }

    retention = std::min({c4_character.retention,
                          c2_character.retention,
                          parity_character.retention});
    if (retention < options.minimum_subspace_retention) return std::nullopt;

    const int dimension = static_cast<int>(group.size());
    const double tolerance = options.character_tolerance;
    if (std::abs(std::abs(parity_character.value) -
                 static_cast<double>(dimension)) > tolerance) {
        return std::nullopt;
    }
    const char parity = parity_character.value >= 0.0 ? 'g' : 'u';

    std::string base;
    if (dimension == 1 &&
        std::abs(c2_character.value - 1.0) <= tolerance) {
        if (std::abs(c4_character.value - 1.0) <= tolerance) base = "A1";
        else if (std::abs(c4_character.value + 1.0) <= tolerance) base = "A2";
    } else if (dimension == 2 &&
               std::abs(c4_character.value) <= tolerance &&
               std::abs(c2_character.value - 2.0) <= tolerance) {
        base = "E";
    } else if (dimension == 3 &&
               std::abs(c2_character.value + 1.0) <= tolerance) {
        if (std::abs(c4_character.value - 1.0) <= tolerance) base = "T1";
        else if (std::abs(c4_character.value + 1.0) <= tolerance) base = "T2";
    }

    if (base.empty()) return std::nullopt;
    return base + parity;
}

} // namespace

OrbitalSymmetryResult derive_orbital_symmetry(
    Wavefunction& wavefunction,
    const OrbitalSymmetryOptions& options) {
    if (wavefunction.orbitals.empty() || wavefunction.basis_count == 0 ||
        wavefunction.ao_overlap.size() !=
            static_cast<std::size_t>(wavefunction.basis_count) *
                wavefunction.basis_count) {
        return {};
    }

    const MolecularSymmetry symmetry = analyse_molecular_symmetry(wavefunction);

    // Oh is handled before the legacy finite-group path so no provisional
    // Derived label can block the deterministic COV-native classifier. Producer
    // labels remain immutable because group_unlabelled() rejects such groups.
    if (symmetry.point_group == "Oh") {
        OrbitalSymmetryResult result;
        result.point_group = symmetry.point_group;
        if (!symmetry.available()) return result;

        for (const auto& group :
             energy_groups(wavefunction, options.degeneracy_tolerance_hartree)) {
            if (!group_unlabelled(wavefunction, group)) continue;
            ++result.groups_examined;

            double retention = 1.0;
            const auto label = classify_oh_native(
                wavefunction, symmetry, group, options, retention);
            result.worst_subspace_retention =
                std::min(result.worst_subspace_retention, retention);
            if (!label) continue;

            ++result.groups_labelled;
            for (const std::size_t index : group) {
                wavefunction.orbitals[index].symmetry = *label;
                wavefunction.orbitals[index].symmetry_provenance =
                    DataProvenance::Derived;
                ++result.orbitals_labelled;
            }
        }
        return result;
    }

    return derive_orbital_symmetry_dispatch_v1(wavefunction, options);
}

} // namespace cov
