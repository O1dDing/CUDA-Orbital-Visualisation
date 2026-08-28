#include "cov/orbital_symmetry.hpp"

// Keep the already-tested finite/linear machinery in one translation unit so
// the final dispatcher can reuse its AO-operation helpers without exposing them
// as public API.  The included implementation is renamed here; CMake compiles
// this file instead of compiling orbital_symmetry_dispatch.cpp directly.
#define derive_orbital_symmetry derive_orbital_symmetry_dispatch_v1
#include "orbital_symmetry_dispatch.cpp"
#undef derive_orbital_symmetry

namespace cov {

OrbitalSymmetryResult derive_orbital_symmetry(
    Wavefunction& wavefunction,
    const OrbitalSymmetryOptions& options) {
    if (wavefunction.orbitals.empty() || wavefunction.basis_count == 0 ||
        wavefunction.ao_overlap.size() !=
            static_cast<std::size_t>(wavefunction.basis_count) * wavefunction.basis_count) {
        return {};
    }

    const MolecularSymmetry symmetry = analyse_molecular_symmetry(wavefunction);

    // Oh is deliberately handled before the legacy finite-group path.  The
    // legacy classifier distinguishes the two octahedral C2 classes by picking
    // a representative operation from an unordered candidate set; that can
    // swap/reject T1/T2 for otherwise valid subspaces.  COV's native Oh path
    // needs only C3, C4 and inversion characters, which uniquely distinguish
    // A1/A2/E/T1/T2 and g/u and are independent of that representative choice.
    if (symmetry.point_group == "Oh") {
        OrbitalSymmetryResult result;
        result.point_group = symmetry.point_group;
        if (!symmetry.available()) return result;

        for (const auto& group :
             energy_groups(wavefunction, options.degeneracy_tolerance_hartree)) {
            if (!group_unlabelled(wavefunction, group)) continue;
            ++result.groups_examined;

            double retention = 1.0;
            const auto label = classify_oh_fallback(
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
