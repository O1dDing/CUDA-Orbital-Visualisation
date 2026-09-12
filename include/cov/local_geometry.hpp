#pragma once

#include "cov/coordination_geometry.hpp"
#include "cov/model.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace cov {

// A local molecular geometry is a structural property of any bonded centre.
// It is deliberately separate from LigandFieldEnvironment: main-group species
// such as PF5, SF6, XeF4 and [SbF6]- have well-defined local geometries but do
// not acquire a transition-metal d-shell splitting diagram merely because
// their coordination number matches a catalogue entry.
struct LocalGeometryEnvironment {
    std::size_t centre_atom = 0;
    GeometryId geometry_id = GeometryId::Unknown;
    std::vector<std::size_t> neighbour_atoms;
    double confidence = 0.0;
    double angular_rms = 0.0;
    double shape_measure = 0.0;
    double radial_cv = 0.0;
    bool ambiguous = false;

    [[nodiscard]] bool available() const noexcept {
        return geometry_id != GeometryId::Unknown && !ambiguous &&
               !neighbour_atoms.empty();
    }

    [[nodiscard]] std::size_t coordination_number() const noexcept {
        return neighbour_atoms.size();
    }
};

// Analyse every atom with a Mayer-supported first shell of size CN 2--10.
// No element, molecule name, atom ordering or file name is consulted.
[[nodiscard]] std::vector<LocalGeometryEnvironment>
analyse_local_molecular_geometries(const Wavefunction& wavefunction);

// Returns a unique, high-confidence principal centre when one exists.  The
// ordering is based on coordination number and fit quality only; a tie remains
// unresolved rather than silently choosing an atom by input order.
[[nodiscard]] std::optional<LocalGeometryEnvironment>
principal_local_molecular_geometry(const Wavefunction& wavefunction);

} // namespace cov
