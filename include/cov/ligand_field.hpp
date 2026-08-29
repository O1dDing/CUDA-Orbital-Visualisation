#pragma once

#include "cov/model.hpp"
#include "cov/orbital_view.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace cov {

// Once a locally degenerate subspace carries this much selected-metal d
// character, its Oh/Td d irrep is determined by its dimension (Eg/E or
// T2g/T2).  Keeping this shared avoids the browser and reduced diagram using
// different cutoffs for ligand-dominated antibonding partners.
inline constexpr double kLocalDIrrepWeightFloor=0.02;

enum class LigandFieldGeometry {
    Unknown,
    Tetrahedral,
    Octahedral,
};

// A local coordination environment is deliberately separate from the full
// molecular point group.  Substituent orientations can lower the latter while
// leaving an unambiguous first-shell ligand field around a transition metal.
struct LigandFieldEnvironment {
    LigandFieldGeometry geometry = LigandFieldGeometry::Unknown;
    std::size_t metal_atom = 0;
    std::vector<std::size_t> ligand_atoms;
    double confidence = 0.0;
    bool equivalent_ligand_elements = false;

    [[nodiscard]] bool available() const noexcept {
        return geometry != LigandFieldGeometry::Unknown &&
               !ligand_atoms.empty();
    }

    [[nodiscard]] std::string local_point_group() const;
};

// Uses density-derived Mayer connectivity to isolate the first coordination
// shell, removes collinear through-ligand contacts, and classifies only
// geometries supported by the ligand-field diagram (currently Td and Oh).
// No molecule name, atom ordering, or canonical-MO number is consulted.
[[nodiscard]] LigandFieldEnvironment analyse_ligand_field_environment(
    const Wavefunction& wavefunction);

// Apply local Oh/Td labels to the same per-MO metadata consumed by the
// orbital browser, diagram tooltips and exports.  Producer labels are never
// overwritten: only unavailable members of a confidently recovered central-
// metal valence subspace are filled.  This keeps the ordinary MO table and
// the reduced ligand-field diagram on one metadata path.
void apply_local_ligand_field_symmetry(
    const Wavefunction& wavefunction,
    std::vector<OrbitalMetadata>& metadata);

} // namespace cov
