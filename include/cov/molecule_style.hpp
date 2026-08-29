#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <vector>

namespace cov {

enum class MoleculeStyle {
    MediumBallAndStick = 0, // UI: Enhanced Ball-and-Stick
    StickDelocalisation,
};

struct MoleculeRenderSettings {
    MoleculeStyle style = MoleculeStyle::MediumBallAndStick;
    // 1.0 is the human-facing reference size. The renderer's base sphere and
    // cylinder radii are calibrated so 1.0 reproduces the proportions that
    // were manually accepted during pre.2 validation.
    float atom_scale = 1.00f;
    float bond_scale = 1.00f;
    float molecule_opacity = 1.00f;
    float orbital_opacity = 0.92f;
    bool show_hydrogens = true;
};

struct BondVisual {
    std::size_t atom_a = 0;
    std::size_t atom_b = 0;
    double distance_bohr = 0.0;
    bool delocalised = false;
    double bond_order = 0.0;
    DataProvenance provenance = DataProvenance::Unavailable;
};

[[nodiscard]] double covalent_radius_angstrom(int atomic_number) noexcept;

// Wavefunction-derived Mayer evidence is used when available, inside a
// conservative structural-neighbour envelope. This keeps non-local electronic
// coupling from being drawn as a bond. Covalent-radius geometry is the fallback
// only when electronic evidence is unavailable. Pairwise Mayer evidence is not
// an electron-counting or multicentre claim.
[[nodiscard]] std::vector<BondVisual> analyse_bonds(const Wavefunction& wavefunction);

} // namespace cov
