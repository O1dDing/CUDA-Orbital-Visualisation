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
    // Manual validation showed the previous 1.35/1.45 defaults still read as
    // pinpoints and hairlines on a large viewport. Start at the top of the
    // current UI ranges so the default silhouette is much closer to the
    // familiar GaussView ball-and-stick weight while the MO remains dominant.
    float atom_scale = 1.80f;
    float bond_scale = 2.00f;
    float molecule_opacity = 1.00f;
    float orbital_opacity = 0.92f;
    bool show_hydrogens = true;
};

struct BondVisual {
    std::size_t atom_a = 0;
    std::size_t atom_b = 0;
    double distance_bohr = 0.0;
    bool delocalised = false;
};

[[nodiscard]] double covalent_radius_angstrom(int atomic_number) noexcept;
[[nodiscard]] std::vector<BondVisual> analyse_bonds(const Wavefunction& wavefunction);

} // namespace cov
