#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <vector>

namespace cov {

enum class MoleculeStyle {
    MediumBallAndStick = 0,
    StickDelocalisation,
};

struct MoleculeRenderSettings {
    MoleculeStyle style = MoleculeStyle::MediumBallAndStick;
    // Deliberately larger than the first MVP defaults: the molecular skeleton
    // should read immediately without becoming the visual subject over the MO.
    float atom_scale = 1.35f;
    float bond_scale = 1.45f;
    float molecule_opacity = 0.96f;
    float orbital_opacity = 0.94f;
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
