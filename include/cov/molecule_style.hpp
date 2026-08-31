#pragma once

#include "cov/interaction_graph.hpp"
#include "cov/model.hpp"

#include <cstddef>
#include <cstdint>
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
    // Strong non-covalent connectivity is visible by default, but it uses
    // dedicated line styles so it cannot be mistaken for an ordinary bond.
    bool show_coordination_contacts = true;
    bool show_multicentre_support = true;
    // Polyhedral cage support is a global skeletal/topological inference, not
    // a 3c electron assignment, and therefore has its own visual layer.
    bool show_polyhedral_cage_support = true;
    // Hydrogen bonds, generic non-covalent contacts and ion-pair contacts are
    // an optional context layer.  Ambiguous contacts remain hidden even when
    // this layer is enabled.
    bool show_weak_interactions = false;
};

enum class InteractionVisualStyle : std::uint8_t {
    Hidden = 0,
    OrdinaryBond,
    CoordinationDash,
    MulticentreDash,
    PolyhedralCageDash,
    HydrogenBondDots,
    NoncovalentDots,
    IonicDash,
};

// A single, testable semantic-to-visual contract used by the renderer.  In
// particular, weak or ambiguous contacts can never fall through to the
// ordinary-bond style.
[[nodiscard]] constexpr InteractionVisualStyle interaction_visual_style(
    const InteractionKind kind,
    const MoleculeRenderSettings& settings) noexcept {
    switch (kind) {
        case InteractionKind::CovalentConnectivity:
            return InteractionVisualStyle::OrdinaryBond;
        case InteractionKind::CoordinationContact:
            return settings.show_coordination_contacts
                       ? InteractionVisualStyle::CoordinationDash
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::MulticentreSupport:
            return settings.show_multicentre_support
                       ? InteractionVisualStyle::MulticentreDash
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::PolyhedralCageSupport:
            return settings.show_polyhedral_cage_support
                       ? InteractionVisualStyle::PolyhedralCageDash
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::HydrogenBond:
            return settings.show_weak_interactions
                       ? InteractionVisualStyle::HydrogenBondDots
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::NoncovalentContact:
            return settings.show_weak_interactions
                       ? InteractionVisualStyle::NoncovalentDots
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::IonicContact:
            return settings.show_weak_interactions
                       ? InteractionVisualStyle::IonicDash
                       : InteractionVisualStyle::Hidden;
        case InteractionKind::AmbiguousContact:
            return InteractionVisualStyle::Hidden;
    }
    return InteractionVisualStyle::Hidden;
}

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
