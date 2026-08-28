#include "cov/molecule_style.hpp"

#include <cmath>
#include <iostream>

int main() {
    const cov::MoleculeRenderSettings defaults;
    if (std::abs(defaults.atom_scale - 1.0f) > 1.0e-6f ||
        std::abs(defaults.bond_scale - 1.0f) > 1.0e-6f ||
        defaults.molecule_opacity < 0.95f) {
        std::cerr << "1.0 molecular-size reference defaults regressed\n";
        return 1;
    }

    cov::Wavefunction wf;
    constexpr double pi = 3.14159265358979323846;
    constexpr double side_angstrom = 1.40;
    const double radius_angstrom = side_angstrom / (2.0 * std::sin(pi / 5.0));

    for (int i = 0; i < 5; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 5.0;
        cov::Atom atom;
        atom.symbol = "C";
        atom.atomic_number = 6;
        atom.x = radius_angstrom * std::cos(angle) * cov::kAngstromToBohr;
        atom.y = radius_angstrom * std::sin(angle) * cov::kAngstromToBohr;
        atom.z = 0.0;
        wf.atoms.push_back(atom);
    }

    const auto bonds = cov::analyse_bonds(wf);
    if (bonds.size() != 5) {
        std::cerr << "expected five ring bonds, got " << bonds.size() << '\n';
        return 2;
    }
    for (const auto& bond : bonds) {
        if (!bond.delocalised) {
            std::cerr << "compact equal five-member ring was not marked delocalised\n";
            return 3;
        }
    }

    cov::Wavefunction saturated = wf;
    // Stretch the ring into a single-bond-like regime. The conservative
    // detector must stop claiming delocalisation.
    for (std::size_t i = 0; i < saturated.atoms.size(); ++i) {
        saturated.atoms[i].x *= 1.10;
        saturated.atoms[i].y *= 1.10;
    }
    const auto saturated_bonds = cov::analyse_bonds(saturated);
    bool any_delocalised = false;
    for (const auto& bond : saturated_bonds) any_delocalised = any_delocalised || bond.delocalised;
    if (any_delocalised) {
        std::cerr << "long single-bond-like ring was incorrectly marked delocalised\n";
        return 4;
    }

    // Once pairwise electronic evidence is available it must outrank the old
    // 1.25*r_cov geometry rule. This synthetic C--C separation lies outside the
    // geometry fallback cutoff but inside the deliberately generous electronic
    // sanity distance; a positive Mayer index should therefore retain the bond.
    cov::Wavefunction electronic;
    electronic.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"C", 6, 2.25 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    electronic.bond_order_provenance = cov::DataProvenance::Derived;
    electronic.bond_orders.push_back({0, 1, 0.18, cov::DataProvenance::Derived});
    const auto electronic_bonds = cov::analyse_bonds(electronic);
    if (electronic_bonds.size() != 1u ||
        electronic_bonds[0].provenance != cov::DataProvenance::Derived ||
        std::abs(electronic_bonds[0].bond_order - 0.18) > 1.0e-12) {
        std::cerr << "Mayer connectivity did not outrank geometry fallback\n";
        return 5;
    }

    // Conversely, electronic availability with no supported pair must not
    // silently re-enable a geometry bond and manufacture connectivity.
    cov::Wavefunction electronic_none;
    electronic_none.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"C", 6, 1.40 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    electronic_none.bond_order_provenance = cov::DataProvenance::Derived;
    if (!cov::analyse_bonds(electronic_none).empty()) {
        std::cerr << "geometry fallback overrode an available electronic analysis\n";
        return 6;
    }

    std::cout << "molecule_style_smoke ok\n";
    return 0;
}
