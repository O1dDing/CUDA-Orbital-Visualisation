#include "cov/molecule_style.hpp"

#include <cmath>
#include <iostream>

int main() {
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
        return 1;
    }
    for (const auto& bond : bonds) {
        if (!bond.delocalised) {
            std::cerr << "compact equal five-member ring was not marked delocalised\n";
            return 2;
        }
    }

    cov::Wavefunction saturated = wf;
    // Stretch one edge/ring geometry into a single-bond-like regime. The
    // conservative detector must stop claiming delocalisation.
    for (std::size_t i = 0; i < saturated.atoms.size(); ++i) {
        saturated.atoms[i].x *= 1.10;
        saturated.atoms[i].y *= 1.10;
    }
    const auto saturated_bonds = cov::analyse_bonds(saturated);
    bool any_delocalised = false;
    for (const auto& bond : saturated_bonds) any_delocalised = any_delocalised || bond.delocalised;
    if (any_delocalised) {
        std::cerr << "long single-bond-like ring was incorrectly marked delocalised\n";
        return 3;
    }

    std::cout << "molecule_style_smoke ok\n";
    return 0;
}
