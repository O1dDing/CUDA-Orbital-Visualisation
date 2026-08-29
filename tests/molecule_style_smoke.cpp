#include "cov/molecule_style.hpp"
#include "cov/wavefunction_io.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <utility>
#include <vector>

namespace {

void add_all_pair_mayer_records(cov::Wavefunction& wf, const double order) {
    wf.bond_order_provenance = cov::DataProvenance::Derived;
    for (std::size_t i = 0; i < wf.atoms.size(); ++i) {
        for (std::size_t j = i + 1; j < wf.atoms.size(); ++j) {
            wf.bond_orders.push_back({
                static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j),
                order, cov::DataProvenance::Derived});
        }
    }
}

bool check_carbon_graph(const cov::Wavefunction& wf,
                        const std::vector<cov::BondVisual>& bonds,
                        const std::size_t expected_edges,
                        const std::vector<int>& expected_degrees) {
    if (bonds.size() != expected_edges || expected_degrees.size() != wf.atoms.size()) {
        return false;
    }
    std::vector<int> degrees(wf.atoms.size(), 0);
    for (const auto& bond : bonds) {
        if (bond.atom_a >= degrees.size() || bond.atom_b >= degrees.size()) return false;
        ++degrees[bond.atom_a];
        ++degrees[bond.atom_b];
    }
    return degrees == expected_degrees;
}

bool check_real_topology(const std::filesystem::path& path,
                         const std::size_t expected_bonds,
                         const std::size_t expected_carbon_bonds) {
    const auto wf = cov::parse_wavefunction(path);
    const auto bonds = cov::analyse_bonds(wf);
    std::size_t carbon_bonds = 0;
    for (const auto& bond : bonds) {
        const auto& a = wf.atoms[bond.atom_a];
        const auto& b = wf.atoms[bond.atom_b];
        if (a.atomic_number == 6 && b.atomic_number == 6) ++carbon_bonds;
        const double radii_bohr =
            (cov::covalent_radius_angstrom(a.atomic_number) +
             cov::covalent_radius_angstrom(b.atomic_number)) * cov::kAngstromToBohr;
        if (bond.distance_bohr > 1.50 * radii_bohr + 1.0e-10) return false;
    }
    return bonds.size() == expected_bonds && carbon_bonds == expected_carbon_bonds;
}

} // namespace

int main(const int argc, char** argv) {
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

    // Delocalised Mayer coupling is not structural adjacency. Even when every
    // carbon pair carries a small positive index, benzene must remain a C6
    // perimeter rather than a complete graph with centre-crossing chords.
    cov::Wavefunction benzene;
    constexpr double benzene_radius_angstrom = side_angstrom;
    for (int i = 0; i < 6; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 6.0;
        benzene.atoms.push_back({
            "C", 6,
            benzene_radius_angstrom * std::cos(angle) * cov::kAngstromToBohr,
            benzene_radius_angstrom * std::sin(angle) * cov::kAngstromToBohr,
            0.0});
    }
    add_all_pair_mayer_records(benzene, 0.06);
    const auto benzene_bonds = cov::analyse_bonds(benzene);
    if (!check_carbon_graph(benzene, benzene_bonds, 6u,
                            std::vector<int>(6u, 2))) {
        std::cerr << "benzene electronic coupling created non-neighbour chords\n";
        return 7;
    }
    for (const auto& bond : benzene_bonds) {
        if (!bond.delocalised) {
            std::cerr << "benzene perimeter was not marked delocalised\n";
            return 8;
        }
    }

    // Two regular fused six-membered rings: the carbon skeleton has exactly
    // eleven edges. Only the two shared-edge atoms have degree three.
    cov::Wavefunction naphthalene;
    const double root3 = std::sqrt(3.0);
    const std::vector<std::pair<double, double>> naphthalene_xy = {
        {0.0, 0.5 * side_angstrom},
        {0.0, -0.5 * side_angstrom},
        {-0.5 * root3 * side_angstrom, side_angstrom},
        {-root3 * side_angstrom, 0.5 * side_angstrom},
        {-root3 * side_angstrom, -0.5 * side_angstrom},
        {-0.5 * root3 * side_angstrom, -side_angstrom},
        {0.5 * root3 * side_angstrom, side_angstrom},
        {root3 * side_angstrom, 0.5 * side_angstrom},
        {root3 * side_angstrom, -0.5 * side_angstrom},
        {0.5 * root3 * side_angstrom, -side_angstrom},
    };
    for (const auto [x_angstrom, y_angstrom] : naphthalene_xy) {
        naphthalene.atoms.push_back({
            "C", 6, x_angstrom * cov::kAngstromToBohr,
            y_angstrom * cov::kAngstromToBohr, 0.0});
    }
    add_all_pair_mayer_records(naphthalene, 0.06);
    const auto naphthalene_bonds = cov::analyse_bonds(naphthalene);
    const std::vector<int> naphthalene_degrees = {3, 3, 2, 2, 2, 2, 2, 2, 2, 2};
    if (!check_carbon_graph(naphthalene, naphthalene_bonds, 11u,
                            naphthalene_degrees)) {
        std::cerr << "naphthalene electronic coupling created non-neighbour chords\n";
        return 9;
    }
    for (const auto& bond : naphthalene_bonds) {
        if (!bond.delocalised) {
            std::cerr << "naphthalene fused-ring perimeter was not marked delocalised\n";
            return 10;
        }
    }

    if (argc == 3) {
        if (!check_real_topology(argv[1], 12u, 6u)) {
            std::cerr << "real benzene topology regression\n";
            return 11;
        }
        if (!check_real_topology(argv[2], 19u, 11u)) {
            std::cerr << "real naphthalene topology regression\n";
            return 12;
        }
    }

    std::cout << "molecule_style_smoke ok\n";
    return 0;
}
