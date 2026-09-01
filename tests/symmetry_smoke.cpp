#include "cov/symmetry.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

cov::Atom atom(const char* symbol, int z, double x, double y, double zc) {
    return {symbol, z, x, y, zc};
}

bool expect_group(const cov::Wavefunction& wf,
                  const std::string& expected,
                  const char* label) {
    const auto result = cov::analyse_molecular_symmetry(wf);
    if (result.point_group != expected) {
        std::cerr << label << ": expected " << expected
                  << ", got " << result.point_group
                  << " (ops=" << result.operations.size()
                  << ", tol=" << result.tolerance_bohr << ")\n";
        return false;
    }
    if (result.operations.empty()) {
        std::cerr << label << ": no validated symmetry operations retained\n";
        return false;
    }
    return true;
}

cov::Wavefunction h2() {
    cov::Wavefunction wf;
    wf.atoms = {
        atom("H", 1, 0.0, 0.0, -0.7),
        atom("H", 1, 0.0, 0.0,  0.7),
    };
    return wf;
}

cov::Wavefunction trigonal_planar() {
    cov::Wavefunction wf;
    const double r = 1.4;
    for (int i = 0; i < 3; ++i) {
        const double a = 2.0 * pi * static_cast<double>(i) / 3.0;
        wf.atoms.push_back(atom("H", 1, r * std::cos(a), r * std::sin(a), 0.0));
    }
    return wf;
}

cov::Wavefunction pentagonal_planar() {
    cov::Wavefunction wf;
    wf.atoms.push_back(atom("Fe", 26, 0.0, 0.0, 0.0));
    const double r = 2.0;
    for (int i = 0; i < 5; ++i) {
        const double a = 2.0 * pi * static_cast<double>(i) / 5.0;
        wf.atoms.push_back(atom("C", 6, r * std::cos(a), r * std::sin(a), 0.0));
    }
    return wf;
}

cov::Wavefunction tetrahedral() {
    cov::Wavefunction wf;
    wf.atoms.push_back(atom("Zn", 30, 0.0, 0.0, 0.0));
    wf.atoms.push_back(atom("Cl", 17,  1.0,  1.0,  1.0));
    wf.atoms.push_back(atom("Cl", 17,  1.0, -1.0, -1.0));
    wf.atoms.push_back(atom("Cl", 17, -1.0,  1.0, -1.0));
    wf.atoms.push_back(atom("Cl", 17, -1.0, -1.0,  1.0));
    return wf;
}

cov::Wavefunction octahedral() {
    cov::Wavefunction wf;
    wf.atoms.push_back(atom("Ti", 22, 0.0, 0.0, 0.0));
    wf.atoms.push_back(atom("F", 9,  1.8, 0.0, 0.0));
    wf.atoms.push_back(atom("F", 9, -1.8, 0.0, 0.0));
    wf.atoms.push_back(atom("F", 9, 0.0,  1.8, 0.0));
    wf.atoms.push_back(atom("F", 9, 0.0, -1.8, 0.0));
    wf.atoms.push_back(atom("F", 9, 0.0, 0.0,  1.8));
    wf.atoms.push_back(atom("F", 9, 0.0, 0.0, -1.8));
    return wf;
}

} // namespace

int main() {
    if (!expect_group(h2(), "Dinfh", "H2")) return EXIT_FAILURE;
    if (!expect_group(trigonal_planar(), "D3h", "D3h planar")) return EXIT_FAILURE;
    if (!expect_group(pentagonal_planar(), "D5h", "D5h planar")) return EXIT_FAILURE;
    if (!expect_group(tetrahedral(), "Td", "Td tetrahedron")) return EXIT_FAILURE;
    if (!expect_group(octahedral(), "Oh", "Oh octahedron")) return EXIT_FAILURE;

    cov::Wavefunction producer = tetrahedral();
    producer.point_group_detected = "PRODUCER-GROUP";
    producer.point_group_used = "PRODUCER-GROUP";
    producer.point_group_provenance = cov::DataProvenance::Producer;
    cov::derive_point_group_from_geometry(producer);
    if (producer.point_group_detected != "PRODUCER-GROUP" ||
        producer.point_group_provenance != cov::DataProvenance::Producer) {
        std::cerr << "producer symmetry provenance was overwritten\n";
        return EXIT_FAILURE;
    }

    cov::Wavefunction derived = tetrahedral();
    cov::derive_point_group_from_geometry(derived);
    if (derived.point_group_detected != "Td" ||
        derived.point_group_provenance != cov::DataProvenance::Derived) {
        std::cerr << "geometry-derived symmetry was not stored with derived provenance\n";
        return EXIT_FAILURE;
    }

    std::cout << "symmetry smoke test passed\n";
    return EXIT_SUCCESS;
}
