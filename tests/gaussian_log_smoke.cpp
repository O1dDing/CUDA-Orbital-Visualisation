#include "cov/gaussian_log.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main() {
    const auto path = std::filesystem::temp_directory_path() /
                      "cov_gaussian_log_enrichment.log";
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) return EXIT_FAILURE;
        out << " Full point group                 C2V     NOp   4\n"
               " Largest Abelian subgroup         C2V     NOp   4\n"
               " Initial guess orbital symmetries:\n"
               "       Occupied  (A1) (A1) (A1) (A1) (A1)\n"
               "       Virtual   (A1) (A1)\n"
               " The electronic state of the initial guess is 1-A1.\n"
               " Orbital symmetries:\n"
               "       Occupied  (A1) (A1) (B2) (A1) (B1)\n"
               "       Virtual   (A1) (B2)\n"
               " The electronic state is 1-A1.\n";
    }

    cov::Wavefunction wf;
    wf.basis_count = 7;
    wf.orbitals.resize(7);
    for (auto& mo : wf.orbitals) mo.spin = cov::Spin::Alpha;

    const auto result = cov::enrich_from_gaussian_log(wf, path);
    std::error_code ec;
    std::filesystem::remove(path, ec);

    const std::vector<std::string> expected = {
        "A1", "A1", "B2", "A1", "B1", "A1", "B2"
    };
    if (!result.opened || result.symmetry_labels_applied != 7u ||
        !result.point_group_applied || wf.point_group_detected != "C2V" ||
        wf.point_group_used != "C2V" ||
        wf.point_group_provenance != cov::DataProvenance::Producer) {
        std::cerr << "Gaussian LOG point-group/symmetry metadata regression\n";
        return EXIT_FAILURE;
    }

    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (wf.orbitals[i].symmetry != expected[i] ||
            wf.orbitals[i].symmetry_provenance != cov::DataProvenance::Producer) {
            std::cerr << "Gaussian LOG orbital symmetry ordering regression at " << i << '\n';
            return EXIT_FAILURE;
        }
    }

    std::cout << "Gaussian LOG enrichment smoke test passed\n";
    return EXIT_SUCCESS;
}
