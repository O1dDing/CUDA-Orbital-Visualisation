#include "cov/molden_parser.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>

namespace {

bool test_no_sym_mo_boundaries() {
    const auto path = std::filesystem::temp_directory_path() /
                      "cov_molden_no_sym_regression.molden";
    {
        std::ofstream out(path, std::ios::binary);
        if (!out) {
            std::cerr << "unable to create temporary no-Sym Molden fixture\n";
            return false;
        }
        out << "[Molden Format]\n"
               "[Atoms] Angs\n"
               "H 1 1 0.0 0.0 0.0\n"
               "H 2 1 0.0 0.0 0.74\n"
               "[GTO]\n"
               "1 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "2 0\n"
               "s 1 1.0\n"
               "1.0 1.0\n"
               "[MO]\n"
               "Ene= -0.500000\n"
               "Spin= Alpha\n"
               "Occup= 2.0\n"
               "1 0.70710678\n"
               "2 0.70710678\n"
               "Ene= 0.200000\n"
               "Spin= Alpha\n"
               "Occup= 0.0\n"
               "1 0.70710678\n"
               "2 -0.70710678\n";
    }

    try {
        const auto wf = cov::parse_molden(path);
        std::error_code ec;
        std::filesystem::remove(path, ec);
        if (wf.orbitals.size() != 2 || wf.basis_count != 2) {
            std::cerr << "no-Sym MO boundary regression: expected 2 MOs / 2 basis, got "
                      << wf.orbitals.size() << " / " << wf.basis_count << '\n';
            return false;
        }
        if (wf.orbitals[0].occupation != 2.0f || wf.orbitals[1].occupation != 0.0f) {
            std::cerr << "no-Sym MO metadata regression\n";
            return false;
        }
        return true;
    } catch (const std::exception& e) {
        std::error_code ec;
        std::filesystem::remove(path, ec);
        std::cerr << "no-Sym MO parsing failed: " << e.what() << '\n';
        return false;
    }
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: cov_parser_smoke <file.molden> [--inspect]\n";
        return EXIT_FAILURE;
    }

    try {
        const auto wf = cov::parse_molden(argv[1]);
        std::cout << "atoms=" << wf.atoms.size()
                  << " shells=" << wf.shells.size()
                  << " basis=" << wf.basis_count
                  << " orbitals=" << wf.orbitals.size() << '\n';
        if (argc == 3 && std::string_view(argv[2]) == "--inspect") {
            return EXIT_SUCCESS;
        }
        if (wf.atoms.size() != 2) {
            std::cerr << "expected 2 atoms, got " << wf.atoms.size() << '\n';
            return EXIT_FAILURE;
        }
        if (wf.basis_count != 2) {
            std::cerr << "expected 2 basis functions, got " << wf.basis_count << '\n';
            return EXIT_FAILURE;
        }
        if (wf.orbitals.size() != 2) {
            std::cerr << "expected 2 orbitals, got " << wf.orbitals.size() << '\n';
            return EXIT_FAILURE;
        }
        if (wf.orbitals[0].coefficients.size() != wf.basis_count) {
            std::cerr << "MO coefficient consistency check failed\n";
            return EXIT_FAILURE;
        }
        if (!test_no_sym_mo_boundaries()) {
            return EXIT_FAILURE;
        }
        std::cout << "parser smoke test passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
