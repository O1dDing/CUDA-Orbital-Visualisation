#include "cov/molden_parser.hpp"

#include <cstdlib>
#include <iostream>
#include <string_view>

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
        std::cout << "parser smoke test passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
