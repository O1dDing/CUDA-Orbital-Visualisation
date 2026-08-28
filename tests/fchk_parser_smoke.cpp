#include "cov/fchk_parser.hpp"
#include "cov/wavefunction_io.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

void field_prefix(std::ofstream& out, const std::string& label, const char type) {
    out << std::left << std::setw(40) << label << type;
}

void scalar_i(std::ofstream& out, const std::string& label, const long long value) {
    field_prefix(out, label, 'I');
    out << std::right << std::setw(16) << value << '\n';
}

void array_i(std::ofstream& out,
             const std::string& label,
             const std::vector<long long>& values) {
    field_prefix(out, label, 'I');
    out << "   N=" << std::right << std::setw(12) << values.size() << '\n';
    for (std::size_t i = 0; i < values.size(); ++i) {
        out << std::setw(12) << values[i];
        if ((i + 1u) % 6u == 0u || i + 1u == values.size()) out << '\n';
    }
}

void array_r(std::ofstream& out,
             const std::string& label,
             const std::vector<double>& values) {
    field_prefix(out, label, 'R');
    out << "   N=" << std::right << std::setw(12) << values.size() << '\n';
    out << std::uppercase << std::scientific << std::setprecision(8);
    for (std::size_t i = 0; i < values.size(); ++i) {
        out << std::setw(16) << values[i];
        if ((i + 1u) % 5u == 0u || i + 1u == values.size()) out << '\n';
    }
}

std::filesystem::path make_restricted_h2_fixture() {
    const auto path = std::filesystem::temp_directory_path() / "cov_h2_restricted.fchk";
    std::ofstream out(path, std::ios::binary);
    out << "H2 FCHK parser regression\n"
           "SP        RHF STO-3G\n";
    scalar_i(out, "Number of atoms", 2);
    scalar_i(out, "Number of alpha electrons", 1);
    scalar_i(out, "Number of beta electrons", 1);
    scalar_i(out, "Number of basis functions", 2);
    scalar_i(out, "Number of independent functions", 2);
    array_i(out, "Atomic numbers", {1, 1});
    array_r(out, "Current cartesian coordinates", {0.0, 0.0, -0.7, 0.0, 0.0, 0.7});
    array_i(out, "Shell types", {0, 0});
    array_i(out, "Number of primitives per shell", {1, 1});
    array_i(out, "Shell to atom map", {1, 2});
    array_r(out, "Primitive exponents", {1.0, 1.0});
    array_r(out, "Contraction coefficients", {1.0, 1.0});
    array_r(out, "Alpha Orbital Energies", {-0.5, 0.2});
    array_r(out, "Alpha MO coefficients", {
        0.7071067811865475, 0.7071067811865475,
        0.7071067811865475, -0.7071067811865475
    });
    return path;
}

std::filesystem::path make_unrestricted_h_fixture() {
    const auto path = std::filesystem::temp_directory_path() / "cov_h_unrestricted.fchk";
    std::ofstream out(path, std::ios::binary);
    out << "H unrestricted FCHK parser regression\n"
           "SP        UHF STO-3G\n";
    scalar_i(out, "Number of atoms", 1);
    scalar_i(out, "Number of alpha electrons", 1);
    scalar_i(out, "Number of beta electrons", 0);
    scalar_i(out, "Number of basis functions", 1);
    scalar_i(out, "Number of independent functions", 1);
    array_i(out, "Atomic numbers", {1});
    array_r(out, "Current cartesian coordinates", {0.0, 0.0, 0.0});
    array_i(out, "Shell types", {0});
    array_i(out, "Number of primitives per shell", {1});
    array_i(out, "Shell to atom map", {1});
    array_r(out, "Primitive exponents", {1.0});
    array_r(out, "Contraction coefficients", {1.0});
    array_r(out, "Alpha Orbital Energies", {-0.4});
    array_r(out, "Alpha MO coefficients", {1.0});
    array_r(out, "Beta Orbital Energies", {-0.3});
    array_r(out, "Beta MO coefficients", {1.0});
    return path;
}

std::filesystem::path make_shell_convention_fixture() {
    const auto path = std::filesystem::temp_directory_path() / "cov_shell_conventions.fchk";
    std::ofstream out(path, std::ios::binary);
    out << "FCHK shell convention regression\n"
           "SP        RHF synthetic\n";
    scalar_i(out, "Number of atoms", 1);
    scalar_i(out, "Number of alpha electrons", 0);
    scalar_i(out, "Number of beta electrons", 0);
    // SP contributes 1+3, then pure 5D+7F+9G => 25 total.
    scalar_i(out, "Number of basis functions", 25);
    scalar_i(out, "Number of independent functions", 25);
    array_i(out, "Atomic numbers", {6});
    array_r(out, "Current cartesian coordinates", {0.0, 0.0, 0.0});
    array_i(out, "Shell types", {-1, -2, -3, -4});
    array_i(out, "Number of primitives per shell", {1, 1, 1, 1});
    array_i(out, "Shell to atom map", {1, 1, 1, 1});
    array_r(out, "Primitive exponents", {4.0, 3.0, 2.0, 1.0});
    array_r(out, "Contraction coefficients", {1.0, 1.0, 1.0, 1.0});
    array_r(out, "P(S=P) Contraction coefficients", {0.5, 0.0, 0.0, 0.0});
    return path;
}

std::filesystem::path make_cartesian_g_fixture() {
    const auto path = std::filesystem::temp_directory_path() / "cov_cartesian_g_order.fchk";
    std::ofstream out(path, std::ios::binary);
    out << "Cartesian g ordering regression\n"
           "SP        ROHF synthetic\n";
    scalar_i(out, "Number of atoms", 1);
    scalar_i(out, "Number of alpha electrons", 1);
    scalar_i(out, "Number of beta electrons", 0);
    scalar_i(out, "Number of basis functions", 15);
    scalar_i(out, "Number of independent functions", 15);
    array_i(out, "Atomic numbers", {6});
    array_r(out, "Current cartesian coordinates", {0.0, 0.0, 0.0});
    array_i(out, "Shell types", {4});
    array_i(out, "Number of primitives per shell", {1});
    array_i(out, "Shell to atom map", {1});
    array_r(out, "Primitive exponents", {1.0});
    array_r(out, "Contraction coefficients", {1.0});
    array_r(out, "Alpha Orbital Energies", {-0.1});

    std::vector<double> coefficients(15);
    for (std::size_t i = 0; i < coefficients.size(); ++i) coefficients[i] = static_cast<double>(i + 1u);
    array_r(out, "Alpha MO coefficients", coefficients);

    std::vector<double> density;
    density.reserve(120);
    for (std::size_t i = 0; i < 15u; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            density.push_back(static_cast<double>(100u * i + j));
        }
    }
    array_r(out, "Total SCF Density", density);
    return path;
}

bool approximately(const double a, const double b, const double eps = 1.0e-6) {
    return std::abs(a - b) <= eps;
}

} // namespace

int main() {
    std::error_code ec;
    const auto h2 = make_restricted_h2_fixture();
    const auto h = make_unrestricted_h_fixture();
    const auto shells = make_shell_convention_fixture();
    const auto cart_g = make_cartesian_g_fixture();

    try {
        const auto restricted = cov::parse_fchk(h2);
        if (restricted.source != cov::WavefunctionSource::Fchk ||
            restricted.atoms.size() != 2u || restricted.basis_count != 2u ||
            restricted.orbitals.size() != 2u ||
            restricted.alpha_electrons != 1u || restricted.beta_electrons != 1u) {
            std::cerr << "restricted FCHK metadata regression\n";
            return EXIT_FAILURE;
        }
        if (restricted.orbitals[0].occupation != 2.0f ||
            restricted.orbitals[1].occupation != 0.0f) {
            std::cerr << "restricted occupation reconstruction regression\n";
            return EXIT_FAILURE;
        }
        if (restricted.total_density_provenance != cov::DataProvenance::Derived ||
            restricted.total_density_packed.size() != 3u ||
            !approximately(restricted.total_density_packed[0], 1.0) ||
            !approximately(restricted.total_density_packed[1], 1.0) ||
            !approximately(restricted.total_density_packed[2], 1.0)) {
            std::cerr << "derived AO density regression\n";
            return EXIT_FAILURE;
        }

        const auto unified = cov::parse_wavefunction(h2);
        if (unified.source != cov::WavefunctionSource::Fchk || unified.basis_count != 2u) {
            std::cerr << "unified FCHK dispatch regression\n";
            return EXIT_FAILURE;
        }

        const auto unrestricted = cov::parse_fchk(h);
        if (unrestricted.orbitals.size() != 2u ||
            unrestricted.orbitals[0].spin != cov::Spin::Alpha ||
            unrestricted.orbitals[1].spin != cov::Spin::Beta ||
            unrestricted.orbitals[0].occupation != 1.0f ||
            unrestricted.orbitals[1].occupation != 0.0f) {
            std::cerr << "unrestricted alpha/beta orbital regression\n";
            return EXIT_FAILURE;
        }

        cov::FchkParseOptions shell_options;
        shell_options.require_orbitals = false;
        shell_options.keep_density = false;
        const auto convention = cov::parse_fchk(shells, shell_options);
        if (convention.basis_count != 25u || convention.shells.size() != 5u ||
            !convention.pure_d || !convention.pure_f || !convention.pure_g) {
            std::cerr << "SP / pure DFG shell convention regression: basis="
                      << convention.basis_count << " shells=" << convention.shells.size() << '\n';
            return EXIT_FAILURE;
        }
        if (convention.shells[0].angular_momentum != 0u ||
            convention.shells[1].angular_momentum != 1u ||
            convention.shells[2].angular_momentum != 2u ||
            convention.shells[3].angular_momentum != 3u ||
            convention.shells[4].angular_momentum != 4u) {
            std::cerr << "FCHK shell ordering regression\n";
            return EXIT_FAILURE;
        }

        const auto g = cov::parse_fchk(cart_g);
        const std::vector<float> expected_g{
            15.0f, 5.0f, 1.0f, 14.0f, 13.0f, 9.0f, 4.0f, 6.0f, 2.0f,
            12.0f, 10.0f, 3.0f, 11.0f, 8.0f, 7.0f
        };
        if (g.orbitals.size() != 1u || g.orbitals[0].coefficients != expected_g) {
            std::cerr << "Cartesian g coefficient permutation regression\n";
            return EXIT_FAILURE;
        }
        if (g.total_density_provenance != cov::DataProvenance::Producer ||
            g.total_density_packed.size() != 120u ||
            !approximately(g.total_density_packed[0], 1414.0) ||
            !approximately(g.total_density_packed[1], 1404.0) ||
            !approximately(g.total_density_packed[2], 404.0)) {
            std::cerr << "Cartesian g packed-density permutation regression\n";
            return EXIT_FAILURE;
        }

        std::filesystem::remove(h2, ec);
        std::filesystem::remove(h, ec);
        std::filesystem::remove(shells, ec);
        std::filesystem::remove(cart_g, ec);
        std::cout << "fchk parser smoke test passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::filesystem::remove(h2, ec);
        std::filesystem::remove(h, ec);
        std::filesystem::remove(shells, ec);
        std::filesystem::remove(cart_g, ec);
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
