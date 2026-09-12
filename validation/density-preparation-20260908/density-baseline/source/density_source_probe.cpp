#include "cov/density.hpp"
#include "cov/numerical_diagnostics.hpp"
#include "cov/wavefunction_io.hpp"

#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

// Instrument the actual production parser/postprocessor and reconstruction.
// There is no alternate density algorithm or expected result in this program.
namespace {
template<class Range> void numbers(const Range& values) {
    std::cout << '[';
    bool first = true;
    for (const auto value : values) {
        if (!first) std::cout << ',';
        first = false;
        cov::numerical_json::number(std::cout, value);
    }
    std::cout << ']';
}

template<class Calculate> void reconstructed(Calculate calculate) {
    try {
        const auto values = calculate();
        std::cout << "{\"returned\":true,\"packed\":";
        numbers(values);
        std::cout << '}';
    } catch (const std::exception& error) {
        std::cout << "{\"returned\":false,\"error\":";
        cov::numerical_json::string(std::cout, error.what());
        std::cout << '}';
    }
}
}

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    try {
        cov::WavefunctionParseOptions options;
        options.max_atoms = 1000;
        options.keep_density = true;
        options.reconstruct_density_if_missing = true;
        const auto wf = cov::parse_wavefunction(argv[1], options);
        std::cout << std::setprecision(17)
                  << "{\"parse_status\":\"returned\",\"nbasis\":" << wf.basis_count
                  << ",\"source\":" << static_cast<int>(wf.source)
                  << ",\"alpha_electrons\":" << wf.alpha_electrons
                  << ",\"beta_electrons\":" << wf.beta_electrons
                  << ",\"electron_counts_provenance\":" << static_cast<int>(wf.electron_counts_provenance)
                  << ",\"total_density_provenance\":" << static_cast<int>(wf.total_density_provenance)
                  << ",\"spin_density_provenance\":" << static_cast<int>(wf.spin_density_provenance)
                  << ",\"total_density_packed\":";
        numbers(wf.total_density_packed);
        std::cout << ",\"spin_density_packed\":";
        numbers(wf.spin_density_packed);
        std::cout << ",\"overlap\":";
        numbers(wf.ao_overlap);
        std::cout << ",\"numerical_diagnostics\":";
        cov::write_numerical_diagnostics_json(std::cout, wf);
        std::cout << ",\"bond_order_provenance\":" << static_cast<int>(wf.bond_order_provenance)
                  << ",\"bond_orders\":[";
        bool first = true;
        for (const auto& bond : wf.bond_orders) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"atoms\":[" << bond.atom_a << ',' << bond.atom_b << "],\"mayer\":";
            cov::numerical_json::number(std::cout, bond.mayer_order);
            std::cout << ",\"provenance\":" << static_cast<int>(bond.provenance) << '}';
        }
        std::cout << "],\"orbitals\":[";
        first = true;
        for (const auto& mo : wf.orbitals) {
            if (!first) std::cout << ',';
            first = false;
            std::cout << "{\"occupation\":";
            cov::numerical_json::number(std::cout, mo.occupation);
            std::cout << ",\"occupation_provenance\":" << static_cast<int>(mo.occupation_provenance)
                      << ",\"spin\":" << static_cast<int>(mo.spin) << ",\"coefficients\":";
            numbers(mo.coefficients);
            std::cout << '}';
        }
        std::cout << "],\"production_reconstructed_total\":";
        reconstructed([&] { return cov::reconstruct_total_density_packed(wf); });
        std::cout << ",\"production_reconstructed_spin\":";
        reconstructed([&] { return cov::reconstruct_spin_density_packed(wf); });
        std::cout << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cout << "{\"parse_status\":\"rejected\",\"error\":";
        cov::numerical_json::string(std::cout, error.what());
        std::cout << "}\n";
        return 2;
    }
}
