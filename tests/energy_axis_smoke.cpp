#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
    // Real crowded valence energies from the uploaded Cp- Cartesian Gaussian
    // reference. Distinct near-neighbours deliberately coexist with exact /
    // tolerance-level degeneracies.
    constexpr std::array<double, 17> cp_energies{
        -0.399340, -0.254871, -0.254868, -0.237168, -0.237168,
        -0.236093, -0.068048, -0.068048, 0.133376, 0.133376,
         0.137838,  0.138589,  0.138590,  0.224213, 0.224222,
         0.225642,  0.236579
    };

    std::vector<double> distinct(cp_energies.begin(), cp_energies.end());
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end(), [](double a, double b) {
        return std::abs(a - b) <= 1.0e-12;
    }), distinct.end());

    const cov::EnergyTransform focus = cov::build_energy_transform(
        distinct, cov::EnergyAxisMode::NonlinearFocus, 0.055);
    const cov::EnergyTransform linear = cov::build_energy_transform(
        distinct, cov::EnergyAxisMode::Linear, 0.055);

    double previous = cov::energy_display_coordinate(distinct.front(), focus);
    double min_gap = 1.0;
    for (std::size_t i = 1; i < distinct.size(); ++i) {
        const double current = cov::energy_display_coordinate(distinct[i], focus);
        if (!(current > previous)) {
            std::cerr << "adaptive transform is not strictly monotone\n";
            return 1;
        }
        min_gap = std::min(min_gap, current - previous);
        previous = current;
    }

    if (min_gap < 0.050) {
        std::cerr << "crowded nondegenerate levels are still too close: " << min_gap << '\n';
        return 2;
    }

    // On the default export's ~720 px drawable axis, the new floor should
    // give at least ~36 px between distinct crowded levels -- roughly double
    // the previous ~20 px reference.
    if (min_gap * 720.0 < 36.0) {
        std::cerr << "pixel-equivalent spacing target failed: " << min_gap * 720.0 << '\n';
        return 3;
    }

    const double degenerate_energy = -0.068048;
    if (cov::energy_display_coordinate(degenerate_energy, focus) !=
        cov::energy_display_coordinate(degenerate_energy, focus)) {
        std::cerr << "degenerate source energies did not retain equal coordinates\n";
        return 4;
    }

    for (const double energy : distinct) {
        const double coordinate = cov::energy_display_coordinate(energy, focus);
        const double inverse = cov::energy_from_display_coordinate(coordinate, focus);
        if (std::abs(inverse - energy) > 1.0e-11 * std::max(1.0, std::abs(energy))) {
            std::cerr << "adaptive inverse/tick consistency failed\n";
            return 5;
        }
        if (std::abs(cov::energy_display_coordinate(energy, linear) - energy) > 1.0e-14) {
            std::cerr << "linear identity transform changed energy\n";
            return 6;
        }
        for (const cov::EnergyUnit unit : {
                 cov::EnergyUnit::Hartree, cov::EnergyUnit::ElectronVolt,
                 cov::EnergyUnit::JoulePerMol, cov::EnergyUnit::KilojoulePerMol,
                 cov::EnergyUnit::CaloriePerMol, cov::EnergyUnit::KilocaloriePerMol}) {
            if (!std::isfinite(cov::convert_hartree(inverse, unit))) {
                std::cerr << "unit-switch stability failed\n";
                return 7;
            }
        }
    }

    if (std::string(cov::energy_transform_name(cov::EnergyAxisMode::NonlinearFocus)) !=
        "adaptive-log-gap-v3") {
        std::cerr << "unexpected adaptive transform name\n";
        return 8;
    }

    std::cout << "energy_axis_smoke ok; minimum normalized gap=" << min_gap
              << "; pixel-equivalent=" << min_gap * 720.0 << " px\n";
    return 0;
}
