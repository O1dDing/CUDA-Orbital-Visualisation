#include "cov/mo_diagram.hpp"

#include <array>
#include <cmath>
#include <iostream>

int main() {
    cov::EnergyTransform focus{cov::EnergyAxisMode::NonlinearFocus, -0.25, 0.01};
    constexpr std::array<double, 7> energies{-2.0, -0.3, -0.25001, -0.25,
                                             -0.24999, 0.1, 3.0};
    double previous = cov::energy_display_coordinate(energies.front(), focus);
    for (std::size_t i = 1; i < energies.size(); ++i) {
        const double current = cov::energy_display_coordinate(energies[i], focus);
        if (!(current > previous)) {
            std::cerr << "nonlinear transform is not strictly monotone\n";
            return 1;
        }
        previous = current;
    }

    const double degenerate = cov::energy_display_coordinate(-0.125, focus);
    if (degenerate != cov::energy_display_coordinate(-0.125, focus)) {
        std::cerr << "equal source energies did not retain equal coordinates\n";
        return 2;
    }

    const double near_gap = cov::energy_display_coordinate(-0.249, focus) -
                            cov::energy_display_coordinate(-0.25, focus);
    cov::EnergyTransform linear{cov::EnergyAxisMode::Linear, -0.25, 0.01};
    const double linear_gap = cov::energy_display_coordinate(-0.249, linear) -
                             cov::energy_display_coordinate(-0.25, linear);
    if (!(near_gap > linear_gap * 50.0)) {
        std::cerr << "focus transform did not separate near nondegenerate levels\n";
        return 3;
    }

    for (const double energy : energies) {
        const double coordinate = cov::energy_display_coordinate(energy, focus);
        const double inverse = cov::energy_from_display_coordinate(coordinate, focus);
        if (std::abs(inverse - energy) > 1.0e-12 * std::max(1.0, std::abs(energy))) {
            std::cerr << "focus inverse/tick consistency failed\n";
            return 4;
        }
        for (const cov::EnergyUnit unit : {cov::EnergyUnit::Hartree,
                                           cov::EnergyUnit::ElectronVolt,
                                           cov::EnergyUnit::JoulePerMol,
                                           cov::EnergyUnit::KilojoulePerMol,
                                           cov::EnergyUnit::CaloriePerMol,
                                           cov::EnergyUnit::KilocaloriePerMol}) {
            if (!std::isfinite(cov::convert_hartree(inverse, unit))) {
                std::cerr << "unit-switch stability failed\n";
                return 5;
            }
        }
    }
    std::cout << "energy_axis_smoke ok\n";
    return 0;
}
