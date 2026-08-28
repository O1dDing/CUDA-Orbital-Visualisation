#pragma once

#include "cov/model.hpp"

#include <vector>

namespace cov {

// Lower-triangular packed AO density, index(i,j) = i*(i+1)/2 + j for i >= j.
// This is the ordinary AO density P = C n C^T reconstructed from the orbitals
// already stored in Wavefunction. It intentionally does not require an overlap
// matrix; S enters population/bond-index analyses later, not this reconstruction.
std::vector<double> reconstruct_total_density_packed(const Wavefunction& wavefunction);

} // namespace cov
