#pragma once

#include "cov/model.hpp"

#include <vector>

namespace cov {

// Lower-triangular packed AO density, index(i,j) = i*(i+1)/2 + j for i >= j.
// This is the ordinary spin-summed AO density P = C n C^T reconstructed from
// the orbitals already stored in Wavefunction. It intentionally does not
// require an overlap matrix; S enters population/bond-index analyses later.
std::vector<double> reconstruct_total_density_packed(const Wavefunction& wavefunction);

// Reconstruct P_alpha - P_beta when the orbital model carries enough spin/
// occupation information. Restricted open-shell files are handled from the
// alpha/beta electron counts and canonical orbital ordering; unrestricted files
// use their explicit alpha/beta orbital blocks.
std::vector<double> reconstruct_spin_density_packed(const Wavefunction& wavefunction);

} // namespace cov
