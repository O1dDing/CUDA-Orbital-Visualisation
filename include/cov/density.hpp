#pragma once

#include "cov/model.hpp"

#include <vector>

namespace cov {

// Common production result for both parsers and subsequent analyses. It
// distinguishes a known zero matrix from missing, invalid or failed data.
DensityReconstruction reconstruct_density(const Wavefunction& wavefunction);

// Preserve valid producer matrices, retain their failures, and optionally use
// the shared MO reconstruction when a producer matrix was absent.
void establish_density(Wavefunction& wavefunction, bool reconstruct_missing = true);

// Update trace(P S) and trace(Q S) after the independent AO metric is available.
// Matrix availability itself does not require an overlap matrix.
void update_density_metric_diagnostics(Wavefunction& wavefunction);

// Molden has no mandatory global counters. Derive integer counters only if its
// explicit spin data or declared shared integer model establishes them; real
// fractional occupation sums remain in the density diagnostics.
void derive_molden_electron_counts(Wavefunction& wavefunction);

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
