#pragma once

#include "cov/model.hpp"

namespace cov {

struct OrbitalChemistryOptions {
    double degeneracy_tolerance_hartree = 1.0e-5;
    double pair_distance_ceiling_bohr = 6.0;
    double pair_mayer_floor = 0.02;
    double pair_atom_weight_floor = 0.08;
    double overlap_bonding_threshold = 0.015;
    double nonbonding_overlap_threshold = 0.004;
    double determined_fraction = 0.82;
    double maximum_undetermined_for_percentages = 0.55;
};

// Builds a molecule-independent, element/ECP-aware minimal chemical-valence
// reference inside the working AO basis, projects every canonical MO onto that
// reference with the AO-overlap metric, selects the canonical valence manifold,
// and derives per-interaction chemistry that FCHK can support.
//
// This is COV's own FCHK-derived reference analysis. It is not called IAO,
// NAO, NBO, Wiberg or producer data.
void derive_orbital_chemistry(
    Wavefunction& wavefunction,
    const OrbitalChemistryOptions& options = {});

[[nodiscard]] const char* chemistry_status_name(ChemistryStatus value) noexcept;
[[nodiscard]] const char* orbital_angular_family_name(
    OrbitalAngularFamily value) noexcept;
[[nodiscard]] const char* orbital_bonding_role_name(
    OrbitalBondingRole value) noexcept;

} // namespace cov
