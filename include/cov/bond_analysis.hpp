#pragma once

#include "cov/model.hpp"

namespace cov {

struct BondAnalysisOptions {
    // Retain weak positive/negative pairwise indices for diagnostics. Rendering
    // applies its own stricter connectivity threshold.
    double record_abs_mayer_floor = 0.02;

    // Multicentre candidates are participation diagnostics only, not a 3c2e or
    // 3c4e assignment. Atomic weights are absolute Mulliken-like MO gross
    // populations normalized over atoms.
    double multicentre_atom_participation = 0.12;
    double multicentre_min_covered_weight = 0.72;
};

// Requires wavefunction.ao_overlap and total density. Populates pairwise Mayer
// bond orders and occupied-MO multicentre participation candidates. No
// molecule-specific chemistry rules are used and no electron-counting class is
// asserted here.
void derive_bond_and_multicentre_analysis(
    Wavefunction& wavefunction,
    const BondAnalysisOptions& options = {});

} // namespace cov
