#pragma once

#include "cov/nbo_salc.hpp"
#include <string>
#include <vector>

namespace cov::ui {
struct NboAomoName {
    std::string label, irrep;
    // Zero means that a complete, unambiguous ordinal could not be established.
    std::size_t ordinal=0;
    // Irrep evidence only: true with ordinal==0 means the symmetry is known
    // while occurrence order / repeated-copy membership remains unresolved.
    bool verified=false;
    std::string detail;
};
struct NboAomoNames {
    std::vector<NboAomoName> canonical, salc;
};
// Presentation-only evidence cache. Call once per immutable attachment, never
// per frame. No coefficients, source labels, orbital ordering or fields change.
NboAomoNames build_nbo_aomo_names(const Wavefunction&,const NboIntegration&,
                                 const NboSalcModel*);
} // namespace cov::ui
