#pragma once

#include "cov/nbo_salc.hpp"
#include <string>
#include <vector>
#include <memory>

namespace cov::ui {
struct NboAomoName {
    std::string label, irrep;
    // Zero means that a complete, unambiguous ordinal could not be established.
    std::size_t ordinal=0;
    // Irrep evidence only: true with ordinal==0 means the symmetry is known
    // while occurrence order / repeated-copy membership remains unresolved.
    bool verified=false;
    std::string detail;
    // A verified single irrep occurrence, independent of whether its position
    // among other occurrences can be numbered. Empty means unproved partners.
    std::string partner_block_id;
    std::size_t partner_block_size=0;
};
struct NboAomoNames {
    std::vector<NboAomoName> canonical, salc;
};
// Presentation-only evidence cache. Call once per immutable attachment, never
// per frame. No coefficients, source labels, orbital ordering or fields change.
NboAomoNames build_nbo_aomo_names(const Wavefunction&,const NboIntegration&,
                                 const NboSalcModel*);
// Canonical naming is available without an NBO attachment. The cache is for
// immutable loaded wavefunctions; it never changes producer data or coefficients.
std::shared_ptr<const NboAomoNames> canonical_mo_names(const Wavefunction&);
void invalidate_canonical_mo_names_cache();
// Source identity is spin-block based; unavailable source indices are explicitly
// identified as list positions. Internal indices remain the selection addresses.
std::string canonical_mo_source_label(const Wavefunction&,std::size_t);
std::string canonical_mo_display_label(const Wavefunction&,std::size_t,
                                     const NboAomoName* = nullptr);
std::string orbital_irrep_display_label(const NboAomoName&);
} // namespace cov::ui
