#pragma once

#include "cov/nbo_salc.hpp"
#include <string>
#include <vector>
#include <memory>

namespace cov::ui {
struct NboIrrepContent {
    std::string irrep;
    std::size_t dimension=0,multiplicity=0;
};
struct NboAomoName {
    std::string label, irrep;
    // Group actually used for the verified character comparison; never a
    // dimension/energy guess or a group inferred from the label spelling.
    std::string point_group;
    // Zero means that a complete, unambiguous ordinal could not be established.
    std::size_t ordinal=0;
    // Irrep evidence only: true with ordinal==0 means the symmetry is known
    // while occurrence order / repeated-copy membership remains unresolved.
    bool verified=false;
    // Number of irreducible copies in the measured containing span.
    // This does not assign an individual copy or occurrence ordinal.
    std::size_t representation_multiplicity=1;
    std::string detail;
    // A verified single irrep occurrence, independent of whether its position
    // among other occurrences can be numbered. Empty means unproved partners.
    std::string partner_block_id;
    std::size_t partner_block_size=0;
    // Scientific containing span is independent of the current display filter.
    std::vector<std::size_t> containing_members;
    std::vector<NboIrrepContent> containing_irreps;
    std::string status,ordinal_scope;
    std::size_t complete_set_ordinal=0;
    std::optional<double> projection_residual;
};
struct NboAomoNames {
    std::vector<NboAomoName> canonical, salc;
};
std::string serialize_orbital_name_json(const NboAomoName&);
std::string serialize_orbital_names_json(const NboAomoNames&);
// Presentation-only evidence cache. Call once per immutable attachment, never
// per frame. No coefficients, source labels, orbital ordering or fields change.
NboAomoNames build_nbo_aomo_names(const Wavefunction&,const NboIntegration&,
                                 const NboSalcModel*);
// Canonical naming is available without an NBO attachment. The cache is for
// immutable loaded wavefunctions; it never changes producer data or coefficients.
std::shared_ptr<const NboAomoNames> canonical_mo_names(const Wavefunction&);
void invalidate_canonical_mo_names_cache();
// Display ordinals count only complete, verified occurrences in this filter.
// Source IDs, complete-set ordinals, evidence and unresolved copies are retained.
NboAomoNames nbo_aomo_names_for_view(const Wavefunction&,const NboAomoNames&,
    const std::vector<std::size_t>& canonical_indices,
    const std::vector<std::size_t>& salc_indices,const NboSalcModel*,
    const std::string& scope);
// Source identity is spin-block based; unavailable source indices are explicitly
// identified as list positions. Internal indices remain the selection addresses.
std::string canonical_mo_source_label(const Wavefunction&,std::size_t);
std::string canonical_mo_display_label(const Wavefunction&,std::size_t,
                                     const NboAomoName* = nullptr);
std::string orbital_irrep_display_label(const NboAomoName&);
std::string canonical_mo_current_irrep(const Wavefunction&,std::size_t,
                                     const NboAomoName* = nullptr);
} // namespace cov::ui
