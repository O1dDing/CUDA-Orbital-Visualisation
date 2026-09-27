#pragma once

#include "cov/ui.hpp"
#include <string>
#include <vector>

namespace cov {
struct Wavefunction;
struct NboIntegration;
struct NboSalcModel;
namespace ui {
struct NboAomoNode;
// Presentation of existing evidence only. No chemistry is derived, and no
// orbital/selection is modified. At most four concise lines, without energies.
std::vector<std::string> nbo_aomo_hover_lines(const NboAomoNode&,
    const Wavefunction&, const NboIntegration&, const NboSalcModel*, Language);
std::string nbo_aomo_hover_glyph_seed(Language);
} // namespace ui
} // namespace cov
