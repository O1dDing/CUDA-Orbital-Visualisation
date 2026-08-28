#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <filesystem>
#include <string>

namespace cov {

struct GaussianLogEnrichmentResult {
    bool opened = false;
    std::size_t symmetry_labels_applied = 0;
    bool point_group_applied = false;
    std::string warning;

    [[nodiscard]] bool applied() const noexcept {
        return symmetry_labels_applied != 0 || point_group_applied;
    }
};

// Apply producer-reported metadata from a Gaussian .log/.out sidecar. At this
// stage the enrichment intentionally covers final orbital symmetry labels and
// point-group metadata only. It never substitutes the log for the wavefunction
// coefficients carried by FCHK/Molden.
[[nodiscard]] GaussianLogEnrichmentResult enrich_from_gaussian_log(
    Wavefunction& wavefunction,
    const std::filesystem::path& path);

// Look for an adjacent producer output with the same stem. Preference is .log,
// then .out. Returns an empty path when neither exists.
[[nodiscard]] std::filesystem::path find_sibling_gaussian_log(
    const std::filesystem::path& wavefunction_path);

} // namespace cov
