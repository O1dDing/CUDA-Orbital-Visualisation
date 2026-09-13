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
    bool spin_squared_applied = false;
    bool scf_diagnostic_applied = false;
    bool stability_diagnostic_applied = false;
    std::string warning;

    [[nodiscard]] bool applied() const noexcept {
        return symmetry_labels_applied != 0 || point_group_applied ||
               spin_squared_applied || scf_diagnostic_applied ||
               stability_diagnostic_applied;
    }
};

// Apply producer-reported metadata from a Gaussian .log/.out sidecar. At this
// stage the enrichment covers final orbital symmetry/point-group metadata and
// explicit producer diagnostics (SCF convergence, stability, and S**2 before
// and after annihilation). It never substitutes the log for the wavefunction
// coefficients carried by FCHK/Molden and never infers a diagnostic merely
// from route keywords.
[[nodiscard]] GaussianLogEnrichmentResult enrich_from_gaussian_log(
    Wavefunction& wavefunction,
    const std::filesystem::path& path);

// Look for a producer output with the same stem. Preference is an adjacent
// .log/.out, followed by an equivalent relative path below a parallel .log
// tree when the wavefunction resides below a .chk directory. Returns an empty
// path when no companion output exists.
[[nodiscard]] std::filesystem::path find_sibling_gaussian_log(
    const std::filesystem::path& wavefunction_path);

} // namespace cov
