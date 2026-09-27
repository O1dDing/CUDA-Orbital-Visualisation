#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <filesystem>

namespace cov {

struct WavefunctionParseOptions {
    std::size_t max_atoms = 100;
    bool require_orbitals = true;
    bool keep_density = true;
    bool reconstruct_density_if_missing = true;

    // FCHK remains the wavefunction authority. A Gaussian output sidecar is
    // optional enrichment for producer-reported point-group/orbital symmetry.
    bool auto_enrich_gaussian_log = true;
    std::filesystem::path gaussian_log_path;
};

Wavefunction parse_wavefunction(const std::filesystem::path& path,
                                const WavefunctionParseOptions& options = {});

} // namespace cov
