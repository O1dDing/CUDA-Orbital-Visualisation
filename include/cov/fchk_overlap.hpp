#pragma once

#include "cov/model.hpp"

#include <filesystem>

namespace cov {

// Read the optional packed AO overlap matrix from a formatted checkpoint and
// transform it into COV's internal AO ordering. Returns true only when a valid
// producer matrix was applied. Missing overlap is not an error because COV can
// fall back to strict MO-derived overlap reconstruction.
bool enrich_fchk_overlap_from_file(
    Wavefunction& wavefunction,
    const std::filesystem::path& fchk_path);

} // namespace cov
