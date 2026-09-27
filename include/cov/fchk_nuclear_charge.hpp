#pragma once

#include "cov/model.hpp"

#include <filesystem>

namespace cov {

// Reads the producer Nuclear charges array from an FCHK/FCH file. Values equal
// Z for all-electron centres and are smaller when an ECP replaces core
// electrons. If the field is absent, zero/unset centres fall back to Z.
bool enrich_fchk_nuclear_charges_from_file(
    Wavefunction& wavefunction,
    const std::filesystem::path& path);

} // namespace cov
