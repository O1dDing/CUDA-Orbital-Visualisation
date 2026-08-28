#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <filesystem>

namespace cov {

struct MoldenParseOptions {
    std::size_t max_atoms = 100;
    bool require_orbitals = true;
};

// Strict Molden-format parser. Most callers should use parse_wavefunction();
// this symbol remains public for format-specific regression tests.
Wavefunction parse_molden_format(const std::filesystem::path& path,
                                 const MoldenParseOptions& options = {});

// Compatibility entry point retained for existing code. It now dispatches
// through the unified input layer, so legacy viewer code can open FCHK without
// losing the existing Molden API surface.
Wavefunction parse_molden(const std::filesystem::path& path,
                          const MoldenParseOptions& options = {});

} // namespace cov
