#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <filesystem>

namespace cov {

struct FchkParseOptions {
    std::size_t max_atoms = 100;
    bool require_orbitals = true;
    bool keep_density = true;
    bool reconstruct_density_if_missing = true;
};

Wavefunction parse_fchk(const std::filesystem::path& path,
                        const FchkParseOptions& options = {});

} // namespace cov
