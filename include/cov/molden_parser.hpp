#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <filesystem>

namespace cov {

struct MoldenParseOptions {
    std::size_t max_atoms = 100;
    bool require_orbitals = true;
};

Wavefunction parse_molden(const std::filesystem::path& path,
                          const MoldenParseOptions& options = {});

} // namespace cov
