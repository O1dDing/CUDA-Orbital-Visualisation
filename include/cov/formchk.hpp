#pragma once

#include "cov/fchk_parser.hpp"

#include <filesystem>

namespace cov {

// Convert a binary Gaussian checkpoint with the locally installed formchk tool,
// parse the temporary FCHK through the same strict path, and remove the
// temporary file. COV_FORMCHK may point to a specific formchk executable;
// otherwise "formchk"/"formchk.exe" is resolved through PATH.
Wavefunction parse_gaussian_chk_via_formchk(
    const std::filesystem::path& chk_path,
    const FchkParseOptions& options = {});

} // namespace cov
