#include "cov/molden_parser.hpp"

#include "cov/wavefunction_io.hpp"

namespace cov {

Wavefunction parse_molden(const std::filesystem::path& path,
                          const MoldenParseOptions& options) {
    WavefunctionParseOptions generic;
    generic.max_atoms = options.max_atoms;
    generic.require_orbitals = options.require_orbitals;
    return parse_wavefunction(path, generic);
}

} // namespace cov
