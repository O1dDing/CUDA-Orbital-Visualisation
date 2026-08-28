#include "cov/wavefunction_io.hpp"

#include "cov/fchk_parser.hpp"
#include "cov/molden_parser.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>
#include <string>

namespace cov {
namespace {

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool looks_like_molden(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    std::string line;
    for (int i = 0; i < 8 && std::getline(in, line); ++i) {
        const auto folded = lower(line);
        if (folded.find("[molden format]") != std::string::npos) return true;
    }
    return false;
}

} // namespace

Wavefunction parse_wavefunction(const std::filesystem::path& path,
                                const WavefunctionParseOptions& options) {
    const std::string filename = lower(path.filename().string());
    const std::string extension = lower(path.extension().string());

    if (extension == ".fchk" || extension == ".fch") {
        FchkParseOptions fchk;
        fchk.max_atoms = options.max_atoms;
        fchk.require_orbitals = options.require_orbitals;
        fchk.keep_density = options.keep_density;
        fchk.reconstruct_density_if_missing = options.reconstruct_density_if_missing;
        return parse_fchk(path, fchk);
    }

    if (extension == ".molden" ||
        filename.ends_with(".molden.input") || filename.ends_with(".molden.inp") ||
        looks_like_molden(path)) {
        MoldenParseOptions molden;
        molden.max_atoms = options.max_atoms;
        molden.require_orbitals = options.require_orbitals;
        Wavefunction wf = parse_molden(path, molden);
        wf.source = WavefunctionSource::Molden;
        return wf;
    }

    throw std::runtime_error(
        "Unsupported wavefunction input. Use Gaussian .fchk/.fch or Molden .molden files.");
}

} // namespace cov
