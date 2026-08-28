#include "cov/wavefunction_io.hpp"

#include "cov/bond_analysis.hpp"
#include "cov/density.hpp"
#include "cov/fchk_parser.hpp"
#include "cov/molden_parser.hpp"
#include "cov/overlap.hpp"

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

void postprocess_wavefunction(Wavefunction& wf,
                              const WavefunctionParseOptions& options) {
    // Normalize provenance at the unified boundary. Molden values are explicit
    // producer fields; FCHK occupations are reconstructed from electron counts.
    for (auto& mo : wf.orbitals) {
        if (wf.source == WavefunctionSource::Molden) {
            if (mo.occupation_provenance == DataProvenance::Unavailable) {
                mo.occupation_provenance = DataProvenance::Producer;
            }
            if (!mo.symmetry.empty() &&
                mo.symmetry_provenance == DataProvenance::Unavailable) {
                mo.symmetry_provenance = DataProvenance::Producer;
            }
        } else if (wf.source == WavefunctionSource::Fchk &&
                   mo.occupation_provenance == DataProvenance::Unavailable) {
            mo.occupation_provenance = DataProvenance::Derived;
        }
    }

    if (options.keep_density && !wf.orbitals.empty()) {
        if (wf.total_density_packed.empty() &&
            options.reconstruct_density_if_missing) {
            wf.total_density_packed = reconstruct_total_density_packed(wf);
            wf.total_density_provenance = DataProvenance::Derived;
        }
        if (wf.spin_density_packed.empty() &&
            options.reconstruct_density_if_missing) {
            wf.spin_density_packed = reconstruct_spin_density_packed(wf);
            wf.spin_density_provenance = DataProvenance::Derived;
        }
    }

    // FCHK and Molden both normally carry complete canonical MO blocks. Recover
    // S from MO orthonormality before introducing a second, convention-sensitive
    // integral engine. If the block is incomplete/ill-conditioned, leave the
    // analysis explicitly unavailable rather than guessing.
    const auto overlap = derive_ao_overlap_from_mos(wf);
    if (overlap.available()) {
        wf.ao_overlap = overlap.matrix;
        wf.ao_overlap_provenance = DataProvenance::Derived;
        wf.ao_overlap_orthonormality_error = overlap.max_orthonormality_error;
    }

    if (!wf.ao_overlap.empty() && !wf.total_density_packed.empty()) {
        derive_bond_and_multicentre_analysis(wf);
    }
}

} // namespace

Wavefunction parse_wavefunction(const std::filesystem::path& path,
                                const WavefunctionParseOptions& options) {
    const std::string filename = lower(path.filename().string());
    const std::string extension = lower(path.extension().string());

    Wavefunction wf;
    if (extension == ".fchk" || extension == ".fch") {
        FchkParseOptions fchk;
        fchk.max_atoms = options.max_atoms;
        fchk.require_orbitals = options.require_orbitals;
        fchk.keep_density = options.keep_density;
        fchk.reconstruct_density_if_missing = options.reconstruct_density_if_missing;
        wf = parse_fchk(path, fchk);
    } else if (extension == ".molden" ||
               filename.ends_with(".molden.input") || filename.ends_with(".molden.inp") ||
               looks_like_molden(path)) {
        MoldenParseOptions molden;
        molden.max_atoms = options.max_atoms;
        molden.require_orbitals = options.require_orbitals;
        wf = parse_molden_format(path, molden);
        wf.source = WavefunctionSource::Molden;
    } else {
        throw std::runtime_error(
            "Unsupported wavefunction input. Use Gaussian .fchk/.fch or Molden .molden files.");
    }

    postprocess_wavefunction(wf, options);
    return wf;
}

} // namespace cov
