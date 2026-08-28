#include "cov/wavefunction_io.hpp"

#include "cov/bond_analysis.hpp"
#include "cov/density.hpp"
#include "cov/fchk_parser.hpp"
#include "cov/formchk.hpp"
#include "cov/gaussian_log.hpp"
#include "cov/molden_parser.hpp"
#include "cov/orbital_symmetry.hpp"
#include "cov/overlap.hpp"
#include "cov/symmetry.hpp"

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

    // Producer-reported Gaussian point-group information wins when a sibling
    // log/out enrichment supplied it. FCHK-only and Molden-only inputs receive a
    // geometry-derived point group from COV's own operation/permutation engine.
    derive_point_group_from_geometry(wf);

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

    const auto overlap = derive_ao_overlap_from_mos(wf);
    if (overlap.available()) {
        wf.ao_overlap = overlap.matrix;
        wf.ao_overlap_provenance = DataProvenance::Derived;
        wf.ao_overlap_orthonormality_error = overlap.max_orthonormality_error;
    }

    // FCHK normally lacks per-MO irreps. Once a validated geometry operation set
    // and AO overlap are available, derive irreps from the transformed AO/MO
    // representation. Producer labels remain immutable and always take priority.
    if (!wf.ao_overlap.empty() && !wf.orbitals.empty()) {
        (void)derive_orbital_symmetry(wf);
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
    if (extension == ".fchk" || extension == ".fch" || extension == ".chk") {
        FchkParseOptions fchk;
        fchk.max_atoms = options.max_atoms;
        fchk.require_orbitals = options.require_orbitals;
        fchk.keep_density = options.keep_density;
        fchk.reconstruct_density_if_missing = options.reconstruct_density_if_missing;
        if (extension == ".chk") {
            wf = parse_gaussian_chk_via_formchk(path, fchk);
        } else {
            wf = parse_fchk(path, fchk);
        }
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
            "Unsupported wavefunction input. Use Gaussian .fchk/.fch/.chk or Molden .molden files.");
    }

    std::filesystem::path enrichment = options.gaussian_log_path;
    if (enrichment.empty() && options.auto_enrich_gaussian_log) {
        enrichment = find_sibling_gaussian_log(path);
    }
    if (!enrichment.empty()) {
        (void)enrich_from_gaussian_log(wf, enrichment);
    }

    postprocess_wavefunction(wf, options);
    return wf;
}

} // namespace cov
