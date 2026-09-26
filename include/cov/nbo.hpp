#pragma once
#include "cov/model.hpp"
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cov {
enum class NboSpin { Total, Alpha, Beta };
struct NboSource {
    std::string path, block, producer_version, raw;
    std::size_t line_begin=0, line_end=0, analysis_segment=0;
};
struct NboPopulation {
    std::size_t atom=0; // producer one-based atom identity
    std::string symbol;
    NboSpin spin=NboSpin::Total;
    double charge=0, core=0, valence=0, rydberg=0, total=0;
    // NBO NPA includes the ECP replacement core in its printed population.
    // These are derived only when the archive reports effective charges.
    std::optional<double> effective_core_electrons, explicit_population;
    std::optional<double> spin_density;
    NboSource source;
};
struct NboNao {
    std::size_t id=0, atom=0;
    NboSpin spin=NboSpin::Total;
    std::string symbol, angular, type;
    double occupation=0;
    std::optional<double> energy_hartree; // NAO diagonal Fock energy, not canonical eigenvalue
    std::optional<double> spin_density; // total-spin NAO table prints Spin instead of Energy
    NboSource source;
};
struct NboLocalComponent {
    std::size_t atom=0;
    double percent=0, coefficient=0;
    std::string hybrid; // literal producer hybrid description, no renormalization
    NboSource source;
};
struct NboOrbital {
    std::size_t id=0, ordinal=0;
    NboSpin spin=NboSpin::Total;
    std::string kind, label;
    std::vector<std::size_t> atoms;
    double occupation=0;
    std::optional<double> diagonal_fock_hartree;
    std::optional<NboSource> energy_source;
    std::vector<NboLocalComponent> components;
    NboSource source;
};
struct NboE2 {
    std::size_t donor=0, acceptor=0;
    NboSpin spin=NboSpin::Total;
    double value=0, energy_gap_hartree=0, fock_hartree=0;
    std::string units="kcal/mol";
    std::optional<double> printing_threshold;
    NboSource source;
};
struct NboE2Section {
    NboSpin spin=NboSpin::Total;
    std::optional<double> printing_threshold;
    std::string units="kcal/mol", missing_reason;
    NboSource source;
};
struct NboWiberg {
    std::size_t atom_a=0, atom_b=0;
    NboSpin spin=NboSpin::Total;
    double value=0;
    NboSource source;
};
struct NboMatrix {
    std::string kind; // AONBO: AO rows / NBO columns; NBOMO: NBO rows / MO columns
    NboSpin spin=NboSpin::Total;
    std::size_t rows=0, columns=0;
    std::vector<double> values; // row-major, always complete or rejected
    NboSource source;
};
struct NboArchive {
    std::vector<Atom> atoms;
    std::size_t basis_count=0;
    bool open_shell=false, density_is_bond_order=false;
    std::vector<int> centers, labels, ncomp, nprim, nptr;
    std::vector<double> exponents, cs, cp, cd, cf, cg;
    std::vector<NboMatrix> matrices;
    NboSource source;
};
struct NboCanonicalEvidence {
    NboSpin spin=NboSpin::Total;
    std::string coefficient_source;
    bool direct_fchk_coefficients=false, density_verified=false;
    std::string detail;
};
struct NboAssociation {
    bool compatible=false;
    std::string status="not_checked", detail;
    double geometry_max_error_bohr=0, overlap_max_error=0, density_max_error=0, canonical_max_error=0;
    // NBO AO row -> literal Gaussian AO row, with C_G = scale * C_NBO.
    std::vector<std::size_t> gaussian_row;
    std::vector<double> coefficient_scale;
    std::vector<NboCanonicalEvidence> canonical_evidence;
};
struct NboDataset {
    std::string producer_version;
    std::vector<NboPopulation> populations;
    std::vector<NboNao> naos;
    std::vector<NboOrbital> orbitals;
    std::vector<NboE2> e2;
    std::vector<NboE2Section> e2_sections;
    std::vector<NboWiberg> wiberg;
    std::vector<NboMatrix> matrices;
    std::vector<NboSource> cmo_summaries; // thresholded text is never a complete matrix
    std::vector<std::string> warnings;
    std::optional<NboArchive> archive;
    NboAssociation association;
    NboSource source;
};
struct NboReadOptions {
    std::optional<std::size_t> analysis_segment; // zero-based; multiple analyses require explicit selection
    std::filesystem::path archive47, aonbo, nbomo;
};
NboDataset read_nbo(const std::filesystem::path& output, const NboReadOptions& options={});
NboArchive read_nbo_archive(const std::filesystem::path& path);
NboAssociation associate_nbo(NboDataset& dataset, const Wavefunction& canonical);
Wavefunction make_nbo_wavefunction(const NboDataset& dataset, const Wavefunction& canonical);
std::string serialize_nbo_json(const NboDataset& dataset);
const char* nbo_spin_name(NboSpin spin) noexcept;
} // namespace cov
