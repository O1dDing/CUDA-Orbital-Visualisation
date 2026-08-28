#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace cov {

constexpr double kAngstromToBohr = 1.8897261254578281;

enum class Spin : std::uint8_t {
    Alpha = 0,
    Beta = 1,
};

enum class WavefunctionSource : std::uint8_t {
    Unknown = 0,
    Molden = 1,
    Fchk = 2,
};

enum class DataProvenance : std::uint8_t {
    Unavailable = 0,
    Producer = 1,
    Derived = 2,
};

struct Atom {
    std::string symbol;
    int atomic_number = 0;
    double x = 0.0; // bohr
    double y = 0.0; // bohr
    double z = 0.0; // bohr

    // Producer-reported effective nuclear charge. This equals atomic_number
    // for all-electron calculations and is smaller when an ECP replaces core
    // electrons. Parsers set it explicitly; zero means "not reported".
    double nuclear_charge = 0.0;
};

struct Primitive {
    float exponent = 0.0f;
    float coefficient = 0.0f;
};

struct Shell {
    std::uint32_t atom_index = 0;
    std::uint32_t primitive_offset = 0;
    std::uint32_t primitive_count = 0;
    std::uint32_t basis_offset = 0;
    std::uint8_t angular_momentum = 0; // s=0, p=1, d=2, ...
    std::uint8_t pure = 0;             // 0 Cartesian, 1 real spherical
};

enum class ChemistryStatus : std::uint8_t {
    Determined = 0,
    Percentages,
    Undetermined,
    NotApplicable,
    Unavailable,
};

enum class OrbitalAngularFamily : std::uint8_t {
    Sigma = 0,
    Pi,
    Delta,
    Phi,
    Mixed,
    NotApplicable,
};

enum class OrbitalBondingRole : std::uint8_t {
    Bonding = 0,
    Antibonding,
    Nonbonding,
    Mixed,
    NotApplicable,
};

struct OrbitalChannelDistribution {
    double sigma = 0.0;
    double pi = 0.0;
    double delta = 0.0;
    double phi = 0.0;
    double undetermined = 1.0;
    OrbitalAngularFamily dominant = OrbitalAngularFamily::Mixed;
    ChemistryStatus status = ChemistryStatus::Unavailable;
};

struct OrbitalBondingDistribution {
    double bonding = 0.0;
    double antibonding = 0.0;
    double nonbonding = 0.0;
    double undetermined = 1.0;
    OrbitalBondingRole dominant = OrbitalBondingRole::Mixed;
    ChemistryStatus status = ChemistryStatus::Unavailable;
};

struct OrbitalAOContribution {
    std::uint32_t atom_index = 0;
    int principal_n = 0;
    int angular_momentum = 0;
    std::string label;
    double weight = 0.0;
};

struct OrbitalPairInteraction {
    std::uint32_t atom_a = 0;
    std::uint32_t atom_b = 0;
    std::string atom_a_label;
    std::string atom_b_label;

    // Total density-level Mayer index for the pair.
    double total_mayer_index = 0.0;

    // Per-unit-occupation Mulliken overlap character for this canonical MO.
    // This is representation/basis dependent and is not a universal bond order.
    double overlap_character = 0.0;
    double occupied_overlap_contribution = 0.0;

    OrbitalChannelDistribution channel;
    OrbitalBondingDistribution bonding;
};

struct OrbitalChemistry {
    bool available = false;
    bool valence_manifold = false;

    double deep_core_weight = 0.0;
    double semicore_weight = 0.0;
    double valence_weight = 0.0;
    double unresolved_weight = 1.0;

    OrbitalChannelDistribution channel;
    OrbitalBondingDistribution bonding;
    std::vector<OrbitalAOContribution> ao_contributions;
    std::vector<OrbitalPairInteraction> interactions;

    std::string multicentre_label;
    std::string family_symbol;
    std::size_t participating_atoms = 0;
    double participating_electrons = 0.0;

    // FCHK-only canonical-MO evidence cannot always define a unique
    // donor/acceptor direction. In that case this remains "UND".
    std::string donor_acceptor = "UND";

    std::string method;
    double confidence = 0.0;
    std::string note;
};

struct MolecularOrbital {
    double energy_hartree = 0.0;
    float occupation = 0.0f;
    Spin spin = Spin::Alpha;
    std::string symmetry;
    std::vector<float> coefficients;

    // Molden can carry occupation/symmetry explicitly; FCHK normally carries
    // electron counts rather than per-orbital occupations and usually does not
    // carry per-MO irreps. Keeping provenance prevents derived values from being
    // presented as producer data later in the UI/export layer.
    DataProvenance occupation_provenance = DataProvenance::Unavailable;
    DataProvenance symmetry_provenance = DataProvenance::Unavailable;

    // Derived after density/overlap/symmetry/bond analysis. Producer data is
    // never fabricated; unsupported conclusions remain percentage-valued or UND.
    OrbitalChemistry chemistry;
};

struct BondOrderRecord {
    std::uint32_t atom_a = 0;
    std::uint32_t atom_b = 0;
    double mayer_order = 0.0;
    DataProvenance provenance = DataProvenance::Unavailable;
};

struct MulticentreCandidate {
    std::uint32_t orbital_index = 0;
    std::vector<std::uint32_t> atoms;
    std::vector<double> participation;
    double occupation = 0.0;
    DataProvenance provenance = DataProvenance::Unavailable;
};

enum class MulticentreKind : std::uint8_t {
    Unclassified = 0,
    ThreeCentreTwoElectron,
    ThreeCentreFourElectron,
};

struct MulticentreAssignment {
    MulticentreKind kind = MulticentreKind::Unclassified;
    std::vector<std::uint32_t> atoms;
    std::vector<std::uint32_t> orbitals;
    double electron_count = 0.0;
    double confidence = 0.0;
    std::string rationale;
    DataProvenance provenance = DataProvenance::Unavailable;
};

struct Wavefunction {
    std::vector<Atom> atoms;
    std::vector<Primitive> primitives;
    std::vector<Shell> shells;
    std::vector<MolecularOrbital> orbitals;
    std::uint32_t basis_count = 0;
    bool pure_d = false;
    bool pure_f = false;
    bool pure_g = false;

    // Input/provenance metadata. These fields are intentionally orthogonal to
    // rendering so Molden and FCHK converge on the same Wavefunction model.
    WavefunctionSource source = WavefunctionSource::Unknown;
    std::uint32_t alpha_electrons = 0;
    std::uint32_t beta_electrons = 0;
    std::string source_title;
    std::string source_route;
    std::string enrichment_source;
    std::string point_group_detected;
    std::string point_group_used;
    DataProvenance point_group_provenance = DataProvenance::Unavailable;

    // Packed lower-triangular AO density matrices when present or safely
    // reconstructable. Size is basis_count*(basis_count+1)/2.
    std::vector<double> total_density_packed;
    std::vector<double> spin_density_packed;
    DataProvenance total_density_provenance = DataProvenance::Unavailable;
    DataProvenance spin_density_provenance = DataProvenance::Unavailable;

    // AO overlap recovered from a complete orthonormal MO block when the input
    // format does not provide S explicitly. Row-major basis_count*basis_count.
    std::vector<double> ao_overlap;
    DataProvenance ao_overlap_provenance = DataProvenance::Unavailable;
    double ao_overlap_orthonormality_error = 0.0;

    // Provenance-aware analyses derived from density + overlap. Mayer bond
    // orders are pairwise evidence only. Multicentre assignments are made only
    // by the separate active-subspace classifier when its stricter evidence
    // requirements are met.
    std::vector<BondOrderRecord> bond_orders;
    DataProvenance bond_order_provenance = DataProvenance::Unavailable;
    std::vector<MulticentreCandidate> multicentre_candidates;
    std::vector<MulticentreAssignment> multicentre_assignments;
};

struct GridBox {
    float min_x = -5.0f;
    float min_y = -5.0f;
    float min_z = -5.0f;
    float max_x = 5.0f;
    float max_y = 5.0f;
    float max_z = 5.0f;
};

inline std::uint32_t cartesian_basis_count(const std::uint8_t l) {
    return static_cast<std::uint32_t>((l + 1u) * (l + 2u) / 2u);
}

inline std::uint32_t spherical_basis_count(const std::uint8_t l) {
    return static_cast<std::uint32_t>(2u * l + 1u);
}

inline std::uint32_t shell_basis_count(const Shell& shell) {
    return shell.pure ? spherical_basis_count(shell.angular_momentum)
                      : cartesian_basis_count(shell.angular_momentum);
}

inline const char* wavefunction_source_name(const WavefunctionSource source) noexcept {
    switch (source) {
        case WavefunctionSource::Molden: return "Molden";
        case WavefunctionSource::Fchk: return "FCHK";
        default: return "Unknown";
    }
}

inline const char* data_provenance_name(const DataProvenance provenance) noexcept {
    switch (provenance) {
        case DataProvenance::Producer: return "producer";
        case DataProvenance::Derived: return "derived";
        default: return "unavailable";
    }
}

inline const char* multicentre_kind_name(const MulticentreKind kind) noexcept {
    switch (kind) {
        case MulticentreKind::ThreeCentreTwoElectron: return "3c2e";
        case MulticentreKind::ThreeCentreFourElectron: return "3c4e";
        default: return "unclassified";
    }
}

} // namespace cov
