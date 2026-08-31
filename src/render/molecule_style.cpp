#include "cov/molecule_style.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <queue>
#include <vector>

namespace cov {
namespace {

constexpr double kMayerRenderFloor = 0.05;
// Mayer indices also contain non-local through-bond and pi-system coupling.
// Keep the electronic value as bond evidence, but do not turn those remote
// couplings into structural chords. The base 1.50 separates the next-neighbour
// shell in ordinary carbon rings. Strong direct electronic evidence can relax
// the envelope smoothly (needed by weakly bound one-electron bonds such as
// equilibrium H2+) without granting the same range to small remote pi Mayer
// terms.
constexpr double kElectronicAdjacencyFactor = 1.50;

double electronic_adjacency_factor(const double mayer_order) noexcept {
    const double strong_fraction=std::clamp(
        (mayer_order-0.20)/0.30,0.0,1.0);
    return kElectronicAdjacencyFactor+0.35*strong_fraction;
}

// Covalent radii are structural data, not element-colour cosmetics: an
// incomplete table can silently remove genuine metal--ligand bonds at the
// electronic-adjacency gate below. Values through Cm follow Cordero et al.,
// Dalton Trans. 2008, 2832 (DOI:10.1039/B801115J), matching the values this
// renderer already used. The remaining super-heavy elements use the
// single-bond radii of Pyykko and Atsumi, Chem. Eur. J. 2009, 15, 186
// (DOI:10.1002/chem.200800987), so every real atomic number accepted by the
// parser has a defined radius. Index zero is a ghost/unknown sentinel.
constexpr std::array<double, 119> kCovalentRadiiAngstrom = {
    0.85,
    // H--Ne
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
    // Na--Ar
    1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
    // K--Kr
    2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24,
    1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
    // Rb--Xe
    2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39,
    1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.40,
    // Cs--Rn
    2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96,
    1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87, 1.75, 1.70, 1.62,
    1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40,
    1.50, 1.50,
    // Fr--Cm
    2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69,
    // Bk--Og (Pyykko--Atsumi single-bond radii)
    1.68, 1.68, 1.65, 1.67, 1.73, 1.76, 1.61, 1.57, 1.49, 1.43,
    1.41, 1.34, 1.29, 1.28, 1.21, 1.22, 1.36, 1.43, 1.62, 1.75,
    1.65, 1.57,
};

bool likely_pi_element(const int z) noexcept {
    switch (z) {
        case 5:  // B
        case 6:  // C
        case 7:  // N
        case 8:  // O
        case 14: // Si
        case 15: // P
        case 16: // S
            return true;
        default:
            return false;
    }
}

double atom_distance_bohr(const Atom& a, const Atom& b) noexcept {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<std::size_t> alternate_path(const std::vector<std::vector<std::size_t>>& adjacency,
                                        const std::size_t source,
                                        const std::size_t target,
                                        const std::size_t blocked_a,
                                        const std::size_t blocked_b,
                                        const std::size_t max_edges) {
    const std::size_t n = adjacency.size();
    std::vector<int> parent(n, -1);
    std::vector<std::size_t> depth(n, max_edges + 1);
    std::queue<std::size_t> queue;
    queue.push(source);
    depth[source] = 0;

    while (!queue.empty()) {
        const std::size_t current = queue.front();
        queue.pop();
        if (depth[current] >= max_edges) continue;

        for (const std::size_t next : adjacency[current]) {
            if ((current == blocked_a && next == blocked_b) ||
                (current == blocked_b && next == blocked_a)) {
                continue;
            }
            if (depth[next] <= depth[current] + 1) continue;
            parent[next] = static_cast<int>(current);
            depth[next] = depth[current] + 1;
            if (next == target) {
                std::vector<std::size_t> path;
                for (std::size_t p = target;; p = static_cast<std::size_t>(parent[p])) {
                    path.push_back(p);
                    if (p == source) break;
                }
                std::reverse(path.begin(), path.end());
                return path;
            }
            queue.push(next);
        }
    }
    return {};
}

std::size_t bond_index_between(const std::vector<BondVisual>& bonds,
                               const std::size_t a,
                               const std::size_t b) {
    for (std::size_t i = 0; i < bonds.size(); ++i) {
        if ((bonds[i].atom_a == a && bonds[i].atom_b == b) ||
            (bonds[i].atom_a == b && bonds[i].atom_b == a)) {
            return i;
        }
    }
    return bonds.size();
}

void mark_conservative_ring_delocalisation(const Wavefunction& wavefunction,
                                           std::vector<BondVisual>& bonds) {
    if (bonds.empty()) return;
    const std::size_t atom_count = wavefunction.atoms.size();

    std::vector<std::vector<std::size_t>> adjacency(atom_count);
    for (const auto& bond : bonds) {
        if (bond.atom_a >= atom_count || bond.atom_b >= atom_count) continue;
        adjacency[bond.atom_a].push_back(bond.atom_b);
        adjacency[bond.atom_b].push_back(bond.atom_a);
    }

    // Delocalisation remains a conservative rendering annotation. Electronic
    // pairwise bond order can improve connectivity, but it is not itself enough
    // to claim aromaticity. Compact 5-7 member pi-capable rings must also have
    // unusually short and near-equal edge lengths.
    for (std::size_t bi = 0; bi < bonds.size(); ++bi) {
        const BondVisual seed = bonds[bi];
        const auto path = alternate_path(adjacency, seed.atom_a, seed.atom_b,
                                         seed.atom_a, seed.atom_b, 6);
        if (path.size() < 5 || path.size() > 7) continue;

        std::vector<std::size_t> cycle_bonds;
        cycle_bonds.push_back(bi);
        bool valid = true;
        double min_ratio = 1.0e9;
        double max_ratio = 0.0;
        double ratio_sum = 0.0;
        std::size_t ratio_count = 0;

        for (const std::size_t atom_index : path) {
            if (atom_index >= atom_count ||
                !likely_pi_element(wavefunction.atoms[atom_index].atomic_number)) {
                valid = false;
                break;
            }
        }
        if (!valid) continue;

        for (std::size_t pi = 0; pi + 1 < path.size(); ++pi) {
            const std::size_t bidx = bond_index_between(bonds, path[pi], path[pi + 1]);
            if (bidx >= bonds.size()) {
                valid = false;
                break;
            }
            cycle_bonds.push_back(bidx);
        }
        if (!valid) continue;

        for (const std::size_t bidx : cycle_bonds) {
            const auto& bond = bonds[bidx];
            const Atom& a = wavefunction.atoms[bond.atom_a];
            const Atom& b = wavefunction.atoms[bond.atom_b];
            const double single_length_bohr =
                (covalent_radius_angstrom(a.atomic_number) +
                 covalent_radius_angstrom(b.atomic_number)) * kAngstromToBohr;
            const double ratio = bond.distance_bohr / std::max(1.0e-12, single_length_bohr);
            min_ratio = std::min(min_ratio, ratio);
            max_ratio = std::max(max_ratio, ratio);
            ratio_sum += ratio;
            ++ratio_count;
        }

        if (!valid || ratio_count == 0) continue;
        const double mean_ratio = ratio_sum / static_cast<double>(ratio_count);
        const bool short_enough = mean_ratio >= 0.82 && mean_ratio <= 0.97;
        const bool near_equal = (max_ratio - min_ratio) <= 0.10;
        if (!short_enough || !near_equal) continue;

        for (const std::size_t bidx : cycle_bonds) bonds[bidx].delocalised = true;
    }
}

} // namespace

double covalent_radius_angstrom(const int z) noexcept {
    if (z > 0 && static_cast<std::size_t>(z) < kCovalentRadiiAngstrom.size()) {
        return kCovalentRadiiAngstrom[static_cast<std::size_t>(z)];
    }
    return kCovalentRadiiAngstrom[0];
}

std::vector<BondVisual> analyse_bonds(const Wavefunction& wavefunction) {
    std::vector<BondVisual> bonds;
    const std::size_t atom_count = wavefunction.atoms.size();

    if (wavefunction.bond_order_provenance != DataProvenance::Unavailable) {
        // Electronic evidence is authoritative when available, while the
        // distance envelope determines whether that pair is a drawable
        // structural neighbour. This prevents delocalised Mayer coupling from
        // becoming chords across aromatic rings.
        for (const auto& record : wavefunction.bond_orders) {
            const std::size_t i = record.atom_a;
            const std::size_t j = record.atom_b;
            if (i >= atom_count || j >= atom_count || i == j ||
                record.mayer_order < kMayerRenderFloor) {
                continue;
            }
            const Atom& a = wavefunction.atoms[i];
            const Atom& b = wavefunction.atoms[j];
            const double distance_bohr = atom_distance_bohr(a, b);
            const double sanity_bohr =
                electronic_adjacency_factor(record.mayer_order) *
                (covalent_radius_angstrom(a.atomic_number) +
                 covalent_radius_angstrom(b.atomic_number)) * kAngstromToBohr;
            if (distance_bohr > sanity_bohr) continue;

            BondVisual bond;
            bond.atom_a = i;
            bond.atom_b = j;
            bond.distance_bohr = distance_bohr;
            bond.bond_order = record.mayer_order;
            bond.provenance = record.provenance;
            bonds.push_back(bond);
        }
    } else {
        // Compatibility fallback for files where a complete overlap/density
        // analysis cannot be obtained. This remains intentionally conservative.
        for (std::size_t i = 0; i < atom_count; ++i) {
            for (std::size_t j = i + 1; j < atom_count; ++j) {
                const Atom& a = wavefunction.atoms[i];
                const Atom& b = wavefunction.atoms[j];
                const double distance_bohr = atom_distance_bohr(a, b);
                const double cutoff_bohr =
                    1.25 * (covalent_radius_angstrom(a.atomic_number) +
                            covalent_radius_angstrom(b.atomic_number)) * kAngstromToBohr;
                if (distance_bohr <= cutoff_bohr) {
                    BondVisual bond;
                    bond.atom_a = i;
                    bond.atom_b = j;
                    bond.distance_bohr = distance_bohr;
                    bond.provenance = DataProvenance::Unavailable;
                    bonds.push_back(bond);
                }
            }
        }
    }

    mark_conservative_ring_delocalisation(wavefunction, bonds);
    return bonds;
}

} // namespace cov
