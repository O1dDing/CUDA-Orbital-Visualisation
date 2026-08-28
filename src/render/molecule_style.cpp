#include "cov/molecule_style.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <vector>

namespace cov {
namespace {

constexpr double kMayerRenderFloor = 0.05;
constexpr double kElectronicDistanceSanityFactor = 2.20;

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
    switch (z) {
        case 1: return 0.31;
        case 5: return 0.84;
        case 6: return 0.76;
        case 7: return 0.71;
        case 8: return 0.66;
        case 9: return 0.57;
        case 14: return 1.11;
        case 15: return 1.07;
        case 16: return 1.05;
        case 17: return 1.02;
        case 21: return 1.70;
        case 22: return 1.60;
        case 23: return 1.53;
        case 24: return 1.39;
        case 25: return 1.39;
        case 26: return 1.32;
        case 27: return 1.26;
        case 28: return 1.24;
        case 29: return 1.32;
        case 30: return 1.22;
        case 35: return 1.20;
        case 53: return 1.39;
        case 54: return 1.40;
        default: return 0.85;
    }
}

std::vector<BondVisual> analyse_bonds(const Wavefunction& wavefunction) {
    std::vector<BondVisual> bonds;
    const std::size_t atom_count = wavefunction.atoms.size();

    if (wavefunction.bond_order_provenance != DataProvenance::Unavailable) {
        // Electronic connectivity is authoritative when available. A generous
        // distance sanity guard rejects remote numerical couplings without
        // reverting to a tight covalent-radius cutoff that would again erase
        // coordination bonds.
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
                kElectronicDistanceSanityFactor *
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
