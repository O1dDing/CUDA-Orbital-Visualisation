#include "cov/bond_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cov {
namespace {

std::size_t packed_index(std::size_t i, std::size_t j) {
    if (j > i) std::swap(i, j);
    return i * (i + 1u) / 2u + j;
}

std::vector<double> unpack_symmetric(const std::vector<double>& packed,
                                     const std::size_t n,
                                     const char* label) {
    const std::size_t expected = n * (n + 1u) / 2u;
    if (packed.size() != expected) {
        throw std::runtime_error(std::string(label) +
                                 " packed dimension does not match AO basis");
    }
    std::vector<double> full(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const double value = packed[packed_index(i, j)];
            full[i * n + j] = value;
            full[j * n + i] = value;
        }
    }
    return full;
}

std::vector<std::uint32_t> basis_atom_map(const Wavefunction& wf) {
    const std::size_t n = wf.basis_count;
    const auto missing = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> map(n, missing);
    for (const auto& shell : wf.shells) {
        const std::size_t begin = shell.basis_offset;
        const std::size_t count = shell_basis_count(shell);
        if (begin + count > n || shell.atom_index >= wf.atoms.size()) {
            throw std::runtime_error("Shell/basis mapping is inconsistent during bond analysis");
        }
        for (std::size_t local = 0; local < count; ++local) {
            const std::size_t basis = begin + local;
            if (map[basis] != missing) {
                throw std::runtime_error("AO basis function is assigned to more than one shell");
            }
            map[basis] = shell.atom_index;
        }
    }
    for (const auto atom : map) {
        if (atom == missing) {
            throw std::runtime_error("AO basis function is not assigned to an atom");
        }
    }
    return map;
}

std::vector<double> multiply(const std::vector<double>& a,
                             const std::vector<double>& b,
                             const std::size_t n) {
    std::vector<double> out(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < n; ++k) {
            const double aik = a[i * n + k];
            if (aik == 0.0) continue;
            for (std::size_t j = 0; j < n; ++j) {
                out[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    return out;
}

} // namespace

void derive_bond_and_multicentre_analysis(Wavefunction& wavefunction,
                                          const BondAnalysisOptions& options) {
    wavefunction.bond_orders.clear();
    wavefunction.bond_order_provenance = DataProvenance::Unavailable;
    wavefunction.multicentre_candidates.clear();

    const std::size_t n = wavefunction.basis_count;
    if (n == 0 || wavefunction.ao_overlap.size() != n * n ||
        wavefunction.total_density_packed.empty()) {
        return;
    }

    const auto basis_atom = basis_atom_map(wavefunction);
    const auto total_density = unpack_symmetric(
        wavefunction.total_density_packed, n, "Total density");
    const auto total_ps = multiply(total_density, wavefunction.ao_overlap, n);

    std::vector<double> spin_ps;
    if (!wavefunction.spin_density_packed.empty()) {
        const auto spin_density = unpack_symmetric(
            wavefunction.spin_density_packed, n, "Spin density");
        spin_ps = multiply(spin_density, wavefunction.ao_overlap, n);
    }

    // Mayer pair index. For an unrestricted wavefunction the total/spin form
    // B_AB = sum[(PS)_mn(PS)_nm + (QS)_mn(QS)_nm] is equivalent to the
    // conventional 2*(alpha + beta) expression and reduces continuously to the
    // ordinary closed-shell formula when Q=0.
    for (std::size_t atom_a = 0; atom_a < wavefunction.atoms.size(); ++atom_a) {
        for (std::size_t atom_b = atom_a + 1; atom_b < wavefunction.atoms.size(); ++atom_b) {
            double order = 0.0;
            for (std::size_t mu = 0; mu < n; ++mu) {
                if (basis_atom[mu] != atom_a) continue;
                for (std::size_t nu = 0; nu < n; ++nu) {
                    if (basis_atom[nu] != atom_b) continue;
                    order += total_ps[mu * n + nu] * total_ps[nu * n + mu];
                    if (!spin_ps.empty()) {
                        order += spin_ps[mu * n + nu] * spin_ps[nu * n + mu];
                    }
                }
            }
            if (std::abs(order) >= options.record_abs_mayer_floor) {
                BondOrderRecord record;
                record.atom_a = static_cast<std::uint32_t>(atom_a);
                record.atom_b = static_cast<std::uint32_t>(atom_b);
                record.mayer_order = order;
                record.provenance = DataProvenance::Derived;
                wavefunction.bond_orders.push_back(record);
            }
        }
    }
    wavefunction.bond_order_provenance = DataProvenance::Derived;

    // Occupied-MO atomic participation. q_A = sum_{mu in A} C_mu (S C)_mu
    // is a Mulliken-like gross MO population and sums to one for an exactly
    // normalized orbital. Absolute q_A values are normalized only for the
    // purpose of flagging spatially shared candidates; the signed values are
    // not reinterpreted as charges or bond orders.
    for (std::size_t orbital_index = 0;
         orbital_index < wavefunction.orbitals.size(); ++orbital_index) {
        const auto& mo = wavefunction.orbitals[orbital_index];
        if (mo.occupation <= 1.0e-4f || mo.coefficients.size() != n) continue;

        std::vector<double> sc(n, 0.0);
        for (std::size_t mu = 0; mu < n; ++mu) {
            double value = 0.0;
            for (std::size_t nu = 0; nu < n; ++nu) {
                value += wavefunction.ao_overlap[mu * n + nu] *
                         static_cast<double>(mo.coefficients[nu]);
            }
            sc[mu] = value;
        }

        std::vector<double> atom_weight(wavefunction.atoms.size(), 0.0);
        for (std::size_t mu = 0; mu < n; ++mu) {
            atom_weight[basis_atom[mu]] +=
                static_cast<double>(mo.coefficients[mu]) * sc[mu];
        }

        double abs_sum = 0.0;
        for (double& value : atom_weight) {
            value = std::abs(value);
            abs_sum += value;
        }
        if (!(abs_sum > 1.0e-12) || !std::isfinite(abs_sum)) continue;
        for (double& value : atom_weight) value /= abs_sum;

        std::vector<std::pair<double, std::uint32_t>> significant;
        double covered = 0.0;
        for (std::size_t atom = 0; atom < atom_weight.size(); ++atom) {
            if (atom_weight[atom] >= options.multicentre_atom_participation) {
                significant.emplace_back(atom_weight[atom],
                                         static_cast<std::uint32_t>(atom));
                covered += atom_weight[atom];
            }
        }
        if (significant.size() < 3 ||
            covered < options.multicentre_min_covered_weight) {
            continue;
        }

        std::sort(significant.begin(), significant.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        MulticentreCandidate candidate;
        candidate.orbital_index = static_cast<std::uint32_t>(orbital_index);
        candidate.occupation = static_cast<double>(mo.occupation);
        candidate.provenance = DataProvenance::Derived;
        candidate.atoms.reserve(significant.size());
        candidate.participation.reserve(significant.size());
        for (const auto& [weight, atom] : significant) {
            candidate.atoms.push_back(atom);
            candidate.participation.push_back(weight);
        }
        wavefunction.multicentre_candidates.push_back(std::move(candidate));
    }
}

} // namespace cov
