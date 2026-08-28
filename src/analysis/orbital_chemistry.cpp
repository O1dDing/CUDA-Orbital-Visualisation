#include "cov/orbital_chemistry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace cov {
namespace {

enum class ReferenceClass : std::uint8_t {
    DeepCore = 0,
    Semicore,
    ChemicalValence,
};

struct ReferenceShell {
    int n = 0;
    int l = 0;
    int electrons = 0;
    ReferenceClass space = ReferenceClass::DeepCore;
    bool removed_by_ecp = false;
};

struct ReferenceColumn {
    std::vector<double> coefficients;
    std::uint32_t atom = 0;
    int n = 0;
    int l = 0;
    int component = 0;
    ReferenceClass space = ReferenceClass::DeepCore;
    std::string label;
};

struct AufbauEntry {
    int n;
    int l;
    int capacity;
};

constexpr std::array<AufbauEntry, 19> kAufbau{{
    {1,0,2}, {2,0,2}, {2,1,6}, {3,0,2}, {3,1,6},
    {4,0,2}, {3,2,10}, {4,1,6}, {5,0,2}, {4,2,10},
    {5,1,6}, {6,0,2}, {4,3,14}, {5,2,10}, {6,1,6},
    {7,0,2}, {5,3,14}, {6,2,10}, {7,1,6},
}};

int period_for_z(const int z) {
    if (z <= 0) return 0;
    if (z <= 2) return 1;
    if (z <= 10) return 2;
    if (z <= 18) return 3;
    if (z <= 36) return 4;
    if (z <= 54) return 5;
    if (z <= 86) return 6;
    return 7;
}

bool transition_metal(const int z) {
    return (z >= 21 && z <= 30) || (z >= 39 && z <= 48) ||
           (z >= 72 && z <= 80) || (z >= 104 && z <= 112);
}

bool f_block(const int z) {
    return (z >= 57 && z <= 71) || (z >= 89 && z <= 103);
}

using ShellKey = std::pair<int,int>;

std::map<ShellKey,int> neutral_configuration(const int z) {
    std::map<ShellKey,int> result;
    int remaining = std::max(0, z);
    for (const auto& shell : kAufbau) {
        if (remaining <= 0) break;
        const int fill = std::min(shell.capacity, remaining);
        result[{shell.n,shell.l}] = fill;
        remaining -= fill;
    }
    return result;
}

std::vector<ReferenceShell> atomic_reference_shells(const Atom& atom) {
    const int z = atom.atomic_number;
    const int period = period_for_z(z);
    auto occupation = neutral_configuration(z);

    int explicit_charge = z;
    if (std::isfinite(atom.nuclear_charge) && atom.nuclear_charge > 0.0) {
        const int rounded = static_cast<int>(std::llround(atom.nuclear_charge));
        if (rounded > 0 && rounded <= z) explicit_charge = rounded;
    }

    int remove = std::max(0, z - explicit_charge);
    std::vector<ShellKey> occupied;
    for (const auto& [key, electrons] : occupation) {
        if (electrons > 0) occupied.push_back(key);
    }
    std::sort(occupied.begin(), occupied.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    });

    std::set<ShellKey> removed;
    for (const auto& key : occupied) {
        if (remove <= 0) break;
        int& electrons = occupation[key];
        const int take = std::min(remove, electrons);
        electrons -= take;
        remove -= take;
        if (electrons == 0) removed.insert(key);
    }

    std::set<ShellKey> valence;
    std::set<ShellKey> semicore;
    if (period == 1) {
        valence.insert({1,0});
    } else if (period > 1) {
        valence.insert({period,0});
        valence.insert({period,1});
        if (transition_metal(z)) {
            valence.insert({period-1,2});
            semicore.insert({period-1,0});
            semicore.insert({period-1,1});
        } else if (f_block(z)) {
            valence.insert({period-2,3});
            valence.insert({period-1,2});
            semicore.insert({period-1,0});
            semicore.insert({period-1,1});
        } else if (period >= 4) {
            semicore.insert({period-1,2});
        }
    }

    std::set<ShellKey> keys;
    for (const auto& [key, electrons] : occupation) {
        if (electrons > 0 || removed.count(key) != 0u) keys.insert(key);
    }
    keys.insert(valence.begin(), valence.end());
    keys.insert(semicore.begin(), semicore.end());

    std::vector<ReferenceShell> result;
    result.reserve(keys.size());
    for (const auto& key : keys) {
        ReferenceShell shell;
        shell.n = key.first;
        shell.l = key.second;
        shell.electrons = occupation[key];
        shell.removed_by_ecp = removed.count(key) != 0u;
        if (valence.count(key) != 0u) shell.space = ReferenceClass::ChemicalValence;
        else if (semicore.count(key) != 0u) shell.space = ReferenceClass::Semicore;
        else shell.space = ReferenceClass::DeepCore;
        result.push_back(shell);
    }
    return result;
}

double shell_effective_exponent(const Wavefunction& wf, const Shell& shell) {
    const std::size_t begin = shell.primitive_offset;
    const std::size_t count = shell.primitive_count;
    if (count == 0u || begin + count > wf.primitives.size()) return 0.0;
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t p = 0; p < count; ++p) {
        const auto& primitive = wf.primitives[begin+p];
        const double weight = std::max(1.0e-12,
            std::abs(static_cast<double>(primitive.coefficient)));
        numerator += weight * std::max(1.0e-20,
            static_cast<double>(primitive.exponent));
        denominator += weight;
    }
    return denominator > 0.0 ? numerator / denominator : 0.0;
}

const char* angular_letter(const int l) {
    static constexpr const char* letters[] = {"s","p","d","f","g","h"};
    return l >= 0 && l < 6 ? letters[l] : "?";
}

std::string reference_label(const Wavefunction& wf,
                            const std::uint32_t atom,
                            const int n,
                            const int l) {
    std::ostringstream out;
    if (atom < wf.atoms.size()) out << wf.atoms[atom].symbol;
    else out << 'X';
    out << (atom + 1u) << ' ' << n << angular_letter(l);
    return out.str();
}

double s_dot(const std::vector<double>& a,
             const std::vector<double>& b,
             const std::vector<double>& overlap,
             const std::size_t n) {
    double result = 0.0;
    for (std::size_t mu = 0; mu < n; ++mu) {
        if (a[mu] == 0.0) continue;
        double sb = 0.0;
        for (std::size_t nu = 0; nu < n; ++nu) {
            sb += overlap[mu*n+nu] * b[nu];
        }
        result += a[mu] * sb;
    }
    return result;
}

double mo_s_dot(const std::vector<double>& a,
                const MolecularOrbital& mo,
                const std::vector<double>& overlap,
                const std::size_t n) {
    double result = 0.0;
    for (std::size_t mu = 0; mu < n; ++mu) {
        if (a[mu] == 0.0) continue;
        double sc = 0.0;
        for (std::size_t nu = 0; nu < n; ++nu) {
            sc += overlap[mu*n+nu] *
                  static_cast<double>(mo.coefficients[nu]);
        }
        result += a[mu] * sc;
    }
    return result;
}

std::vector<std::size_t> largest_gap_boundaries(
    const std::vector<std::pair<std::size_t,double>>& shells,
    const std::size_t group_count) {
    std::vector<std::pair<double,std::size_t>> gaps;
    for (std::size_t i = 0; i + 1u < shells.size(); ++i) {
        const double a = std::log(std::max(1.0e-30, shells[i].second));
        const double b = std::log(std::max(1.0e-30, shells[i+1u].second));
        gaps.push_back({a-b,i+1u});
    }
    std::sort(gaps.begin(), gaps.end(), [](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second < b.second;
    });
    std::vector<std::size_t> boundaries;
    const std::size_t wanted = group_count > 0u ? group_count-1u : 0u;
    for (std::size_t i=0;i<std::min(wanted,gaps.size());++i) {
        boundaries.push_back(gaps[i].second);
    }
    std::sort(boundaries.begin(),boundaries.end());
    return boundaries;
}

std::vector<double> radial_weights(
    const Wavefunction& wf,
    const std::vector<std::size_t>& shell_indices) {
    std::vector<double> result(shell_indices.size(),0.0);
    if (shell_indices.empty()) return result;
    double max_log = -std::numeric_limits<double>::infinity();
    for (const auto s : shell_indices) {
        max_log = std::max(max_log,
            std::log(std::max(1.0e-30,
                shell_effective_exponent(wf,wf.shells[s]))));
    }
    double norm = 0.0;
    for (std::size_t i=0;i<shell_indices.size();++i) {
        const double loga = std::log(std::max(1.0e-30,
            shell_effective_exponent(wf,wf.shells[shell_indices[i]])));
        result[i] = std::exp(0.65 * (loga-max_log));
        norm += result[i]*result[i];
    }
    if (norm > 0.0) {
        const double inv = 1.0/std::sqrt(norm);
        for (double& v : result) v *= inv;
    }
    return result;
}

void append_standard_components(
    const Wavefunction& wf,
    const std::vector<std::size_t>& shell_indices,
    const std::vector<double>& radial,
    const ReferenceShell& ref,
    const std::uint32_t atom,
    std::vector<ReferenceColumn>& raw) {
    if (shell_indices.empty() || radial.size()!=shell_indices.size()) return;
    const Shell& first = wf.shells[shell_indices.front()];
    const int l = static_cast<int>(first.angular_momentum);
    for (const auto s : shell_indices) {
        if (wf.shells[s].angular_momentum != first.angular_momentum ||
            wf.shells[s].pure != first.pure ||
            shell_basis_count(wf.shells[s]) != shell_basis_count(first)) {
            return;
        }
    }

    const std::size_t nbf = wf.basis_count;
    const std::size_t local_count = shell_basis_count(first);

    auto make_direct = [&](const int standard_component,
                           const std::vector<std::pair<int,double>>& local) {
        ReferenceColumn column;
        column.coefficients.assign(nbf,0.0);
        column.atom=atom;
        column.n=ref.n;
        column.l=l;
        column.component=standard_component;
        column.space=ref.space;
        column.label=reference_label(wf,atom,ref.n,l);
        for (std::size_t r=0;r<shell_indices.size();++r) {
            const Shell& shell=wf.shells[shell_indices[r]];
            for (const auto& [component,coefficient] : local) {
                if (component < 0 ||
                    static_cast<std::size_t>(component)>=local_count) continue;
                column.coefficients[shell.basis_offset+
                    static_cast<std::size_t>(component)] += radial[r]*coefficient;
            }
        }
        raw.push_back(std::move(column));
    };

    if (l <= 1 || first.pure != 0u) {
        const std::size_t chemical_count =
            l >= 0 ? static_cast<std::size_t>(2*l+1) : 0u;
        const std::size_t count = std::min(local_count,chemical_count);
        for (std::size_t c=0;c<count;++c) {
            make_direct(static_cast<int>(c),{{static_cast<int>(c),1.0}});
        }
        return;
    }

    if (l == 2 && local_count == 6u) {
        make_direct(0,{{0,-1.0/std::sqrt(6.0)},
                       {1,-1.0/std::sqrt(6.0)},
                       {2, 2.0/std::sqrt(6.0)}});
        make_direct(1,{{4,1.0}});
        make_direct(2,{{5,1.0}});
        make_direct(3,{{0, 1.0/std::sqrt(2.0)},
                       {1,-1.0/std::sqrt(2.0)}});
        make_direct(4,{{3,1.0}});
    }
}

std::vector<ReferenceColumn> build_raw_reference(const Wavefunction& wf) {
    std::vector<ReferenceColumn> raw;
    for (std::uint32_t atom=0;atom<wf.atoms.size();++atom) {
        const auto refs_all=atomic_reference_shells(wf.atoms[atom]);
        for (int l=0;l<=4;++l) {
            std::vector<ReferenceShell> refs;
            for (const auto& ref:refs_all) {
                if (!ref.removed_by_ecp && ref.l==l) refs.push_back(ref);
            }
            if (refs.empty()) continue;
            std::sort(refs.begin(),refs.end(),[](const auto& a,const auto& b) {
                return a.n<b.n;
            });

            std::vector<std::pair<std::size_t,double>> working;
            for (std::size_t s=0;s<wf.shells.size();++s) {
                const auto& shell=wf.shells[s];
                if (shell.atom_index==atom &&
                    static_cast<int>(shell.angular_momentum)==l) {
                    working.push_back({s,shell_effective_exponent(wf,shell)});
                }
            }
            if (working.empty()) continue;
            std::sort(working.begin(),working.end(),[](const auto& a,const auto& b) {
                if (a.second!=b.second) return a.second>b.second;
                return a.first<b.first;
            });

            const std::size_t groups=std::min(refs.size(),working.size());
            if (groups==0u) continue;
            const auto boundaries=largest_gap_boundaries(working,groups);
            std::size_t begin=0;
            for (std::size_t group=0;group<groups;++group) {
                const std::size_t end =
                    group<boundaries.size()?boundaries[group]:working.size();
                std::vector<std::size_t> shell_indices;
                for (std::size_t i=begin;i<end;++i) {
                    shell_indices.push_back(working[i].first);
                }
                if (shell_indices.empty() && begin<working.size()) {
                    shell_indices.push_back(working[begin].first);
                }
                const auto radial=radial_weights(wf,shell_indices);
                const std::size_t ref_offset=refs.size()-groups;
                append_standard_components(
                    wf,shell_indices,radial,refs[ref_offset+group],atom,raw);
                begin=end;
            }
        }
    }
    return raw;
}

int reference_class_order(const ReferenceClass value) {
    switch (value) {
        case ReferenceClass::DeepCore: return 0;
        case ReferenceClass::Semicore: return 1;
        default: return 2;
    }
}

std::vector<ReferenceColumn> orthonormalise_reference(
    std::vector<ReferenceColumn> raw,
    const Wavefunction& wf) {
    const std::size_t n=wf.basis_count;
    std::stable_sort(raw.begin(),raw.end(),[](const auto& a,const auto& b) {
        const int ca=reference_class_order(a.space);
        const int cb=reference_class_order(b.space);
        if (ca!=cb) return ca<cb;
        if (a.atom!=b.atom) return a.atom<b.atom;
        if (a.n!=b.n) return a.n<b.n;
        if (a.l!=b.l) return a.l<b.l;
        return a.component<b.component;
    });

    std::vector<ReferenceColumn> q;
    for (auto& candidate:raw) {
        if (candidate.coefficients.size()!=n) continue;
        for (const auto& previous:q) {
            const double projection=s_dot(
                previous.coefficients,candidate.coefficients,wf.ao_overlap,n);
            for (std::size_t mu=0;mu<n;++mu) {
                candidate.coefficients[mu]-=
                    projection*previous.coefficients[mu];
            }
        }
        const double norm2=s_dot(
            candidate.coefficients,candidate.coefficients,wf.ao_overlap,n);
        if (!(norm2>1.0e-10) || !std::isfinite(norm2)) continue;
        const double inv=1.0/std::sqrt(norm2);
        for (double& v:candidate.coefficients) v*=inv;
        q.push_back(std::move(candidate));
    }
    return q;
}

std::vector<std::uint32_t> basis_atom_map(const Wavefunction& wf) {
    const auto missing=std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> result(wf.basis_count,missing);
    for (const auto& shell:wf.shells) {
        const std::size_t count=shell_basis_count(shell);
        if (shell.basis_offset+count>wf.basis_count) return {};
        for (std::size_t c=0;c<count;++c) {
            result[shell.basis_offset+c]=shell.atom_index;
        }
    }
    if (std::find(result.begin(),result.end(),missing)!=result.end()) return {};
    return result;
}

std::vector<double> atomic_weights(
    const Wavefunction& wf,
    const MolecularOrbital& mo,
    const std::vector<std::uint32_t>& basis_atom) {
    const std::size_t n=wf.basis_count;
    std::vector<double> result(wf.atoms.size(),0.0);
    if (mo.coefficients.size()!=n ||
        wf.ao_overlap.size()!=n*n ||
        basis_atom.size()!=n) return result;
    std::vector<double> sc(n,0.0);
    for (std::size_t mu=0;mu<n;++mu) {
        for (std::size_t nu=0;nu<n;++nu) {
            sc[mu]+=wf.ao_overlap[mu*n+nu]*
                static_cast<double>(mo.coefficients[nu]);
        }
    }
    double total=0.0;
    for (std::size_t mu=0;mu<n;++mu) {
        result[basis_atom[mu]]+=
            static_cast<double>(mo.coefficients[mu])*sc[mu];
    }
    for (double& v:result) {
        v=std::abs(v);
        total+=v;
    }
    if (total>1.0e-12) for (double& v:result) v/=total;
    return result;
}

struct DegenerateGroup {
    Spin spin=Spin::Alpha;
    std::vector<std::size_t> orbitals;
    double score=0.0;
};

std::vector<DegenerateGroup> degeneracy_groups(
    const Wavefunction& wf,
    const double tolerance) {
    std::vector<DegenerateGroup> groups;
    std::size_t begin=0;
    while (begin<wf.orbitals.size()) {
        DegenerateGroup group;
        group.spin=wf.orbitals[begin].spin;
        group.orbitals.push_back(begin);
        std::size_t end=begin+1u;
        while (end<wf.orbitals.size()) {
            const auto& first=wf.orbitals[begin];
            const auto& current=wf.orbitals[end];
            if (current.spin!=first.spin) break;
            if (std::abs(current.energy_hartree-first.energy_hartree)>tolerance) break;
            if (!first.symmetry.empty() && !current.symmetry.empty() &&
                first.symmetry!=current.symmetry) break;
            group.orbitals.push_back(end);
            ++end;
        }
        groups.push_back(std::move(group));
        begin=end;
    }
    return groups;
}

void select_valence_manifold(
    Wavefunction& wf,
    const std::size_t reference_rank,
    const double tolerance) {
    if (reference_rank==0u || wf.orbitals.empty()) return;
    auto groups=degeneracy_groups(wf,tolerance);
    for (auto& group:groups) {
        for (const auto i:group.orbitals) {
            group.score+=wf.orbitals[i].chemistry.valence_weight;
        }
    }

    for (const Spin spin : {Spin::Alpha,Spin::Beta}) {
        std::vector<std::size_t> candidates;
        std::size_t orbital_count=0;
        for (std::size_t g=0;g<groups.size();++g) {
            if (groups[g].spin==spin) {
                candidates.push_back(g);
                orbital_count+=groups[g].orbitals.size();
            }
        }
        if (candidates.empty()) continue;
        const std::size_t target=std::min(reference_rank,orbital_count);
        const double neg=-std::numeric_limits<double>::infinity();
        std::vector<double> dp(orbital_count+1u,neg);
        std::vector<std::vector<std::size_t>> chosen(orbital_count+1u);
        dp[0]=0.0;
        for (const auto gi:candidates) {
            const std::size_t size=groups[gi].orbitals.size();
            const double value=groups[gi].score;
            for (std::size_t used=orbital_count+1u;used-- > size;) {
                if (!std::isfinite(dp[used-size])) continue;
                const double trial=dp[used-size]+value;
                if (trial>dp[used]+1.0e-12) {
                    dp[used]=trial;
                    chosen[used]=chosen[used-size];
                    chosen[used].push_back(gi);
                }
            }
        }

        std::size_t best=target;
        if (!std::isfinite(dp[best])) {
            std::size_t distance=orbital_count+1u;
            double best_score=neg;
            for (std::size_t count=0;count<=orbital_count;++count) {
                if (!std::isfinite(dp[count])) continue;
                const std::size_t d=count>target?count-target:target-count;
                if (d<distance || (d==distance && dp[count]>best_score)) {
                    distance=d;
                    best_score=dp[count];
                    best=count;
                }
            }
        }
        for (const auto gi:chosen[best]) {
            for (const auto i:groups[gi].orbitals) {
                wf.orbitals[i].chemistry.valence_manifold=true;
            }
        }
    }
}

double pair_mayer(const Wavefunction& wf,
                  const std::uint32_t a,
                  const std::uint32_t b) {
    for (const auto& record:wf.bond_orders) {
        if ((record.atom_a==a && record.atom_b==b) ||
            (record.atom_a==b && record.atom_b==a)) {
            return record.mayer_order;
        }
    }
    return 0.0;
}

double distance_bohr(const Atom& a,const Atom& b) {
    const double dx=a.x-b.x;
    const double dy=a.y-b.y;
    const double dz=a.z-b.z;
    return std::sqrt(dx*dx+dy*dy+dz*dz);
}

std::array<double,3> unit_axis(const Atom& a,const Atom& b) {
    std::array<double,3> z{b.x-a.x,b.y-a.y,b.z-a.z};
    const double n=std::sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
    if (n>1.0e-14) for (double& v:z) v/=n;
    return z;
}

std::array<double,3> cross3(
    const std::array<double,3>& a,
    const std::array<double,3>& b) {
    return {a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
}

double dot3(const std::array<double,3>& a,
            const std::array<double,3>& b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

struct RawChannels {
    double sigma=0.0;
    double pi=0.0;
    double delta=0.0;
    double phi=0.0;
    double undetermined=0.0;
};

using Mat3=std::array<double,9>;

Mat3 mat_add_scaled(Mat3 a,const Mat3& b,const double scale) {
    for (std::size_t i=0;i<9u;++i) a[i]+=scale*b[i];
    return a;
}

double frobenius(const Mat3& a,const Mat3& b) {
    double result=0.0;
    for (std::size_t i=0;i<9u;++i) result+=a[i]*b[i];
    return result;
}

Mat3 transform_tensor(const Mat3& t,
                      const std::array<double,3>& ex,
                      const std::array<double,3>& ey,
                      const std::array<double,3>& ez) {
    const std::array<std::array<double,3>,3> e{ex,ey,ez};
    Mat3 out{};
    for (int i=0;i<3;++i) {
        for (int j=0;j<3;++j) {
            double value=0.0;
            for (int a=0;a<3;++a) {
                for (int b=0;b<3;++b) {
                    value+=e[i][a]*t[3*a+b]*e[j][b];
                }
            }
            out[3*i+j]=value;
        }
    }
    return out;
}

void add_atom_channels(
    const std::uint32_t atom,
    const std::array<double,3>& axis,
    const std::vector<ReferenceColumn>& reference,
    const std::vector<double>& amplitudes,
    RawChannels& out) {
    std::map<std::pair<int,int>,std::vector<double>> by_shell;
    for (std::size_t r=0;r<reference.size();++r) {
        const auto& ref=reference[r];
        if (ref.atom!=atom || ref.space!=ReferenceClass::ChemicalValence) continue;
        auto& values=by_shell[{ref.n,ref.l}];
        const std::size_t needed=static_cast<std::size_t>(
            std::max(0,2*ref.l+1));
        if (values.size()<needed) values.resize(needed,0.0);
        if (ref.component>=0 &&
            static_cast<std::size_t>(ref.component)<values.size()) {
            values[static_cast<std::size_t>(ref.component)]=amplitudes[r];
        }
    }

    std::array<double,3> ez=axis;
    std::array<double,3> seed=
        std::abs(ez[2])<0.85?std::array<double,3>{0.0,0.0,1.0}
                            :std::array<double,3>{1.0,0.0,0.0};
    auto ex=cross3(seed,ez);
    double exn=std::sqrt(dot3(ex,ex));
    if (exn<1.0e-12) {
        out.undetermined+=1.0;
        return;
    }
    for (double& v:ex) v/=exn;
    auto ey=cross3(ez,ex);

    for (const auto& [key,values]:by_shell) {
        const int l=key.second;
        double norm=0.0;
        for (const double v:values) norm+=v*v;
        if (norm<=1.0e-14) continue;
        if (l==0) {
            out.sigma+=norm;
        } else if (l==1 && values.size()>=3u) {
            const std::array<double,3> p{values[0],values[1],values[2]};
            const double along=dot3(p,ez);
            const double sigma=along*along;
            out.sigma+=sigma;
            out.pi+=std::max(0.0,norm-sigma);
        } else if (l==2 && values.size()>=5u) {
            const double inv_sqrt6=1.0/std::sqrt(6.0);
            const double inv_sqrt2=1.0/std::sqrt(2.0);
            const Mat3 b0{-inv_sqrt6,0,0, 0,-inv_sqrt6,0, 0,0,2*inv_sqrt6};
            const Mat3 b1{0,0,inv_sqrt2, 0,0,0, inv_sqrt2,0,0};
            const Mat3 b2{0,0,0, 0,0,inv_sqrt2, 0,inv_sqrt2,0};
            const Mat3 b3{inv_sqrt2,0,0, 0,-inv_sqrt2,0, 0,0,0};
            const Mat3 b4{0,inv_sqrt2,0, inv_sqrt2,0,0, 0,0,0};
            Mat3 tensor{};
            tensor=mat_add_scaled(tensor,b0,values[0]);
            tensor=mat_add_scaled(tensor,b1,values[1]);
            tensor=mat_add_scaled(tensor,b2,values[2]);
            tensor=mat_add_scaled(tensor,b3,values[3]);
            tensor=mat_add_scaled(tensor,b4,values[4]);
            const Mat3 local=transform_tensor(tensor,ex,ey,ez);
            const double c0=frobenius(local,b0);
            const double c1=frobenius(local,b1);
            const double c2=frobenius(local,b2);
            const double c3=frobenius(local,b3);
            const double c4=frobenius(local,b4);
            out.sigma+=c0*c0;
            out.pi+=c1*c1+c2*c2;
            out.delta+=c3*c3+c4*c4;
        } else {
            if (std::abs(std::abs(ez[2])-1.0)<1.0e-6) {
                if (!values.empty()) out.sigma+=values[0]*values[0];
                if (values.size()>2u) out.pi+=values[1]*values[1]+values[2]*values[2];
                if (values.size()>4u) out.delta+=values[3]*values[3]+values[4]*values[4];
                if (values.size()>6u) out.phi+=values[5]*values[5]+values[6]*values[6];
            } else {
                out.undetermined+=norm;
            }
        }
    }
}

OrbitalChannelDistribution finish_channels(
    const RawChannels& raw,
    const OrbitalChemistryOptions& options) {
    OrbitalChannelDistribution result;
    const double total=raw.sigma+raw.pi+raw.delta+raw.phi+raw.undetermined;
    if (!(total>1.0e-14) || !std::isfinite(total)) {
        result.status=ChemistryStatus::Undetermined;
        return result;
    }
    result.sigma=raw.sigma/total;
    result.pi=raw.pi/total;
    result.delta=raw.delta/total;
    result.phi=raw.phi/total;
    result.undetermined=raw.undetermined/total;

    const std::array<std::pair<double,OrbitalAngularFamily>,4> values{{
        {result.sigma,OrbitalAngularFamily::Sigma},
        {result.pi,OrbitalAngularFamily::Pi},
        {result.delta,OrbitalAngularFamily::Delta},
        {result.phi,OrbitalAngularFamily::Phi},
    }};
    const auto best=std::max_element(values.begin(),values.end(),
        [](const auto& a,const auto& b){return a.first<b.first;});
    result.dominant=best->second;
    if (best->first>=options.determined_fraction &&
        result.undetermined<=0.10) {
        result.status=ChemistryStatus::Determined;
    } else if (result.undetermined<=
               options.maximum_undetermined_for_percentages) {
        result.status=ChemistryStatus::Percentages;
    } else {
        result.status=ChemistryStatus::Undetermined;
        result.dominant=OrbitalAngularFamily::Mixed;
    }
    return result;
}

OrbitalBondingDistribution pair_bonding(
    const double overlap_character,
    const bool applicable,
    const OrbitalChemistryOptions& options) {
    OrbitalBondingDistribution result;
    if (!applicable) {
        result.status=ChemistryStatus::NotApplicable;
        result.dominant=OrbitalBondingRole::NotApplicable;
        result.undetermined=0.0;
        return result;
    }
    if (overlap_character>=options.overlap_bonding_threshold) {
        result.bonding=1.0;
        result.undetermined=0.0;
        result.dominant=OrbitalBondingRole::Bonding;
        result.status=ChemistryStatus::Determined;
    } else if (overlap_character<=-options.overlap_bonding_threshold) {
        result.antibonding=1.0;
        result.undetermined=0.0;
        result.dominant=OrbitalBondingRole::Antibonding;
        result.status=ChemistryStatus::Determined;
    } else if (std::abs(overlap_character)<=
               options.nonbonding_overlap_threshold) {
        result.nonbonding=1.0;
        result.undetermined=0.0;
        result.dominant=OrbitalBondingRole::Nonbonding;
        result.status=ChemistryStatus::Determined;
    } else {
        const double magnitude=std::abs(overlap_character);
        const double fraction=std::min(0.80,
            magnitude/options.overlap_bonding_threshold);
        if (overlap_character>0.0) result.bonding=fraction;
        else result.antibonding=fraction;
        result.nonbonding=1.0-fraction;
        result.undetermined=0.0;
        result.dominant=overlap_character>0.0
            ?OrbitalBondingRole::Bonding
            :OrbitalBondingRole::Antibonding;
        result.status=ChemistryStatus::Percentages;
    }
    return result;
}

OrbitalBondingDistribution aggregate_bonding(
    const std::vector<OrbitalPairInteraction>& interactions) {
    OrbitalBondingDistribution result;
    double total=0.0;
    for (const auto& interaction:interactions) {
        if (interaction.bonding.status==ChemistryStatus::NotApplicable) continue;
        const double weight=std::max(0.01,
            std::abs(interaction.overlap_character));
        result.bonding+=weight*interaction.bonding.bonding;
        result.antibonding+=weight*interaction.bonding.antibonding;
        result.nonbonding+=weight*interaction.bonding.nonbonding;
        result.undetermined+=weight*interaction.bonding.undetermined;
        total+=weight;
    }
    if (!(total>0.0)) {
        result.status=ChemistryStatus::Undetermined;
        result.undetermined=1.0;
        return result;
    }
    result.bonding/=total;
    result.antibonding/=total;
    result.nonbonding/=total;
    result.undetermined/=total;
    const std::array<std::pair<double,OrbitalBondingRole>,3> values{{
        {result.bonding,OrbitalBondingRole::Bonding},
        {result.antibonding,OrbitalBondingRole::Antibonding},
        {result.nonbonding,OrbitalBondingRole::Nonbonding},
    }};
    const auto best=std::max_element(values.begin(),values.end(),
        [](const auto& a,const auto& b){return a.first<b.first;});
    result.dominant=best->second;
    result.status=best->first>=0.82
        ?ChemistryStatus::Determined
        :ChemistryStatus::Percentages;
    return result;
}

OrbitalChannelDistribution aggregate_channels(
    const std::vector<OrbitalPairInteraction>& interactions,
    const OrbitalChemistryOptions& options) {
    RawChannels raw;
    for (const auto& interaction:interactions) {
        const double weight=std::max(0.01,
            std::abs(interaction.overlap_character));
        raw.sigma+=weight*interaction.channel.sigma;
        raw.pi+=weight*interaction.channel.pi;
        raw.delta+=weight*interaction.channel.delta;
        raw.phi+=weight*interaction.channel.phi;
        raw.undetermined+=weight*interaction.channel.undetermined;
    }
    return finish_channels(raw,options);
}

std::string family_symbol(const OrbitalAngularFamily family) {
    switch (family) {
        case OrbitalAngularFamily::Sigma: return "sigma";
        case OrbitalAngularFamily::Pi: return "pi";
        case OrbitalAngularFamily::Delta: return "delta";
        case OrbitalAngularFamily::Phi: return "phi";
        default: return "mixed";
    }
}

bool near(const double value,const double target,const double tolerance) {
    return std::abs(value-target)<=tolerance;
}

void attach_multicentre_assignments(Wavefunction& wf) {
    for (const auto& assignment:wf.multicentre_assignments) {
        std::string label;
        if (assignment.kind==MulticentreKind::ThreeCentreTwoElectron) label="3c2e";
        else if (assignment.kind==MulticentreKind::ThreeCentreFourElectron) label="3c4e";
        else continue;
        for (const auto index:assignment.orbitals) {
            if (index>=wf.orbitals.size()) continue;
            auto& chemistry=wf.orbitals[index].chemistry;
            chemistry.multicentre_label=label;
            chemistry.participating_atoms=assignment.atoms.size();
            chemistry.participating_electrons=assignment.electron_count;
            chemistry.family_symbol=family_symbol(chemistry.channel.dominant);
        }
    }
}

void generic_three_centre_fallback(Wavefunction& wf) {
    std::vector<std::size_t> selected;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        if (wf.orbitals[i].spin==Spin::Alpha &&
            wf.orbitals[i].chemistry.valence_manifold) selected.push_back(i);
    }
    if (selected.size()!=3u) return;

    std::set<std::uint32_t> atoms;
    double electrons=0.0;
    bool sigma=true;
    for (const auto i:selected) {
        electrons+=static_cast<double>(wf.orbitals[i].occupation);
        const auto& chemistry=wf.orbitals[i].chemistry;
        sigma=sigma &&
            chemistry.channel.sigma>=0.70 &&
            chemistry.channel.undetermined<=0.25;
        for (const auto& contribution:chemistry.ao_contributions) {
            if (contribution.weight>=0.05) atoms.insert(contribution.atom_index);
        }
    }
    if (atoms.size()!=3u || !sigma) return;
    int supported=0;
    std::vector<std::uint32_t> centres(atoms.begin(),atoms.end());
    for (std::size_t a=0;a<centres.size();++a) {
        for (std::size_t b=a+1u;b<centres.size();++b) {
            if (std::abs(pair_mayer(wf,centres[a],centres[b]))>=0.02) ++supported;
        }
    }
    if (supported<2) return;

    std::string label;
    if (near(electrons,2.0,0.20)) label="3c2e";
    else if (near(electrons,4.0,0.20)) label="3c4e";
    else return;
    for (const auto i:selected) {
        auto& chemistry=wf.orbitals[i].chemistry;
        if (!chemistry.multicentre_label.empty()) continue;
        chemistry.multicentre_label=label;
        chemistry.participating_atoms=3u;
        chemistry.participating_electrons=electrons;
        chemistry.family_symbol="sigma";
    }
}

} // namespace

const char* chemistry_status_name(const ChemistryStatus value) noexcept {
    switch (value) {
        case ChemistryStatus::Determined: return "determined";
        case ChemistryStatus::Percentages: return "percentages";
        case ChemistryStatus::Undetermined: return "undetermined";
        case ChemistryStatus::NotApplicable: return "not-applicable";
        default: return "unavailable";
    }
}

const char* orbital_angular_family_name(
    const OrbitalAngularFamily value) noexcept {
    switch (value) {
        case OrbitalAngularFamily::Sigma: return "sigma";
        case OrbitalAngularFamily::Pi: return "pi";
        case OrbitalAngularFamily::Delta: return "delta";
        case OrbitalAngularFamily::Phi: return "phi";
        case OrbitalAngularFamily::NotApplicable: return "not-applicable";
        default: return "mixed";
    }
}

const char* orbital_bonding_role_name(
    const OrbitalBondingRole value) noexcept {
    switch (value) {
        case OrbitalBondingRole::Bonding: return "bonding";
        case OrbitalBondingRole::Antibonding: return "antibonding";
        case OrbitalBondingRole::Nonbonding: return "nonbonding";
        case OrbitalBondingRole::NotApplicable: return "not-applicable";
        default: return "mixed";
    }
}

void derive_orbital_chemistry(
    Wavefunction& wf,
    const OrbitalChemistryOptions& options) {
    const std::size_t n=wf.basis_count;
    if (n==0u || wf.orbitals.empty() ||
        wf.ao_overlap.size()!=n*n) {
        for (auto& mo:wf.orbitals) {
            mo.chemistry.available=false;
            mo.chemistry.note="AO overlap metric unavailable";
        }
        return;
    }

    const auto raw=build_raw_reference(wf);
    const auto reference=orthonormalise_reference(raw,wf);
    if (reference.empty()) {
        for (auto& mo:wf.orbitals) {
            mo.chemistry.available=false;
            mo.chemistry.note="COV minimal atomic reference could not be built";
        }
        return;
    }

    std::size_t valence_rank=0;
    for (const auto& ref:reference) {
        if (ref.space==ReferenceClass::ChemicalValence) ++valence_rank;
    }
    const auto basis_atom=basis_atom_map(wf);
    if (basis_atom.size()!=n || valence_rank==0u) {
        for (auto& mo:wf.orbitals) {
            mo.chemistry.available=false;
            mo.chemistry.note="AO-to-atom map or chemical-valence rank unavailable";
        }
        return;
    }

    std::vector<std::vector<double>> amplitudes(
        wf.orbitals.size(),std::vector<double>(reference.size(),0.0));
    std::vector<std::vector<double>> atom_weight(wf.orbitals.size());

    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        auto& mo=wf.orbitals[i];
        auto& chemistry=mo.chemistry;
        chemistry=OrbitalChemistry{};
        chemistry.available=true;
        chemistry.method="COV FCHK S-metric minimal atomic-reference projection";
        atom_weight[i]=atomic_weights(wf,mo,basis_atom);

        std::map<std::tuple<std::uint32_t,int,int>,double> grouped;
        for (std::size_t r=0;r<reference.size();++r) {
            const double amplitude=mo_s_dot(
                reference[r].coefficients,mo,wf.ao_overlap,n);
            amplitudes[i][r]=amplitude;
            const double weight=amplitude*amplitude;
            switch (reference[r].space) {
                case ReferenceClass::DeepCore:
                    chemistry.deep_core_weight+=weight; break;
                case ReferenceClass::Semicore:
                    chemistry.semicore_weight+=weight; break;
                case ReferenceClass::ChemicalValence:
                    chemistry.valence_weight+=weight;
                    grouped[{reference[r].atom,reference[r].n,reference[r].l}]+=weight;
                    break;
            }
        }
        const double assigned=chemistry.deep_core_weight+
                              chemistry.semicore_weight+
                              chemistry.valence_weight;
        chemistry.unresolved_weight=std::clamp(1.0-assigned,0.0,1.0);
        for (const auto& [key,weight]:grouped) {
            OrbitalAOContribution contribution;
            contribution.atom_index=std::get<0>(key);
            contribution.principal_n=std::get<1>(key);
            contribution.angular_momentum=std::get<2>(key);
            contribution.label=reference_label(
                wf,contribution.atom_index,
                contribution.principal_n,
                contribution.angular_momentum);
            contribution.weight=weight;
            chemistry.ao_contributions.push_back(std::move(contribution));
        }
        std::sort(chemistry.ao_contributions.begin(),
                  chemistry.ao_contributions.end(),
                  [](const auto& a,const auto& b){return a.weight>b.weight;});
        chemistry.confidence=std::clamp(
            chemistry.valence_weight+
            0.5*(chemistry.deep_core_weight+chemistry.semicore_weight),
            0.0,1.0);
    }

    select_valence_manifold(
        wf,valence_rank,options.degeneracy_tolerance_hartree);

    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        auto& mo=wf.orbitals[i];
        auto& chemistry=mo.chemistry;
        for (std::uint32_t a=0;a<wf.atoms.size();++a) {
            for (std::uint32_t b=a+1u;b<wf.atoms.size();++b) {
                const double mayer=pair_mayer(wf,a,b);
                const bool geometric=
                    distance_bohr(wf.atoms[a],wf.atoms[b])<=
                    options.pair_distance_ceiling_bohr;
                const bool atoms_present=
                    a<atom_weight[i].size() && b<atom_weight[i].size() &&
                    atom_weight[i][a]>=options.pair_atom_weight_floor &&
                    atom_weight[i][b]>=options.pair_atom_weight_floor;
                if (std::abs(mayer)<options.pair_mayer_floor &&
                    !(geometric && atoms_present)) continue;

                double overlap_character=0.0;
                for (std::size_t mu=0;mu<n;++mu) {
                    if (basis_atom[mu]!=a) continue;
                    for (std::size_t nu=0;nu<n;++nu) {
                        if (basis_atom[nu]!=b) continue;
                        overlap_character+=
                            2.0*static_cast<double>(mo.coefficients[mu])*
                            wf.ao_overlap[mu*n+nu]*
                            static_cast<double>(mo.coefficients[nu]);
                    }
                }

                OrbitalPairInteraction interaction;
                interaction.atom_a=a;
                interaction.atom_b=b;
                interaction.atom_a_label=
                    wf.atoms[a].symbol+std::to_string(a+1u);
                interaction.atom_b_label=
                    wf.atoms[b].symbol+std::to_string(b+1u);
                interaction.total_mayer_index=mayer;
                interaction.overlap_character=overlap_character;
                interaction.occupied_overlap_contribution=
                    static_cast<double>(mo.occupation)*overlap_character;

                RawChannels channels;
                const auto axis=unit_axis(wf.atoms[a],wf.atoms[b]);
                add_atom_channels(a,axis,reference,amplitudes[i],channels);
                add_atom_channels(b,axis,reference,amplitudes[i],channels);
                interaction.channel=finish_channels(channels,options);
                interaction.bonding=pair_bonding(
                    overlap_character,atoms_present,options);
                chemistry.interactions.push_back(std::move(interaction));
            }
        }

        chemistry.channel=aggregate_channels(
            chemistry.interactions,options);
        chemistry.bonding=aggregate_bonding(
            chemistry.interactions);

        if (chemistry.interactions.empty()) {
            if (chemistry.deep_core_weight>=0.85) {
                chemistry.channel.status=ChemistryStatus::NotApplicable;
                chemistry.channel.dominant=OrbitalAngularFamily::NotApplicable;
                chemistry.channel.undetermined=0.0;
                chemistry.bonding.status=ChemistryStatus::NotApplicable;
                chemistry.bonding.dominant=OrbitalBondingRole::NotApplicable;
                chemistry.bonding.undetermined=0.0;
            } else {
                chemistry.note=
                    "No stable atom-pair interaction frame; chemistry remains UND";
            }
        }

        if (!chemistry.valence_manifold &&
            chemistry.valence_weight<0.10 &&
            chemistry.deep_core_weight<0.50) {
            chemistry.note=
                "Outside the selected minimal chemical-valence canonical manifold";
        }
    }

    attach_multicentre_assignments(wf);
    generic_three_centre_fallback(wf);
}

} // namespace cov
