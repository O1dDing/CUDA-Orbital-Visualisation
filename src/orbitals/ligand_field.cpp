#include "cov/ligand_field.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace cov {
namespace {

bool transition_metal(const int z) noexcept {
    return (z>=21 && z<=30) || (z>=39 && z<=48) ||
           (z>=72 && z<=80) || (z>=104 && z<=112);
}

struct Neighbour {
    std::size_t atom = 0;
    double mayer = 0.0;
    double distance = 0.0;
    std::array<double,3> direction{};
};

double dot(const std::array<double,3>& a,
           const std::array<double,3>& b) noexcept {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

Neighbour make_neighbour(const Wavefunction& wavefunction,
                         const std::size_t metal,
                         const std::size_t atom,
                         const double mayer) {
    const auto& m=wavefunction.atoms[metal];
    const auto& a=wavefunction.atoms[atom];
    const std::array<double,3> displacement{a.x-m.x,a.y-m.y,a.z-m.z};
    const double distance=std::sqrt(dot(displacement,displacement));
    Neighbour result;
    result.atom=atom;
    result.mayer=std::abs(mayer);
    result.distance=distance;
    if (distance>1.0e-12) {
        result.direction={displacement[0]/distance,
                          displacement[1]/distance,
                          displacement[2]/distance};
    }
    return result;
}

double tetrahedral_confidence(const std::vector<Neighbour>& shell);
double octahedral_confidence(const std::vector<Neighbour>& shell);

std::vector<Neighbour> first_shell(const Wavefunction& wavefunction,
                                   const std::size_t metal) {
    std::vector<Neighbour> candidates;
    for (const auto& bond:wavefunction.bond_orders) {
        std::size_t other=wavefunction.atoms.size();
        if (bond.atom_a==metal) other=bond.atom_b;
        else if (bond.atom_b==metal) other=bond.atom_a;
        if (other>=wavefunction.atoms.size() ||
            transition_metal(wavefunction.atoms[other].atomic_number)) continue;
        const double magnitude=std::abs(bond.mayer_order);
        if (magnitude<0.02) continue;
        candidates.push_back(make_neighbour(
            wavefunction,metal,other,bond.mayer_order));
    }
    if (candidates.empty()) return {};
    candidates.erase(std::remove_if(
        candidates.begin(),candidates.end(),[](const auto& item) {
            return item.distance<=1.0e-12;
        }),candidates.end());
    std::sort(candidates.begin(),candidates.end(),[](const auto& a,const auto& b) {
        if (a.distance!=b.distance) return a.distance<b.distance;
        return a.mayer>b.mayer;
    });

    // A carbonyl oxygen, for example, can carry a small density-level Mayer
    // coupling to the metal but lies behind the directly coordinated carbon.
    // Retain only the nearest atom along each ligand ray.
    std::vector<Neighbour> radial;
    for (const auto& candidate:candidates) {
        const bool shadowed=std::any_of(
            radial.begin(),radial.end(),[&](const auto& retained) {
                return dot(candidate.direction,retained.direction)>0.94;
            });
        if (!shadowed) radial.push_back(candidate);
    }
    if (radial.size()>6u) {
        std::stable_sort(radial.begin(),radial.end(),[](const auto& a,const auto& b) {
            const double sa=a.mayer/std::max(1.0e-9,a.distance);
            const double sb=b.mayer/std::max(1.0e-9,b.distance);
            return sa>sb;
        });
        radial.resize(6u);
    }

    // Set the electronic cutoff only after collinear through-ligand atoms
    // have been removed.  Otherwise a remote M...O/N density coupling can
    // raise the floor and delete the nearer direct donor before radial
    // shadowing has a chance to act.
    const double strongest=std::accumulate(
        radial.begin(),radial.end(),0.0,[](const double value,const auto& item) {
            return std::max(value,item.mayer);
        });
    const double floor=std::max(0.035,0.08*strongest);
    auto electronically_strong=radial;
    electronically_strong.erase(std::remove_if(
        electronically_strong.begin(),electronically_strong.end(),
        [&](const auto& item) {return item.mayer<floor;}),
        electronically_strong.end());
    const auto supported=[](const auto& shell) {
        return (shell.size()==4u && tetrahedral_confidence(shell)>=0.55) ||
               (shell.size()==6u && octahedral_confidence(shell)>=0.55);
    };
    if (supported(electronically_strong)) return electronically_strong;
    // Ionic first-shell bonds can all lie below 0.035.  The low-Mayer retry
    // is accepted only when the complete radial-nearest set independently
    // forms a supported Td/Oh geometry; this does not lower the global floor.
    if (supported(radial)) return radial;
    return electronically_strong;
}

double tetrahedral_confidence(const std::vector<Neighbour>& shell) {
    if (shell.size()!=4u) return 0.0;
    double squared=0.0;
    std::size_t count=0u;
    for (std::size_t i=0;i<shell.size();++i) {
        for (std::size_t j=i+1u;j<shell.size();++j) {
            const double error=dot(shell[i].direction,shell[j].direction)+1.0/3.0;
            squared+=error*error;
            ++count;
        }
    }
    const double rms=std::sqrt(squared/static_cast<double>(count));
    return std::clamp(1.0-rms/0.30,0.0,1.0);
}

double octahedral_confidence(const std::vector<Neighbour>& shell) {
    if (shell.size()!=6u) return 0.0;
    struct Pair { std::size_t a=0; std::size_t b=0; double value=0.0; };
    std::vector<Pair> pairs;
    for (std::size_t i=0;i<shell.size();++i) {
        for (std::size_t j=i+1u;j<shell.size();++j) {
            pairs.push_back({i,j,dot(shell[i].direction,shell[j].direction)});
        }
    }
    std::sort(pairs.begin(),pairs.end(),[](const auto& a,const auto& b) {
        return a.value<b.value;
    });
    std::set<std::size_t> used;
    std::vector<Pair> opposite;
    for (const auto& pair:pairs) {
        if (used.count(pair.a) || used.count(pair.b)) continue;
        used.insert(pair.a);
        used.insert(pair.b);
        opposite.push_back(pair);
        if (opposite.size()==3u) break;
    }
    if (opposite.size()!=3u ||
        std::any_of(opposite.begin(),opposite.end(),[](const auto& pair) {
            return pair.value>-0.72;
        })) return 0.0;

    double opposite_squared=0.0;
    for (const auto& pair:opposite) {
        const double error=pair.value+1.0;
        opposite_squared+=error*error;
    }
    double right_squared=0.0;
    std::size_t right_count=0u;
    for (const auto& pair:pairs) {
        const bool is_opposite=std::any_of(
            opposite.begin(),opposite.end(),[&](const auto& selected) {
                return selected.a==pair.a && selected.b==pair.b;
            });
        if (is_opposite) continue;
        right_squared+=pair.value*pair.value;
        ++right_count;
    }
    const double opposite_rms=std::sqrt(opposite_squared/3.0);
    const double right_rms=std::sqrt(
        right_squared/static_cast<double>(std::max<std::size_t>(1u,right_count)));
    return std::clamp(1.0-0.45*opposite_rms/0.25-
                          0.55*right_rms/0.35,0.0,1.0);
}

std::string normalised_symmetry(std::string value) {
    value.erase(std::remove_if(value.begin(),value.end(),[](unsigned char c) {
        return c==' ' || c=='\t' || c=='\r' || c=='\n';
    }),value.end());
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char c) {
        if (c>='A' && c<='Z') return static_cast<char>(c-'A'+'a');
        return static_cast<char>(c);
    });
    if (value=="?" || value=="n/a" || value=="na" ||
        value=="none" || value=="-") return {};
    return value;
}

double central_metal_weight(const MolecularOrbital& orbital,
                            const std::size_t metal_atom,
                            const int angular_momentum) {
    double result=0.0;
    for (const auto& contribution:orbital.chemistry.ao_contributions) {
        if (contribution.atom_index==metal_atom &&
            contribution.angular_momentum==angular_momentum) {
            result+=contribution.weight;
        }
    }
    return result;
}

struct LocalOrbitalGroup {
    std::vector<std::size_t> members;
    Spin spin=Spin::Alpha;
    double energy=0.0;
    std::array<double,3> metal{};
    std::string symmetry;
};

int dominant_family(const LocalOrbitalGroup& group) {
    const auto found=std::max_element(group.metal.begin(),group.metal.end());
    if (found==group.metal.end()) return -1;
    const int family=static_cast<int>(std::distance(group.metal.begin(),found));
    const double floor=family==2?kLocalDIrrepWeightFloor:0.08;
    return *found>=floor?family:-1;
}

LocalOrbitalGroup merge_groups(const LocalOrbitalGroup& left,
                               const LocalOrbitalGroup& right) {
    LocalOrbitalGroup result=left;
    const double left_count=static_cast<double>(left.members.size());
    const double right_count=static_cast<double>(right.members.size());
    const double total=std::max(1.0,left_count+right_count);
    result.energy=(left.energy*left_count+right.energy*right_count)/total;
    for (std::size_t family=0;family<result.metal.size();++family) {
        result.metal[family]=(
            left.metal[family]*left_count+right.metal[family]*right_count)/total;
    }
    result.members.insert(
        result.members.end(),right.members.begin(),right.members.end());
    if (result.symmetry.empty()) result.symmetry=right.symmetry;
    return result;
}

bool merge_resolved_five_d_run(std::vector<LocalOrbitalGroup>& groups) {
    constexpr double d_tolerance_hartree=0.005;
    const auto try_run=[&](const std::vector<std::size_t>& run) {
        if (run.size()!=5u) return false;
        std::array<double,4> gaps{};
        for (std::size_t i=0;i<gaps.size();++i) {
            gaps[i]=groups[run[i+1u]].energy-groups[run[i]].energy;
        }
        const auto largest=std::max_element(gaps.begin(),gaps.end());
        const std::size_t split=static_cast<std::size_t>(
            std::distance(gaps.begin(),largest))+1u;
        if (split!=2u && split!=3u) return false;
        double second=0.0;
        for (std::size_t i=0;i<gaps.size();++i) {
            if (i+1u==split) continue;
            second=std::max(second,gaps[i]);
        }
        if (*largest<1.0e-6 || *largest<1.25*second) return false;

        LocalOrbitalGroup lower=groups[run.front()];
        for (std::size_t i=1u;i<split;++i) {
            lower=merge_groups(lower,groups[run[i]]);
        }
        LocalOrbitalGroup upper=groups[run[split]];
        for (std::size_t i=split+1u;i<run.size();++i) {
            upper=merge_groups(upper,groups[run[i]]);
        }
        const std::set<std::size_t> consumed(run.begin(),run.end());
        std::vector<LocalOrbitalGroup> rebuilt;
        rebuilt.reserve(groups.size()-3u);
        for (std::size_t i=0;i<groups.size();++i) {
            if (i==run.front()) rebuilt.push_back(lower);
            else if (i==run[split]) rebuilt.push_back(upper);
            else if (consumed.count(i)==0u) rebuilt.push_back(groups[i]);
        }
        groups=std::move(rebuilt);
        std::stable_sort(groups.begin(),groups.end(),[](const auto& a,const auto& b) {
            if (a.spin!=b.spin) return a.spin==Spin::Alpha;
            return a.energy<b.energy;
        });
        return true;
    };

    for (const auto spin:{Spin::Alpha,Spin::Beta}) {
        std::vector<std::size_t> run;
        auto flush=[&]() {
            const bool merged=try_run(run);
            run.clear();
            return merged;
        };
        for (std::size_t i=0;i<groups.size();++i) {
            if (groups[i].spin!=spin) continue;
            const int family=dominant_family(groups[i]);
            const bool candidate=family==2 && groups[i].members.size()==1u &&
                                 normalised_symmetry(groups[i].symmetry).empty();
            if (!candidate) {
                if (family==2 && flush()) return true;
                continue;
            }
            if (!run.empty() && groups[i].energy-groups[run.back()].energy>
                                    d_tolerance_hartree) {
                if (flush()) return true;
            }
            run.push_back(i);
        }
        if (flush()) return true;
    }
    return false;
}

std::vector<LocalOrbitalGroup> initial_orbital_groups(
    const Wavefunction& wavefunction,
    const std::vector<OrbitalMetadata>& metadata,
    const std::size_t metal_atom) {
    std::vector<LocalOrbitalGroup> result;
    for (std::size_t base=0;base<wavefunction.orbitals.size();) {
        const Spin spin=wavefunction.orbitals[base].spin;
        const std::size_t requested=base<metadata.size()
            ?std::max<std::size_t>(1u,metadata[base].degeneracy_size):1u;
        std::size_t end=std::min(wavefunction.orbitals.size(),base+requested);
        for (std::size_t index=base+1u;index<end;++index) {
            if (wavefunction.orbitals[index].spin!=spin) {
                end=index;
                break;
            }
        }
        LocalOrbitalGroup group;
        group.spin=spin;
        if (base<metadata.size()) {
            const std::string current=normalised_symmetry(metadata[base].symmetry);
            if (!current.empty() && current!="?" && current!="n/a") {
                group.symmetry=metadata[base].symmetry;
            }
        }
        for (std::size_t index=base;index<end;++index) {
            const auto& orbital=wavefunction.orbitals[index];
            group.members.push_back(index);
            group.energy+=orbital.energy_hartree;
            for (int family=0;family<3;++family) {
                group.metal[static_cast<std::size_t>(family)]+=
                    central_metal_weight(orbital,metal_atom,family);
            }
        }
        const double count=static_cast<double>(
            std::max<std::size_t>(1u,group.members.size()));
        group.energy/=count;
        for (auto& weight:group.metal) weight/=count;
        result.push_back(std::move(group));
        base=std::max(base+1u,end);
    }
    return result;
}

void merge_local_pseudodegenerate_groups(
    std::vector<LocalOrbitalGroup>& groups) {
    while (merge_resolved_five_d_run(groups)) {
        // Restart because non-adjacent ligand-centred rows may have been
        // crossed while forming the two resolved d subspaces.
    }
    for (std::size_t index=0;index+1u<groups.size();) {
        const auto& left=groups[index];
        const auto& right=groups[index+1u];
        const std::string left_symmetry=normalised_symmetry(left.symmetry);
        const std::string right_symmetry=normalised_symmetry(right.symmetry);
        const bool labels_compatible=
            (left_symmetry.empty() && right_symmetry.empty()) ||
            (!left_symmetry.empty() && left_symmetry==right_symmetry);
        const int family=dominant_family(left);
        if (left.spin!=right.spin || !labels_compatible || family<1 ||
            family!=dominant_family(right)) {
            ++index;
            continue;
        }
        const std::size_t left_size=left.members.size();
        const std::size_t right_size=right.members.size();
        const std::size_t combined=left_size+right_size;
        // A C1/Cs producer commonly emits a local metal-p T2/T1u block as
        // three single rows.  Permit the intermediate 1+1 -> 2 merge so the
        // following pass can complete 2+1 -> 3; requiring combined==3 here
        // made that recovery impossible.
        const bool expected_p=family==1 && combined<=3u &&
            (left_size<3u || right_size<3u);
        const bool expected_d=family==2 &&
            ((combined==3u && (left_size==1u || right_size==1u)) ||
             (combined==2u && left_size==1u && right_size==1u));
        const double tolerance=family==1?0.015:0.005;
        if ((!expected_p && !expected_d) ||
            std::abs(right.energy-left.energy)>tolerance) {
            ++index;
            continue;
        }
        groups[index]=merge_groups(left,right);
        groups.erase(groups.begin()+static_cast<std::ptrdiff_t>(index+1u));
    }
}

std::string local_group_label(const std::string& point_group,
                              const LocalOrbitalGroup& group) {
    if (!group.symmetry.empty()) return group.symmetry;
    const std::size_t degeneracy=group.members.size();
    if (point_group=="Oh") {
        if (degeneracy==3u &&
            group.metal[2]>=kLocalDIrrepWeightFloor) return "T2g";
        if (degeneracy==2u &&
            group.metal[2]>=kLocalDIrrepWeightFloor) return "Eg";
        if (degeneracy==3u && group.metal[1]>=0.08) return "T1u";
        if (degeneracy==1u && group.metal[0]>=0.08) return "A1g";
    } else if (point_group=="Td") {
        if (degeneracy==2u &&
            group.metal[2]>=kLocalDIrrepWeightFloor) return "E";
        if (degeneracy==3u &&
            (group.metal[2]>=kLocalDIrrepWeightFloor ||
             group.metal[1]>=0.08)) return "T2";
        if (degeneracy==1u && group.metal[0]>=0.08) return "A1";
    }
    return {};
}

} // namespace

std::string LigandFieldEnvironment::local_point_group() const {
    switch (geometry) {
        case LigandFieldGeometry::Tetrahedral: return "Td";
        case LigandFieldGeometry::Octahedral: return "Oh";
        default: return {};
    }
}

LigandFieldEnvironment analyse_ligand_field_environment(
    const Wavefunction& wavefunction) {
    LigandFieldEnvironment best;
    for (std::size_t metal=0;metal<wavefunction.atoms.size();++metal) {
        if (!transition_metal(wavefunction.atoms[metal].atomic_number)) continue;
        const auto shell=first_shell(wavefunction,metal);
        LigandFieldEnvironment candidate;
        candidate.metal_atom=metal;
        const double oh=octahedral_confidence(shell);
        const double td=tetrahedral_confidence(shell);
        if (oh>=0.55 && oh>=td) {
            candidate.geometry=LigandFieldGeometry::Octahedral;
            candidate.confidence=oh;
        } else if (td>=0.55) {
            candidate.geometry=LigandFieldGeometry::Tetrahedral;
            candidate.confidence=td;
        } else {
            continue;
        }
        for (const auto& item:shell) candidate.ligand_atoms.push_back(item.atom);
        if (!candidate.ligand_atoms.empty()) {
            const int first_z=wavefunction.atoms[
                candidate.ligand_atoms.front()].atomic_number;
            candidate.equivalent_ligand_elements=std::all_of(
                candidate.ligand_atoms.begin(),candidate.ligand_atoms.end(),
                [&](const auto atom) {
                    return wavefunction.atoms[atom].atomic_number==first_z;
                });
        }
        if (!best.available() || candidate.confidence>best.confidence) {
            best=std::move(candidate);
        }
    }
    return best;
}

void apply_local_ligand_field_symmetry(
    const Wavefunction& wavefunction,
    std::vector<OrbitalMetadata>& metadata) {
    if (metadata.size()!=wavefunction.orbitals.size()) return;
    const LigandFieldEnvironment environment=
        analyse_ligand_field_environment(wavefunction);
    if (!environment.available() ||
        !environment.equivalent_ligand_elements) return;

    auto groups=initial_orbital_groups(
        wavefunction,metadata,environment.metal_atom);
    merge_local_pseudodegenerate_groups(groups);
    const std::string point_group=environment.local_point_group();
    for (const auto& group:groups) {
        const std::string symmetry=local_group_label(point_group,group);
        if (symmetry.empty()) continue;
        for (const auto member:group.members) {
            if (member>=metadata.size()) continue;
            const std::string current=normalised_symmetry(
                metadata[member].symmetry);
            if (current.empty() || current=="?" || current=="n/a") {
                metadata[member].symmetry=symmetry;
            }
        }
    }
}

} // namespace cov
