#include "cov/mo_diagram.hpp"
#include "cov/ligand_field.hpp"
#include "cov/local_orbital_symmetry.hpp"
#include "cov/point_group_catalog.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace cov {

OrbitalAnnotation annotate_orbital_legacy(const MolecularOrbital& orbital);
DiagramSelectionPlan build_valence_selection_plan_legacy(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::vector<OrbitalMetadata>& metadata);
MODiagramData build_mo_diagram_data_legacy(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options);

namespace {

bool chemistry_available(const Wavefunction& wavefunction) {
    return std::any_of(
        wavefunction.orbitals.begin(),wavefunction.orbitals.end(),
        [](const MolecularOrbital& orbital) {
            return orbital.chemistry.available;
        });
}

bool occupied(const MolecularOrbital& orbital,const double threshold) {
    return static_cast<double>(orbital.occupation)>threshold;
}

std::string family_name(const OrbitalAngularFamily family) {
    switch (family) {
        case OrbitalAngularFamily::Sigma: return "sigma";
        case OrbitalAngularFamily::Pi: return "pi";
        case OrbitalAngularFamily::Delta: return "delta";
        case OrbitalAngularFamily::Phi: return "phi";
        default: return "unavailable";
    }
}

std::string script_number(const std::size_t value,const bool superscript) {
    static constexpr const char* superscript_digits[]={
        "⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹"};
    static constexpr const char* subscript_digits[]={
        "₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"};
    const std::string plain=std::to_string(value);
    std::string result;
    for (const char digit:plain) {
        const auto index=static_cast<std::size_t>(digit-'0');
        result+=superscript?superscript_digits[index]:subscript_digits[index];
    }
    return result;
}

std::string delocalised_pi_label(const std::size_t atoms,
                                 const double electrons) {
    const auto rounded=static_cast<std::size_t>(
        std::max(0LL,std::llround(electrons)));
    return std::string("Π")+script_number(atoms,true)+
           script_number(rounded,false);
}

OrbitalAnnotation chemistry_annotation(const MolecularOrbital& orbital) {
    const auto& chemistry=orbital.chemistry;
    OrbitalAnnotation result;
    if (!chemistry.available) return result;

    if (chemistry.channel.status==ChemistryStatus::Determined ||
        chemistry.channel.status==ChemistryStatus::Percentages) {
        result.family=family_name(chemistry.channel.dominant);
        result.family_source=AnnotationSource::Derived;
        result.family_confidence=chemistry.confidence;
    }

    if (chemistry.bonding.status==ChemistryStatus::Determined ||
        chemistry.bonding.status==ChemistryStatus::Percentages) {
        switch (chemistry.bonding.dominant) {
            case OrbitalBondingRole::Bonding:
                result.bonding_class=BondingClass::Bonding; break;
            case OrbitalBondingRole::Antibonding:
                result.bonding_class=BondingClass::Antibonding; break;
            case OrbitalBondingRole::Nonbonding:
                result.bonding_class=BondingClass::Nonbonding; break;
            default: break;
        }
        result.bonding_source=AnnotationSource::Derived;
        result.bonding_confidence=chemistry.confidence;
    }

    if (!chemistry.multicentre_label.empty()) {
        result.multicentre.available=true;
        result.multicentre.centres=chemistry.participating_atoms;
        result.multicentre.electrons=chemistry.participating_electrons;
        result.multicentre.atom_indices.assign(
            chemistry.participating_atom_indices.begin(),
            chemistry.participating_atom_indices.end());
        result.multicentre.label=chemistry.multicentre_label;
        result.multicentre.source=AnnotationSource::Derived;
        result.multicentre.confidence=chemistry.confidence;
        result.multicentre.heuristic=false;
    }

    if (chemistry.channel.dominant==OrbitalAngularFamily::Pi &&
        chemistry.participating_atoms>2u &&
        (!chemistry.delocalised_family_id.empty() ||
         chemistry.multicentre_label=="delocalised-pi")) {
        result.delocalised_pi.available=true;
        result.delocalised_pi.participating_atoms=
            chemistry.participating_atoms;
        result.delocalised_pi.participating_electrons=
            chemistry.participating_electrons;
        result.delocalised_pi.atom_indices.assign(
            chemistry.participating_atom_indices.begin(),
            chemistry.participating_atom_indices.end());
        result.delocalised_pi.orbital_indices.assign(
            chemistry.delocalised_family_orbitals.begin(),
            chemistry.delocalised_family_orbitals.end());
        result.delocalised_pi.family_id=
            chemistry.delocalised_family_id;
        result.delocalised_pi.label=delocalised_pi_label(
            chemistry.participating_atoms,
            chemistry.participating_electrons);
        result.delocalised_pi.source=AnnotationSource::Derived;
        result.delocalised_pi.confidence=chemistry.confidence;
        result.delocalised_pi.heuristic=false;
    }
    result.heuristic=false;
    return result;
}

std::vector<std::size_t> expand_degenerate(
    std::vector<std::size_t> selected,
    const std::vector<OrbitalLabel>& labels) {
    std::set<std::size_t> expanded(selected.begin(),selected.end());
    for (const auto index:selected) {
        if (index>=labels.size()) continue;
        const auto& label=labels[index];
        if (label.group_size<=1u) continue;
        const std::size_t begin=label.group_base_number>0u
            ?label.group_base_number-1u:index;
        for (std::size_t member=0;member<label.group_size;++member) {
            if (begin+member<labels.size()) expanded.insert(begin+member);
        }
    }
    return {expanded.begin(),expanded.end()};
}

double group_layout_energy(
    const std::vector<OrbitalMetadata>& metadata,
    const std::vector<OrbitalLabel>& labels,
    const std::size_t index) {
    if (index>=metadata.size() || index>=labels.size()) return 0.0;
    const auto& label=labels[index];
    if (label.group_size<=1u) return metadata[index].energy_hartree;
    const std::size_t begin=label.group_base_number>0u
        ?label.group_base_number-1u:index;
    double sum=0.0;
    std::size_t count=0;
    for (std::size_t member=0;member<label.group_size;++member) {
        if (begin+member>=metadata.size()) break;
        sum+=metadata[begin+member].energy_hartree;
        ++count;
    }
    return count>0u?sum/static_cast<double>(count)
                   :metadata[index].energy_hartree;
}

bool transition_metal(const int z) noexcept {
    return (z>=21 && z<=30) || (z>=39 && z<=48) ||
           (z>=72 && z<=80) || (z>=104 && z<=112);
}

enum class LocalLigandPiRole {
    Unresolved,
    SigmaOnly,
    Donor,
    Acceptor,
    Ambiguous,
};

struct LigandScope {
    bool available=false;
    std::size_t metal=0;
    std::set<std::size_t> direct_donors;
    std::set<std::size_t> fragment_atoms;
    LocalLigandPiRole pi_role=LocalLigandPiRole::Unresolved;
};

bool any_metal_ligand_pair(const Wavefunction& wavefunction,
                           const OrbitalPairInteraction& interaction) {
    if (interaction.atom_a>=wavefunction.atoms.size() ||
        interaction.atom_b>=wavefunction.atoms.size()) return false;
    const bool a=transition_metal(
        wavefunction.atoms[interaction.atom_a].atomic_number);
    const bool b=transition_metal(
        wavefunction.atoms[interaction.atom_b].atomic_number);
    return a!=b;
}

LigandScope make_ligand_scope(
    const Wavefunction& wavefunction,
    const LigandFieldEnvironment& environment) {
    LigandScope scope;
    if (!environment.available() ||
        environment.metal_atom>=wavefunction.atoms.size()) return scope;
    scope.available=true;
    scope.metal=environment.metal_atom;
    scope.direct_donors.insert(
        environment.ligand_atoms.begin(),environment.ligand_atoms.end());

    // Grow each ligand fragment without crossing a transition-metal centre.
    // This includes O in CO and N in CN for ligand-character analysis while
    // keeping the direct M-donor interaction restricted to the first shell.
    std::vector<std::size_t> pending(
        scope.direct_donors.begin(),scope.direct_donors.end());
    scope.fragment_atoms=scope.direct_donors;
    while (!pending.empty()) {
        const std::size_t atom=pending.back();
        pending.pop_back();
        for (const auto& bond:wavefunction.bond_orders) {
            if (std::abs(bond.mayer_order)<0.05) continue;
            std::size_t other=wavefunction.atoms.size();
            if (bond.atom_a==atom) other=bond.atom_b;
            else if (bond.atom_b==atom) other=bond.atom_a;
            if (other>=wavefunction.atoms.size() || other==scope.metal ||
                transition_metal(wavefunction.atoms[other].atomic_number)) {
                continue;
            }
            if (scope.fragment_atoms.insert(other).second) {
                pending.push_back(other);
            }
        }
    }

    std::size_t sigma_only=0u;
    std::size_t donor=0u;
    std::size_t acceptor=0u;
    std::size_t ambiguous=0u;
    for (const auto atom:scope.direct_donors) {
        if (atom>=wavefunction.atoms.size()) continue;
        const int z=wavefunction.atoms[atom].atomic_number;
        std::size_t external_neighbour_count=0u;
        bool multiple_heavy_bond=false;
        for (const auto& bond:wavefunction.bond_orders) {
            std::size_t other=wavefunction.atoms.size();
            if (bond.atom_a==atom) other=bond.atom_b;
            else if (bond.atom_b==atom) other=bond.atom_a;
            if (other>=wavefunction.atoms.size() || other==scope.metal ||
                transition_metal(wavefunction.atoms[other].atomic_number)) {
                continue;
            }
            const int other_z=wavefunction.atoms[other].atomic_number;
            if (std::abs(bond.mayer_order)>=0.10) {
                ++external_neighbour_count;
            }
            if (other_z>1 && std::abs(bond.mayer_order)>=0.10) {
                multiple_heavy_bond=multiple_heavy_bond ||
                    std::abs(bond.mayer_order)>=1.15;
            }
        }
        if (z==6 && multiple_heavy_bond) {
            ++acceptor; // CO, CN and related unsaturated carbon donors.
        } else if ((z==7 || z==15) && !multiple_heavy_bond &&
                   external_neighbour_count>=3u) {
            ++sigma_only; // Saturated NH3/amine/phosphine-like donors.
        } else if ((z==7 || z==15) && multiple_heavy_bond) {
            ++acceptor;
        } else if (z==9 || z==17 || z==35 || z==53 ||
                   z==8 || z==16 || z==34) {
            ++donor;
        } else {
            ++ambiguous;
        }
    }
    const std::size_t classified=sigma_only+donor+acceptor+ambiguous;
    if (classified==0u) scope.pi_role=LocalLigandPiRole::Unresolved;
    else if (sigma_only==classified) scope.pi_role=LocalLigandPiRole::SigmaOnly;
    else if (acceptor>0u && donor==0u && ambiguous==0u) {
        scope.pi_role=LocalLigandPiRole::Acceptor;
    } else if (donor>0u && acceptor==0u && ambiguous==0u) {
        scope.pi_role=LocalLigandPiRole::Donor;
    } else {
        scope.pi_role=LocalLigandPiRole::Ambiguous;
    }
    return scope;
}

bool scoped_metal_ligand_pair(
    const Wavefunction& wavefunction,
    const OrbitalPairInteraction& interaction,
    const LigandScope* scope) {
    if (scope==nullptr || !scope->available) {
        return any_metal_ligand_pair(wavefunction,interaction);
    }
    const std::size_t a=interaction.atom_a;
    const std::size_t b=interaction.atom_b;
    return (a==scope->metal && scope->direct_donors.count(b)>0u) ||
           (b==scope->metal && scope->direct_donors.count(a)>0u);
}

double bounded_overlap_character(const double value) noexcept {
    // Individual canonical-MO Mulliken terms are basis dependent and may be
    // arbitrarily large for diffuse virtuals.  tanh preserves sign and the
    // linear weak-coupling regime while preventing one outlier from taking
    // over a group score.
    return std::tanh(value);
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

struct GroupCandidate {
    MODiagramLevel level;
    std::size_t base_index=0;
    bool selected_by_reference=false;
    bool selected_by_raw=false;
    bool include=false;
    bool locally_grouped=false;
    bool suppressed_spin_counterpart=false;
    std::uint8_t local_irrep_copy=0;
    bool locally_classified=false;
};

Spin group_spin(const Wavefunction& wavefunction,
                const GroupCandidate& group) {
    if (group.level.member_indices.empty() ||
        group.level.member_indices.front()>=wavefunction.orbitals.size()) {
        return Spin::Alpha;
    }
    return wavefunction.orbitals[group.level.member_indices.front()].spin;
}

int dominant_metal_family(const MODiagramLevel& level) {
    const std::array<double,3> weights{
        level.metal_s_weight,level.metal_p_weight,level.metal_d_weight};
    const auto found=std::max_element(weights.begin(),weights.end());
    if (found==weights.end()) return -1;
    const int family=static_cast<int>(std::distance(weights.begin(),found));
    const double floor=family==2?kLocalDIrrepWeightFloor:0.08;
    return *found>=floor?family:-1;
}

GroupCandidate merge_local_groups(const Wavefunction& wavefunction,
                                  const GroupCandidate& left,
                                  const GroupCandidate& right) {
    GroupCandidate result=left;
    auto& level=result.level;
    const auto left_count=std::max<std::size_t>(
        1u,left.level.member_indices.size());
    const auto right_count=std::max<std::size_t>(
        1u,right.level.member_indices.size());
    const double total=static_cast<double>(left_count+right_count);
    const auto weighted=[&](const double a,const double b) {
        return (a*static_cast<double>(left_count)+
                b*static_cast<double>(right_count))/total;
    };
    level.layout_energy_hartree=weighted(
        left.level.layout_energy_hartree,right.level.layout_energy_hartree);
    level.metadata.energy_hartree=level.layout_energy_hartree;
    level.total_occupation=left.level.total_occupation+
                           right.level.total_occupation;
    level.metadata.occupation=static_cast<float>(
        level.total_occupation/total);
    level.metal_s_weight=weighted(
        left.level.metal_s_weight,right.level.metal_s_weight);
    level.metal_p_weight=weighted(
        left.level.metal_p_weight,right.level.metal_p_weight);
    level.metal_d_weight=weighted(
        left.level.metal_d_weight,right.level.metal_d_weight);
    level.direct_ligand_p_weight=weighted(
        left.level.direct_ligand_p_weight,
        right.level.direct_ligand_p_weight);
    level.ligand_p_weight=weighted(
        left.level.ligand_p_weight,right.level.ligand_p_weight);
    level.sigma_fraction=weighted(
        left.level.sigma_fraction,right.level.sigma_fraction);
    level.pi_fraction=weighted(
        left.level.pi_fraction,right.level.pi_fraction);
    level.metal_ligand_overlap=weighted(
        left.level.metal_ligand_overlap,right.level.metal_ligand_overlap);
    level.homo=left.level.homo || right.level.homo;
    level.lumo=left.level.lumo || right.level.lumo;
    level.metadata.selected=left.level.metadata.selected ||
                            right.level.metadata.selected;
    level.raw_data_fallback=left.level.raw_data_fallback ||
                            right.level.raw_data_fallback;
    level.member_indices.insert(level.member_indices.end(),
        right.level.member_indices.begin(),right.level.member_indices.end());
    level.member_electrons.insert(level.member_electrons.end(),
        right.level.member_electrons.begin(),right.level.member_electrons.end());
    level.metadata.degeneracy_size=level.member_indices.size();
    if (!level.member_electrons.empty()) level.electrons=level.member_electrons.front();

    double minimum=std::numeric_limits<double>::infinity();
    double maximum=-std::numeric_limits<double>::infinity();
    for (const auto index:level.member_indices) {
        if (index>=wavefunction.orbitals.size()) continue;
        minimum=std::min(minimum,wavefunction.orbitals[index].energy_hartree);
        maximum=std::max(maximum,wavefunction.orbitals[index].energy_hartree);
    }
    level.energy_spread_hartree=
        std::isfinite(minimum) && std::isfinite(maximum)
            ?std::max(0.0,maximum-minimum):0.0;
    if (level.annotation.family=="unavailable" &&
        right.level.annotation.family!="unavailable") {
        level.annotation=right.level.annotation;
    }
    result.selected_by_reference=left.selected_by_reference ||
                                 right.selected_by_reference;
    result.selected_by_raw=left.selected_by_raw || right.selected_by_raw;
    result.include=left.include || right.include;
    result.locally_grouped=true;
    if (left.local_irrep_copy!=right.local_irrep_copy) {
        result.local_irrep_copy=0;
    }
    result.locally_classified=left.locally_classified &&
                              right.locally_classified;
    return result;
}

bool merge_resolved_five_d_run(const Wavefunction& wavefunction,
                               std::vector<GroupCandidate>& groups) {
    constexpr double d_tolerance_hartree=0.005;
    const auto try_run=[&](const std::vector<std::size_t>& run) {
        if (run.size()!=5u) return false;
        std::array<double,4> gaps{};
        for (std::size_t i=0;i<gaps.size();++i) {
            gaps[i]=groups[run[i+1u]].level.layout_energy_hartree-
                    groups[run[i]].level.layout_energy_hartree;
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

        GroupCandidate lower=groups[run.front()];
        for (std::size_t i=1u;i<split;++i) {
            lower=merge_local_groups(wavefunction,lower,groups[run[i]]);
        }
        GroupCandidate upper=groups[run[split]];
        for (std::size_t i=split+1u;i<run.size();++i) {
            upper=merge_local_groups(wavefunction,upper,groups[run[i]]);
        }
        const std::set<std::size_t> consumed(run.begin(),run.end());
        std::vector<GroupCandidate> rebuilt;
        rebuilt.reserve(groups.size()-3u);
        for (std::size_t i=0;i<groups.size();++i) {
            if (i==run.front()) rebuilt.push_back(lower);
            else if (i==run[split]) rebuilt.push_back(upper);
            else if (consumed.count(i)==0u) rebuilt.push_back(groups[i]);
        }
        groups=std::move(rebuilt);
        std::stable_sort(groups.begin(),groups.end(),[&](const auto& a,const auto& b) {
            const Spin a_spin=group_spin(wavefunction,a);
            const Spin b_spin=group_spin(wavefunction,b);
            if (a_spin!=b_spin) return a_spin==Spin::Alpha;
            return a.level.layout_energy_hartree<b.level.layout_energy_hartree;
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
            if (group_spin(wavefunction,groups[i])!=spin) continue;
            const int family=dominant_metal_family(groups[i].level);
            const bool candidate=family==2 &&
                groups[i].level.member_indices.size()==1u &&
                normalised_symmetry(groups[i].level.metadata.symmetry).empty();
            if (!candidate) {
                if (family==2 && flush()) return true;
                continue;
            }
            if (!run.empty() &&
                groups[i].level.layout_energy_hartree-
                    groups[run.back()].level.layout_energy_hartree>
                        d_tolerance_hartree) {
                if (flush()) return true;
            }
            run.push_back(i);
        }
        if (flush()) return true;
    }
    return false;
}

void merge_local_pseudodegenerate_groups(
    const Wavefunction& wavefunction,
    const LigandFieldEnvironment& environment,
    std::vector<GroupCandidate>& groups) {
    if (!environment.available() ||
        !environment.equivalent_ligand_elements) return;
    const std::string point_group=environment.local_point_group();
    const bool unique_two_three_d=
        classify_local_irrep_by_dimension(
            point_group,MetalAOShell::D,2u).has_value() &&
        classify_local_irrep_by_dimension(
            point_group,MetalAOShell::D,3u).has_value();
    if (unique_two_three_d) {
        while (merge_resolved_five_d_run(wavefunction,groups)) {
            // Re-evaluate indices after forming a non-adjacent
            // five-dimensional d run into two resolved subspaces.
        }
    }
    for (std::size_t i=0;i+1u<groups.size();) {
        auto& left=groups[i];
        auto& right=groups[i+1u];
        const std::string left_symmetry=normalised_symmetry(
            left.level.metadata.symmetry);
        const std::string right_symmetry=normalised_symmetry(
            right.level.metadata.symmetry);
        const bool local_copies_compatible=
            (!left.locally_classified && !right.locally_classified) ||
            (left.locally_classified && right.locally_classified &&
             left.local_irrep_copy!=0u &&
             left.local_irrep_copy==right.local_irrep_copy);
        const bool labels_compatible=
            (left_symmetry.empty() && right_symmetry.empty()) ||
            (!left_symmetry.empty() && left_symmetry==right_symmetry &&
             local_copies_compatible);
        if (group_spin(wavefunction,left)!=group_spin(wavefunction,right) ||
            !labels_compatible) {
            ++i;
            continue;
        }
        const int family=dominant_metal_family(left.level);
        if (family<1 || family!=dominant_metal_family(right.level)) {
            ++i;
            continue;
        }
        const std::size_t left_size=left.level.member_indices.size();
        const std::size_t right_size=right.level.member_indices.size();
        const std::size_t combined=left_size+right_size;
        const bool unlabeled=left_symmetry.empty() && right_symmetry.empty();
        const bool dimension_unambiguous=!unlabeled ||
            classify_local_irrep_by_dimension(
                point_group,static_cast<MetalAOShell>(family),combined)
                .has_value();
        const bool expected_p=family==1 && combined<=3u &&
            (left_size<3u || right_size<3u);
        const bool expected_d=family==2 &&
            ((combined==3u && (left_size==1u || right_size==1u)) ||
             (combined==2u && left_size==1u && right_size==1u));
        const double gap=std::abs(
            right.level.layout_energy_hartree-left.level.layout_energy_hartree);
        const double local_tolerance_hartree=family==1?0.015:0.005;
        if (!dimension_unambiguous || (!expected_p && !expected_d) ||
            gap>local_tolerance_hartree) {
            ++i;
            continue;
        }
        groups[i]=merge_local_groups(wavefunction,left,right);
        groups.erase(groups.begin()+static_cast<std::ptrdiff_t>(i+1u));
    }
}

void recover_local_ligand_field_symmetry(
    const Wavefunction& wavefunction,
    const LigandFieldEnvironment& environment,
    std::vector<GroupCandidate>& groups) {
    if (!environment.available()) return;
    const std::string point_group=environment.local_point_group();
    for (auto& group:groups) {
        auto& level=group.level;
        const int family=dominant_metal_family(level);
        // Do not promote numerical AO noise on a ligand-only SALC into a
        // central-metal irrep.  Besides overwriting a valid producer label,
        // that can create dimensionally impossible rows (for example a
        // three-member Eg group) and then prevent alpha/beta spatial pairing.
        // The same chemically meaningful metal-family floor used by the
        // compact selector is therefore also the admission gate here.
        std::optional<LocalIrrepAssignment> assignment;
        if (family>=0) {
            assignment=classify_local_metal_irrep(
                wavefunction,level.member_indices,environment.metal_atom,
                point_group,environment.rotation_reference_to_input);
        }
        if (!assignment) {
            if (family>=0) {
                assignment=classify_local_irrep_by_dimension(
                    point_group,static_cast<MetalAOShell>(family),
                    level.member_indices.size());
            }
        }
        if (!assignment) continue;
        level.metadata.symmetry=assignment->label;
        group.local_irrep_copy=assignment->copy_index;
        group.locally_classified=true;
    }
}

class SpinOverlapWorkspace {
public:
    explicit SpinOverlapWorkspace(const Wavefunction& wavefunction)
        : wavefunction_(wavefunction),
          basis_count_(wavefunction.basis_count),
          alpha_position_(wavefunction.orbitals.size(),invalid_index()),
          beta_position_(wavefunction.orbitals.size(),invalid_index()) {
        if (basis_count_==0u ||
            wavefunction_.ao_overlap.size()!=basis_count_*basis_count_) {
            return;
        }
        for (std::size_t orbital=0;orbital<wavefunction_.orbitals.size();++orbital) {
            const auto& item=wavefunction_.orbitals[orbital];
            if (item.coefficients.size()!=basis_count_) continue;
            if (item.spin==Spin::Beta) {
                beta_position_[orbital]=beta_orbitals_.size();
                beta_orbitals_.push_back(orbital);
            } else {
                alpha_position_[orbital]=alpha_orbitals_.size();
                alpha_orbitals_.push_back(orbital);
            }
        }
        if (alpha_orbitals_.empty() || beta_orbitals_.empty()) return;

        // The old pairwise implementation recomputed S*C_beta for every
        // alpha/beta group candidate.  Open-shell FCHK files contain a full
        // alpha block followed by a full beta block, so that repeated matrix
        // vector product dominated every UI frame.  Transform each beta MO
        // once, then reuse ordinary dot products for every subspace score and
        // member match.  The mathematical S-metric overlap is unchanged.
        transformed_beta_.assign(
            beta_orbitals_.size()*basis_count_,0.0);
        for (std::size_t beta=0;beta<beta_orbitals_.size();++beta) {
            const auto& coefficients=
                wavefunction_.orbitals[beta_orbitals_[beta]].coefficients;
            double* transformed=transformed_beta_.data()+beta*basis_count_;
            for (std::size_t mu=0;mu<basis_count_;++mu) {
                double value=0.0;
                for (std::size_t nu=0;nu<basis_count_;++nu) {
                    value+=wavefunction_.ao_overlap[mu*basis_count_+nu]*
                           static_cast<double>(coefficients[nu]);
                }
                transformed[mu]=value;
            }
        }

        overlaps_.assign(alpha_orbitals_.size()*beta_orbitals_.size(),0.0);
        for (std::size_t alpha=0;alpha<alpha_orbitals_.size();++alpha) {
            const auto& coefficients=
                wavefunction_.orbitals[alpha_orbitals_[alpha]].coefficients;
            for (std::size_t beta=0;beta<beta_orbitals_.size();++beta) {
                const double* transformed=
                    transformed_beta_.data()+beta*basis_count_;
                double value=0.0;
                for (std::size_t mu=0;mu<basis_count_;++mu) {
                    value+=static_cast<double>(coefficients[mu])*
                           transformed[mu];
                }
                overlaps_[alpha*beta_orbitals_.size()+beta]=value;
            }
        }
    }

    [[nodiscard]] double overlap(const std::size_t left,
                                 const std::size_t right) const noexcept {
        if (left>=alpha_position_.size() || right>=beta_position_.size()) {
            return 0.0;
        }
        std::size_t alpha=alpha_position_[left];
        std::size_t beta=beta_position_[right];
        if (alpha==invalid_index() || beta==invalid_index()) {
            if (right>=alpha_position_.size() || left>=beta_position_.size()) {
                return 0.0;
            }
            alpha=alpha_position_[right];
            beta=beta_position_[left];
        }
        if (alpha==invalid_index() || beta==invalid_index() ||
            beta_orbitals_.empty()) return 0.0;
        return overlaps_[alpha*beta_orbitals_.size()+beta];
    }

private:
    [[nodiscard]] static constexpr std::size_t invalid_index() noexcept {
        return std::numeric_limits<std::size_t>::max();
    }

    const Wavefunction& wavefunction_;
    std::size_t basis_count_=0u;
    std::vector<std::size_t> alpha_orbitals_;
    std::vector<std::size_t> beta_orbitals_;
    std::vector<std::size_t> alpha_position_;
    std::vector<std::size_t> beta_position_;
    std::vector<double> transformed_beta_;
    std::vector<double> overlaps_;
};

double subspace_overlap(const SpinOverlapWorkspace& workspace,
                        const std::vector<std::size_t>& alpha_members,
                        const std::vector<std::size_t>& beta_members) {
    if (alpha_members.empty() ||
        alpha_members.size()!=beta_members.size()) {
        return 0.0;
    }
    double squared=0.0;
    for (const auto a:alpha_members) {
        for (const auto b:beta_members) {
            const double overlap=workspace.overlap(a,b);
            squared+=overlap*overlap;
        }
    }
    return squared/static_cast<double>(alpha_members.size());
}

double subspace_overlap(const SpinOverlapWorkspace& workspace,
                        const GroupCandidate& alpha,
                        const GroupCandidate& beta) {
    return subspace_overlap(workspace,alpha.level.member_indices,
                            beta.level.member_indices);
}

void combine_spin_occupations(const Wavefunction& wavefunction,
                              const SpinOverlapWorkspace& workspace,
                              GroupCandidate& alpha,
                              const GroupCandidate& beta,
                              const double occupation_threshold) {
    const std::string alpha_symmetry=normalised_symmetry(
        alpha.level.metadata.symmetry);
    const std::string beta_symmetry=normalised_symmetry(
        beta.level.metadata.symmetry);
    if ((alpha_symmetry.empty() || alpha_symmetry=="?" ||
         alpha_symmetry=="n/a") &&
        !beta_symmetry.empty() && beta_symmetry!="?" &&
        beta_symmetry!="n/a") {
        alpha.level.metadata.symmetry=beta.level.metadata.symmetry;
    }
    struct Match { std::size_t a=0; std::size_t b=0; double score=0.0; };
    std::vector<Match> matches;
    for (std::size_t a=0;a<alpha.level.member_indices.size();++a) {
        for (std::size_t b=0;b<beta.level.member_indices.size();++b) {
            matches.push_back({a,b,std::abs(workspace.overlap(
                alpha.level.member_indices[a],
                beta.level.member_indices[b]))});
        }
    }
    std::sort(matches.begin(),matches.end(),[](const auto& a,const auto& b) {
        return a.score>b.score;
    });
    std::vector<std::size_t> beta_for_alpha(
        alpha.level.member_indices.size(),std::numeric_limits<std::size_t>::max());
    std::set<std::size_t> used_beta;
    for (const auto& match:matches) {
        if (beta_for_alpha[match.a]!=std::numeric_limits<std::size_t>::max() ||
            used_beta.count(match.b)) continue;
        beta_for_alpha[match.a]=match.b;
        used_beta.insert(match.b);
    }
    alpha.level.member_electrons.assign(
        alpha.level.member_indices.size(),ElectronGlyphs{});
    alpha.level.member_spin_counterparts.assign(
        alpha.level.member_indices.size(),
        std::numeric_limits<std::size_t>::max());
    for (std::size_t a=0;a<alpha.level.member_indices.size();++a) {
        const auto alpha_index=alpha.level.member_indices[a];
        if (alpha_index<wavefunction.orbitals.size() &&
            wavefunction.orbitals[alpha_index].occupation>
                occupation_threshold) {
            alpha.level.member_electrons[a].alpha=1;
        }
        const auto matched=beta_for_alpha[a];
        if (matched<beta.level.member_indices.size()) {
            const auto beta_index=beta.level.member_indices[matched];
            alpha.level.member_spin_counterparts[a]=beta_index;
            if (beta_index<wavefunction.orbitals.size() &&
                wavefunction.orbitals[beta_index].occupation>
                    occupation_threshold) {
                alpha.level.member_electrons[a].beta=1;
            }
        }
    }
    if (!alpha.level.member_electrons.empty()) {
        alpha.level.electrons=alpha.level.member_electrons.front();
    }
    alpha.level.total_occupation+=beta.level.total_occupation;
    alpha.level.metadata.occupation=static_cast<float>(
        alpha.level.total_occupation/
        static_cast<double>(std::max<std::size_t>(
            1u,alpha.level.member_indices.size())));
    alpha.level.metadata.selected=alpha.level.metadata.selected ||
                                  beta.level.metadata.selected;
    alpha.selected_by_reference=alpha.selected_by_reference ||
                                beta.selected_by_reference;
    alpha.selected_by_raw=alpha.selected_by_raw || beta.selected_by_raw;
    alpha.include=alpha.include || beta.include;
}

struct SpinCollapseResult {
    std::size_t paired_groups=0u;
    std::size_t paired_members=0u;
};

struct SpinPairCandidate {
    std::size_t alpha=0u;
    std::size_t beta=0u;
    double score=0.0;
};

std::vector<SpinPairCandidate> maximum_cardinality_spin_matching(
    const std::vector<SpinPairCandidate>& candidates,
    const std::size_t group_count) {
    if (candidates.empty()) return {};

    std::vector<std::size_t> alpha_groups;
    std::vector<std::size_t> beta_groups;
    for (const auto& candidate:candidates) {
        alpha_groups.push_back(candidate.alpha);
        beta_groups.push_back(candidate.beta);
    }
    std::sort(alpha_groups.begin(),alpha_groups.end());
    alpha_groups.erase(
        std::unique(alpha_groups.begin(),alpha_groups.end()),alpha_groups.end());
    std::sort(beta_groups.begin(),beta_groups.end());
    beta_groups.erase(
        std::unique(beta_groups.begin(),beta_groups.end()),beta_groups.end());

    const int source=0;
    const int alpha_begin=1;
    const int beta_begin=alpha_begin+static_cast<int>(alpha_groups.size());
    const int sink=beta_begin+static_cast<int>(beta_groups.size());
    const int node_count=sink+1;
    std::vector<int> alpha_node(group_count,-1);
    std::vector<int> beta_node(group_count,-1);
    for (std::size_t i=0;i<alpha_groups.size();++i) {
        if (alpha_groups[i]<group_count) {
            alpha_node[alpha_groups[i]]=alpha_begin+static_cast<int>(i);
        }
    }
    for (std::size_t i=0;i<beta_groups.size();++i) {
        if (beta_groups[i]<group_count) {
            beta_node[beta_groups[i]]=beta_begin+static_cast<int>(i);
        }
    }

    struct FlowEdge {
        int to=0;
        std::size_t reverse=0u;
        int capacity=0;
        double cost=0.0;
    };
    std::vector<std::vector<FlowEdge>> graph(
        static_cast<std::size_t>(node_count));
    const auto add_edge=[&](const int from,const int to,const double cost) {
        const std::size_t forward=graph[static_cast<std::size_t>(from)].size();
        const std::size_t reverse=graph[static_cast<std::size_t>(to)].size();
        graph[static_cast<std::size_t>(from)].push_back(
            {to,reverse,1,cost});
        graph[static_cast<std::size_t>(to)].push_back(
            {from,forward,0,-cost});
        return forward;
    };
    for (const auto group:alpha_groups) {
        add_edge(source,alpha_node[group],0.0);
    }
    for (const auto group:beta_groups) {
        add_edge(beta_node[group],sink,0.0);
    }
    struct PairArc {
        std::size_t candidate=0u;
        int from=0;
        std::size_t edge=0u;
    };
    std::vector<PairArc> pair_arcs;
    pair_arcs.reserve(candidates.size());
    for (std::size_t i=0;i<candidates.size();++i) {
        const int from=alpha_node[candidates[i].alpha];
        const int to=beta_node[candidates[i].beta];
        pair_arcs.push_back({i,from,add_edge(from,to,-candidates[i].score)});
    }

    // Successive shortest augmenting paths continue until no source/sink
    // path remains.  Consequently cardinality is the primary objective; the
    // negative overlap cost makes total S-metric overlap the secondary one.
    // Residual reverse arcs permit an earlier choice to be replaced, which is
    // precisely the reassignment the former sorted-edge greedy pass lacked.
    constexpr double infinity=std::numeric_limits<double>::infinity();
    for (;;) {
        std::vector<double> distance(
            static_cast<std::size_t>(node_count),infinity);
        std::vector<int> previous_node(
            static_cast<std::size_t>(node_count),-1);
        std::vector<std::size_t> previous_edge(
            static_cast<std::size_t>(node_count),0u);
        std::vector<bool> queued(static_cast<std::size_t>(node_count),false);
        std::queue<int> pending;
        distance[static_cast<std::size_t>(source)]=0.0;
        pending.push(source);
        queued[static_cast<std::size_t>(source)]=true;
        while (!pending.empty()) {
            const int from=pending.front();
            pending.pop();
            queued[static_cast<std::size_t>(from)]=false;
            const auto& edges=graph[static_cast<std::size_t>(from)];
            for (std::size_t edge_index=0;edge_index<edges.size();++edge_index) {
                const auto& edge=edges[edge_index];
                if (edge.capacity<=0) continue;
                const double proposed=distance[static_cast<std::size_t>(from)]+
                                      edge.cost;
                if (proposed>=distance[static_cast<std::size_t>(edge.to)]-
                                 1.0e-12) continue;
                distance[static_cast<std::size_t>(edge.to)]=proposed;
                previous_node[static_cast<std::size_t>(edge.to)]=from;
                previous_edge[static_cast<std::size_t>(edge.to)]=edge_index;
                if (!queued[static_cast<std::size_t>(edge.to)]) {
                    pending.push(edge.to);
                    queued[static_cast<std::size_t>(edge.to)]=true;
                }
            }
        }
        if (previous_node[static_cast<std::size_t>(sink)]<0) break;
        for (int node=sink;node!=source;
             node=previous_node[static_cast<std::size_t>(node)]) {
            const int from=previous_node[static_cast<std::size_t>(node)];
            const std::size_t edge_index=
                previous_edge[static_cast<std::size_t>(node)];
            auto& edge=graph[static_cast<std::size_t>(from)][edge_index];
            --edge.capacity;
            ++graph[static_cast<std::size_t>(node)][edge.reverse].capacity;
        }
    }

    std::vector<SpinPairCandidate> result;
    for (const auto& arc:pair_arcs) {
        if (graph[static_cast<std::size_t>(arc.from)][arc.edge].capacity==0) {
            result.push_back(candidates[arc.candidate]);
        }
    }
    std::sort(result.begin(),result.end(),[](const auto& left,const auto& right) {
        if (left.alpha!=right.alpha) return left.alpha<right.alpha;
        return left.beta<right.beta;
    });
    return result;
}

std::optional<std::size_t> formal_irrep_dimension(
    const std::string& point_group,const std::string& symmetry) {
    const auto* definition=find_point_group(point_group);
    if (definition==nullptr || symmetry.empty()) return std::nullopt;
    const auto found=std::find_if(
        definition->irreps.begin(),definition->irreps.end(),[&](const auto& irrep) {
            return normalised_symmetry(std::string(irrep.label))==symmetry;
        });
    if (found==definition->irreps.end() || found->dimension==0u) {
        return std::nullopt;
    }
    return static_cast<std::size_t>(found->dimension);
}

std::vector<std::vector<std::size_t>> group_subsets_with_member_count(
    const std::vector<GroupCandidate>& groups,
    const std::vector<std::size_t>& candidates,
    const std::size_t target) {
    std::vector<std::vector<std::size_t>> result;
    std::vector<std::size_t> current;
    const auto visit=[&](const auto& self,const std::size_t begin,
                         const std::size_t count)->void {
        if (count==target) {
            if (current.size()>=2u) result.push_back(current);
            return;
        }
        for (std::size_t i=begin;i<candidates.size();++i) {
            const std::size_t group=candidates[i];
            const std::size_t size=groups[group].level.member_indices.size();
            if (size==0u || count+size>target) continue;
            current.push_back(group);
            self(self,i+1u,count+size);
            current.pop_back();
        }
    };
    visit(visit,0u,0u);
    return result;
}

GroupCandidate merge_spin_partition(
    const Wavefunction& wavefunction,
    const std::vector<GroupCandidate>& groups,
    const std::vector<std::size_t>& partition) {
    GroupCandidate merged=groups[partition.front()];
    for (std::size_t i=1u;i<partition.size();++i) {
        merged=merge_local_groups(wavefunction,merged,groups[partition[i]]);
    }
    return merged;
}

void collapse_split_spin_partitions(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::string& point_group,
    const SpinOverlapWorkspace& overlap_workspace,
    std::vector<GroupCandidate>& groups,
    std::set<std::size_t>& used_alpha,
    std::set<std::size_t>& used_beta,
    SpinCollapseResult& result) {
    struct ResidualCandidate {
        std::vector<std::size_t> alpha;
        std::vector<std::size_t> beta;
        double score=0.0;
    };
    std::vector<ResidualCandidate> candidates;
    const auto add_full_against_split=[&](const Spin full_spin) {
        const Spin split_spin=full_spin==Spin::Alpha?Spin::Beta:Spin::Alpha;
        const auto& used_full=full_spin==Spin::Alpha?used_alpha:used_beta;
        const auto& used_split=split_spin==Spin::Alpha?used_alpha:used_beta;
        for (std::size_t full=0;full<groups.size();++full) {
            if (group_spin(wavefunction,groups[full])!=full_spin ||
                used_full.count(full)>0u) continue;
            const std::string symmetry=normalised_symmetry(
                groups[full].level.metadata.symmetry);
            const auto dimension=formal_irrep_dimension(point_group,symmetry);
            if (!dimension || *dimension<2u ||
                groups[full].level.member_indices.size()!=*dimension) continue;
            std::vector<std::size_t> pieces;
            for (std::size_t split=0;split<groups.size();++split) {
                if (group_spin(wavefunction,groups[split])!=split_spin ||
                    used_split.count(split)>0u ||
                    groups[split].level.member_indices.empty() ||
                    groups[split].level.member_indices.size()>=*dimension ||
                    normalised_symmetry(
                        groups[split].level.metadata.symmetry)!=symmetry) continue;
                pieces.push_back(split);
            }
            for (auto partition:group_subsets_with_member_count(
                     groups,pieces,*dimension)) {
                std::vector<std::size_t> alpha;
                std::vector<std::size_t> beta;
                if (full_spin==Spin::Alpha) {
                    alpha={full};
                    beta=std::move(partition);
                } else {
                    alpha=std::move(partition);
                    beta={full};
                }
                std::vector<std::size_t> alpha_members;
                std::vector<std::size_t> beta_members;
                for (const auto group:alpha) {
                    alpha_members.insert(alpha_members.end(),
                        groups[group].level.member_indices.begin(),
                        groups[group].level.member_indices.end());
                }
                for (const auto group:beta) {
                    beta_members.insert(beta_members.end(),
                        groups[group].level.member_indices.begin(),
                        groups[group].level.member_indices.end());
                }
                const double score=subspace_overlap(
                    overlap_workspace,alpha_members,beta_members);
                if (score>=0.80) {
                    candidates.push_back(
                        {std::move(alpha),std::move(beta),score});
                }
            }
        }
    };
    add_full_against_split(Spin::Alpha);
    add_full_against_split(Spin::Beta);
    std::sort(candidates.begin(),candidates.end(),[](const auto& left,
                                                     const auto& right) {
        if (left.score!=right.score) return left.score>right.score;
        if (left.alpha!=right.alpha) return left.alpha<right.alpha;
        return left.beta<right.beta;
    });

    for (const auto& candidate:candidates) {
        if (std::any_of(candidate.alpha.begin(),candidate.alpha.end(),
                [&](const auto group) {return used_alpha.count(group)>0u;}) ||
            std::any_of(candidate.beta.begin(),candidate.beta.end(),
                [&](const auto group) {return used_beta.count(group)>0u;})) {
            continue;
        }
        GroupCandidate alpha=merge_spin_partition(
            wavefunction,groups,candidate.alpha);
        const GroupCandidate beta=merge_spin_partition(
            wavefunction,groups,candidate.beta);
        combine_spin_occupations(
            wavefunction,overlap_workspace,alpha,beta,
            options.filter.occupation_threshold);

        const std::size_t alpha_primary=candidate.alpha.front();
        groups[alpha_primary]=std::move(alpha);
        groups[alpha_primary].suppressed_spin_counterpart=false;
        for (std::size_t i=1u;i<candidate.alpha.size();++i) {
            auto& piece=groups[candidate.alpha[i]];
            piece.include=false;
            piece.suppressed_spin_counterpart=true;
        }
        for (const auto group:candidate.beta) {
            groups[group].include=false;
            groups[group].suppressed_spin_counterpart=true;
        }
        used_alpha.insert(candidate.alpha.begin(),candidate.alpha.end());
        used_beta.insert(candidate.beta.begin(),candidate.beta.end());
        ++result.paired_groups;
        result.paired_members+=
            groups[alpha_primary].level.member_indices.size();
    }
}

SpinCollapseResult collapse_spin_counterparts(
                                const Wavefunction& wavefunction,
                                const MODiagramOptions& options,
                                const std::string& point_group,
    std::vector<GroupCandidate>& groups) {
    SpinCollapseResult result;
    const bool has_alpha_orbitals=std::any_of(
        wavefunction.orbitals.begin(),wavefunction.orbitals.end(),
        [](const auto& orbital) {return orbital.spin==Spin::Alpha;});
    const bool has_beta_orbitals=std::any_of(
        wavefunction.orbitals.begin(),wavefunction.orbitals.end(),
        [](const auto& orbital) {return orbital.spin==Spin::Beta;});
    if (!has_alpha_orbitals || !has_beta_orbitals ||
        wavefunction.ao_overlap.empty()) return result;
    const SpinOverlapWorkspace overlap_workspace(wavefunction);
    std::vector<SpinPairCandidate> candidates;
    for (std::size_t a=0;a<groups.size();++a) {
        if (group_spin(wavefunction,groups[a])!=Spin::Alpha) continue;
        for (std::size_t b=0;b<groups.size();++b) {
            if (group_spin(wavefunction,groups[b])!=Spin::Beta ||
                groups[a].level.member_indices.size()!=
                    groups[b].level.member_indices.size()) continue;
            const std::string sa=normalised_symmetry(
                groups[a].level.metadata.symmetry);
            const std::string sb=normalised_symmetry(
                groups[b].level.metadata.symmetry);
            if (!sa.empty() && !sb.empty() && sa!=sb) continue;
            const int alpha_family=dominant_metal_family(groups[a].level);
            const int beta_family=dominant_metal_family(groups[b].level);
            const bool same_valid_irrep=!sa.empty() && !sb.empty() && sa==sb;
            if (!same_valid_irrep && alpha_family>=0 && beta_family>=0 &&
                alpha_family!=beta_family) continue;
            const double score=subspace_overlap(
                overlap_workspace,groups[a],groups[b]);
            // The squared S-metric subspace score is invariant to rotations
            // inside a degenerate block.  UDFT spin polarisation can lower it
            // well below the old 0.50 cutoff even for the same spatial d
            // block; metal-family and symmetry compatibility guard the more
            // permissive floor against unrelated matches.
            if (score>=0.30) candidates.push_back({a,b,score});
        }
    }
    const auto pairs=maximum_cardinality_spin_matching(
        candidates,groups.size());
    std::set<std::size_t> used_alpha;
    std::set<std::size_t> used_beta;
    for (const auto& pair:pairs) {
        used_alpha.insert(pair.alpha);
        used_beta.insert(pair.beta);
        combine_spin_occupations(
            wavefunction,overlap_workspace,
            groups[pair.alpha],groups[pair.beta],
            options.filter.occupation_threshold);
        groups[pair.beta].include=false;
        groups[pair.beta].suppressed_spin_counterpart=true;
        ++result.paired_groups;
        result.paired_members+=groups[pair.alpha].level.member_indices.size();
    }
    collapse_split_spin_partitions(
        wavefunction,options,point_group,overlap_workspace,groups,
        used_alpha,used_beta,result);
    return result;
}

double member_metal_weight(const Wavefunction& wavefunction,
                           const MolecularOrbital& orbital,
                           const int angular_momentum,
                           const LigandScope* scope) {
    double result=0.0;
    for (const auto& contribution:orbital.chemistry.ao_contributions) {
        if (contribution.atom_index>=wavefunction.atoms.size()) continue;
        if (scope!=nullptr && scope->available) {
            if (contribution.atom_index!=scope->metal) continue;
        } else if (!transition_metal(
                       wavefunction.atoms[contribution.atom_index].atomic_number)) {
            continue;
        }
        if (contribution.angular_momentum==angular_momentum) {
            result+=contribution.weight;
        }
    }
    return result;
}

double member_ligand_p_weight(const Wavefunction& wavefunction,
                              const MolecularOrbital& orbital,
                              const LigandScope* scope,
                              const bool direct_only) {
    double result=0.0;
    for (const auto& contribution:orbital.chemistry.ao_contributions) {
        if (contribution.atom_index>=wavefunction.atoms.size()) continue;
        if (scope!=nullptr && scope->available) {
            const auto& allowed=direct_only?scope->direct_donors
                                           :scope->fragment_atoms;
            if (allowed.count(contribution.atom_index)==0u) continue;
        } else if (transition_metal(
                       wavefunction.atoms[contribution.atom_index].atomic_number)) {
            continue;
        }
        if (contribution.angular_momentum==1) result+=contribution.weight;
    }
    return result;
}

GroupCandidate make_group_candidate(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const FrontierOrbitals& frontier,
    const std::vector<OrbitalMetadata>& metadata,
    const std::vector<OrbitalAnnotation>& annotations,
    const std::vector<OrbitalLabel>& labels,
    const std::size_t base,
    const LigandScope* scope) {
    GroupCandidate candidate;
    candidate.base_index=base;
    if (base>=wavefunction.orbitals.size() || base>=metadata.size() ||
        base>=labels.size()) return candidate;

    const std::size_t requested=std::max<std::size_t>(1u,labels[base].group_size);
    const std::size_t end=std::min(wavefunction.orbitals.size(),base+requested);
    auto& level=candidate.level;
    level.metadata=metadata[base];
    level.annotation=annotations[base];
    level.chemistry=wavefunction.orbitals[base].chemistry;
    level.metadata.display_label=std::to_string(base+1u);
    level.metadata.degeneracy_size=end-base;
    level.metadata.selected=false;

    double energy_sum=0.0;
    double energy_min=std::numeric_limits<double>::infinity();
    double energy_max=-std::numeric_limits<double>::infinity();
    double sigma_weighted=0.0;
    double pi_weighted=0.0;
    double channel_weight=0.0;
    double overlap=0.0;
    double valence_weight=0.0;
    bool all_valence=true;

    for (std::size_t index=base;index<end;++index) {
        const auto& orbital=wavefunction.orbitals[index];
        level.member_indices.push_back(index);
        level.member_electrons.push_back(electron_glyphs_for_orbital(
            orbital,frontier.separate_spin_sets));
        level.total_occupation+=static_cast<double>(orbital.occupation);
        energy_sum+=orbital.energy_hartree;
        energy_min=std::min(energy_min,orbital.energy_hartree);
        energy_max=std::max(energy_max,orbital.energy_hartree);
        level.homo=level.homo || (frontier.homo && *frontier.homo==index);
        level.lumo=level.lumo || (frontier.lumo && *frontier.lumo==index);
        level.metadata.selected=level.metadata.selected ||
                                options.selected_index==index;
        candidate.selected_by_reference=candidate.selected_by_reference ||
            (orbital.chemistry.available && orbital.chemistry.valence_manifold);
        all_valence=all_valence && orbital.chemistry.available &&
                    orbital.chemistry.valence_manifold;
        valence_weight+=orbital.chemistry.valence_weight;

        level.metal_s_weight+=member_metal_weight(wavefunction,orbital,0,scope);
        level.metal_p_weight+=member_metal_weight(wavefunction,orbital,1,scope);
        level.metal_d_weight+=member_metal_weight(wavefunction,orbital,2,scope);
        level.direct_ligand_p_weight+=member_ligand_p_weight(
            wavefunction,orbital,scope,true);
        level.ligand_p_weight+=member_ligand_p_weight(
            wavefunction,orbital,scope,false);

        for (const auto& interaction:orbital.chemistry.interactions) {
            if (!scoped_metal_ligand_pair(wavefunction,interaction,scope)) continue;
            const double bounded=bounded_overlap_character(
                interaction.overlap_character);
            const double magnitude=std::abs(bounded);
            const double determined=interaction.channel.sigma+
                                    interaction.channel.pi+
                                    interaction.channel.delta+
                                    interaction.channel.phi;
            sigma_weighted+=magnitude*interaction.channel.sigma;
            pi_weighted+=magnitude*interaction.channel.pi;
            channel_weight+=magnitude*determined;
            overlap+=bounded;
        }
    }

    const double count=static_cast<double>(std::max<std::size_t>(1u,end-base));
    level.layout_energy_hartree=energy_sum/count;
    level.metadata.energy_hartree=level.layout_energy_hartree;
    level.energy_spread_hartree=std::max(0.0,energy_max-energy_min);
    level.metadata.occupation=static_cast<float>(level.total_occupation/count);
    level.electrons=level.member_electrons.empty()?ElectronGlyphs{}
                                                   :level.member_electrons.front();
    level.metal_s_weight/=count;
    level.metal_p_weight/=count;
    level.metal_d_weight/=count;
    level.direct_ligand_p_weight/=count;
    level.ligand_p_weight/=count;
    const double donor_count=scope!=nullptr && scope->available
        ?static_cast<double>(std::max<std::size_t>(1u,scope->direct_donors.size()))
        :1.0;
    level.metal_ligand_overlap=std::clamp(
        overlap/(count*donor_count),-1.0,1.0);
    valence_weight/=count;
    if (channel_weight>1.0e-12) {
        level.sigma_fraction=std::clamp(sigma_weighted/channel_weight,0.0,1.0);
        level.pi_fraction=std::clamp(pi_weighted/channel_weight,0.0,1.0);
    }

    // Producer FCHK commonly omits individual irreps.  When the full
    // symmetry projector cannot retain a label for a numerically mixed set,
    // recover only the unambiguous central-metal valence cases.  This is a
    // raw-MO subspace assignment, not a molecule-name template.
    if (normalised_symmetry(level.metadata.symmetry).empty()) {
        const std::string& point_group=wavefunction.point_group_detected;
        const std::size_t degeneracy=level.metadata.degeneracy_size;
        if (point_group=="Oh") {
            if (degeneracy==3u && level.metal_d_weight>=0.08) {
                level.metadata.symmetry="T2g";
            } else if (degeneracy==2u && level.metal_d_weight>=0.08) {
                level.metadata.symmetry="Eg";
            } else if (degeneracy==3u && level.metal_p_weight>=0.08) {
                level.metadata.symmetry="T1u";
            } else if (degeneracy==1u && level.metal_s_weight>=0.08) {
                level.metadata.symmetry="A1g";
            }
        } else if (point_group=="Td") {
            if (degeneracy==2u && level.metal_d_weight>=0.08) {
                level.metadata.symmetry="E";
            } else if (degeneracy==3u &&
                       (level.metal_d_weight>=0.08 ||
                        level.metal_p_weight>=0.08)) {
                level.metadata.symmetry="T2";
            } else if (degeneracy==1u && level.metal_s_weight>=0.08) {
                level.metadata.symmetry="A1";
            }
        }
    }

    const double metal_valence=level.metal_s_weight+
                               level.metal_p_weight+
                               level.metal_d_weight;
    const bool occupied_group=level.total_occupation>options.filter.occupation_threshold;
    const double virtual_ceiling=frontier.lumo && *frontier.lumo<wavefunction.orbitals.size()
        ?wavefunction.orbitals[*frontier.lumo].energy_hartree+
             std::min(0.75,options.filter.virtual_window_hartree)
        :std::numeric_limits<double>::infinity();
    const bool energy_relevant=occupied_group ||
        level.layout_energy_hartree<=virtual_ceiling;
    const bool raw_relevant=energy_relevant &&
        level.layout_energy_hartree>=
            options.filter.core_energy_cutoff_hartree &&
        (metal_valence>=0.08 ||
         (valence_weight>=0.20 &&
          std::abs(level.metal_ligand_overlap)>=
              options.weak_metal_ligand_overlap));
    candidate.selected_by_raw=!candidate.selected_by_reference && raw_relevant;
    level.raw_data_fallback=candidate.selected_by_raw;
    const bool selected_inspection_only=
        options.hide_ligand_centred_intermediates &&
        level.metadata.selected;
    candidate.include=(candidate.selected_by_reference &&
                       level.layout_energy_hartree>=
                           options.filter.core_energy_cutoff_hartree &&
                       energy_relevant) ||
                      candidate.selected_by_raw ||
                      (level.metadata.selected && !selected_inspection_only);

    // A sigma-framework group that the minimal canonical manifold missed is
    // recovered from the complete MO block rather than silently disappearing.
    if (!candidate.include && energy_relevant &&
        level.layout_energy_hartree>=
            options.filter.core_energy_cutoff_hartree &&
        level.sigma_fraction>=0.55 &&
        std::abs(level.metal_ligand_overlap)>=options.weak_metal_ligand_overlap) {
        candidate.include=true;
        candidate.selected_by_raw=true;
        level.raw_data_fallback=true;
    }

    if (all_valence && level.annotation.family=="unavailable") {
        if (level.pi_fraction>=0.55) level.annotation.family="pi";
        else if (level.sigma_fraction>=0.55) level.annotation.family="sigma";
        if (level.annotation.family!="unavailable") {
            level.annotation.family_source=AnnotationSource::Derived;
            level.annotation.family_confidence=std::max(
                level.pi_fraction,level.sigma_fraction);
        }
    }

    if (std::abs(level.metal_ligand_overlap)<
        options.weak_metal_ligand_overlap) {
        level.annotation.bonding_class=BondingClass::Nonbonding;
    } else if (level.metal_ligand_overlap>0.0) {
        level.annotation.bonding_class=BondingClass::Bonding;
    } else {
        level.annotation.bonding_class=BondingClass::Antibonding;
    }
    level.annotation.bonding_source=AnnotationSource::Derived;
    level.annotation.bonding_confidence=std::clamp(
        std::abs(level.metal_ligand_overlap)/0.10,0.0,1.0);
    return candidate;
}

struct RawPiPair {
    std::size_t lower=0;
    std::size_t upper=0;
    PiInteractionKind kind=PiInteractionKind::Coupled;
    double split=0.0;
    double confidence=0.0;
    std::size_t retained=0;
    std::string symmetry;
};

bool local_pi_irrep(const std::string& point_group,
                    const std::string& symmetry) {
    if (symmetry.empty() || symmetry=="?" || symmetry=="n/a") return false;
    const auto decomposition=decompose_metal_ao_shell(
        point_group,MetalAOShell::D);
    if (!decomposition) return false;
    return std::any_of(
        decomposition->begin(),decomposition->end(),[&](const auto& copy) {
            return normalised_symmetry(std::string(copy.label))==symmetry;
        });
}

std::vector<RawPiPair> find_pi_pairs(
    const Wavefunction& wavefunction,
    std::vector<GroupCandidate>& groups,
    const MODiagramOptions& options,
    const std::string& point_group,
    const LigandScope* scope) {
    struct ScoredPair { RawPiPair pair; double score=0.0; };
    std::vector<ScoredPair> scored;
    const LocalLigandPiRole role=scope!=nullptr
        ?scope->pi_role:LocalLigandPiRole::Unresolved;
    for (std::size_t lower=0;lower<groups.size();++lower) {
        if (groups[lower].suppressed_spin_counterpart) continue;
        const auto& a=groups[lower].level;
        if (a.pi_fraction<0.60 ||
            a.metal_d_weight+a.ligand_p_weight<0.18) continue;
        const std::string symmetry_a=normalised_symmetry(a.metadata.symmetry);
        for (std::size_t upper=lower+1u;upper<groups.size();++upper) {
            if (groups[upper].suppressed_spin_counterpart) continue;
            const auto& b=groups[upper].level;
            if (b.layout_energy_hartree<=a.layout_energy_hartree ||
                b.metadata.degeneracy_size!=a.metadata.degeneracy_size ||
                group_spin(wavefunction,groups[upper])!=
                    group_spin(wavefunction,groups[lower]) ||
                b.pi_fraction<0.60 ||
                b.metal_d_weight+b.ligand_p_weight<0.18) continue;
            const std::string symmetry_b=normalised_symmetry(
                b.metadata.symmetry);
            const bool a_known=local_pi_irrep(point_group,symmetry_a);
            const bool b_known=local_pi_irrep(point_group,symmetry_b);
            std::string symmetry;
            if (a_known && b_known && symmetry_a==symmetry_b) {
                symmetry=symmetry_a;
            } else if (a_known &&
                       (symmetry_b.empty() || symmetry_b=="?" ||
                        symmetry_b=="n/a")) {
                symmetry=symmetry_a;
            } else if (b_known &&
                       (symmetry_a.empty() || symmetry_a=="?" ||
                        symmetry_a=="n/a")) {
                symmetry=symmetry_b;
            } else {
                continue;
            }
            if (!local_pi_irrep(point_group,symmetry)) continue;
            if (role==LocalLigandPiRole::SigmaOnly) continue;
            const double split=b.layout_energy_hartree-a.layout_energy_hartree;
            if (split>1.50) continue;
            const double lower_ligand_minus_metal=
                a.ligand_p_weight-a.metal_d_weight;
            const double upper_ligand_minus_metal=
                b.ligand_p_weight-b.metal_d_weight;
            const double donor=lower_ligand_minus_metal-
                               upper_ligand_minus_metal;
            const double acceptor=-donor;
            const double opposite=std::max(0.0,
                -a.metal_ligand_overlap*b.metal_ligand_overlap);
            const bool weak=split<=options.weak_pi_split_hartree &&
                std::abs(a.metal_ligand_overlap)<=
                    options.weak_metal_ligand_overlap &&
                std::abs(b.metal_ligand_overlap)<=
                    options.weak_metal_ligand_overlap;
            const bool complementary_composition=
                lower_ligand_minus_metal*upper_ligand_minus_metal<0.0 &&
                std::abs(lower_ligand_minus_metal)>=0.05 &&
                std::abs(upper_ligand_minus_metal)>=0.05;
            if (!weak && !complementary_composition) continue;
            const bool directed=a.metal_ligand_overlap*
                                    b.metal_ligand_overlap<0.0 &&
                std::abs(a.metal_ligand_overlap)>=
                    options.weak_metal_ligand_overlap &&
                std::abs(b.metal_ligand_overlap)>=
                    options.weak_metal_ligand_overlap;
            const double expected_contrast=
                role==LocalLigandPiRole::Donor?donor:
                (role==LocalLigandPiRole::Acceptor?acceptor:
                 std::max(donor,acceptor));
            const bool composition_directed=expected_contrast>=0.18;
            if (!weak && !directed && !composition_directed) continue;
            if (weak && std::max(a.metal_d_weight,b.metal_d_weight)<0.08) {
                continue;
            }
            if ((role==LocalLigandPiRole::Donor ||
                 role==LocalLigandPiRole::Acceptor) &&
                expected_contrast<0.18) continue;
            const double complement=std::abs(donor);
            const double pi_quality=std::min(a.pi_fraction,b.pi_fraction);
            const double score=2.0*pi_quality+2.0*std::max(0.0,expected_contrast)+
                               2.0*std::sqrt(opposite)-0.08*split;
            if (score<0.75) continue;

            RawPiPair pair;
            pair.lower=lower;
            pair.upper=upper;
            pair.split=split;
            pair.symmetry=a_known?a.metadata.symmetry:b.metadata.symmetry;
            if (weak) {
                pair.kind=PiInteractionKind::WeakNearNonbonding;
            } else if (role==LocalLigandPiRole::Donor) {
                pair.kind=PiInteractionKind::Donor;
            } else if (role==LocalLigandPiRole::Acceptor) {
                pair.kind=PiInteractionKind::Acceptor;
            } else if (donor-acceptor>=0.15) {
                pair.kind=PiInteractionKind::Donor;
            } else if (acceptor-donor>=0.15) {
                pair.kind=PiInteractionKind::Acceptor;
            } else {
                pair.kind=PiInteractionKind::Coupled;
            }
            pair.confidence=std::clamp(
                0.45*pi_quality+0.35*std::min(1.0,complement)+
                0.20*std::min(1.0,std::sqrt(opposite)/0.05),0.0,1.0);
            const double a_metal=a.metal_s_weight+a.metal_p_weight+a.metal_d_weight;
            const double b_metal=b.metal_s_weight+b.metal_p_weight+b.metal_d_weight;
            pair.retained=b_metal>a_metal?upper:lower;
            scored.push_back({pair,score});
        }
    }

    // In a weak-field d manifold the two crystal-field components have
    // different irreps, so same-symmetry pi pairing cannot find them.  Treat
    // an E/T2 (Td) or Eg/T2g (Oh) pair as one unresolved, approximately
    // nonbonding split only when both the energy gap and the actual
    // metal-ligand mixing are below the user-facing thresholds.
    if (role!=LocalLigandPiRole::SigmaOnly) {
    for (std::size_t first=0;first<groups.size();++first) {
        if (groups[first].suppressed_spin_counterpart) continue;
        const auto& a=groups[first].level;
        if (a.metal_d_weight<0.60) continue;
        for (std::size_t second=first+1u;second<groups.size();++second) {
            if (groups[second].suppressed_spin_counterpart) continue;
            const auto& b=groups[second].level;
            if (b.metal_d_weight<0.60 ||
                group_spin(wavefunction,groups[second])!=
                    group_spin(wavefunction,groups[first])) continue;
            const std::string sa=normalised_symmetry(a.metadata.symmetry);
            const std::string sb=normalised_symmetry(b.metadata.symmetry);
            const bool tetrahedral=(sa=="e" && sb=="t2") ||
                                   (sa=="t2" && sb=="e");
            const bool octahedral=(sa=="eg" && sb=="t2g") ||
                                  (sa=="t2g" && sb=="eg");
            if (!tetrahedral && !octahedral) continue;
            const double split=std::abs(
                b.layout_energy_hartree-a.layout_energy_hartree);
            if (split>options.weak_pi_split_hartree ||
                std::abs(a.metal_ligand_overlap)>
                    options.weak_metal_ligand_overlap ||
                std::abs(b.metal_ligand_overlap)>
                    options.weak_metal_ligand_overlap) continue;
            RawPiPair pair;
            if (a.layout_energy_hartree<=b.layout_energy_hartree) {
                pair.lower=first;
                pair.upper=second;
            } else {
                pair.lower=second;
                pair.upper=first;
            }
            pair.kind=PiInteractionKind::WeakNearNonbonding;
            pair.split=split;
            pair.confidence=std::clamp(
                1.0-split/options.weak_pi_split_hartree,0.0,1.0);
            pair.retained=a.metal_d_weight>=b.metal_d_weight?first:second;
            pair.symmetry=tetrahedral?"E/T2":"Eg/T2g";
            scored.push_back({pair,10.0+pair.confidence});
        }
    }
    }

    std::sort(scored.begin(),scored.end(),[](const auto& a,const auto& b) {
        return a.score>b.score;
    });
    std::set<std::size_t> used;
    std::vector<RawPiPair> result;
    for (const auto& item:scored) {
        if (used.count(item.pair.lower) || used.count(item.pair.upper)) continue;
        used.insert(item.pair.lower);
        used.insert(item.pair.upper);
        groups[item.pair.lower].include=true;
        groups[item.pair.upper].include=true;
        for (const auto group_index:{item.pair.lower,item.pair.upper}) {
            auto& symmetry=groups[group_index].level.metadata.symmetry;
            const std::string current=normalised_symmetry(symmetry);
            if ((current.empty() || current=="?" || current=="n/a") &&
                !item.pair.symmetry.empty()) {
                symmetry=item.pair.symmetry;
            }
        }
        if (item.pair.kind==PiInteractionKind::WeakNearNonbonding) {
            auto& lower=groups[item.pair.lower].level;
            auto& upper=groups[item.pair.upper].level;
            lower.approximate_nonbonding=true;
            upper.approximate_nonbonding=true;
            lower.annotation.bonding_class=BondingClass::Nonbonding;
            upper.annotation.bonding_class=BondingClass::Nonbonding;
            const std::size_t dropped=item.pair.retained==item.pair.lower
                ?item.pair.upper:item.pair.lower;
            if (options.hide_ligand_centred_intermediates ||
                !groups[dropped].level.metadata.selected) {
                groups[dropped].include=false;
            }
        }
        result.push_back(item.pair);
    }
    return result;
}

} // namespace

const char* pi_interaction_kind_name(const PiInteractionKind kind) noexcept {
    switch (kind) {
        case PiInteractionKind::Donor: return "pi-donor splitting";
        case PiInteractionKind::Acceptor: return "pi-acceptor splitting";
        case PiInteractionKind::WeakNearNonbonding:
            return "weak-field split; approximately nonbonding";
        default: return "pi-coupled splitting";
    }
}

OrbitalAnnotation annotate_orbital(const MolecularOrbital& orbital) {
    if (orbital.chemistry.available) {
        return chemistry_annotation(orbital);
    }
    return annotate_orbital_legacy(orbital);
}

DiagramSelectionPlan build_valence_selection_plan(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options,
    const std::vector<OrbitalMetadata>& metadata) {
    if (!chemistry_available(wavefunction)) {
        return build_valence_selection_plan_legacy(
            wavefunction,options,metadata);
    }

    DiagramSelectionPlan plan;
    const LigandFieldEnvironment ligand_field=
        analyse_ligand_field_environment(wavefunction);
    const LigandScope ligand_scope=make_ligand_scope(
        wavefunction,ligand_field);
    const LigandScope* scope=ligand_scope.available?&ligand_scope:nullptr;
    std::size_t raw_supplements=0u;
    for (std::size_t i=0;i<wavefunction.orbitals.size();++i) {
        const auto& orbital=wavefunction.orbitals[i];
        if (orbital.chemistry.available &&
            orbital.chemistry.valence_manifold) {
            plan.included_indices.push_back(i);
            continue;
        }

        double metal=0.0;
        for (const auto& contribution:orbital.chemistry.ao_contributions) {
            const bool selected_metal=scope!=nullptr
                ?contribution.atom_index==scope->metal
                :(contribution.atom_index<wavefunction.atoms.size() &&
                  transition_metal(wavefunction.atoms[
                      contribution.atom_index].atomic_number));
            if (selected_metal && contribution.angular_momentum<=2) {
                metal+=contribution.weight;
            }
        }
        double interaction=0.0;
        for (const auto& pair:orbital.chemistry.interactions) {
            if (scoped_metal_ligand_pair(wavefunction,pair,scope)) {
                interaction+=std::abs(bounded_overlap_character(
                    pair.overlap_character));
            }
        }
        if (metal>=0.08 || interaction>=options.weak_metal_ligand_overlap ||
            (i==options.selected_index &&
             !options.hide_ligand_centred_intermediates)) {
            plan.included_indices.push_back(i);
            ++raw_supplements;
        }
    }
    const auto labels=build_orbital_labels(
        wavefunction.orbitals,options.degeneracy);
    plan.included_indices=expand_degenerate(
        std::move(plan.included_indices),labels);
    plan.hidden_count=metadata.size()>plan.included_indices.size()
        ?metadata.size()-plan.included_indices.size():0u;

    for (const auto index:plan.included_indices) {
        if (index>=wavefunction.orbitals.size()) continue;
        if (occupied(wavefunction.orbitals[index],
                     options.filter.occupation_threshold)) {
            ++plan.valence_occupied_count;
        } else {
            ++plan.frontier_virtual_count;
        }
    }

    std::ostringstream summary;
    summary<<"chemical-valence reference: "
           <<plan.included_indices.size()<<'/'<<metadata.size()
           <<" canonical MOs; occupied="<<plan.valence_occupied_count
           <<"; virtual="<<plan.frontier_virtual_count
           <<"; raw-MO supplements="<<raw_supplements
           <<"; core/polarisation/Rydberg hidden";
    plan.summary=summary.str();
    return plan;
}

MODiagramData build_mo_diagram_data(
    const Wavefunction& wavefunction,
    const MODiagramOptions& options) {
    if (!chemistry_available(wavefunction)) {
        return build_mo_diagram_data_legacy(wavefunction,options);
    }

    MODiagramData data;
    data.mode=MODiagramMode::ValenceCentral;
    data.plan=choose_diagram_plan(wavefunction);
    data.plan.machine_reason=
        "COV S-metric minimal atomic chemical-valence reference";
    data.frontier=find_frontier_orbitals(
        wavefunction.orbitals,options.filter.occupation_threshold);
    data.metadata=build_orbital_metadata(
        wavefunction,options.selected_index,
        options.degeneracy,options.filter);

    const LigandFieldEnvironment ligand_field=
        analyse_ligand_field_environment(wavefunction);
    const LigandScope ligand_scope=make_ligand_scope(
        wavefunction,ligand_field);
    if (ligand_field.available()) {
        data.ligand_field_point_group=ligand_field.local_point_group();
        data.ligand_field_geometry_id=ligand_field.geometry_machine_id();
        data.ligand_field_geometry_name=ligand_field.geometry_name();
        data.ligand_field_coordination_number=
            ligand_field.coordination_number();
        data.ligand_field_metal_atom=ligand_field.metal_atom;
        data.ligand_field_ligand_atoms=ligand_field.ligand_atoms;
        data.ligand_field_confidence=ligand_field.confidence;
        data.ligand_field_angular_rms=ligand_field.angular_rms;
        data.ligand_field_shape_measure=ligand_field.shape_measure;
        data.ligand_field_radial_cv=ligand_field.radial_cv;
    }

    data.annotations.reserve(wavefunction.orbitals.size());
    for (const auto& orbital:wavefunction.orbitals) {
        data.annotations.push_back(annotate_orbital(orbital));
    }

    data.selection=build_valence_selection_plan(
        wavefunction,options,data.metadata);
    const auto labels=build_orbital_labels(
        wavefunction.orbitals,options.degeneracy);
    std::vector<GroupCandidate> groups;
    groups.reserve(wavefunction.orbitals.size());
    for (std::size_t base=0;base<wavefunction.orbitals.size();) {
        groups.push_back(make_group_candidate(
            wavefunction,options,data.frontier,data.metadata,
            data.annotations,labels,base,
            ligand_scope.available?&ligand_scope:nullptr));
        const std::size_t step=base<labels.size()
            ?std::max<std::size_t>(1u,labels[base].group_size):1u;
        base+=step;
    }

    if (ligand_field.available()) {
        recover_local_ligand_field_symmetry(
            wavefunction,ligand_field,groups);
        merge_local_pseudodegenerate_groups(
            wavefunction,ligand_field,groups);
        recover_local_ligand_field_symmetry(
            wavefunction,ligand_field,groups);
    }

    // Some FCHK producers omit member-level irreps even though the complete
    // exactly-degenerate subspace has an unambiguous ligand-field label.  A
    // recovered Eg/T2g/E/T2 label belongs to every canonical member of that
    // subspace; copy it back to the member metadata so per-line selection and
    // tooltips do not regress to N/A.
    for (const auto& group:groups) {
        const std::string recovered=group.level.metadata.symmetry;
        const std::string normalised=normalised_symmetry(recovered);
        if (normalised.empty() || normalised=="?" || normalised=="n/a") continue;
        for (const auto member:group.level.member_indices) {
            if (member>=data.metadata.size()) continue;
            const std::string current=normalised_symmetry(
                data.metadata[member].symmetry);
            if (current.empty() || current=="?" || current=="n/a") {
                data.metadata[member].symmetry=recovered;
            }
        }
    }

    const SpinCollapseResult spin_collapse=collapse_spin_counterparts(
        wavefunction,options,
        data.ligand_field_point_group.empty()
            ?wavefunction.point_group_detected
            :data.ligand_field_point_group,
        groups);
    data.spin_counterpart_pair_count=spin_collapse.paired_groups;

    // Matching may recover a label from either spin channel.  Apply it to
    // both canonical members so the browser, hover text and exported member
    // metadata agree with the spatial-row diagram.
    for (const auto& group:groups) {
        const std::string recovered=group.level.metadata.symmetry;
        const std::string normalised=normalised_symmetry(recovered);
        if (normalised.empty() || normalised=="?" || normalised=="n/a") continue;
        for (const auto member:group.level.member_indices) {
            if (member<data.metadata.size()) {
                const std::string current=normalised_symmetry(
                    data.metadata[member].symmetry);
                if (current.empty() || current=="?" || current=="n/a") {
                    data.metadata[member].symmetry=recovered;
                }
            }
        }
        for (const auto counterpart:group.level.member_spin_counterparts) {
            if (counterpart>=data.metadata.size()) continue;
            const std::string current=normalised_symmetry(
                data.metadata[counterpart].symmetry);
            if (current.empty() || current=="?" || current=="n/a") {
                data.metadata[counterpart].symmetry=recovered;
            }
        }
    }

    const auto raw_pairs=find_pi_pairs(
        wavefunction,groups,options,data.ligand_field_point_group,
        ligand_scope.available?&ligand_scope:nullptr);
    for (const auto& pair:raw_pairs) {
        for (const auto group_index:{pair.lower,pair.upper}) {
            if (group_index>=groups.size()) continue;
            const auto& group=groups[group_index].level;
            const std::string recovered=group.metadata.symmetry;
            const std::string normalised=normalised_symmetry(recovered);
            if (normalised.empty() || normalised=="?" ||
                normalised=="n/a") continue;
            for (const auto member:group.member_indices) {
                if (member>=data.metadata.size()) continue;
                const std::string current=normalised_symmetry(
                    data.metadata[member].symmetry);
                if (current.empty() || current=="?" || current=="n/a") {
                    data.metadata[member].symmetry=recovered;
                }
            }
        }
    }

    // Keep the reduced diagram readable.  The budget applies to degenerate
    // rows rather than canonical MOs.  In the compact ligand-field view the
    // canonical row set is purely chemical and therefore independent of the
    // MO currently inspected in the browser or 3-D viewport.  Both members
    // of a resolved donor/acceptor pair are never removed; the selected row
    // is pinned only in the expanded view.
    std::size_t row_budget=options.max_levels>0u
        ?std::max<std::size_t>(4u,options.max_levels)
        :std::clamp<std::size_t>(2u*options.neighbourhood+1u,10u,48u);
    std::vector<bool> essential(groups.size(),false);
    for (std::size_t group=0;group<groups.size();++group) {
        essential[group]=!options.hide_ligand_centred_intermediates &&
                         groups[group].level.metadata.selected;
    }
    for (const auto& pair:raw_pairs) {
        essential[pair.retained]=true;
        if (pair.kind!=PiInteractionKind::WeakNearNonbonding) {
            essential[pair.lower]=true;
            essential[pair.upper]=true;
        }
    }
    std::size_t hidden_intermediate_count=0u;
    if (options.hide_ligand_centred_intermediates) {
        const std::string local_group=data.ligand_field_point_group.empty()
            ?wavefunction.point_group_detected
            :data.ligand_field_point_group;
        std::map<std::string,std::size_t> d_irrep_multiplicity;
        if (const auto decomposition=decompose_metal_ao_shell(
                local_group,MetalAOShell::D)) {
            for (const auto& copy:*decomposition) {
                ++d_irrep_multiplicity[normalised_symmetry(
                    std::string(copy.label))];
            }
        }
        const bool supported_ligand_field=!d_irrep_multiplicity.empty();
        std::map<std::string,std::vector<std::size_t>> candidates;
        for (std::size_t group=0;group<groups.size();++group) {
            if (!groups[group].include) continue;
            const std::string symmetry=normalised_symmetry(
                groups[group].level.metadata.symmetry);
            if (d_irrep_multiplicity.count(symmetry)) {
                candidates[symmetry].push_back(group);
            }
        }
        std::set<std::size_t> anchors;
        for (auto& [symmetry,indices]:candidates) {
            std::sort(indices.begin(),indices.end(),[&](const auto a,const auto b) {
                const auto& left=groups[a].level;
                const auto& right=groups[b].level;
                const double left_score=left.metal_d_weight+
                    0.25*left.direct_ligand_p_weight;
                const double right_score=right.metal_d_weight+
                    0.25*right.direct_ligand_p_weight;
                if (left_score!=right_score) return left_score>right_score;
                return a<b;
            });
            const std::size_t needed=std::max<std::size_t>(
                1u,d_irrep_multiplicity[symmetry]);
            for (std::size_t i=0;i<std::min(needed,indices.size());++i) {
                anchors.insert(indices[i]);
            }
        }
        for (const auto group:anchors) essential[group]=true;

        // Preserve one bonding and one antibonding representative for each
        // central-metal s/p/d sigma-framework symmetry.  The measured sigma
        // channel is the gate; the point-group catalogue supplies the allowed
        // local labels without hard-coding Td/Oh names.
        using SigmaKey=std::pair<std::string,bool>;
        std::map<SigmaKey,std::pair<std::size_t,double>> sigma_representatives;
        std::set<std::string> local_spd_labels;
        if (const auto decomposition=decompose_metal_spd(local_group)) {
            for (const auto* block:{&decomposition->s,&decomposition->p,
                                   &decomposition->d}) {
                for (const auto& copy:*block) {
                    local_spd_labels.insert(normalised_symmetry(
                        std::string(copy.label)));
                }
            }
        }
        for (std::size_t group=0;group<groups.size();++group) {
            if (!groups[group].include) continue;
            const auto& level=groups[group].level;
            const double metal_sp=level.metal_s_weight+level.metal_p_weight;
            const double metal_spd=metal_sp+level.metal_d_weight;
            const std::string symmetry=normalised_symmetry(
                level.metadata.symmetry);
            if (symmetry.empty() || symmetry=="?" || symmetry=="n/a") continue;
            const bool local_sigma_label=local_spd_labels.count(symmetry)>0u;
            const double sigma_floor=local_sigma_label?0.50:0.55;
            if (level.sigma_fraction<sigma_floor || metal_spd<0.025) continue;
            const SigmaKey key{
                symmetry,
                level.metal_ligand_overlap<0.0};
            const double balanced_mixing=2.0*std::min(
                metal_spd,level.direct_ligand_p_weight);
            const double representative_score=balanced_mixing+
                0.20*std::min(1.0,std::abs(level.metal_ligand_overlap))+
                0.05*metal_spd;
            const auto current=sigma_representatives.find(key);
            if (current==sigma_representatives.end() ||
                representative_score>current->second.second) {
                sigma_representatives[key]={group,representative_score};
            }
        }
        for (const auto& [key,representative]:sigma_representatives) {
            (void)key;
            essential[representative.first]=true;
        }

        // One recovered row for every copy in the formal d decomposition is a
        // complete minimal ligand-field manifold.  Repeated irreps (for
        // example 2A1 in C2v or 2E in C3v) are counted explicitly.
        const bool complete_d_manifold=supported_ligand_field &&
            std::all_of(d_irrep_multiplicity.begin(),
                d_irrep_multiplicity.end(),
                [&](const auto& expected) {
                    const auto found=candidates.find(expected.first);
                    return found!=candidates.end() &&
                           found->second.size()>=expected.second;
                });
        if (complete_d_manifold) {
            for (std::size_t group=0;group<groups.size();++group) {
                if (!groups[group].include) continue;
                if (groups[group].level.metal_d_weight>=0.15) {
                    essential[group]=true;
                }
                if (!essential[group]) {
                    groups[group].include=false;
                    ++hidden_intermediate_count;
                }
            }
        }
    }
    std::size_t included_count=static_cast<std::size_t>(std::count_if(
        groups.begin(),groups.end(),[](const auto& group) {
            return group.include;
        }));
    if (included_count>row_budget) {
        const double frontier_energy=data.frontier.homo && data.frontier.lumo
            ?0.5*(wavefunction.orbitals[*data.frontier.homo].energy_hartree+
                  wavefunction.orbitals[*data.frontier.lumo].energy_hartree)
            :0.0;
        struct RankedGroup { std::size_t index=0; double score=0.0; };
        std::vector<RankedGroup> ranked;
        std::size_t essential_count=0u;
        for (std::size_t group=0;group<groups.size();++group) {
            if (!groups[group].include) continue;
            if (essential[group]) {
                ++essential_count;
                continue;
            }
            const auto& level=groups[group].level;
            const double metal=level.metal_s_weight+level.metal_p_weight+
                               level.metal_d_weight;
            const double score=4.0*metal+1.5*level.ligand_p_weight+
                0.35*std::max(level.sigma_fraction,level.pi_fraction)+
                1.5*std::min(0.20,std::abs(level.metal_ligand_overlap))-
                0.20*std::abs(level.layout_energy_hartree-frontier_energy)+
                (level.total_occupation>
                     options.filter.occupation_threshold?0.35:0.0);
            ranked.push_back({group,score});
            groups[group].include=false;
        }
        std::sort(ranked.begin(),ranked.end(),[](const auto& a,const auto& b) {
            if (a.score!=b.score) return a.score>b.score;
            return a.index<b.index;
        });
        const std::size_t available=row_budget>essential_count
            ?row_budget-essential_count:0u;
        for (std::size_t i=0;i<std::min(available,ranked.size());++i) {
            groups[ranked[i].index].include=true;
        }
    }

    std::vector<std::size_t> group_to_level(
        groups.size(),std::numeric_limits<std::size_t>::max());
    std::vector<double> axis_energies;
    for (std::size_t group=0;group<groups.size();++group) {
        if (!groups[group].include) continue;
        group_to_level[group]=data.levels.size();
        axis_energies.push_back(groups[group].level.layout_energy_hartree);
        data.levels.push_back(groups[group].level);
    }

    for (const auto& pair:raw_pairs) {
        const std::size_t lower=group_to_level[pair.lower];
        const std::size_t upper=group_to_level[pair.upper];
        const std::size_t retained=group_to_level[pair.retained];
        if (retained==std::numeric_limits<std::size_t>::max()) continue;
        PiInteractionDescriptor descriptor;
        descriptor.lower_level=lower==std::numeric_limits<std::size_t>::max()
            ?retained:lower;
        descriptor.upper_level=upper==std::numeric_limits<std::size_t>::max()
            ?retained:upper;
        descriptor.lower_orbitals=groups[pair.lower].level.member_indices;
        descriptor.upper_orbitals=groups[pair.upper].level.member_indices;
        descriptor.symmetry=pair.symmetry.empty()
            ?groups[pair.lower].level.metadata.symmetry:pair.symmetry;
        descriptor.kind=pair.kind;
        descriptor.splitting_hartree=pair.split;
        descriptor.confidence=pair.confidence;
        descriptor.lower_visible=lower!=std::numeric_limits<std::size_t>::max();
        descriptor.upper_visible=upper!=std::numeric_limits<std::size_t>::max();
        descriptor.retained_level=retained;
        data.pi_interactions.push_back(std::move(descriptor));
    }

    std::size_t raw_groups=0u;
    for (const auto& level:data.levels) {
        if (level.raw_data_fallback) ++raw_groups;
    }
    data.spin_counterpart_unmatched_visible=static_cast<std::size_t>(
        std::count_if(groups.begin(),groups.end(),[&](const auto& group) {
            return group.include && !group.suppressed_spin_counterpart &&
                   group_spin(wavefunction,group)==Spin::Beta;
        }));
    data.spin_counterparts_collapsed=
        data.spin_counterpart_pair_count>0u &&
        data.spin_counterpart_unmatched_visible==0u;
    data.spin_counterparts_partial=
        data.spin_counterpart_pair_count>0u &&
        data.spin_counterpart_unmatched_visible>0u;
    std::ostringstream summary;
    summary<<"ligand-field valence groups: "<<data.levels.size()
           <<"; local field="<<(data.ligand_field_point_group.empty()
                 ?"unresolved":data.ligand_field_point_group)
           <<"; geometry="<<(data.ligand_field_geometry_id.empty()
                ?"unresolved":data.ligand_field_geometry_id)
           <<"; CN="<<data.ligand_field_coordination_number
           <<"; spin counterparts="<<(data.spin_counterparts_collapsed
                ?"collapsed":(data.spin_counterparts_partial
                    ?"partial":"separate"))
           <<" (pairs="<<data.spin_counterpart_pair_count
           <<", unmatched-visible="<<data.spin_counterpart_unmatched_visible
           <<')'
           <<"; pi pairs="<<data.pi_interactions.size()
           <<"; raw-MO recovered groups="<<raw_groups
           <<"; intermediate groups hidden="<<hidden_intermediate_count
           <<"; degeneracies collapsed";
    data.selection.summary=summary.str();

    data.energy_transform=build_energy_transform(
        axis_energies,options.energy_axis_mode,
        options.nonlinear_minimum_gap_weight);
    return data;
}

} // namespace cov
