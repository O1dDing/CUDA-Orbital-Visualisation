#include "cov/mo_diagram.hpp"
#include "cov/ligand_field.hpp"
#include "cov/wavefunction_io.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

cov::Wavefunction synthetic_oh_environment() {
    cov::Wavefunction wavefunction;
    wavefunction.atoms.resize(7u);
    wavefunction.atoms[0].symbol="Cr";
    wavefunction.atoms[0].atomic_number=24;
    const std::array<std::array<double,3>,6> ligand_positions{{
        {{4.0,0.0,0.0}},{{-4.0,0.0,0.0}},
        {{0.0,4.0,0.0}},{{0.0,-4.0,0.0}},
        {{0.0,0.0,4.0}},{{0.0,0.0,-4.0}},
    }};
    for (std::size_t i=0;i<ligand_positions.size();++i) {
        auto& atom=wavefunction.atoms[i+1u];
        atom.symbol="F";
        atom.atomic_number=9;
        atom.x=ligand_positions[i][0];
        atom.y=ligand_positions[i][1];
        atom.z=ligand_positions[i][2];
        cov::BondOrderRecord bond;
        bond.atom_a=0u;
        bond.atom_b=static_cast<std::uint32_t>(i+1u);
        bond.mayer_order=0.20;
        wavefunction.bond_orders.push_back(bond);
    }
    return wavefunction;
}

cov::Wavefunction synthetic_td_environment() {
    cov::Wavefunction wavefunction;
    wavefunction.atoms.resize(5u);
    wavefunction.atoms[0].symbol="Cr";
    wavefunction.atoms[0].atomic_number=24;
    const std::array<std::array<double,3>,4> ligand_positions{{
        {{2.3,2.3,2.3}},{{2.3,-2.3,-2.3}},
        {{-2.3,2.3,-2.3}},{{-2.3,-2.3,2.3}},
    }};
    for (std::size_t i=0;i<ligand_positions.size();++i) {
        auto& atom=wavefunction.atoms[i+1u];
        atom.symbol="F";
        atom.atomic_number=9;
        atom.x=ligand_positions[i][0];
        atom.y=ligand_positions[i][1];
        atom.z=ligand_positions[i][2];
        cov::BondOrderRecord bond;
        bond.atom_a=0u;
        bond.atom_b=static_cast<std::uint32_t>(i+1u);
        bond.mayer_order=0.20;
        wavefunction.bond_orders.push_back(bond);
    }
    return wavefunction;
}

bool validate_radial_first_shell_retry() {
    auto wavefunction=synthetic_oh_environment();
    for (auto& bond:wavefunction.bond_orders) bond.mayer_order=0.025;
    const std::size_t near_count=wavefunction.atoms.size();
    wavefunction.atoms.resize(13u);
    for (std::size_t i=0;i<6u;++i) {
        const auto& near=wavefunction.atoms[i+1u];
        auto& remote=wavefunction.atoms[i+7u];
        remote.symbol="O";
        remote.atomic_number=8;
        remote.x=2.0*near.x;
        remote.y=2.0*near.y;
        remote.z=2.0*near.z;
        cov::BondOrderRecord bond;
        bond.atom_a=0u;
        bond.atom_b=static_cast<std::uint32_t>(i+7u);
        bond.mayer_order=0.08;
        wavefunction.bond_orders.push_back(bond);
    }
    const auto environment=cov::analyse_ligand_field_environment(wavefunction);
    if (near_count!=7u || environment.local_point_group()!="Oh" ||
        environment.ligand_atoms.size()!=6u ||
        std::any_of(environment.ligand_atoms.begin(),
                    environment.ligand_atoms.end(),[&](const auto atom) {
            return atom==0u || atom>=near_count ||
                   wavefunction.atoms[atom].atomic_number!=9;
        })) {
        std::cerr<<"radial-first low-Mayer Oh retry selected a remote atom\n";
        return false;
    }
    return true;
}

bool validate_local_symmetry_missing_markers() {
    auto make_block=[](const int angular_momentum,const double weight) {
        auto wavefunction=synthetic_oh_environment();
        for (std::size_t i=0;i<3u;++i) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=-0.10+0.001*static_cast<double>(i);
            orbital.occupation=1.0f;
            orbital.chemistry.available=true;
            cov::OrbitalAOContribution contribution;
            contribution.atom_index=0u;
            contribution.angular_momentum=angular_momentum;
            contribution.weight=weight;
            orbital.chemistry.ao_contributions.push_back(contribution);
            wavefunction.orbitals.push_back(std::move(orbital));
        }
        return wavefunction;
    };
    const std::array<std::string,3> missing{{"","?","N/A"}};
    auto metadata_for=[&](const auto& wavefunction) {
        std::vector<cov::OrbitalMetadata> metadata(wavefunction.orbitals.size());
        for (std::size_t i=0;i<metadata.size();++i) {
            metadata[i].energy_hartree=wavefunction.orbitals[i].energy_hartree;
            metadata[i].degeneracy_size=1u;
            metadata[i].symmetry=missing[i%missing.size()];
        }
        return metadata;
    };

    auto above=make_block(2,0.021);
    auto above_metadata=metadata_for(above);
    cov::apply_local_ligand_field_symmetry(above,above_metadata);
    if (!std::all_of(above_metadata.begin(),above_metadata.end(),[](const auto& item) {
            return item.symmetry=="T2g";
        })) {
        std::cerr<<"missing local-symmetry markers were not recovered\n";
        return false;
    }

    auto below=make_block(2,0.019);
    auto below_metadata=metadata_for(below);
    cov::apply_local_ligand_field_symmetry(below,below_metadata);
    if (std::any_of(below_metadata.begin(),below_metadata.end(),[](const auto& item) {
            return item.symmetry=="T2g";
        })) {
        std::cerr<<"sub-threshold metal-d noise was grouped\n";
        return false;
    }

    auto low_p=make_block(1,0.021);
    auto low_p_metadata=metadata_for(low_p);
    cov::apply_local_ligand_field_symmetry(low_p,low_p_metadata);
    if (std::any_of(low_p_metadata.begin(),low_p_metadata.end(),[](const auto& item) {
            return item.symmetry=="T1u";
        })) {
        std::cerr<<"low metal-p noise used the metal-d grouping floor\n";
        return false;
    }
    return true;
}

bool validate_resolved_five_d_runs() {
    const auto make_case=[](cov::Wavefunction wavefunction,
                            const std::array<double,5>& energies) {
        std::vector<cov::OrbitalMetadata> metadata(energies.size());
        for (std::size_t i=0;i<energies.size();++i) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=energies[i];
            orbital.occupation=1.0f;
            orbital.chemistry.available=true;
            cov::OrbitalAOContribution contribution;
            contribution.atom_index=0u;
            contribution.angular_momentum=2;
            contribution.weight=0.30;
            orbital.chemistry.ao_contributions.push_back(contribution);
            wavefunction.orbitals.push_back(std::move(orbital));
            metadata[i].energy_hartree=energies[i];
            metadata[i].degeneracy_size=1u;
            metadata[i].symmetry=i%2u==0u?"N/A":"?";
        }
        cov::apply_local_ligand_field_symmetry(wavefunction,metadata);
        return metadata;
    };
    const auto td=make_case(
        synthetic_td_environment(),{{-0.1000,-0.0998,-0.0968,-0.0966,-0.0964}});
    const bool td_ok=std::all_of(td.begin(),td.begin()+2u,[](const auto& item) {
            return item.symmetry=="E";
        }) && std::all_of(td.begin()+2u,td.end(),[](const auto& item) {
            return item.symmetry=="T2";
        });
    const auto oh=make_case(
        synthetic_oh_environment(),{{-0.1000,-0.0998,-0.0996,-0.0966,-0.0964}});
    const bool oh_ok=std::all_of(oh.begin(),oh.begin()+3u,[](const auto& item) {
            return item.symmetry=="T2g";
        }) && std::all_of(oh.begin()+3u,oh.end(),[](const auto& item) {
            return item.symmetry=="Eg";
        });
    if (!td_ok || !oh_ok) {
        std::cerr<<"resolved weak-field five-d run was greedily mispartitioned\n";
        return false;
    }
    return true;
}

bool validate_minimal_compact_ligand_field_framework() {
    const auto add_group=[](cov::Wavefunction& wavefunction,
                            const std::string& symmetry,
                            const std::size_t degeneracy,
                            const double energy,
                            const int metal_angular_momentum,
                            const double metal_weight,
                            const double ligand_p_weight,
                            const bool sigma) {
        const std::size_t ligand_count=wavefunction.atoms.size()-1u;
        for (std::size_t member=0;member<degeneracy;++member) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=energy;
            orbital.occupation=2.0f;
            orbital.symmetry=symmetry;
            orbital.chemistry.available=true;
            orbital.chemistry.valence_manifold=true;
            orbital.chemistry.valence_weight=1.0;
            orbital.chemistry.unresolved_weight=0.0;
            if (metal_weight>0.0) {
                cov::OrbitalAOContribution metal;
                metal.atom_index=0u;
                metal.angular_momentum=metal_angular_momentum;
                metal.weight=metal_weight;
                orbital.chemistry.ao_contributions.push_back(metal);
            }
            cov::OrbitalAOContribution ligand;
            ligand.atom_index=static_cast<std::uint32_t>(
                1u+member%ligand_count);
            ligand.angular_momentum=1;
            ligand.weight=ligand_p_weight;
            orbital.chemistry.ao_contributions.push_back(ligand);

            cov::OrbitalPairInteraction interaction;
            interaction.atom_a=0u;
            interaction.atom_b=ligand.atom_index;
            interaction.overlap_character=0.06;
            interaction.channel.sigma=sigma?1.0:0.0;
            interaction.channel.pi=sigma?0.0:1.0;
            interaction.channel.undetermined=0.0;
            orbital.chemistry.interactions.push_back(interaction);
            wavefunction.orbitals.push_back(std::move(orbital));
        }
    };
    const auto diagram_symmetries=[](const cov::Wavefunction& wavefunction,
                                     const bool compact) {
        cov::MODiagramOptions options;
        options.hide_ligand_centred_intermediates=compact;
        options.selected_index=wavefunction.orbitals.size();
        const auto data=cov::build_mo_diagram_data(wavefunction,options);
        std::set<std::string> result;
        for (const auto& level:data.levels) {
            result.insert(level.metadata.symmetry);
        }
        return std::pair{data.ligand_field_point_group,std::move(result)};
    };

    auto oh=synthetic_oh_environment();
    oh.point_group_detected="Oh";
    add_group(oh,"Eg",2u,-0.20,2,0.70,0.20,true);
    add_group(oh,"T2g",3u,-0.10,2,0.70,0.20,false);
    add_group(oh,"T1g",3u,-0.05,2,0.00,0.80,false);
    add_group(oh,"T1u",3u,0.05,1,0.25,0.30,true);
    const auto [oh_expanded_group,oh_expanded]=diagram_symmetries(oh,false);
    const auto [oh_group,oh_symmetries]=diagram_symmetries(oh,true);
    if (oh_expanded_group!="Oh" || oh_expanded.count("T1g")==0u ||
        oh_group!="Oh" || oh_symmetries.count("Eg")==0u ||
        oh_symmetries.count("T2g")==0u ||
        oh_symmetries.count("T1u")==0u ||
        oh_symmetries.count("T1g")!=0u) {
        std::cerr<<"minimal compact Oh framework pruning failed\n";
        return false;
    }

    auto td=synthetic_td_environment();
    td.point_group_detected="Td";
    add_group(td,"E",2u,-0.20,2,0.70,0.20,false);
    add_group(td,"T2",3u,-0.10,2,0.70,0.20,true);
    add_group(td,"T1",3u,-0.05,2,0.00,0.80,false);
    const auto [td_expanded_group,td_expanded]=diagram_symmetries(td,false);
    const auto [td_group,td_symmetries]=diagram_symmetries(td,true);
    if (td_expanded_group!="Td" || td_expanded.count("T1")==0u ||
        td_group!="Td" || td_symmetries.count("E")==0u ||
        td_symmetries.count("T2")==0u ||
        td_symmetries.count("T1")!=0u) {
        std::cerr<<"minimal compact Td framework pruning failed\n";
        return false;
    }
    return true;
}

cov::Wavefunction synthetic_oh_pi_pair(const int donor_z,
                                       const int external_z,
                                       const double external_mayer,
                                       const double lower_metal_d,
                                       const double lower_ligand_p,
                                       const double upper_metal_d,
                                       const double upper_ligand_p) {
    auto wavefunction=synthetic_oh_environment();
    wavefunction.point_group_detected="Oh";
    for (std::size_t i=1u;i<wavefunction.atoms.size();++i) {
        wavefunction.atoms[i].atomic_number=donor_z;
        wavefunction.atoms[i].symbol=donor_z==6?"C":(donor_z==7?"N":"X");
    }
    if (external_z>0) {
        wavefunction.atoms.resize(13u);
        for (std::size_t i=0;i<6u;++i) {
            const auto& direct=wavefunction.atoms[i+1u];
            auto& external=wavefunction.atoms[i+7u];
            external.atomic_number=external_z;
            external.symbol=external_z==6?"C":(external_z==7?"N":"X");
            external.x=1.7*direct.x;
            external.y=1.7*direct.y;
            external.z=1.7*direct.z;
            cov::BondOrderRecord bond;
            bond.atom_a=static_cast<std::uint32_t>(i+1u);
            bond.atom_b=static_cast<std::uint32_t>(i+7u);
            bond.mayer_order=external_mayer;
            wavefunction.bond_orders.push_back(bond);
        }
    }
    const auto add_group=[&](const double energy,
                             const float occupation,
                             const double metal_d,
                             const double ligand_p,
                             const double overlap) {
        for (std::size_t member=0;member<3u;++member) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=energy;
            orbital.occupation=occupation;
            orbital.symmetry="T2g";
            orbital.chemistry.available=true;
            orbital.chemistry.valence_manifold=true;
            orbital.chemistry.confidence=1.0;
            cov::OrbitalAOContribution metal;
            metal.atom_index=0u;
            metal.angular_momentum=2;
            metal.weight=metal_d;
            orbital.chemistry.ao_contributions.push_back(metal);
            cov::OrbitalAOContribution ligand;
            ligand.atom_index=static_cast<std::uint32_t>(member+1u);
            ligand.angular_momentum=1;
            ligand.weight=ligand_p;
            orbital.chemistry.ao_contributions.push_back(ligand);
            cov::OrbitalPairInteraction interaction;
            interaction.atom_a=0u;
            interaction.atom_b=static_cast<std::uint32_t>(member+1u);
            interaction.overlap_character=overlap;
            interaction.channel.pi=1.0;
            interaction.channel.undetermined=0.0;
            interaction.channel.dominant=cov::OrbitalAngularFamily::Pi;
            interaction.channel.status=cov::ChemistryStatus::Determined;
            orbital.chemistry.interactions.push_back(interaction);
            wavefunction.orbitals.push_back(std::move(orbital));
        }
    };
    add_group(-0.20,2.0f,lower_metal_d,lower_ligand_p,1.0);
    add_group(0.10,0.0f,upper_metal_d,upper_ligand_p,-1.0);
    return wavefunction;
}

bool validate_pi_topology_and_two_sided_composition() {
    cov::MODiagramOptions options;
    options.selected_index=0u;
    options.max_levels=0u;
    const auto one_sided=cov::build_mo_diagram_data(
        synthetic_oh_pi_pair(6,7,1.50,0.70,0.10,0.10,0.12),options);
    if (!one_sided.pi_interactions.empty()) {
        std::cerr<<"one-sided pi composition produced a strong pair\n";
        return false;
    }
    const auto alkyl=cov::build_mo_diagram_data(
        synthetic_oh_pi_pair(6,6,0.80,0.10,0.70,0.70,0.10),options);
    const bool alkyl_donor=alkyl.pi_interactions.size()==1u &&
        alkyl.pi_interactions.front().kind==cov::PiInteractionKind::Donor;
    if (!alkyl_donor) {
        std::cerr<<"single-bonded carbon donor was hard-coded as acceptor\n";
        return false;
    }
    const auto terminal_n=cov::build_mo_diagram_data(
        synthetic_oh_pi_pair(7,0,0.0,0.10,0.70,0.70,0.10),options);
    const bool terminal_n_donor=terminal_n.pi_interactions.size()==1u &&
        terminal_n.pi_interactions.front().kind==cov::PiInteractionKind::Donor;
    if (!terminal_n_donor) {
        std::cerr<<"low-coordinate terminal nitrogen was forced sigma-only\n";
        return false;
    }
    return true;
}

bool validate_equal_count_unrestricted_collapse() {
    auto unrestricted=synthetic_oh_environment();
    unrestricted.basis_count=3u;
    unrestricted.alpha_electrons=3u;
    unrestricted.beta_electrons=3u;
    unrestricted.ao_overlap={1.0,0.0,0.0,
                             0.0,1.0,0.0,
                             0.0,0.0,1.0};
    for (const auto spin:{cov::Spin::Alpha,cov::Spin::Beta}) {
        for (std::size_t i=0;i<3u;++i) {
            cov::MolecularOrbital orbital;
            orbital.energy_hartree=-0.10+0.001*static_cast<double>(i);
            orbital.occupation=1.0f;
            orbital.spin=spin;
            orbital.coefficients.assign(3u,0.0f);
            orbital.coefficients[i]=1.0f;
            orbital.chemistry.available=true;
            orbital.chemistry.valence_manifold=true;
            orbital.chemistry.confidence=1.0;
            cov::OrbitalAOContribution contribution;
            contribution.atom_index=0u;
            contribution.angular_momentum=2;
            contribution.weight=0.60;
            orbital.chemistry.ao_contributions.push_back(contribution);
            unrestricted.orbitals.push_back(std::move(orbital));
        }
    }
    cov::MODiagramOptions options;
    options.selected_index=1u;
    options.max_levels=0u;
    const auto unrestricted_data=cov::build_mo_diagram_data(
        unrestricted,options);
    if (!unrestricted_data.spin_counterparts_collapsed ||
        unrestricted_data.spin_counterparts_partial ||
        unrestricted_data.spin_counterpart_pair_count==0u ||
        unrestricted_data.spin_counterpart_unmatched_visible!=0u) {
        std::cerr<<"equal-count unrestricted spin sets were not collapsed\n";
        return false;
    }

    auto restricted=unrestricted;
    restricted.orbitals.resize(3u);
    for (auto& orbital:restricted.orbitals) orbital.occupation=2.0f;
    const auto restricted_data=cov::build_mo_diagram_data(restricted,options);
    if (restricted_data.spin_counterpart_pair_count!=0u ||
        restricted_data.spin_counterparts_collapsed ||
        restricted_data.spin_counterparts_partial) {
        std::cerr<<"alpha-only restricted orbitals were treated as UHF\n";
        return false;
    }
    return true;
}

bool same_compact_structure(const cov::MODiagramData& left,
                            const cov::MODiagramData& right) {
    if (left.ligand_field_point_group!=right.ligand_field_point_group ||
        left.levels.size()!=right.levels.size() ||
        left.pi_interactions.size()!=right.pi_interactions.size()) {
        return false;
    }
    for (std::size_t i=0;i<left.levels.size();++i) {
        const auto& a=left.levels[i];
        const auto& b=right.levels[i];
        if (a.member_indices!=b.member_indices ||
            a.member_spin_counterparts!=b.member_spin_counterparts ||
            a.metadata.symmetry!=b.metadata.symmetry ||
            a.approximate_nonbonding!=b.approximate_nonbonding ||
            std::abs(a.layout_energy_hartree-b.layout_energy_hartree)>1.0e-12 ||
            std::abs(a.energy_spread_hartree-b.energy_spread_hartree)>1.0e-12 ||
            std::abs(cov::energy_display_coordinate(
                a.layout_energy_hartree,left.energy_transform)-
                cov::energy_display_coordinate(
                    b.layout_energy_hartree,right.energy_transform))>1.0e-12) {
            return false;
        }
    }
    for (std::size_t i=0;i<left.pi_interactions.size();++i) {
        const auto& a=left.pi_interactions[i];
        const auto& b=right.pi_interactions[i];
        if (a.lower_level!=b.lower_level || a.upper_level!=b.upper_level ||
            a.lower_orbitals!=b.lower_orbitals ||
            a.upper_orbitals!=b.upper_orbitals ||
            a.symmetry!=b.symmetry || a.kind!=b.kind ||
            std::abs(a.splitting_hartree-b.splitting_hartree)>1.0e-12 ||
            a.lower_visible!=b.lower_visible ||
            a.upper_visible!=b.upper_visible ||
            a.retained_level!=b.retained_level) {
            return false;
        }
    }
    return true;
}

bool inspect_real_fchk(const std::filesystem::path& path) {
    const auto wf=cov::parse_wavefunction(path);
    cov::MODiagramOptions options;
    options.selected_index=0u;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        if (wf.orbitals[i].occupation>1.0e-4f) options.selected_index=i;
    }
    options.max_levels=0u;
    const auto data=cov::build_mo_diagram_data(wf,options);
    const auto browser_metadata=cov::build_orbital_metadata(
        wf,options.selected_index,options.degeneracy,options.filter);
    std::cout<<"DIAGNOSTIC "<<path.filename().string()
             <<" point-group="<<wf.point_group_detected
             <<" alpha="<<wf.alpha_electrons
             <<" beta="<<wf.beta_electrons
             <<" orbitals="<<wf.orbitals.size()
             <<" rows="<<data.levels.size()
             <<" selected="<<data.selection.included_indices.size()
             <<" local-field="<<data.ligand_field_point_group
             <<" first-shell="<<data.ligand_field_ligand_atoms.size()
             <<" spin-collapsed="<<data.spin_counterparts_collapsed
             <<" spin-partial="<<data.spin_counterparts_partial
             <<" spin-pairs="<<data.spin_counterpart_pair_count
             <<" spin-unmatched="<<data.spin_counterpart_unmatched_visible
             <<" pi-pairs="<<data.pi_interactions.size()<<'\n';
    std::cout<<"  diagnostic-first-shell metal="<<data.ligand_field_metal_atom;
    for (const auto atom:data.ligand_field_ligand_atoms) {
        if (atom<wf.atoms.size()) {
            std::cout<<' '<<atom<<':'<<wf.atoms[atom].symbol;
        }
    }
    std::cout<<'\n';
    for (const auto& bond:wf.bond_orders) {
        if (bond.atom_a>=wf.atoms.size() || bond.atom_b>=wf.atoms.size()) continue;
        const int za=wf.atoms[bond.atom_a].atomic_number;
        const int zb=wf.atoms[bond.atom_b].atomic_number;
        const auto tm=[](const int z) {
            return (z>=21 && z<=30) || (z>=39 && z<=48) ||
                   (z>=72 && z<=80) || (z>=104 && z<=112);
        };
        if (!tm(za) && !tm(zb)) continue;
        std::cout<<"  diagnostic-bond "<<bond.atom_a<<'-'<<bond.atom_b
                 <<' '<<wf.atoms[bond.atom_a].symbol<<'-'
                 <<wf.atoms[bond.atom_b].symbol
                 <<" Mayer="<<bond.mayer_order<<'\n';
    }
    for (std::size_t i=0;i<data.levels.size();++i) {
        const auto& level=data.levels[i];
        std::cout<<"  diagnostic-row "<<i
                 <<" MO"<<level.metadata.raw_mo_number
                 <<" spin="<<static_cast<int>(
                        wf.orbitals[level.member_indices.front()].spin)
                 <<" sym="<<level.metadata.symmetry
                 <<" deg="<<level.member_indices.size()
                 <<" E="<<level.layout_energy_hartree
                 <<" occ="<<level.total_occupation
                 <<" s="<<level.metal_s_weight
                 <<" p="<<level.metal_p_weight
                 <<" d="<<level.metal_d_weight
                 <<" direct-Lp="<<level.direct_ligand_p_weight
                 <<" Lp="<<level.ligand_p_weight
                 <<" sigma="<<level.sigma_fraction
                 <<" pi="<<level.pi_fraction
                 <<" ML="<<level.metal_ligand_overlap
                 <<" raw="<<level.raw_data_fallback<<'\n';
    }
    if (data.levels.empty() || data.selection.included_indices.empty()) {
        std::cerr<<path.filename().string()<<": empty ligand-field diagram\n";
        return false;
    }
    for (const auto& level:data.levels) {
        const double allowed_spread=data.ligand_field_point_group.empty()
            ?options.degeneracy.tolerance_hartree:0.01501;
        if (level.member_indices.empty() ||
            level.member_indices.size()!=level.metadata.degeneracy_size ||
            level.energy_spread_hartree>allowed_spread) {
            std::cerr<<path.filename().string()<<": broken degenerate row\n";
            return false;
        }
    }
    const std::size_t row_budget=std::clamp<std::size_t>(
        2u*options.neighbourhood+1u,10u,48u);
    if (data.levels.size()>row_budget ||
        data.levels.size()>=data.selection.included_indices.size()) {
        std::cerr<<path.filename().string()
                 <<": valence reduction did not reduce the raw MO set\n";
        return false;
    }
    // Producer-derived ligand-field manifolds can be complete without using
    // the generic raw-MO supplement.  Raw fallback is therefore an optional
    // recovery path, not a correctness requirement for every fixture.
    auto compact_options=options;
    compact_options.hide_ligand_centred_intermediates=true;
    const auto compact=cov::build_mo_diagram_data(wf,compact_options);
    const std::string filename=path.filename().string();
    for (std::size_t i=0;i<compact.levels.size();++i) {
        const auto& level=compact.levels[i];
        std::cout<<"  compact-row "<<i
                 <<" MO"<<level.metadata.raw_mo_number
                 <<" sym="<<level.metadata.symmetry
                 <<" deg="<<level.member_indices.size()
                 <<" E="<<level.layout_energy_hartree
                 <<" coord="<<cov::energy_display_coordinate(
                        level.layout_energy_hartree,compact.energy_transform)
                 <<" occ="<<level.total_occupation
                 <<" d="<<level.metal_d_weight
                 <<" direct-Lp="<<level.direct_ligand_p_weight
                 <<" sigma="<<level.sigma_fraction
                 <<" ML="<<level.metal_ligand_overlap<<'\n';
    }
    if (compact.levels.size()>=data.levels.size() ||
        compact.pi_interactions.size()!=data.pi_interactions.size()) {
        std::cerr<<path.filename().string()
                 <<": intermediate-orbital toggle did not reduce safely ("
                 <<data.levels.size()<<" -> "<<compact.levels.size()
                 <<" rows; "<<data.pi_interactions.size()<<" -> "
                 <<compact.pi_interactions.size()<<" pairs)\n";
        return false;
    }
    const bool selection_invariance_fixture=
        filename.find("TiF6")!=std::string::npos ||
        filename.find("ZnCl4")!=std::string::npos ||
        filename.find("19_FeF6")!=std::string::npos;
    if (selection_invariance_fixture) {
        auto canonical_options=compact_options;
        canonical_options.selected_index=wf.orbitals.size();
        const auto canonical=cov::build_mo_diagram_data(
            wf,canonical_options);
        if (!same_compact_structure(compact,canonical)) {
            std::cerr<<filename
                     <<": compact row set changed with the initial selection\n";
            return false;
        }
        const bool octahedral=canonical.ligand_field_point_group=="Oh";
        const std::string hidden_symmetry=octahedral?"T1g":"T1";
        const std::set<std::string> probe_symmetries=octahedral
            ?std::set<std::string>{"T1g","T1u","T2g"}
            :std::set<std::string>{"T1","T2","E"};
        std::vector<std::size_t> probes;
        bool source_contains_hidden=false;
        for (std::size_t i=0;i<data.metadata.size();++i) {
            const auto& symmetry=data.metadata[i].symmetry;
            if (symmetry==hidden_symmetry) source_contains_hidden=true;
            if (probe_symmetries.count(symmetry)) probes.push_back(i);
        }
        std::sort(probes.begin(),probes.end());
        probes.erase(std::unique(probes.begin(),probes.end()),probes.end());
        for (const auto selected:probes) {
            auto selected_options=canonical_options;
            selected_options.selected_index=selected;
            const auto selected_data=cov::build_mo_diagram_data(
                wf,selected_options);
            if (!same_compact_structure(canonical,selected_data)) {
                std::cerr<<filename
                         <<": compact row set changed after selecting MO"
                         <<selected+1u<<'\n';
                return false;
            }
        }
        const auto compact_has=[&](const std::string& symmetry) {
            return std::any_of(
                canonical.levels.begin(),canonical.levels.end(),
                [&](const auto& level) {
                    return level.metadata.symmetry==symmetry;
                });
        };
        const bool required_framework=octahedral
            ?compact_has("T1u") && compact_has("T2g")
            :compact_has("T2");
        if (!source_contains_hidden || compact_has(hidden_symmetry) ||
            !required_framework) {
            std::cerr<<filename
                     <<": compact T1/T1g exclusion damaged the metal framework\n";
            return false;
        }
    }
    const auto valid_symmetry=[](const std::string& symmetry) {
        return !symmetry.empty() && symmetry!="?" &&
               symmetry!="N/A" && symmetry!="n/a";
    };
    const auto level_metrics_are_bounded=[](const auto& levels) {
        return std::all_of(levels.begin(),levels.end(),[](const auto& level) {
            return std::isfinite(level.metal_ligand_overlap) &&
                   std::abs(level.metal_ligand_overlap)<=1.000001;
        });
    };
    if (!level_metrics_are_bounded(data.levels) ||
        !level_metrics_are_bounded(compact.levels)) {
        std::cerr<<path.filename().string()
                 <<": first-shell metal-ligand metric is not finite/bounded\n";
        return false;
    }
    if (!data.ligand_field_point_group.empty()) {
        const bool bad_shell=data.ligand_field_metal_atom>=wf.atoms.size() ||
            data.ligand_field_ligand_atoms.empty() ||
            std::any_of(data.ligand_field_ligand_atoms.begin(),
                        data.ligand_field_ligand_atoms.end(),[&](const auto atom) {
                return atom>=wf.atoms.size() ||
                       atom==data.ligand_field_metal_atom ||
                       std::count(data.ligand_field_ligand_atoms.begin(),
                                  data.ligand_field_ligand_atoms.end(),atom)!=1;
            });
        const bool bad_compact_label=std::any_of(
            compact.levels.begin(),compact.levels.end(),[&](const auto& level) {
                if (!valid_symmetry(level.metadata.symmetry)) return true;
                return std::any_of(
                    level.member_indices.begin(),level.member_indices.end(),
                    [&](const auto member) {
                        return member>=wf.orbitals.size() ||
                               member>=data.metadata.size() ||
                               member>=browser_metadata.size() ||
                               data.metadata[member].symmetry!=
                                   level.metadata.symmetry ||
                               browser_metadata[member].symmetry!=
                                   level.metadata.symmetry ||
                               std::count(level.member_indices.begin(),
                                          level.member_indices.end(),member)!=1;
                    });
            });
        if (bad_shell || bad_compact_label) {
            for (const auto& level:compact.levels) {
                for (const auto member:level.member_indices) {
                    if (member>=data.metadata.size() ||
                        member>=browser_metadata.size() ||
                        data.metadata[member].symmetry!=level.metadata.symmetry ||
                        browser_metadata[member].symmetry!=level.metadata.symmetry) {
                        std::cerr<<"  member "<<member
                                 <<" level="<<level.metadata.symmetry
                                 <<" diagram="<<(member<data.metadata.size()
                                     ?data.metadata[member].symmetry:"<out>")
                                 <<" browser="<<(member<browser_metadata.size()
                                     ?browser_metadata[member].symmetry:"<out>")
                                 <<'\n';
                    }
                }
            }
            std::cerr<<path.filename().string()
                     <<": local-field shell/compact member metadata regression\n";
            return false;
        }
    }
    for (const auto& interaction:data.pi_interactions) {
        const auto validate_members=[&](const auto& members) {
            return std::all_of(members.begin(),members.end(),[&](const auto member) {
                if (member>=data.metadata.size()) return false;
                return valid_symmetry(data.metadata[member].symmetry);
            });
        };
        if (!validate_members(interaction.lower_orbitals) ||
            !validate_members(interaction.upper_orbitals)) {
            std::cerr<<path.filename().string()
                     <<": split-orbital member symmetry fell back to N/A\n";
            return false;
        }
    }

    const auto expected_pair=[&](const cov::PiInteractionKind kind,
                                 const std::string& symmetry) {
        return std::find_if(
            data.pi_interactions.begin(),data.pi_interactions.end(),
            [&](const auto& pair) {
                return pair.kind==kind && pair.symmetry==symmetry;
            });
    };
    const auto all_members_have=[&](const auto& members,
                                     const std::string& symmetry) {
        return std::all_of(members.begin(),members.end(),[&](const auto member) {
            return member<data.metadata.size() &&
                   data.metadata[member].symmetry==symmetry;
        });
    };
    const auto compact_group=[&](const std::string& symmetry,
                                 const std::size_t degeneracy) {
        return std::find_if(
            compact.levels.begin(),compact.levels.end(),[&](const auto& level) {
                return level.metadata.symmetry==symmetry &&
                       level.member_indices.size()==degeneracy;
            });
    };
    const auto expected_shell=[&](const std::string& point_group,
                                  const std::size_t ligand_count,
                                  const int metal_z,
                                  const int ligand_z) {
        return data.ligand_field_point_group==point_group &&
            data.ligand_field_metal_atom<wf.atoms.size() &&
            wf.atoms[data.ligand_field_metal_atom].atomic_number==metal_z &&
            data.ligand_field_ligand_atoms.size()==ligand_count &&
            std::all_of(data.ligand_field_ligand_atoms.begin(),
                        data.ligand_field_ligand_atoms.end(),[&](const auto atom) {
                return atom<wf.atoms.size() &&
                       wf.atoms[atom].atomic_number==ligand_z;
            });
    };
    const auto fully_collapsed_open_shell=[&]() {
        return data.spin_counterparts_collapsed &&
               !data.spin_counterparts_partial &&
               data.spin_counterpart_pair_count>0u &&
               data.spin_counterpart_unmatched_visible==0u;
    };
    if (filename.find("TiF6")!=std::string::npos) {
        const auto pair=expected_pair(cov::PiInteractionKind::Donor,"T2g");
        if (!expected_shell("Oh",6u,22,9) ||
            data.pi_interactions.size()!=1u || pair==data.pi_interactions.end() ||
            !pair->lower_visible || !pair->upper_visible ||
            pair->splitting_hartree<0.35 || pair->splitting_hartree>0.40 ||
            !all_members_have(pair->lower_orbitals,"T2g") ||
            !all_members_have(pair->upper_orbitals,"T2g")) {
            std::cerr<<filename<<": expected resolved T2g pi-donor pair\n";
            return false;
        }
    } else if (filename.find("ZnCl4")!=std::string::npos) {
        const auto pair=expected_pair(
            cov::PiInteractionKind::WeakNearNonbonding,"E/T2");
        if (!expected_shell("Td",4u,30,17) ||
            data.pi_interactions.size()!=1u || pair==data.pi_interactions.end() ||
            pair->lower_visible==pair->upper_visible ||
            pair->splitting_hartree>options.weak_pi_split_hartree ||
            pair->retained_level>=data.levels.size() ||
            !data.levels[pair->retained_level].approximate_nonbonding ||
            data.levels[pair->retained_level].metadata.symmetry!="E" ||
            !all_members_have(pair->lower_orbitals,"T2") ||
            !all_members_have(pair->upper_orbitals,"E")) {
            std::cerr<<filename
                     <<": expected reduced approximately-nonbonding E/T2 split\n";
            return false;
        }
    } else if (filename.find("CrCO6")!=std::string::npos) {
        const auto pair=expected_pair(cov::PiInteractionKind::Acceptor,"T2g");
        const bool compact_pair=std::any_of(
            compact.pi_interactions.begin(),compact.pi_interactions.end(),
            [](const auto& item) {
                return item.kind==cov::PiInteractionKind::Acceptor &&
                       item.symmetry=="T2g";
            });
        const bool sigma_bonding=std::any_of(
            compact.levels.begin(),compact.levels.end(),[](const auto& level) {
                return level.sigma_fraction>=0.55 &&
                       level.metal_s_weight+level.metal_p_weight>=0.03 &&
                       level.metal_ligand_overlap>0.0;
            });
        const bool sigma_antibonding=std::any_of(
            compact.levels.begin(),compact.levels.end(),[](const auto& level) {
                return level.sigma_fraction>=0.55 &&
                       level.metal_s_weight+level.metal_p_weight>=0.03 &&
                       level.metal_ligand_overlap<0.0;
            });
        if (!expected_shell("Oh",6u,24,6) ||
            data.pi_interactions.size()!=1u || pair==data.pi_interactions.end() ||
            !pair->lower_visible || !pair->upper_visible ||
            pair->splitting_hartree<0.24 || pair->splitting_hartree>0.28 ||
            compact.levels.size()>10u || !compact_pair ||
            !sigma_bonding || !sigma_antibonding ||
            !all_members_have(pair->lower_orbitals,"T2g") ||
            !all_members_have(pair->upper_orbitals,"T2g")) {
            std::cerr<<filename<<": expected resolved T2g pi-acceptor pair\n";
            return false;
        }
    } else if (filename.find("CrNH3_6")!=std::string::npos) {
        const auto find_compact=[&](const std::string& symmetry,
                                    const std::size_t degeneracy,
                                    const bool antibonding) {
            return std::find_if(
                compact.levels.begin(),compact.levels.end(),
                [&](const auto& level) {
                    return level.metadata.symmetry==symmetry &&
                           level.member_indices.size()==degeneracy &&
                           (level.metal_ligand_overlap<0.0)==antibonding;
                });
        };
        const auto t2g=std::find_if(
            compact.levels.begin(),compact.levels.end(),[](const auto& level) {
                return level.metadata.symmetry=="T2g" &&
                       level.member_indices.size()==3u &&
                       level.total_occupation>2.9 &&
                       level.total_occupation<3.1;
            });
        const bool t2g_arrows=t2g!=compact.levels.end() &&
            t2g->member_electrons.size()==3u &&
            std::all_of(t2g->member_electrons.begin(),
                        t2g->member_electrons.end(),[](const auto& electrons) {
                return electrons.alpha==1 && electrons.beta==0;
            }) && t2g->member_spin_counterparts.size()==3u &&
            std::all_of(t2g->member_spin_counterparts.begin(),
                        t2g->member_spin_counterparts.end(),[&](const auto index) {
                return index<wf.orbitals.size() &&
                       wf.orbitals[index].spin==cov::Spin::Beta &&
                       index<data.metadata.size() &&
                       data.metadata[index].symmetry=="T2g";
            });
        const auto labelled_count=[&](const std::string& symmetry) {
            return static_cast<std::size_t>(std::count_if(
                data.metadata.begin(),data.metadata.end(),
                [&](const auto& item) {return item.symmetry==symmetry;}));
        };
        const auto browser_labelled_count=[&](const std::string& symmetry) {
            return static_cast<std::size_t>(std::count_if(
                browser_metadata.begin(),browser_metadata.end(),
                [&](const auto& item) {return item.symmetry==symmetry;}));
        };
        if (wf.point_group_detected=="Oh" ||
            !expected_shell("Oh",6u,24,7) ||
            data.ligand_field_confidence<0.90 ||
            !data.spin_counterparts_collapsed ||
            !data.pi_interactions.empty() || compact.levels.size()>8u ||
            find_compact("A1g",1u,false)==compact.levels.end() ||
            find_compact("A1g",1u,true)==compact.levels.end() ||
            find_compact("T1u",3u,false)==compact.levels.end() ||
            find_compact("T1u",3u,true)==compact.levels.end() ||
            find_compact("Eg",2u,false)==compact.levels.end() ||
            find_compact("Eg",2u,true)==compact.levels.end() ||
            !t2g_arrows || labelled_count("T2g")<6u ||
            labelled_count("Eg")<8u ||
            browser_labelled_count("T2g")<6u ||
            browser_labelled_count("Eg")<8u) {
            std::cerr<<filename
                     <<": local-Oh open-shell browser/diagram metadata regression\n";
            return false;
        }
    } else if (filename.find("19_FeF6")!=std::string::npos) {
        if (!expected_shell("Oh",6u,26,9) ||
            !fully_collapsed_open_shell() ||
            !data.pi_interactions.empty() || compact.levels.size()>10u ||
            compact_group("Eg",2u)==compact.levels.end() ||
            compact_group("T2g",3u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": weak-field open-shell Oh reduction regression\n";
            return false;
        }
    } else if (filename.find("20_CoCl4")!=std::string::npos) {
        if (!expected_shell("Td",4u,27,17) ||
            !fully_collapsed_open_shell() ||
            !data.pi_interactions.empty() || compact.levels.size()>8u ||
            compact_group("E",2u)==compact.levels.end() ||
            compact_group("T2",3u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": weak-field open-shell Td reduction regression\n";
            return false;
        }
    } else if (filename.find("21_CoNH3_6")!=std::string::npos) {
        if (wf.point_group_detected=="Oh" ||
            !expected_shell("Oh",6u,27,7) ||
            !data.pi_interactions.empty() || compact.levels.size()>8u ||
            compact_group("A1g",1u)==compact.levels.end() ||
            compact_group("T1u",3u)==compact.levels.end() ||
            compact_group("Eg",2u)==compact.levels.end() ||
            compact_group("T2g",3u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": sigma-only local-Oh framework regression\n";
            return false;
        }
    } else if (filename.find("22_ZnNH3_4")!=std::string::npos) {
        if (wf.point_group_detected=="Td" ||
            !expected_shell("Td",4u,30,7) ||
            !data.pi_interactions.empty() || compact.levels.size()>9u ||
            compact_group("A1",1u)==compact.levels.end() ||
            compact_group("E",2u)==compact.levels.end() ||
            compact_group("T2",3u)==compact.levels.end()) {
            std::cerr<<filename
                     <<": sigma-only local-Td framework regression\n";
            return false;
        }
    } else if (filename.find("23_NiCO4")!=std::string::npos) {
        const auto pair=expected_pair(cov::PiInteractionKind::Acceptor,"E");
        if (!expected_shell("Td",4u,28,6) ||
            data.pi_interactions.size()!=1u || pair==data.pi_interactions.end() ||
            !pair->lower_visible || !pair->upper_visible ||
            pair->splitting_hartree<=options.weak_pi_split_hartree ||
            !all_members_have(pair->lower_orbitals,"E") ||
            !all_members_have(pair->upper_orbitals,"E") ||
            compact.levels.size()>8u) {
            std::cerr<<filename
                     <<": strong-field Td pi-acceptor regression\n";
            return false;
        }
    } else if (filename.find("24_CrCN6")!=std::string::npos) {
        const auto pair=expected_pair(cov::PiInteractionKind::Acceptor,"T2g");
        if (!expected_shell("Oh",6u,24,6) ||
            !fully_collapsed_open_shell() ||
            data.pi_interactions.size()!=1u || pair==data.pi_interactions.end() ||
            !pair->lower_visible || !pair->upper_visible ||
            pair->splitting_hartree<=0.10 ||
            !all_members_have(pair->lower_orbitals,"T2g") ||
            !all_members_have(pair->upper_orbitals,"T2g") ||
            compact_group("Eg",2u)==compact.levels.end() ||
            compact_group("T2g",3u)==compact.levels.end() ||
            compact.levels.size()>11u) {
            std::cerr<<filename
                     <<": strong-field open-shell Oh pi-acceptor regression\n";
            return false;
        }
    }

    std::cout<<"REAL "<<path.filename().string()
             <<" point-group="<<wf.point_group_detected
             <<" rows="<<data.levels.size()
             <<" compact-rows="<<compact.levels.size()
             <<" raw-MOs="<<data.selection.included_indices.size()
             <<" pi-pairs="<<data.pi_interactions.size()<<'\n';
    for (std::size_t i=0;i<data.levels.size();++i) {
        const auto& level=data.levels[i];
        std::cout<<"  row "<<i<<" MO"<<level.metadata.raw_mo_number
                 <<" sym="<<level.metadata.symmetry
                 <<" deg="<<level.member_indices.size()
                 <<" E="<<level.layout_energy_hartree
                 <<" d="<<level.metal_d_weight
                 <<" Lp="<<level.ligand_p_weight
                 <<" sigma="<<level.sigma_fraction
                 <<" pi="<<level.pi_fraction
                 <<" ML="<<level.metal_ligand_overlap
                 <<" raw="<<level.raw_data_fallback
                 <<" approxNB="<<level.approximate_nonbonding<<'\n';
    }
    for (const auto& pair:data.pi_interactions) {
        std::cout<<"  PI "<<cov::pi_interaction_kind_name(pair.kind)
                 <<" sym="<<pair.symmetry
                 <<" split="<<pair.splitting_hartree
                 <<" confidence="<<pair.confidence
                 <<" lower-visible="<<pair.lower_visible
                 <<" upper-visible="<<pair.upper_visible<<'\n';
    }
    return true;
}

} // namespace

int main(int argc,char** argv) {
    if (!validate_radial_first_shell_retry()) return 14;
    if (!validate_local_symmetry_missing_markers()) return 12;
    if (!validate_resolved_five_d_runs()) return 15;
    if (!validate_minimal_compact_ligand_field_framework()) return 17;
    if (!validate_pi_topology_and_two_sided_composition()) return 16;
    if (!validate_equal_count_unrestricted_collapse()) return 13;

    cov::Wavefunction wf;
    wf.atoms.resize(4);

    cov::MolecularOrbital core;
    core.energy_hartree = -3.0;
    core.occupation = 2.0f;
    core.symmetry = "core";
    wf.orbitals.push_back(core);

    cov::MolecularOrbital occupied;
    occupied.energy_hartree = -0.5;
    occupied.occupation = 2.0f;
    occupied.symmetry = "sigma_g bonding";
    wf.orbitals.push_back(occupied);

    cov::MolecularOrbital virtual_a;
    virtual_a.energy_hartree = 0.2;
    virtual_a.occupation = 0.0f;
    virtual_a.symmetry = "pi_u";
    wf.orbitals.push_back(virtual_a);

    cov::MolecularOrbital virtual_b = virtual_a;
    virtual_b.energy_hartree = 0.200000004;
    wf.orbitals.push_back(virtual_b);

    cov::MolecularOrbital producer_multicentre;
    producer_multicentre.energy_hartree = 0.45;
    producer_multicentre.occupation = 0.0f;
    producer_multicentre.symmetry = "sigma 3c2e antibonding";
    wf.orbitals.push_back(producer_multicentre);

    cov::MolecularOrbital very_high;
    very_high.energy_hartree = 7.0;
    very_high.occupation = 0.0f;
    very_high.symmetry = "?";
    wf.orbitals.push_back(very_high);

    cov::MODiagramOptions options;
    options.selected_index = 1;
    options.energy_unit = cov::EnergyUnit::ElectronVolt;
    options.width = 760;
    options.height = 720;
    options.filter.core_energy_cutoff_hartree = -1.5;
    options.filter.virtual_window_hartree = 1.5;
    options.degeneracy.tolerance_hartree = 1.0e-5;
    // Keep all three low-lying virtuals so the explicit producer 3c2e fixture
    // participates in the human and machine export checks.
    options.max_levels = 7;

    const auto data = cov::build_mo_diagram_data(wf, options);
    if (!data.frontier.homo || !data.frontier.lumo) return 1;
    if (data.mode != cov::MODiagramMode::ValenceCentral || data.levels.size() != 4) {
        std::cerr << "valence-central selection failed\n";
        return 2;
    }
    if (data.levels[1].metadata.display_label != "3-a" ||
        data.levels[2].metadata.display_label != "3-b") {
        std::cerr << "degenerate labels failed\n";
        return 3;
    }
    if (data.levels[0].electrons.alpha != 1 || data.levels[0].electrons.beta != 1) {
        std::cerr << "electron population failed\n";
        return 4;
    }
    if (data.levels[0].annotation.family != "sigma" ||
        data.levels[0].annotation.bonding_class != cov::BondingClass::Bonding ||
        data.levels[1].annotation.family != "pi") {
        std::cerr << "explicit family/bonding annotation failed\n";
        return 5;
    }
    if (!data.levels[3].annotation.multicentre.available ||
        data.levels[3].annotation.multicentre.label != "3c2e" ||
        data.levels[3].annotation.multicentre.source != cov::AnnotationSource::ParsedLabel) {
        std::cerr << "explicit 3c2e producer annotation failed\n";
        return 6;
    }

    const auto temp = std::filesystem::temp_directory_path() / "cov_mo_diagram_smoke";
    const auto result = cov::export_mo_diagram_bundle(wf, options, temp);
    if (!result.svg || !result.png || !result.json || !result.csv) {
        std::cerr << "diagram export failed: " << result.error << '\n';
        return 7;
    }

    const auto png_path = std::filesystem::path(temp.string() + ".mo.png");
    std::ifstream png(png_path, std::ios::binary);
    std::array<std::uint8_t, 8> signature{};
    png.read(reinterpret_cast<char*>(signature.data()), static_cast<std::streamsize>(signature.size()));
    const std::array<std::uint8_t, 8> expected{137,80,78,71,13,10,26,10};
    if (signature != expected) return 8;

    const auto svg_path = std::filesystem::path(temp.string() + ".mo.svg");
    std::ifstream svg_file(svg_path, std::ios::binary);
    std::ostringstream svg_buffer;
    svg_buffer << svg_file.rdbuf();
    const std::string svg = svg_buffer.str();
    if (svg.find("Valence MO diagram") == std::string::npos ||
        svg.find("MO numbering is intentionally omitted") == std::string::npos ||
        svg.find(">3-a<") != std::string::npos) {
        std::cerr << "human export policy failed\n";
        return 9;
    }

    const auto json_path = std::filesystem::path(temp.string() + ".mo.json");
    std::ifstream json_file(json_path, std::ios::binary);
    std::ostringstream json_buffer;
    json_buffer << json_file.rdbuf();
    const std::string json = json_buffer.str();
    if (json.find("\"mode\": \"valence-central\"") == std::string::npos ||
        json.find("\"label\": \"3-a\"") == std::string::npos ||
        json.find("\"multicentre_label\": \"3c2e\"") == std::string::npos ||
        json.find("\"strict_salc_claimed\": false") == std::string::npos) {
        std::cerr << "machine-readable markers missing\n";
        return 10;
    }

    std::error_code ec;
    std::filesystem::remove(svg_path, ec);
    std::filesystem::remove(png_path, ec);
    std::filesystem::remove(json_path, ec);
    std::filesystem::remove(std::filesystem::path(temp.string() + ".mo.csv"), ec);

    for (int i=1;i<argc;++i) {
        if (!inspect_real_fchk(argv[i])) return 11;
    }

    std::cout << "mo_diagram_smoke ok\n";
    return 0;
}
