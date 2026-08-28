#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
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
        result.multicentre.label=chemistry.multicentre_label;
        result.multicentre.source=AnnotationSource::Derived;
        result.multicentre.confidence=chemistry.confidence;
        result.multicentre.heuristic=false;
    }

    if (chemistry.channel.dominant==OrbitalAngularFamily::Pi &&
        chemistry.participating_atoms>2u) {
        result.delocalised_pi.available=true;
        result.delocalised_pi.participating_atoms=
            chemistry.participating_atoms;
        result.delocalised_pi.participating_electrons=
            chemistry.participating_electrons;
        result.delocalised_pi.label="delocalised-pi";
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

} // namespace

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
    for (std::size_t i=0;i<wavefunction.orbitals.size();++i) {
        if (wavefunction.orbitals[i].chemistry.available &&
            wavefunction.orbitals[i].chemistry.valence_manifold) {
            plan.included_indices.push_back(i);
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

    data.annotations.reserve(wavefunction.orbitals.size());
    for (const auto& orbital:wavefunction.orbitals) {
        data.annotations.push_back(annotate_orbital(orbital));
    }

    data.selection=build_valence_selection_plan(
        wavefunction,options,data.metadata);
    const auto labels=build_orbital_labels(
        wavefunction.orbitals,options.degeneracy);
    std::vector<double> axis_energies;
    axis_energies.reserve(data.selection.included_indices.size());
    data.levels.reserve(data.selection.included_indices.size());

    for (const auto index:data.selection.included_indices) {
        if (index>=data.metadata.size()) continue;
        MODiagramLevel level;
        level.metadata=data.metadata[index];
        level.annotation=data.annotations[index];
        level.chemistry=wavefunction.orbitals[index].chemistry;
        level.electrons=electron_glyphs_for_orbital(
            wavefunction.orbitals[index],
            data.frontier.separate_spin_sets);
        level.homo=data.frontier.homo &&
                   *data.frontier.homo==index;
        level.lumo=data.frontier.lumo &&
                   *data.frontier.lumo==index;
        level.layout_energy_hartree=group_layout_energy(
            data.metadata,labels,index);
        axis_energies.push_back(level.layout_energy_hartree);
        data.levels.push_back(std::move(level));
    }

    data.energy_transform=build_energy_transform(
        axis_energies,options.energy_axis_mode,
        options.nonlinear_minimum_gap_weight);
    return data;
}

} // namespace cov
