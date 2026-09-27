#include "cov/nbo_molecular_overlay.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <iomanip>
namespace cov {
MoleculeOverlay make_nbo_molecule_overlay(const NboIntegration& data,const InteractionGraph& graph,
    std::size_t atom_count,const std::vector<std::size_t>& selected_atoms,
    std::optional<std::size_t> selected_structure,AtomScalarMode mode,bool show_indices,bool show_e2,
    const RoutedAnalysis* routed) {
    MoleculeOverlay out;out.colour_mode=mode;out.selected_atoms=selected_atoms;
    out.show_bond_indices=show_indices;out.atom_values.resize(atom_count);out.scalar_range=0;
    // Scalar values come from the same complete, capability-checked result as
    // the UI. In particular alpha/beta NPA charges cannot overwrite total NPA.
    if(routed) {
        const auto& values=mode==AtomScalarMode::NaturalCharge?
            routed->total_atomic_charge:routed->atomic_spin;
        if(mode!=AtomScalarMode::Element && values.available() &&
           values.value->size()==atom_count)
            for(std::size_t atom=0;atom<atom_count;++atom) {
                out.atom_values[atom]=(*values.value)[atom];
                out.scalar_range=std::max(out.scalar_range,std::abs(*out.atom_values[atom]));
            }
    }
    std::map<std::pair<std::size_t,std::size_t>,InteractionKind> connectivity;
    for(const auto& e:graph.edges)if(e.strength==InteractionStrength::StrongConnectivity)
        connectivity[std::minmax(e.atom_a,e.atom_b)]=e.kind;
    for(std::size_t i=0;i<data.structure.size();++i){
        const auto& e=data.structure[i];const bool selected=selected_structure && *selected_structure==i;
        if(!routed &&
           ((e.kind=="atomic_charge" && e.spin==NboSpin::Total &&
             mode==AtomScalarMode::NaturalCharge && nbo_capability(data,"charges") &&
             nbo_capability(data,"charges")->available()) ||
            (e.kind=="atomic_spin" && mode==AtomScalarMode::SpinPopulation &&
             nbo_capability(data,"spin") && nbo_capability(data,"spin")->available()))){
            if(e.value && e.atoms.size()==1 && e.atoms[0]<atom_count){
                out.atom_values[e.atoms[0]]=e.value;
                out.scalar_range=std::max(out.scalar_range,std::abs(*e.value));
            }
        }else if(e.kind=="bond" && e.atoms.size()==2){
            const auto key=std::minmax(e.atoms[0],e.atoms[1]);const auto found=connectivity.find(key);
            // A continuous index between distant atoms is not a new covalent bond.
            if(!e.lewis_bond_count && found==connectivity.end())continue;
            MoleculeOverlayBond b;b.atom_a=e.atoms[0];b.atom_b=e.atoms[1];b.evidence_index=i;
            b.continuous_index=e.wiberg;b.selected=selected;
            b.multiplicity=static_cast<int>(e.lewis_bond_count.value_or(1));
            if(found!=connectivity.end() && found->second==InteractionKind::CoordinationContact &&
               (!e.lewis_bond_count || *e.lewis_bond_count==1))b.style=OverlayBondStyle::Coordination;
            else if(found!=connectivity.end() && found->second==InteractionKind::MulticentreSupport)
                b.style=OverlayBondStyle::Multicentre;
            else if(!e.lewis_bond_count)b.style=OverlayBondStyle::Unresolved;
            out.bonds.push_back(b);
        }else if(e.kind=="multicentre" && e.atoms.size()>2){
            out.multicentre.push_back({i,e.atoms,selected});
        }else if(e.kind=="donor_acceptor" && (show_e2 || selected) && e.orbitals.size()>=2){
            const auto* donor=nbo_orbital(data,e.orbitals[0]);const auto* acceptor=nbo_orbital(data,e.orbitals[1]);
            if(donor && acceptor && !donor->atoms.empty() && !acceptor->atoms.empty())
                out.relations.push_back({i,donor->atoms,acceptor->atoms,selected});
        }
    }
    out.scalar_range=std::max(out.scalar_range,1e-12);return out;
}
std::string serialize_molecule_overlay_scalars_json(const MoleculeOverlay& overlay,
                                                     const RoutedAnalysis* routed) {
    std::ostringstream out;out<<std::setprecision(17)
        <<"{\"schema\":\"cov.overlay.scalars.v1\",\"mode\":\"";
    const RoutedResult<std::vector<double>>* result=nullptr;
    switch(overlay.colour_mode) {
        case AtomScalarMode::Element:out<<"element";break;
        case AtomScalarMode::NaturalCharge:out<<"natural_charge";
            if(routed)result=&routed->total_atomic_charge;break;
        case AtomScalarMode::SpinPopulation:out<<"spin_population";
            if(routed)result=&routed->atomic_spin;break;
    }
    out<<"\",\"status\":\""<<(result?routed_status_name(result->status):"not_analysed")
       <<"\",\"provider\":\""<<(result?routed_provider_name(result->provider):"none")
       <<"\",\"values\":[";
    for(std::size_t i=0;i<overlay.atom_values.size();++i) {
        if(i)out<<',';
        if(overlay.atom_values[i])out<<*overlay.atom_values[i];else out<<"null";
    }
    out<<"]}";return out.str();
}
}
