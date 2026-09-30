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
    struct Connectivity { bool covalent=false,coordination=false; };
    std::map<std::pair<std::size_t,std::size_t>,Connectivity> connectivity;
    for(const auto& e:graph.edges)if(e.strength==InteractionStrength::StrongConnectivity) {
        auto& pair=connectivity[std::minmax(e.atom_a,e.atom_b)];
        pair.covalent|=e.kind==InteractionKind::CovalentConnectivity;
        pair.coordination|=e.kind==InteractionKind::CoordinationContact;
    }
    std::map<std::pair<NboSpin,std::size_t>,const std::string*> raw_types;
    for(const auto& row:data.dataset.orbitals)if(row.id)raw_types[{row.spin,row.id-1}]=&row.kind;
    const auto has_type=[&](const NboStructureEvidence& e,const char* kind) {
        return std::any_of(e.orbitals.begin(),e.orbitals.end(),[&](const auto& ref) {
            if(ref.kind!=NboOrbitalKind::NBO)return false;
            const auto found=raw_types.find({ref.spin,ref.index});
            return found!=raw_types.end() && *found->second==kind;
        });
    };
    std::map<std::vector<std::size_t>,std::size_t> hyperedges;
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
            // A continuous index or membership in a multicentre group does not
            // create a two-centre bond. Conversely, absent total BD counts do
            // not invalidate existing connectivity or validated spin BD records.
            const bool lewis=(e.lewis_bond_count && *e.lewis_bond_count>0) || has_type(e,"BD");
            const bool covalent=found!=connectivity.end() && found->second.covalent;
            const bool coordination=found!=connectivity.end() && found->second.coordination;
            if(!lewis && !covalent && !coordination)continue;
            MoleculeOverlayBond b;b.atom_a=e.atoms[0];b.atom_b=e.atoms[1];b.evidence_index=i;
            b.continuous_index=e.wiberg;b.selected=selected;
            b.multiplicity=static_cast<int>(e.lewis_bond_count.value_or(1));
            if(!covalent && coordination &&
               (!e.lewis_bond_count || *e.lewis_bond_count==1))b.style=OverlayBondStyle::Coordination;
            out.bonds.push_back(b);
        }else if(e.kind=="multicentre" && e.atoms.size()>2){
            // 3Cn and 3C* remain inspectable orbitals, not extra chemical bonds.
            if(!has_type(e,"3C"))continue;
            auto atoms=e.atoms;std::sort(atoms.begin(),atoms.end());
            atoms.erase(std::unique(atoms.begin(),atoms.end()),atoms.end());
            if(atoms.size()<3 || atoms.back()>=atom_count)continue;
            const auto [it,inserted]=hyperedges.emplace(atoms,out.multicentre.size());
            if(inserted)out.multicentre.push_back({i,atoms,selected,{i}});
            else {auto& group=out.multicentre[it->second];group.selected|=selected;group.evidence_indices.push_back(i);}
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
