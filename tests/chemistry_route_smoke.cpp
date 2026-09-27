#include "cov/chemistry_route.hpp"
#include "cov/nbo_molecular_overlay.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace {
void require(bool condition,const char* message){
    if(!condition)throw std::runtime_error(message);
}
cov::Wavefunction canonical(){
    cov::Wavefunction wf;
    wf.atoms={{"O",8,0,0,0,8},{"O",8,0,0,2.4,8}};
    return wf;
}
cov::NboStructureEvidence value(std::string kind,cov::NboSpin spin,
                                std::size_t atom,double number){
    cov::NboStructureEvidence row;
    row.kind=std::move(kind);row.spin=spin;row.atoms={atom};row.value=number;
    return row;
}
void verify(const cov::Wavefunction& wf,const cov::NboIntegration& data){
    const auto routed=cov::route_chemistry(wf,&data);
    require(routed.total_atomic_charge.available() &&
            routed.total_atomic_charge.provider==cov::RoutedProvider::Nbo,
            "complete independently validated total NPA charges must route");
    require(routed.atomic_spin.available(),"complete NPA spin must route");
    for(auto atom:std::array<std::size_t,2>{0,1}){
        require(std::abs((*routed.total_atomic_charge.value)[atom])<1e-12,
                "alpha/beta charge rows must never overwrite total NPA charge");
        require(std::abs((*routed.atomic_spin.value)[atom]-1.0)<1e-12,
                "spin must retain its independent validated value");
    }
    const auto charge=cov::make_nbo_molecule_overlay(data,*routed.interaction_graph.value,
        wf.atoms.size(),{},std::nullopt,cov::AtomScalarMode::NaturalCharge,false,false,&routed);
    const auto spin=cov::make_nbo_molecule_overlay(data,*routed.interaction_graph.value,
        wf.atoms.size(),{},std::nullopt,cov::AtomScalarMode::SpinPopulation,false,false,&routed);
    for(auto atom:std::array<std::size_t,2>{0,1}){
        require(charge.atom_values[atom] && std::abs(*charge.atom_values[atom])<1e-12,
                "overlay must show routed total charge including actual zero");
        require(spin.atom_values[atom] && std::abs(*spin.atom_values[atom]-1)<1e-12,
                "overlay must show routed spin independently");
    }
}
} // namespace

int main(){
    const auto wf=canonical();
    cov::NboIntegration data;
    data.canonical_fingerprint=cov::nbo_canonical_fingerprint(wf);
    data.dataset.association.compatible=true;
    data.capabilities={{"charges",cov::NboCapabilityState::Available,"complete",{}},
                       {"spin",cov::NboCapabilityState::Available,"complete",{}}};
    for(std::size_t atom=0;atom<2;++atom){
        data.structure.push_back(value("atomic_charge",cov::NboSpin::Total,atom,0));
        data.structure.push_back(value("atomic_charge",cov::NboSpin::Alpha,atom,-0.5));
        data.structure.push_back(value("atomic_charge",cov::NboSpin::Beta,atom,0.5));
        data.structure.push_back(value("atomic_spin",cov::NboSpin::Total,atom,1));
    }
    verify(wf,data);
    std::reverse(data.structure.begin(),data.structure.end());verify(wf,data);
    std::rotate(data.structure.begin(),data.structure.begin()+3,data.structure.end());verify(wf,data);
    data.capabilities[1].state=cov::NboCapabilityState::Rejected;
    auto rejected=cov::route_chemistry(wf,&data);
    require(rejected.total_atomic_charge.available(),"rejected spin must not reject total charge");
    require(!rejected.atomic_spin.available() &&
            rejected.atomic_spin.status==cov::RoutedStatus::Rejected,
            "rejected spin must remain rejected");
    const auto overlay=cov::make_nbo_molecule_overlay(data,*rejected.interaction_graph.value,
        wf.atoms.size(),{},std::nullopt,cov::AtomScalarMode::SpinPopulation,false,false,&rejected);
    require(!overlay.atom_values[0]&&!overlay.atom_values[1],
            "rejected spin must not render a pseudo-value");
    data.structure.erase(std::find_if(data.structure.begin(),data.structure.end(),[](const auto& row){
        return row.kind=="atomic_charge"&&row.spin==cov::NboSpin::Total&&row.atoms==std::vector<std::size_t>{0};
    }));
    const auto missing=cov::route_chemistry(wf,&data);
    require(!missing.total_atomic_charge.available() &&
            !missing.total_atomic_charge.value,
            "missing one total record is incomplete, never zero");
}
