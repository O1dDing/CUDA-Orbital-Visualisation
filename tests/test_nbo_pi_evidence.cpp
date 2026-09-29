#include "cov/pi_coupling.hpp"
#include <iostream>
#include <stdexcept>

namespace {
void require(bool ok,const char* reason){if(!ok)throw std::runtime_error(reason);}
}
int main(){try{
    // Deliberate missing-descriptor fault after valid transform/operator gates.
    // This fixture tests failure semantics, not the chemistry of a molecule.
    cov::Wavefunction w;w.basis_count=6;
    w.atoms={{"C",6,0,0,0,6},{"C",6,0,0,2,6}};
    w.ao_overlap.assign(36,0);
    cov::NboIntegration d;d.dataset.association.compatible=true;
    cov::NboMatrix eye;eye.rows=eye.columns=6;eye.values.assign(36,0);
    for(std::size_t k=0;k<6;++k){
        eye.values[k*6+k]=1;w.ao_overlap[k*6+k]=1;
        cov::MolecularOrbital mo;mo.source_orbital_index=k;mo.energy_hartree=1;
        mo.coefficients.assign(6,0);mo.coefficients[k]=1;w.orbitals.push_back(mo);
        cov::NboNao row;row.id=k+1;row.atom=k/3+1;row.type="Val(2p)";
        row.angular=k%3==0?"px":k%3==1?"py":"pz";d.dataset.naos.push_back(row);
    }
    cov::NboArchive archive;
    for(const char* kind:{"OVERLAP","FOCK"}){eye.kind=kind;archive.matrices.push_back(eye);}
    d.dataset.archive=archive;
    for(const char* kind:{"AONAO","NAOMO"}){eye.kind=kind;d.dataset.matrices.push_back(eye);}
    cov::NboNaoValidation v;v.available=v.direct_fchk_coefficients=true;v.effective_mo_columns=6;
    d.dataset.nao_validation.push_back(v);d.canonical_fingerprint=cov::nbo_canonical_fingerprint(w);
    const auto absent=cov::analyse_nbo_pi_couplings(w,d,{});
    require(absent.status=="not_applicable"&&!absent.evidence.empty(),"no-candidate control lost operator evidence");
    const auto failed=cov::analyse_nbo_pi_couplings(w,d,{{0,1}});
    require(!failed.evidence.empty(),"fault did not reach the angular gate");
    require(failed.status=="insufficient_evidence","missing angular evidence misreported as not applicable");
    require(failed.reason.find("descriptor")!=std::string::npos,"specific angular failure reason lost");
    std::cout<<"NBO pi evidence: no candidate and failed angular projection remain distinct\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
