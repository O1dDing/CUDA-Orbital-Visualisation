#include "cov/nbo_aomo_labels.hpp"
#include "cov/orbital_symmetry.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using Mat=std::array<double,9>;
void require(bool value,const char* why){if(!value)throw std::runtime_error(why);}
Mat mul(const Mat& a,const Mat& b){Mat c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];return c;}
Mat transpose(const Mat& a){Mat b{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)b[3*i+j]=a[3*j+i];return b;}
void shell(cov::Wavefunction& w,std::size_t atom,unsigned l,double exponent){cov::Shell s;s.atom_index=static_cast<std::uint32_t>(atom);s.angular_momentum=l;s.primitive_offset=static_cast<std::uint32_t>(w.primitives.size());s.primitive_count=1;s.basis_offset=w.basis_count;s.pure=0;w.primitives.push_back({float(exponent),1.f});w.shells.push_back(s);w.basis_count+=l==1?3:1;}
struct Fixture {cov::Wavefunction w;cov::NboIntegration data;cov::NboSalcModel salc;};
Fixture c2v(){
    Fixture f;auto& w=f.w;
    // Orthogonal model basis, on a bent triatomic in yz; every supplied
    // character follows from the actual matrices and atom permutations.
    w.atoms={{"O",8,0,0,.2},{"H",1,0,1,-.8},{"H",1,0,-1,-.8}};
    shell(w,0,0,2);shell(w,0,0,1);shell(w,0,1,1);shell(w,1,0,1);shell(w,2,0,1);
    w.ao_overlap.assign(49,0);for(int i=0;i<7;++i)w.ao_overlap[i*7+i]=1;
    const double q=1/std::sqrt(2.);
    const std::vector<std::vector<double>> coeff={
        {1,0,0,0,0,0,0}, {0,q,0,0,0,.5,.5},
        {0,0,0,q,0,.5,-.5}, {0,0,0,0,1,0,0},
        {0,0,1,0,0,0,0}, {0,q,0,0,0,-.5,-.5},
        {0,0,0,q,0,-.5,.5}};
    const double energy[]={-20,-1,-.5,-.3,-.2,.1,.2};
    for(int i=0;i<7;++i){cov::MolecularOrbital mo;mo.coefficients=coeff[i];mo.energy_hartree=energy[i];mo.occupation=i<5?2:0;w.orbitals.push_back(mo);}
    f.salc.available=true;f.salc.group_verified=true;f.salc.point_group=f.salc.used_group="C2v";
    const Mat ops[]={Mat{1,0,0,0,1,0,0,0,1},Mat{-1,0,0,0,1,0,0,0,1},Mat{1,0,0,0,-1,0,0,0,1},Mat{-1,0,0,0,-1,0,0,0,1}};
    for(int i=0;i<4;++i){cov::SymmetryOperation op;op.matrix=ops[i];op.atom_permutation=i<2?std::vector<std::size_t>{0,1,2}:std::vector<std::size_t>{0,2,1};f.salc.operations.push_back(op);}
    // Intentionally reverse source ordering versus energy ordering.
    for(int i=0;i<2;++i){cov::NboSalcOrbital o;o.id=i?"plus":"minus";o.fragment_id="pair";o.subspace_id=o.id;o.energy_hartree=i?.11:.26;o.atoms={1,2};f.salc.orbitals.push_back(o);cov::NboSalcSubspace s;s.id=o.id;s.fragment_id="pair";s.orbital_indices={std::size_t(i)};s.dimension=s.irrep_dimension=s.multiplicity=1;s.symmetry_verified=true;s.characters=i?std::vector<double>{1,1,1,1}:std::vector<double>{1,1,-1,-1};f.salc.subspaces.push_back(s);}
    return f;
}
void rotate(Fixture& f,const Mat& r){
    for(auto& a:f.w.atoms){const std::array<double,3> p{a.x,a.y,a.z};a.x=r[0]*p[0]+r[1]*p[1]+r[2]*p[2];a.y=r[3]*p[0]+r[4]*p[1]+r[5]*p[2];a.z=r[6]*p[0]+r[7]*p[1]+r[8]*p[2];}
    for(auto& mo:f.w.orbitals){const std::array<double,3> p{mo.coefficients[2],mo.coefficients[3],mo.coefficients[4]};for(int a=0;a<3;++a)mo.coefficients[2+a]=r[3*a]*p[0]+r[3*a+1]*p[1]+r[3*a+2]*p[2];}
    for(auto& op:f.salc.operations)op.matrix=mul(mul(r,op.matrix),transpose(r));
}
void compare(const cov::ui::NboAomoNames& a,const cov::ui::NboAomoNames& b){require(a.canonical.size()==b.canonical.size(),"canonical count changed");for(std::size_t i=0;i<a.canonical.size();++i)require(a.canonical[i].label==b.canonical[i].label&&a.canonical[i].verified==b.canonical[i].verified,"rotation/operation-order changed canonical name");for(std::size_t i=0;i<a.salc.size();++i)require(a.salc[i].label==b.salc[i].label,"rotation/operation-order changed SALC name");}
Fixture linear(){
    Fixture f;auto& w=f.w;w.atoms={{"C",6,0,0,0},{"O",8,0,0,-2},{"O",8,0,0,2}};
    shell(w,0,0,3);shell(w,0,0,1);shell(w,0,1,1);shell(w,0,0,.4);shell(w,1,0,1);shell(w,2,0,1);
    w.ao_overlap.assign(64,0);for(int i=0;i<8;++i)w.ao_overlap[i*8+i]=1;
    const double q=1/std::sqrt(2.);const std::vector<std::vector<double>> coeff={{1,0,0,0,0,0,0,0},{0,0,0,0,0,0,q,q},{0,0,0,0,1,0,0,0},{0,0,1,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,1,0,0,0,0,0,0},{0,0,0,0,0,1,0,0},{0,0,0,0,0,0,q,-q}};
    const double energies[]={-20,-1,-.6,-.4,-.4,-.4,.2,.3};
    for(int i=0;i<8;++i){cov::MolecularOrbital mo;mo.coefficients=coeff[i];mo.energy_hartree=energies[i];w.orbitals.push_back(mo);cov::NboOrbitalDescriptor descriptor;descriptor.ref={cov::NboOrbitalKind::NAO,cov::NboSpin::Total,std::size_t(i)};descriptor.coefficients=coeff[i];f.data.orbitals.push_back(descriptor);cov::NboSalcOrbital o;o.fragment_id="fixed";o.energy_hartree=energies[i];o.terms.push_back({descriptor.ref,1});f.salc.orbitals.push_back(o);}
    f.salc.available=f.salc.group_verified=true;f.salc.point_group="Dinfh";f.salc.used_group="finite sampling subgroup of Dinfh";
    for(int inv=0;inv<2;++inv)for(int reflection=0;reflection<2;++reflection)for(int k=0;k<3;++k){const double a=2*3.14159265358979323846*k/3,c=std::cos(a),s=std::sin(a);cov::SymmetryOperation op;op.matrix={c,-s,0,s,c,0,0,0,1};if(reflection)for(int row=0;row<3;++row)op.matrix[3*row+1]*=-1;if(inv)for(auto& x:op.matrix)x=-x;op.atom_permutation=inv?std::vector<std::size_t>{0,2,1}:std::vector<std::size_t>{0,1,2};f.salc.operations.push_back(op);}
    // A reducible central p family and repeated Sigma copies are intentional.
    for(const auto& indices:std::vector<std::vector<std::size_t>>{{0},{1},{2,3,4},{5,6},{7}}){cov::NboSalcSubspace sub;sub.fragment_id="fixed";sub.dimension=indices.size();sub.orbital_indices=indices;sub.symmetry_verified=true;f.salc.subspaces.push_back(sub);}
    return f;
}
void rotate_general(Fixture& f,const Mat& rotation){
    cov::SymmetryOperation op;op.matrix=rotation;for(std::size_t i=0;i<f.w.atoms.size();++i)op.atom_permutation.push_back(i);
    for(auto& mo:f.w.orbitals)mo.coefficients=cov::apply_orbital_symmetry_operation(f.w,op,mo.coefficients,1);
    for(auto& descriptor:f.data.orbitals)descriptor.coefficients=cov::apply_orbital_symmetry_operation(f.w,op,descriptor.coefficients,1);
    for(auto& atom:f.w.atoms){const std::array<double,3> p{atom.x,atom.y,atom.z};atom.x=rotation[0]*p[0]+rotation[1]*p[1]+rotation[2]*p[2];atom.y=rotation[3]*p[0]+rotation[4]*p[1]+rotation[5]*p[2];atom.z=rotation[6]*p[0]+rotation[7]*p[1]+rotation[8]*p[2];}
    for(auto& operation:f.salc.operations)operation.matrix=mul(mul(rotation,operation.matrix),transpose(rotation));
}
Fixture axial(unsigned order,bool horizontal){
    Fixture f;auto& w=f.w;w.atoms.push_back({"C",6,0,0,horizontal?0.:.4});
    for(unsigned k=0;k<order;++k){const double angle=2*3.14159265358979323846*k/order;w.atoms.push_back({"H",1,std::cos(angle),std::sin(angle),horizontal?0.:-.5});}
    shell(w,0,1,1);w.ao_overlap={1,0,0,0,1,0,0,0,1};for(unsigned i=0;i<3;++i){cov::MolecularOrbital mo;mo.coefficients.assign(3,0);mo.coefficients[i]=1;mo.energy_hartree=i==2?-.8:-.4;w.orbitals.push_back(mo);}
    f.salc.available=f.salc.group_verified=true;f.salc.point_group=f.salc.used_group=std::string(horizontal?"D":"C")+std::to_string(order)+(horizontal?"h":"v");
    for(int z=0;z<(horizontal?2:1);++z)for(int y=0;y<2;++y)for(unsigned k=0;k<order;++k){const double angle=2*3.14159265358979323846*k/order,c=std::cos(angle),s=std::sin(angle);cov::SymmetryOperation op;op.matrix=mul(Mat{c,-s,0,s,c,0,0,0,1},Mat{1,0,0,0,y?-1.:1.,0,0,0,z?-1.:1.});
        for(const auto& atom:w.atoms){const std::array<double,3> p{atom.x,atom.y,atom.z};std::array<double,3> target{};for(int r=0;r<3;++r)for(int j=0;j<3;++j)target[r]+=op.matrix[3*r+j]*p[j];std::size_t found=w.atoms.size();for(std::size_t j=0;j<w.atoms.size();++j){const auto& other=w.atoms[j];if(other.atomic_number==atom.atomic_number&&std::abs(other.x-target[0])+std::abs(other.y-target[1])+std::abs(other.z-target[2])<1e-9)found=j;}require(found<w.atoms.size(),"axial group atom permutation invalid");op.atom_permutation.push_back(found);}f.salc.operations.push_back(op);
    }return f;
}
Fixture spectral_copies(double first,double second){
    auto f=linear();shell(f.w,0,1,.3);shell(f.w,0,1,.2);const auto n=f.w.basis_count;
    for(auto& mo:f.w.orbitals)mo.coefficients.resize(n,0);f.w.ao_overlap.assign(n*n,0);for(std::size_t i=0;i<n;++i)f.w.ao_overlap[i*n+i]=1;
    for(std::size_t i=8;i<14;++i){cov::MolecularOrbital mo;mo.coefficients.assign(n,0);mo.coefficients[i]=1;mo.energy_hartree=i<11?first:second;if(i==10||i==13)mo.energy_hartree+=.7;f.w.orbitals.push_back(mo);}
    f.data.orbitals.clear();f.salc.orbitals.clear();f.salc.subspaces.clear();f.salc.links.clear();
    const double c=std::cos(.4),s=std::sin(.4);
    const std::vector<std::vector<std::pair<std::size_t,double>>> entries={{{8,c},{11,s}},{{9,1}},{{8,-s},{11,c}},{{12,1}},{{2,1}},{{3,1}}};
    for(std::size_t i=0;i<entries.size();++i){cov::NboOrbitalDescriptor descriptor;descriptor.ref={cov::NboOrbitalKind::NAO,cov::NboSpin::Total,i};descriptor.coefficients.assign(n,0);for(const auto& [a,value]:entries[i])descriptor.coefficients[a]=value;f.data.orbitals.push_back(descriptor);cov::NboSalcOrbital o;o.fragment_id="copies";o.terms={{descriptor.ref,1}};double energy=0;for(std::size_t j=0;j<f.w.orbitals.size();++j){double overlap=0;for(std::size_t a=0;a<n;++a)overlap+=descriptor.coefficients[a]*f.w.orbitals[j].coefficients[a];energy+=overlap*overlap*f.w.orbitals[j].energy_hartree;f.salc.links.push_back({i,j,overlap,overlap*overlap});}o.energy_hartree=energy;f.salc.orbitals.push_back(o);}
    for(const auto& members:std::vector<std::vector<std::size_t>>{{0,1,2,3},{4,5}}){cov::NboSalcSubspace sub;sub.fragment_id="copies";sub.dimension=members.size();sub.orbital_indices=members;sub.symmetry_verified=true;f.salc.subspaces.push_back(sub);}
    cov::NboSalcEnergyEvidence proof;proof.available=proof.electronic_symmetry_verified=true;proof.status="verified_same_operator";proof.canonical_columns_checked=f.w.orbitals.size();f.salc.energies={proof};return f;
}
Fixture actual_mixed_salc_span(double mixing){
    auto f=axial(3,false);const double p=std::sqrt(1-mixing*mixing);
    const std::vector<std::vector<double>> columns{{1,0,0},{0,p,mixing},{0,-mixing,p}};
    for(std::size_t i=0;i<columns.size();++i){cov::NboOrbitalDescriptor descriptor;descriptor.ref={cov::NboOrbitalKind::NAO,cov::NboSpin::Total,i};descriptor.coefficients=columns[i];f.data.orbitals.push_back(descriptor);
        cov::NboSalcOrbital o;o.fragment_id="actual mixed";o.energy_hartree=-.8+.2*double(i);o.terms={{descriptor.ref,1}};f.salc.orbitals.push_back(o);}
    cov::NboSalcSubspace sub;sub.fragment_id="actual mixed";sub.dimension=3;sub.orbital_indices={0,1,2};sub.symmetry_verified=true;
    for(const auto& op:f.salc.operations)sub.characters.push_back(op.matrix[0]+op.matrix[4]+op.matrix[8]);
    f.salc.subspaces={sub};return f;
}
Fixture slightly_coupled_copies(double second_mixing){
    auto f=axial(3,false);shell(f.w,0,1,.4);shell(f.w,0,0,.3);shell(f.w,0,1,.2);
    const std::size_t n=f.w.basis_count;f.w.ao_overlap.assign(n*n,0);for(std::size_t i=0;i<n;++i)f.w.ao_overlap[i*n+i]=1;
    f.w.orbitals.clear();const double energies[]={-.4,-.4,-.8,.2,.2,.9,1.1,2,2,3};
    for(std::size_t i=0;i<n;++i){cov::MolecularOrbital mo;mo.coefficients.assign(n,0);mo.coefficients[i]=1;mo.energy_hartree=energies[i];f.w.orbitals.push_back(mo);}
    // Two E copies are connected through a slightly contaminated A1 source
    // column in the full operation graph. Their actual source partner spaces
    // may still pass the unchanged strict complete-operation closure test.
    for(const auto& [member,angle]:std::vector<std::pair<std::size_t,double>>{{1,1e-4},{4,second_mixing}}){
        const auto a=f.w.orbitals[member].coefficients,b=f.w.orbitals[6].coefficients;
        for(std::size_t j=0;j<n;++j){f.w.orbitals[member].coefficients[j]=std::cos(angle)*a[j]+std::sin(angle)*b[j];f.w.orbitals[6].coefficients[j]=-std::sin(angle)*a[j]+std::cos(angle)*b[j];}}
    f.data.orbitals.clear();f.salc.orbitals.clear();f.salc.subspaces.clear();f.salc.links.clear();
    for(std::size_t i=0;i<n;++i){cov::NboOrbitalDescriptor descriptor;descriptor.ref={cov::NboOrbitalKind::NAO,cov::NboSpin::Total,i};descriptor.coefficients=f.w.orbitals[i].coefficients;f.data.orbitals.push_back(descriptor);
        cov::NboSalcOrbital o;o.fragment_id="unchanged source copies";o.energy_hartree=energies[i];o.terms={{descriptor.ref,1}};f.salc.orbitals.push_back(o);
        for(std::size_t j=0;j<n;++j)f.salc.links.push_back({i,j,i==j?1.:0.,i==j?1.:0.});}
    cov::NboSalcSubspace sub;sub.fragment_id="unchanged source copies";sub.dimension=n;sub.symmetry_verified=true;for(std::size_t i=0;i<n;++i)sub.orbital_indices.push_back(i);f.salc.subspaces={sub};
    cov::NboSalcEnergyEvidence proof;proof.available=proof.electronic_symmetry_verified=true;proof.status="verified_same_operator";proof.canonical_columns_checked=n;f.salc.energies={proof};return f;
}
}
int main(){try{
    auto f=c2v();const auto before=f.w.orbitals;const auto a=cov::ui::build_nbo_aomo_names(f.w,f.data,&f.salc);
    const char* irreps[]={"A1","A1","B2","A1","B1","A1","B2"};const std::size_t ordinals[]={1,2,1,3,1,4,2};
    for(std::size_t i=0;i<7;++i){require(a.canonical[i].verified&&a.canonical[i].irrep==irreps[i]&&a.canonical[i].ordinal==ordinals[i],"full-set canonical ordinal/irrep wrong");require(f.w.orbitals[i].coefficients==before[i].coefficients,"naming modified coefficients");}
    require(a.canonical[1].label=="2a\xE2\x82\x81","ordinal must precede Unicode irrep subscript");
    const auto standalone=cov::ui::canonical_mo_names(f.w);
    require(standalone->salc.empty(),"standalone canonical naming must not fabricate SALCs");
    for(std::size_t i=0;i<7;++i){
        require(standalone->canonical[i].verified&&standalone->canonical[i].irrep==irreps[i]&&
                standalone->canonical[i].ordinal==ordinals[i],"canonical naming incorrectly depends on NBO capability");
        require(cov::ui::canonical_mo_display_label(f.w,i)==a.canonical[i].label,"standalone main label differs from attached name");
    }
    auto asymmetric=f.w;asymmetric.atoms.push_back({"He",2,.37,.19,.63});
    const auto c1=cov::ui::canonical_mo_names(asymmetric);
    require(c1->canonical[0].verified&&c1->canonical[0].irrep=="A"&&c1->canonical[0].ordinal==1,
            "C1 canonical labels must be available without an NBO frame");
    auto plane=f.w;plane.atoms[1].y=1.35;plane.atoms[1].z=-.83;
    const auto cs=cov::ui::canonical_mo_names(plane);
    require(cs->canonical[0].verified&&cs->canonical[0].irrep=="A'"&&
            cs->canonical[4].verified&&cs->canonical[4].irrep=="A''", "Cs canonical labels must follow mirror characters without NBO");
    cov::Wavefunction fallback;fallback.orbital_occupation_model=cov::OrbitalOccupationModel::ExplicitSpin;
    fallback.orbitals.resize(3);fallback.orbitals[0].source_orbital_index=6;
    fallback.orbitals[1].source_orbital_index=2;fallback.orbitals[1].spin=cov::Spin::Beta;
    cov::ui::NboAomoName unknown;
    fallback.orbitals[0].symmetry="E1g";
    require(cov::ui::canonical_mo_current_irrep(fallback,0,&unknown)=="?",
            "unverified source aliases must not become current full-MO symmetry");
    cov::ui::NboAomoName verified_pi;verified_pi.verified=true;verified_pi.irrep="Pi_g";
    require(cov::ui::canonical_mo_current_irrep(fallback,0,&verified_pi)=="πg" &&
            fallback.orbitals[0].symmetry=="E1g",
            "current irrep display must retain original source spelling separately");
    require(cov::ui::canonical_mo_display_label(fallback,0,&unknown)=="MO 7 [alpha]"&&
            cov::ui::canonical_mo_display_label(fallback,1,&unknown)=="MO 3 [beta]", "fallback must use the actual source spin block rather than list position");
    require(cov::ui::canonical_mo_display_label(fallback,2,&unknown)=="MO [list] 3 [alpha]",
            "unavailable source index must be explicitly distinguished from a source number");
    auto alpha_name=a.canonical[1];
    auto alpha_only=fallback;alpha_only.orbitals.resize(1);
    require(cov::ui::canonical_mo_display_label(alpha_only,0,&alpha_name)=="2a\xE2\x82\x81 [alpha]",
            "an explicit alpha-only orbital set must retain its spin on scientific labels");
    auto known_without_order=a.canonical[1];known_without_order.ordinal=0;known_without_order.label="?a1";
    require(cov::ui::canonical_mo_display_label(fallback,0,&known_without_order)=="a\xE2\x82\x81 [MO 7] [alpha]"&&
            cov::ui::orbital_irrep_display_label(known_without_order)=="a\xE2\x82\x81",
            "a known irrep must lead the stable source identity without inventing an occurrence ordinal");
    auto reloaded=f.w;
    require(cov::ui::canonical_mo_names(reloaded)->canonical[0].ordinal==1,"initial reload fixture label missing");
    reloaded.orbitals[1].energy_hartree=reloaded.orbitals[0].energy_hartree;
    cov::ui::invalidate_canonical_mo_names_cache();
    const auto reload_names=cov::ui::canonical_mo_names(reloaded);
    require(reload_names->canonical[0].verified&&reload_names->canonical[0].ordinal==0&&
            reload_names->canonical[1].verified&&reload_names->canonical[1].ordinal==0,
            "input reload must invalidate labels even when object and buffers retain their addresses");
    require(a.salc[0].irrep=="B2"&&a.salc[1].irrep=="A1"&&a.salc[0].ordinal==1&&a.salc[1].ordinal==1,"pair labels must follow measured characters, not source ordering or sign text");
    auto tied=f;tied.w.orbitals[1].energy_hartree=tied.w.orbitals[0].energy_hartree;const auto tied_names=cov::ui::build_nbo_aomo_names(tied.w,tied.data,&tied.salc);require(tied_names.canonical[0].verified&&tied_names.canonical[1].verified&&tied_names.canonical[0].ordinal==0&&tied_names.canonical[1].ordinal==0&&tied_names.canonical[3].ordinal==3,"source order cannot break an unresolved equal-energy same-irrep copy tie, but later total counts remain known");
    auto rotated=f;const double c=std::cos(.713),s=std::sin(.713);rotate(rotated,Mat{c,-s,0,s,c,0,0,0,1});compare(a,cov::ui::build_nbo_aomo_names(rotated.w,rotated.data,&rotated.salc));
    auto swapped=f;rotate(swapped,Mat{0,-1,0,1,0,0,0,0,1});std::swap(swapped.salc.operations[1],swapped.salc.operations[2]);for(auto& sub:swapped.salc.subspaces)std::swap(sub.characters[1],sub.characters[2]);compare(a,cov::ui::build_nbo_aomo_names(swapped.w,swapped.data,&swapped.salc));
    auto mixed=f;const auto x=mixed.w.orbitals[0].coefficients,y=mixed.w.orbitals[4].coefficients;const double q=.02,p=std::sqrt(1-q*q);for(std::size_t j=0;j<7;++j){mixed.w.orbitals[0].coefficients[j]=p*x[j]+q*y[j];mixed.w.orbitals[4].coefficients[j]=-q*x[j]+p*y[j];}
    const auto m=cov::ui::build_nbo_aomo_names(mixed.w,mixed.data,&mixed.salc);require(!m.canonical[0].verified&&!m.canonical[4].verified,"99.96 percent dominant weight must not certify a mixed eigenfunction");require(m.canonical[1].verified&&m.canonical[1].ordinal==0,"unknown earlier state must not silently corrupt symmetry ordinal");
    const auto mixed_standalone=cov::ui::canonical_mo_names(mixed.w);
    require(!mixed_standalone->canonical[0].verified&&!mixed_standalone->canonical[4].verified&&
            mixed_standalone->canonical[1].verified&&mixed_standalone->canonical[1].ordinal==0,
            "standalone naming must preserve mixed-orbital and incomplete ordinal gates");
    auto low_mixed=mixed;low_mixed.w.orbitals[4].energy_hartree=-19;const double mild=.0003,mild_p=std::sqrt(1-mild*mild);for(std::size_t j=0;j<7;++j){low_mixed.w.orbitals[0].coefficients[j]=mild_p*x[j]+mild*y[j];low_mixed.w.orbitals[4].coefficients[j]=-mild*x[j]+mild_p*y[j];}const auto counted_core=cov::ui::build_nbo_aomo_names(low_mixed.w,low_mixed.data,&low_mixed.salc);require(!counted_core.canonical[0].verified&&!counted_core.canonical[4].verified&&counted_core.canonical[1].ordinal==2&&counted_core.canonical[2].ordinal==1,"complete invariant mixed core must supply integer counts without relabelling its individual members");
    auto reused_derived=low_mixed;
    for(const auto index:{std::size_t(0),std::size_t(4)}){
        const std::string label=index==0?"A1":"B1";
        reused_derived.w.orbitals[index].symmetry=label;
        reused_derived.w.orbitals[index].symmetry_provenance=cov::DataProvenance::Derived;
        cov::DerivedOrbitalSymmetryAssignment proof;proof.point_group="C2v";proof.label=label;proof.orbital_indices={index};proof.subspace_retention=.999;
        reused_derived.w.derived_orbital_symmetry_assignments.push_back(proof);
    }
    const auto recounted=cov::ui::build_nbo_aomo_names(reused_derived.w,reused_derived.data,&reused_derived.salc);
    require(!recounted.canonical[0].verified&&!recounted.canonical[4].verified&&recounted.canonical[1].ordinal==2&&recounted.canonical[2].ordinal==1,
            "looser native derived labels must not remove mixed partners from complete-set counting");
    require(reused_derived.w.orbitals[0].symmetry=="A1"&&reused_derived.w.orbitals[4].symmetry=="B1","strict naming changed literal derived records");
    auto reducible=f;reducible.salc.subspaces.resize(1);auto& red=reducible.salc.subspaces[0];red.orbital_indices={0,1};red.dimension=2;red.irrep_dimension=0;red.multiplicity=0;red.characters={2,2,0,0};const auto bad=cov::ui::build_nbo_aomo_names(reducible.w,reducible.data,&reducible.salc);require(!bad.salc[0].verified&&!bad.salc[1].verified,"reducible SALC must not get a dimension-derived irrep");
    auto actual_mixed=actual_mixed_salc_span(.0003);
    const auto actual_coefficients=actual_mixed.data.orbitals;const auto actual_source_orbitals=actual_mixed.salc.orbitals;
    const auto actual_names=cov::ui::build_nbo_aomo_names(actual_mixed.w,actual_mixed.data,&actual_mixed.salc);
    require(actual_names.salc[0].verified&&actual_names.salc[0].irrep=="E"&&actual_names.salc[0].ordinal==0&&
        actual_names.salc[0].status=="verified_isotypic_member"&&actual_names.salc[0].partner_block_id.empty(),
        "pure source SALC in a verified mixed span needs its measured irrep without fabricated partners");
    for(std::size_t i=1;i<3;++i)require(!actual_names.salc[i].verified&&actual_names.salc[i].status=="mixed_source_member"&&
        actual_names.salc[i].projection_residual&&*actual_names.salc[i].projection_residual>4e-8&&
        actual_names.salc[i].partner_block_id.empty(),"dominant SALC weight must not certify an actually mixed source column");
    for(std::size_t i=0;i<3;++i){require(actual_names.salc[i].containing_members==std::vector<std::size_t>({0,1,2})&&
        actual_names.salc[i].containing_irreps.size()==2,"mixed SALC containing span/content lost");
        require(actual_mixed.data.orbitals[i].coefficients==actual_coefficients[i].coefficients&&
            actual_mixed.salc.orbitals[i].energy_hartree==actual_source_orbitals[i].energy_hartree&&
            actual_mixed.salc.orbitals[i].terms.front().coefficient==actual_source_orbitals[i].terms.front().coefficient,
            "source SALC coefficients, terms or energies changed during naming");}
    require(actual_names.salc[0].projection_residual&&*actual_names.salc[0].projection_residual<1e-12,
        "pure SALC projector residual missing");
    const auto actual_json=cov::ui::serialize_orbital_names_json(actual_names);
    require(actual_json.find("verified_isotypic_member")!=std::string::npos&&actual_json.find("mixed_source_member")!=std::string::npos&&
        actual_json.find("projection_residual_squared")!=std::string::npos,"mixed SALC JSON evidence incomplete");
    auto actual_rotated=actual_mixed;const double ar=std::cos(.437),as=std::sin(.437);
    rotate_general(actual_rotated,Mat{ar,0,as,0,1,0,-as,0,ar});compare(actual_names,cov::ui::build_nbo_aomo_names(actual_rotated.w,actual_rotated.data,&actual_rotated.salc));
    auto omitted_actual=actual_mixed;omitted_actual.salc.orbitals[0].terms.clear();
    const auto omitted_names=cov::ui::build_nbo_aomo_names(omitted_actual.w,omitted_actual.data,&omitted_actual.salc);
    for(const auto& name:omitted_names.salc)require(!name.verified&&!name.projection_residual&&name.status=="stored_span_unclassified",
        "stored reducible characters cannot identify any individual SALC when actual columns are omitted");
    auto failed_actual=actual_mixed;failed_actual.salc.subspaces[0].orbital_indices={0};failed_actual.salc.subspaces[0].dimension=1;
    failed_actual.salc.subspaces[0].characters.assign(failed_actual.salc.operations.size(),1);
    const auto failed_names=cov::ui::build_nbo_aomo_names(failed_actual.w,failed_actual.data,&failed_actual.salc);
    require(!failed_names.salc[0].verified&&failed_names.salc[0].status=="subspace_not_closed",
        "stored A1 signature must not override a present actual E column whose containing span is not closed");
    auto open=f;for(auto mo:before){mo.spin=cov::Spin::Beta;mo.energy_hartree+=.01;open.w.orbitals.push_back(mo);}const auto spin=cov::ui::build_nbo_aomo_names(open.w,open.data,&open.salc);require(spin.canonical[1].ordinal==2&&spin.canonical[8].ordinal==2&&spin.canonical[1].label.find("[alpha]")!=std::string::npos&&spin.canonical[8].label.find("[beta]")!=std::string::npos,"spin counters must be independent");
    auto unsupported=f;unsupported.salc.used_group="finite sampling subgroup of Cinfv";const auto u=cov::ui::build_nbo_aomo_names(unsupported.w,unsupported.data,&unsupported.salc);require(!u.canonical[0].verified&&!u.salc[0].verified&&u.canonical[0].label=="?","unsupported subgroup must remain unknown without inventing an irrep ordinal");
    auto ambiguous=f;ambiguous.w.atoms.push_back({"H",1,1,0,1});ambiguous.w.atoms.push_back({"H",1,-1,0,1});ambiguous.salc.operations[0].atom_permutation={0,1,2,3,4};ambiguous.salc.operations[1].atom_permutation={0,1,2,4,3};ambiguous.salc.operations[2].atom_permutation={0,2,1,3,4};ambiguous.salc.operations[3].atom_permutation={0,2,1,4,3};const auto axes=cov::ui::build_nbo_aomo_names(ambiguous.w,ambiguous.data,&ambiguous.salc);require(axes.canonical[0].verified&&axes.canonical[2].verified&&axes.salc[0].verified&&axes.canonical[2].irrep=="B2"&&!axes.canonical[2].detail.empty(),"tied molecular mirrors must use a recorded coordinate convention with x=B1/y=B2");
    // Octahedral central p triplet: an existing whole-orbital assignment is
    // one irrep occurrence, never three consecutive T1u ordinals.
    cov::Wavefunction oct;oct.atoms={{"M",26,0,0,0},{"H",1,1,0,0},{"H",1,-1,0,0},{"H",1,0,1,0},{"H",1,0,-1,0},{"H",1,0,0,1},{"H",1,0,0,-1}};shell(oct,0,1,1);oct.ao_overlap={1,0,0,0,1,0,0,0,1};
    for(int i=0;i<3;++i){cov::MolecularOrbital mo;mo.energy_hartree=-.5;mo.coefficients.assign(3,0);mo.coefficients[i]=1;mo.symmetry="T1u";mo.symmetry_provenance=cov::DataProvenance::Derived;oct.orbitals.push_back(mo);}cov::DerivedOrbitalSymmetryAssignment assignment;assignment.point_group="Oh";assignment.label="T1u";assignment.orbital_indices={0,1,2};assignment.subspace_retention=1;oct.derived_orbital_symmetry_assignments.push_back(assignment);const auto triplet=cov::ui::build_nbo_aomo_names(oct,f.data,nullptr);for(const auto& name:triplet.canonical)require(name.verified&&name.ordinal==1&&name.irrep=="T1u","multiplet members must share one ordinal");
    cov::NboSalcModel oct_salc;oct_salc.available=oct_salc.group_verified=true;oct_salc.point_group=oct_salc.used_group="Oh";
    std::array<int,3> perm{0,1,2};do{for(int signs=0;signs<8;++signs){cov::SymmetryOperation op;op.matrix={};for(int r=0;r<3;++r)op.matrix[3*r+perm[r]]=(signs&(1<<r))?-1:1;for(const auto& atom:oct.atoms){const std::array<double,3> xyz{atom.x,atom.y,atom.z};std::array<double,3> transformed{};for(int r=0;r<3;++r)for(int c=0;c<3;++c)transformed[r]+=op.matrix[3*r+c]*xyz[c];std::size_t target=oct.atoms.size();for(std::size_t j=0;j<oct.atoms.size();++j){const auto& b=oct.atoms[j];if(b.atomic_number==atom.atomic_number&&std::abs(b.x-transformed[0])+std::abs(b.y-transformed[1])+std::abs(b.z-transformed[2])<1e-10)target=j;}require(target<oct.atoms.size(),"octahedral fixture operation mapping failed");op.atom_permutation.push_back(target);}oct_salc.operations.push_back(op);}}while(std::next_permutation(perm.begin(),perm.end()));
    cov::NboSalcSubspace sub;sub.id="p-triplet";sub.fragment_id="centre";sub.dimension=sub.irrep_dimension=3;sub.multiplicity=1;sub.symmetry_verified=true;sub.orbital_indices={0,1,2};for(const auto& op:oct_salc.operations)sub.characters.push_back(op.matrix[0]+op.matrix[4]+op.matrix[8]);oct_salc.subspaces.push_back(sub);for(int i=0;i<3;++i){cov::NboSalcOrbital o;o.fragment_id="centre";o.energy_hartree=-.4;oct_salc.orbitals.push_back(o);}const auto multi=cov::ui::build_nbo_aomo_names(oct,f.data,&oct_salc);for(const auto& name:multi.salc)require(name.verified&&name.irrep=="T1u"&&name.ordinal==1,"true multidimensional SALC needs measured named signature and shared ordinal");
    auto repeated=oct_salc;repeated.orbitals.insert(repeated.orbitals.end(),oct_salc.orbitals.begin(),oct_salc.orbitals.end());auto& repeated_sub=repeated.subspaces[0];repeated_sub.dimension=6;repeated_sub.multiplicity=2;repeated_sub.orbital_indices={0,1,2,3,4,5};for(auto& chi:repeated_sub.characters)chi*=2;const auto copies=cov::ui::build_nbo_aomo_names(oct,f.data,&repeated);for(const auto& name:copies.salc)require(name.verified&&name.irrep=="T1u"&&name.ordinal==0,"repeated irreps retain known symmetry without inventing copy numbering");
    for(const auto& name:multi.salc)require(!name.partner_block_id.empty()&&name.partner_block_size==3&&name.point_group=="Oh"&&name.partner_block_id==multi.salc.front().partner_block_id,"verified triplet lost partner identity");
    for(const auto& name:copies.salc)require(name.partner_block_id.empty()&&name.partner_block_size==0,"unresolved repeated span invented a partner block");
    auto no_source=oct;no_source.derived_orbital_symmetry_assignments.clear();for(auto& mo:no_source.orbitals){mo.symmetry.clear();mo.symmetry_provenance=cov::DataProvenance::Unavailable;}const auto calculated=cov::ui::build_nbo_aomo_names(no_source,f.data,&oct_salc);for(const auto& name:calculated.canonical)require(name.verified&&name.irrep=="T1u"&&name.ordinal==1,"missing producer label must be calculated from actual finite-group coefficients");
    const auto independent_octahedral=cov::ui::canonical_mo_names(no_source);for(const auto& name:independent_octahedral->canonical)require(name.verified&&name.irrep=="T1u"&&name.ordinal==1,"standalone finite-group naming must close geometry generators without an attached SALC frame");
    auto lin=linear();const auto linear_before=lin.w.orbitals;const auto linear_names=cov::ui::build_nbo_aomo_names(lin.w,lin.data,&lin.salc);
    const auto standalone_linear=cov::ui::canonical_mo_names(lin.w);
    for(std::size_t i=0;i<lin.w.orbitals.size();++i)
        require(standalone_linear->canonical[i].verified&&
                standalone_linear->canonical[i].irrep==linear_names.canonical[i].irrep&&
                standalone_linear->canonical[i].ordinal==linear_names.canonical[i].ordinal,
                "standalone linear naming must resolve actual angular momenta without an NBO attachment");
    for(const auto& name:linear_names.canonical)require(name.verified,"unlabelled linear canonical must be classified");
    require(linear_names.canonical[0].irrep=="Sigma_g+"&&linear_names.canonical[0].ordinal==1&&linear_names.canonical[1].ordinal==2&&linear_names.canonical[5].ordinal==3,"hidden core must contribute to full-set Sigma numbering");
    require(linear_names.canonical[3].irrep=="Pi_u"&&linear_names.canonical[3].ordinal==1&&linear_names.canonical[4].ordinal==1&&linear_names.canonical[5].irrep=="Sigma_g+","accidental Pi plus Sigma energy coincidence must not merge irreps");
    require(linear_names.canonical[0].label.find("σ")!=std::string::npos&&linear_names.canonical[3].label.find("π")!=std::string::npos,"one-electron linear orbital labels need lowercase sigma/pi");
    require(linear_names.salc[2].irrep=="Sigma_u+"&&linear_names.salc[3].irrep=="Pi_u"&&linear_names.salc[4].ordinal==linear_names.salc[3].ordinal,"actual coefficients must resolve reducible central p into Sigma plus Pi");
    require(linear_names.salc[5].irrep=="Sigma_g+"&&linear_names.salc[6].irrep=="Sigma_g+"&&linear_names.salc[5].ordinal!=linear_names.salc[6].ordinal,"separate invariant radial copies need separate occurrences");
    auto linear_rotated=lin;const double c2=std::cos(.531),s2=std::sin(.531);rotate_general(linear_rotated,Mat{c2,0,s2,0,1,0,-s2,0,c2});compare(linear_names,cov::ui::build_nbo_aomo_names(linear_rotated.w,linear_rotated.data,&linear_rotated.salc));
    for(std::size_t i=0;i<lin.w.orbitals.size();++i)require(lin.w.orbitals[i].coefficients==linear_before[i].coefficients&&lin.w.orbitals[i].energy_hartree==linear_before[i].energy_hartree,"classification must not mutate producer coefficients or physical energies");
    auto broken=lin;broken.w.orbitals[4].energy_hartree+=.05;const auto split=cov::ui::build_nbo_aomo_names(broken.w,broken.data,&broken.salc);require(split.canonical[3].verified&&split.canonical[4].verified&&split.canonical[3].irrep=="Pi_u"&&broken.w.orbitals[4].energy_hartree==lin.w.orbitals[4].energy_hartree+.05,"spatial irrep evidence must not be confused with energy degeneracy or overwrite source eigenvalues");
    auto near_partner=lin;near_partner.w.orbitals[4].energy_hartree+=1.2669e-5;
    const auto near_names=cov::ui::build_nbo_aomo_names(near_partner.w,near_partner.data,&near_partner.salc);
    require(near_names.canonical[3].verified&&near_names.canonical[4].verified&&near_names.canonical[3].partner_block_id==near_names.canonical[4].partner_block_id,"measured symmetry partners must survive the old fixed energy boundary");
    const auto filtered=cov::ui::nbo_aomo_names_for_view(lin.w,linear_names,{3,4,5},{},nullptr,"test valence filter");
    require(filtered.canonical[5].ordinal==1&&filtered.canonical[5].complete_set_ordinal==3&&filtered.canonical[3].ordinal==filtered.canonical[4].ordinal&&filtered.canonical[3].ordinal==1,"display counters must start within each visible irrep, preserving complete-set identities and true partners");
    const auto partial=cov::ui::nbo_aomo_names_for_view(lin.w,linear_names,{3,5},{},nullptr,"incomplete filter");
    require(partial.canonical[3].verified&&partial.canonical[3].ordinal==0,"filtering one partner cannot invent a complete visible occurrence");
    auto subgroup=broken;subgroup.w.orbitals[3].symmetry="B1u";subgroup.w.orbitals[3].symmetry_provenance=cov::DataProvenance::Producer;cov::OrbitalSymmetrySourceRecord record;record.source_path="fixture.log";record.detected_group_context="Dinfh";record.abelian_group_context="D2h";record.orbital_indices={3};record.labels={"B1u"};subgroup.w.orbital_symmetry_source_records.push_back(record);const auto subgroup_names=cov::ui::build_nbo_aomo_names(subgroup.w,subgroup.data,&subgroup.salc);require(subgroup_names.canonical[3].verified&&subgroup_names.canonical[3].irrep=="Pi_u"&&subgroup.w.orbitals[3].symmetry=="B1u","derived full-group identity must preserve the producer subgroup literal independently");
    auto unresolved=lin.salc;unresolved.orbitals.clear();unresolved.subspaces.clear();for(int i=0;i<7;++i){cov::NboSalcOrbital o;o.fragment_id="copies";o.energy_hartree=i<4?-.4:i<6?.4:.6;unresolved.orbitals.push_back(o);}
    for(int block=0;block<3;++block){cov::NboSalcSubspace sub;sub.fragment_id="copies";sub.symmetry_verified=true;sub.orbital_indices=block==0?std::vector<std::size_t>{0,1,2,3}:block==1?std::vector<std::size_t>{4,5}:std::vector<std::size_t>{6};sub.dimension=sub.orbital_indices.size();for(const auto& op:unresolved.operations)sub.characters.push_back(block==2?1:(block==0?2:1)*(op.matrix[0]+op.matrix[4]));unresolved.subspaces.push_back(sub);}
    const auto unresolved_names=cov::ui::build_nbo_aomo_names(lin.w,lin.data,&unresolved);require(unresolved_names.salc[0].verified&&unresolved_names.salc[0].irrep=="Pi_u"&&unresolved_names.salc[0].ordinal==0&&unresolved_names.salc[4].verified&&unresolved_names.salc[4].ordinal==0&&unresolved_names.salc[6].irrep=="Sigma_g+"&&unresolved_names.salc[6].ordinal==1,"unresolved copies block only their own irrep numbering, retaining other known counts");
    for(std::size_t i=0;i<4;++i)require(unresolved_names.salc[i].partner_block_id.empty(),"unresolved copies gained false partner identity");
    require(unresolved_names.salc[4].ordinal==0&&unresolved_names.salc[4].partner_block_size==2&&unresolved_names.salc[4].point_group=="Dinfh"&&
        !unresolved_names.salc[4].partner_block_id.empty()&&unresolved_names.salc[4].partner_block_id==unresolved_names.salc[5].partner_block_id,
        "unknown global ordinal must not erase a verified pair");
    const auto high_copies=spectral_copies(2,3);const auto high_names=cov::ui::build_nbo_aomo_names(high_copies.w,high_copies.data,&high_copies.salc);require(high_names.salc[0].verified&&high_names.salc[0].ordinal==0&&high_names.salc[4].irrep=="Pi_u"&&high_names.salc[4].ordinal==1&&high_names.salc[0].detail.find("spectral bounds")!=std::string::npos,"certified high repeated spectrum must not block a low single occurrence");
    for(std::size_t i=0;i<4;++i)require(high_names.salc[i].status=="verified_isotypic_member"&&
        high_names.salc[i].representation_multiplicity==2&&high_names.salc[i].projection_residual&&
        high_names.salc[i].partner_block_id.empty(),"actual repeated isotypic source SALCs must retain projection evidence without copy ordinals or partners");
    const auto early_copies=spectral_copies(-3,-2);const auto early_names=cov::ui::build_nbo_aomo_names(early_copies.w,early_copies.data,&early_copies.salc);require(early_names.salc[4].ordinal==3,"known number of copies with certified lower spectrum must be counted");
    const auto crossing_copies=spectral_copies(-1,1);const auto crossing_names=cov::ui::build_nbo_aomo_names(crossing_copies.w,crossing_copies.data,&crossing_copies.salc);require(crossing_names.salc[4].verified&&crossing_names.salc[4].ordinal==0,"interleaved unresolved spectrum must still withhold the occurrence number");
    auto coupled=slightly_coupled_copies(1e-4);const auto coupled_before=coupled.w.orbitals;const auto recovered=cov::ui::build_nbo_aomo_names(coupled.w,coupled.data,&coupled.salc);
    for(const auto* rows:{&recovered.canonical,&recovered.salc}){
        require((*rows)[0].containing_irreps.size()>1&&(*rows)[0].containing_members.size()==5,"source copies fixture must exercise an actual connected reducible span");
        require((*rows)[0].ordinal==1&&(*rows)[1].ordinal==1&&(*rows)[3].ordinal==2&&(*rows)[4].ordinal==2&&(*rows)[7].ordinal==3&&(*rows)[8].ordinal==3,"verified independent copies must be counted once each, including later source occurrences");
        require((*rows)[0].partner_block_size==2&&(*rows)[0].partner_block_id==(*rows)[1].partner_block_id&&(*rows)[3].partner_block_id==(*rows)[4].partner_block_id&&(*rows)[0].partner_block_id!=(*rows)[3].partner_block_id,"two certified source E copies need distinct true partner blocks");
    }
    auto partially_mixed=slightly_coupled_copies(.02);const auto partial_copies=cov::ui::build_nbo_aomo_names(partially_mixed.w,partially_mixed.data,&partially_mixed.salc);
    for(const auto* rows:{&partial_copies.canonical,&partial_copies.salc}){
        require((*rows)[0].verified&&(*rows)[0].ordinal==1&&(*rows)[1].ordinal==1&&(*rows)[0].partner_block_size==2,"one genuine source copy must survive a mixed second copy in its containing span");
        require((*rows)[3].verified&&(*rows)[3].ordinal==0&&(*rows)[3].partner_block_id.empty()&&!(*rows)[4].verified&&!(*rows)[6].verified,"incomplete pure partners and mixed members must retain uncertainty");
        require((*rows)[7].ordinal==3&&(*rows)[8].ordinal==3,"remaining copy count must exclude the already named occurrence and retain exactly one unresolved earlier copy");
    }
    for(std::size_t i=0;i<coupled_before.size();++i)require(coupled.w.orbitals[i].coefficients==coupled_before[i].coefficients&&coupled.w.orbitals[i].energy_hartree==coupled_before[i].energy_hartree,"copy identification must not rotate source columns or change their energies");
    auto equal_copies=coupled;for(auto i:{0,1,3,4})equal_copies.w.orbitals[i].energy_hartree=-.4;
    const auto equal_names=cov::ui::build_nbo_aomo_names(equal_copies.w,equal_copies.data,&equal_copies.salc);
    require(equal_names.canonical[0].ordinal==0&&equal_names.canonical[3].ordinal==0&&equal_names.canonical[0].partner_block_size==2&&equal_names.canonical[7].ordinal==3,"independent source partner blocks do not resolve an equal-energy ordinal tie, but later total counts remain exact");
    auto incomplete_copies=high_copies;incomplete_copies.salc.links.pop_back();incomplete_copies.salc.links.erase(incomplete_copies.salc.links.begin());const auto incomplete_names=cov::ui::build_nbo_aomo_names(incomplete_copies.w,incomplete_copies.data,&incomplete_copies.salc);require(incomplete_names.salc[4].ordinal==0,"missing canonical projection evidence must not certify a spectral bound");
    auto missing_metric=lin;missing_metric.w.ao_overlap.clear();const auto absent=cov::ui::build_nbo_aomo_names(missing_metric.w,missing_metric.data,&missing_metric.salc);require(!absent.canonical[3].verified,"missing S must not be guessed");
    const auto hex=axial(6,true);const auto hex_names=cov::ui::build_nbo_aomo_names(hex.w,hex.data,&hex.salc);require(hex_names.canonical[0].verified&&hex_names.canonical[0].irrep=="E1u"&&hex_names.canonical[1].ordinal==1&&hex_names.canonical[2].irrep=="A2u","Dnh must not depend on the smaller central-metal group catalogue");
    const auto pyramid=axial(3,false);const auto pyramid_names=cov::ui::build_nbo_aomo_names(pyramid.w,pyramid.data,&pyramid.salc);require(pyramid_names.canonical[0].verified&&pyramid_names.canonical[0].irrep=="E"&&pyramid_names.canonical[1].ordinal==1&&pyramid_names.canonical[2].irrep=="A1","source-free Cnv A1 and E must follow actual rotation/reflection characters");
    auto nonisometric=pyramid;nonisometric.w.ao_overlap[4]=1.0001;
    for(auto& mo:nonisometric.w.orbitals)mo.coefficients[1]/=std::sqrt(1.0001);
    const auto rejected_metric=cov::ui::build_nbo_aomo_names(nonisometric.w,nonisometric.data,&nonisometric.salc);
    require(!rejected_metric.canonical[0].verified&&rejected_metric.canonical[0].status=="symmetry_action_not_isometric"&&
        rejected_metric.canonical[0].detail.find("transformed norm error=")!=std::string::npos,
        "nonisometric action must be rejected with its measured reason, not called subspace leakage");
    Fixture orthorhombic;
    orthorhombic.w.atoms={{"C",6,0,0,0},{"H",1,3,0,0},{"H",1,-3,0,0},{"H",1,0,1,0},{"H",1,0,-1,0},{"H",1,0,0,2},{"H",1,0,0,-2}};
    shell(orthorhombic.w,0,1,1);orthorhombic.w.ao_overlap={1,0,0,0,1,0,0,0,1};
    for(std::size_t i=0;i<3;++i){cov::MolecularOrbital mo;mo.coefficients.assign(3,0);mo.coefficients[i]=1;mo.energy_hartree=-.8+.2*double(i);orthorhombic.w.orbitals.push_back(mo);}
    const auto native_d2h=cov::derive_orbital_symmetry(orthorhombic.w);
    require(native_d2h.point_group=="D2h"&&native_d2h.orbitals_labelled==3,"native D2h classifier incomplete");
    const auto d2h_names=cov::ui::canonical_mo_names(orthorhombic.w);const char* d2h_labels[]={"B3u","B2u","B1u"};
    for(std::size_t i=0;i<3;++i)require(d2h_names->canonical[i].verified&&d2h_names->canonical[i].ordinal==1&&d2h_names->canonical[i].irrep==d2h_labels[i]&&orthorhombic.w.orbitals[i].symmetry==d2h_labels[i],"native/standalone D2h axis convention mismatch");
    std::cout<<"AO/MO names evidence and invariance smoke passed\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
