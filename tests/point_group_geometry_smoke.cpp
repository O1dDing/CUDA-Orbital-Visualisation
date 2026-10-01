#include "cov/molecular_point_group_frame.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

namespace {
cov::Wavefunction read(const std::string& file){
    std::ifstream input(std::string(COV_SOURCE_DIR)+"/tests/fixtures/point_groups/"+file+".xyz");
    std::size_t n=0;input>>n;std::string line;std::getline(input,line);std::getline(input,line);
    const std::map<std::string,unsigned> elements{{"H",1},{"B",5},{"C",6},{"O",8},{"Mg",12},{"Si",14}};
    cov::Wavefunction w;for(std::size_t i=0;i<n;++i){cov::Atom atom;input>>atom.symbol>>atom.x>>atom.y>>atom.z;
        if(!input||!elements.contains(atom.symbol))throw std::runtime_error("invalid public fixture "+file);
        atom.atomic_number=elements.at(atom.symbol);atom.x*=1.889726124626;atom.y*=1.889726124626;atom.z*=1.889726124626;w.atoms.push_back(atom);}
    if(!n)throw std::runtime_error("missing public fixture "+file);return w;
}
bool check(const cov::Wavefunction& w,const std::string& name,const std::string& expected){
    const auto geometry=cov::analyse_molecular_symmetry(w);const auto full=cov::complete_molecular_point_group(w,geometry);
    const auto table=cov::molecular_point_group_irreps(w,full);
    std::cout<<name<<": group="<<geometry.point_group<<" operations="<<full.operations.size()<<" expected="<<expected<<" table="<<table.valid<<" "<<table.reason<<'\n';
    return geometry.point_group==expected&&full.operations.size()==cov::finite_point_group_order(expected)&&table.valid;
}
}
int main(){try{
    bool okay=true;
    const std::pair<const char*,const char*> cases[]={{"allene-D2d","D2d"},{"ethane-staggered-D3d","D3d"},{"ethane-eclipsed-D3h","D3h"},{"neopentane-T","T"},{"neopentane-Td","Td"},{"mg-aqua-Th","Th"},{"octamethyl-POSS-O","O"},{"tropylium-D7h","D7h"},{"dodecaborate-Ih","Ih"},{"cyclopentadienyl-D5h","D5h"},{"corannulene-C5v","C5v"}};
    for(const auto& [file,group]:cases){auto w=read(file);okay=check(w,file,group)&&okay;
        const double c=std::cos(.371),s=std::sin(.371);for(auto& a:w.atoms){const double x=a.x,y=a.y;a.x=c*x-s*y+3.1;a.y=s*x+c*y-1.2;a.z+=.9;}
        std::reverse(w.atoms.begin(),w.atoms.end());okay=check(w,std::string(file)+" rotated/reordered",group)&&okay;}
    for(const auto& [file,group]:std::vector<std::pair<std::string,std::string>>{{"neopentane-T","Td"},{"mg-aqua-Th","Oh"},{"octamethyl-POSS-O","Oh"}}){auto w=read(file);
        w.atoms.erase(std::remove_if(w.atoms.begin(),w.atoms.end(),[](const auto& a){return a.atomic_number==1;}),w.atoms.end());
        okay=check(w,file+" explicit non-H analysis",group)&&okay;}
    auto perturbed=read("ethane-staggered-D3d");perturbed.atoms[0].x+=1e-6;okay=check(perturbed,"within recorded geometry tolerance","D3d")&&okay;
    perturbed.atoms[0].x+=.03;if(cov::analyse_molecular_symmetry(perturbed).point_group=="D3d"){std::cerr<<"above-tolerance perturbation silently symmetrized\n";okay=false;}
    return okay?0:1;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
