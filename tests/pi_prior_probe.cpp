#include "cov/pi_pair_evidence.hpp"
#include <iostream>
#include <limits>
int main(){
    using namespace cov;
    for(int scenario=0;scenario<7;++scenario)for(int prior=0;prior<5;++prior){
        PiPartnerComponents a{-.5,.1,.8,.9,.1},b{0,.8,.1,.9,-.1};
        if(scenario==1){a.metal_d=.8;a.ligand_p=.1;b.metal_d=.1;b.ligand_p=.8;}
        if(scenario==2){b.energy_hartree=-.495;a.metal_ligand_overlap=.001;b.metal_ligand_overlap=.002;}
        if(scenario==3)a.pi_fraction=.59;
        if(scenario==4){b.metal_d=.1;b.ligand_p=.8;}
        if(scenario==5)b.energy_hartree=2;
        if(scenario==6)a.pi_fraction=std::numeric_limits<double>::quiet_NaN();
        std::cout<<"{\"scenario\":"<<scenario<<",\"prior\":"<<prior<<",\"assessment\":"
            <<pi_partner_assessment_json(assess_pi_partner(a,b,static_cast<LigandPiPrior>(prior)))<<"}\n";
    }
}
