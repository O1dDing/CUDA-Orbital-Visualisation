#include "cov/pi_pair_evidence.hpp"
#include <iostream>
#include <stdexcept>
int main(){using namespace cov;auto require=[](bool ok,const char* why){if(!ok)throw std::runtime_error(why);};
 // The original screenshot's scalar observations cannot establish direction.
 PiPartnerComponents deep{-.5003894956667,.0455684785,.9203470735,.9995669512,.0057445997};
 PiPartnerComponents frontier{-.275306531,.6091618954,.3248724218,.9981487152,.0272649687};
 require(!assess_pi_partner(deep,frontier,LigandPiPrior::Acceptor).accepted,"composition falsely certified");
 PiPartnerComponents upper{-.0182660992,.22,.71,.9997,-.0552};
 PiPartnerChannelEvidence e;e.channel_id="sample-independent:pi-antibonding";e.canonical_fingerprint="dataset-x";
 e.spin="alpha";e.lower_members={1,2,3};e.upper_members={6,7,8};e.same_operator_verified=true;e.complete_membership_verified=true;
 e.operator_kind="canonical-same-operator";e.occupations_verified=true;e.direction="centre_to_ligand";
 e.lower_character="bonding_mixing";e.upper_character="antibonding_mixing";
 e.lower_cross_fock_max_hartree=-.1;e.upper_cross_fock_min_hartree=.09;e.operator_error_hartree=1e-8;
 e.operator_tolerance_hartree=2e-5;
 const auto accepted=assess_pi_channel_partner(frontier,upper,LigandPiPrior::Donor,e);
 require(accepted.accepted&&accepted.direction==PiPairDirection::Acceptor,"channel evidence must outrank catalogue");
 auto missing=e;missing.same_operator_verified=false;require(!assess_pi_channel_partner(frontier,upper,LigandPiPrior::Acceptor,missing).accepted,"unverified operator");
 missing=e;missing.operator_error_hartree=1;require(!assess_pi_channel_partner(frontier,upper,LigandPiPrior::Acceptor,missing).accepted,"operator residual exceeds tolerance");
 missing=e;missing.operator_tolerance_hartree=0;require(!assess_pi_channel_partner(frontier,upper,LigandPiPrior::Acceptor,missing).accepted,"missing operator tolerance");
 missing=e;missing.lower_members={1,1,3};require(!assess_pi_channel_partner(frontier,upper,LigandPiPrior::Acceptor,missing).accepted,"duplicate members");
 missing=e;missing.upper_members={1,7,8};require(!assess_pi_channel_partner(frontier,upper,LigandPiPrior::Acceptor,missing).accepted,"overlap members");
 missing=e;missing.occupations_verified=false;const auto unresolved=assess_pi_channel_partner(frontier,upper,LigandPiPrior::Acceptor,missing);
 require(unresolved.accepted&&unresolved.direction==PiPairDirection::Coupled,"coupling without density must remain nondirectional");
 missing=e;missing.operator_kind="independently-verified-physical-spin-fock";
 require(assess_pi_channel_partner(frontier,upper,LigandPiPrior::Acceptor,missing).accepted,"physical operator relationship distinct from eigenvalue operator");
 const auto json=pi_partner_assessment_json(accepted);
 require(json.find("lower_cross_fock_max_hartree")!=std::string::npos&&json.find("canonical-group-energy-separation")!=std::string::npos,"missing reproducible evidence");
 require(json.find("\"ranking_unit\":\"hartree\"")!=std::string::npos,"channel ranking unit");
 std::cout<<"pi channel assessment checks passed\n";
}
