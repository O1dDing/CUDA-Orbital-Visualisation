#include "cov/pi_pair_evidence.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace cov {
const char* ligand_pi_prior_name(LigandPiPrior p) noexcept {
    switch(p) {
        case LigandPiPrior::SigmaOnly:return "sigma-only";
        case LigandPiPrior::Donor:return "donor";
        case LigandPiPrior::Acceptor:return "acceptor";
        case LigandPiPrior::Ambiguous:return "ambiguous";
        default:return "unresolved";
    }
}
const char* pi_pair_direction_name(PiPairDirection p) noexcept {
    switch(p) {
        case PiPairDirection::Donor:return "donor";
        case PiPairDirection::Acceptor:return "acceptor";
        case PiPairDirection::Coupled:return "coupled";
        case PiPairDirection::WeakNearNonbonding:return "weak-near-nonbonding";
        default:return "unresolved";
    }
}
PiPartnerAssessment assess_pi_partner(const PiPartnerComponents& a,const PiPartnerComponents& b,
    LigandPiPrior prior,double weak_split,double weak_overlap) {
    PiPartnerAssessment result;result.prior=prior;result.lower=a;result.upper=b;
    const auto valid=[](const auto& c) {
        return std::isfinite(c.energy_hartree) && std::isfinite(c.metal_d) && c.metal_d>=0 &&
               std::isfinite(c.ligand_p) && c.ligand_p>=0 && std::isfinite(c.pi_fraction) &&
               c.pi_fraction>=0 && std::isfinite(c.metal_ligand_overlap);
    };
    result.input_valid=valid(a)&&valid(b)&&std::isfinite(weak_split)&&weak_split>=0&&
        std::isfinite(weak_overlap)&&weak_overlap>=0&&b.energy_hartree>=a.energy_hartree;
    if(!result.input_valid){result.detail="invalid-input";return result;}
    result.splitting_hartree=b.energy_hartree-a.energy_hartree;
    const double low=a.ligand_p-a.metal_d,high=b.ligand_p-b.metal_d;
    const double donor=low-high,contrast=std::abs(donor);
    result.composition_contrast=donor;
    result.complementary_composition=low*high<0&&std::abs(low)>=.05&&std::abs(high)>=.05;
    result.opposite_overlap_signs=a.metal_ligand_overlap*b.metal_ligand_overlap<0&&
        std::abs(a.metal_ligand_overlap)>=weak_overlap&&std::abs(b.metal_ligand_overlap)>=weak_overlap;
    result.weak=result.splitting_hartree<=weak_split&&std::abs(a.metal_ligand_overlap)<=weak_overlap&&
        std::abs(b.metal_ligand_overlap)<=weak_overlap;
    if(!result.weak) {
        result.detail="composition-only-candidate-needs-verified-local-channel";
        return result;
    }
    if(a.pi_fraction<.60||b.pi_fraction<.60||a.metal_d+a.ligand_p<.18||b.metal_d+b.ligand_p<.18||
       result.splitting_hartree>1.50||(!result.weak&&!result.complementary_composition)||
       (!result.weak&&!result.opposite_overlap_signs&&contrast<.18)||
       (result.weak&&std::max(a.metal_d,b.metal_d)<.08)) {
        result.detail="orbital-evidence-insufficient";return result;
    }
    const double opposite=std::max(0.0,-a.metal_ligand_overlap*b.metal_ligand_overlap);
    const double quality=std::min(a.pi_fraction,b.pi_fraction);
    result.ranking_score=2*quality+2*contrast+2*std::sqrt(opposite)-.08*result.splitting_hartree;
    if(result.ranking_score<.75){result.detail="pair-ranking-support-insufficient";return result;}
    result.accepted=true;
    if(result.weak)result.direction=PiPairDirection::WeakNearNonbonding;
    else if(2*donor>=.15)result.direction=PiPairDirection::Donor;
    else if(-2*donor>=.15)result.direction=PiPairDirection::Acceptor;
    else result.direction=PiPairDirection::Coupled;
    result.support_score=std::clamp(.45*quality+.35*std::min(1.0,contrast)+
        .20*std::min(1.0,std::sqrt(opposite)/.05),0.0,1.0);
    // The catalogue is context. It cannot veto the actual canonical-pair evidence
    // or select a direction contradicted by the observed composition.
    if(prior==LigandPiPrior::Unresolved||prior==LigandPiPrior::Ambiguous||result.weak)
        result.prior_relation="undetermined";
    else if((prior==LigandPiPrior::Donor&&result.direction==PiPairDirection::Donor)||
            (prior==LigandPiPrior::Acceptor&&result.direction==PiPairDirection::Acceptor))
        result.prior_relation="consistent";
    else result.prior_relation="contradicted";
    result.detail="canonical-pair-interpretation-from-orbital-evidence";
    return result;
}
PiPartnerAssessment assess_pi_channel_partner(const PiPartnerComponents& a,
    const PiPartnerComponents& b,LigandPiPrior prior,const PiPartnerChannelEvidence& channel) {
    // Reuse finite-input validation and descriptive quantities, not its
    // composition-only acceptance decision or heuristic direction.
    auto result=assess_pi_partner(a,b,prior);
    result.channel=channel;result.accepted=false;result.direction=PiPairDirection::Unresolved;
    result.prior_relation="undetermined";
    if(!result.input_valid)return result;
    if(channel.channel_id.empty()||channel.canonical_fingerprint.empty()||
       channel.spin.empty()||channel.lower_members.empty()||channel.upper_members.empty()||
       !channel.same_operator_verified||!channel.complete_membership_verified||channel.operator_kind.empty()||
       !std::isfinite(channel.operator_error_hartree)||
       channel.operator_error_hartree<0||!std::isfinite(channel.operator_tolerance_hartree)||
       channel.operator_tolerance_hartree<=0||channel.operator_error_hartree>channel.operator_tolerance_hartree||
       !std::isfinite(channel.lower_cross_fock_max_hartree)||
       !std::isfinite(channel.upper_cross_fock_min_hartree)) {
        result.detail="verified-channel-identity-or-operator-unavailable";return result;
    }
    for(const auto* members:{&channel.lower_members,&channel.upper_members}) {
        auto copy=*members;std::sort(copy.begin(),copy.end());
        if(std::adjacent_find(copy.begin(),copy.end())!=copy.end()) {
            result.detail="duplicate-counterpart-member-identity";return result;
        }
    }
    for(const auto i:channel.lower_members)
        if(std::find(channel.upper_members.begin(),channel.upper_members.end(),i)!=channel.upper_members.end()) {
            result.detail="overlapping-counterpart-members";return result;
        }
    if(channel.lower_character!="bonding_mixing"||channel.upper_character!="antibonding_mixing"||
       !(channel.lower_cross_fock_max_hartree<0)||!(channel.upper_cross_fock_min_hartree>0)) {
        result.detail="same-channel-counterpart-cross-fock-signs-unresolved";return result;
    }
    result.accepted=true;result.weak=false;result.direction=PiPairDirection::Coupled;
    if(channel.occupations_verified&&channel.direction=="centre_to_ligand")result.direction=PiPairDirection::Acceptor;
    if(channel.occupations_verified&&channel.direction=="ligand_to_centre")result.direction=PiPairDirection::Donor;
    result.ranking_score=std::min(-channel.lower_cross_fock_max_hartree,channel.upper_cross_fock_min_hartree);
    result.support_score=std::numeric_limits<double>::quiet_NaN();
    if(result.direction==PiPairDirection::Acceptor||result.direction==PiPairDirection::Donor) {
        if(prior==LigandPiPrior::Acceptor||prior==LigandPiPrior::Donor)
            result.prior_relation=(prior==LigandPiPrior::Acceptor)==(result.direction==PiPairDirection::Acceptor)
                ?"consistent":"contradicted";
    }
    result.detail="shared-verified-local-channel-canonical-energy-separation-not-isolated-two-level-splitting";
    return result;
}
namespace {
void number(std::ostream& out,double v){if(std::isfinite(v))out<<v;else out<<"null";}
void json_string(std::ostream& out,const std::string& value) {
    static constexpr char hex[]="0123456789abcdef";
    out << '"';
    for(unsigned char c:value) {
        if(c=='"'||c=='\\')out << '\\' << static_cast<char>(c);
        else if(c<0x20)out << "\\u00" << hex[c>>4] << hex[c&15];
        else out << static_cast<char>(c);
    }
    out << '"';
}
void components(std::ostream& out,const PiPartnerComponents& v){
    out<<"{\"energy_hartree\":";number(out,v.energy_hartree);
    out<<",\"metal_d\":";number(out,v.metal_d);out<<",\"ligand_p\":";number(out,v.ligand_p);
    out<<",\"pi_fraction\":";number(out,v.pi_fraction);out<<",\"metal_ligand_overlap\":";number(out,v.metal_ligand_overlap);out<<'}';
}
}
std::string pi_partner_assessment_json(const PiPartnerAssessment& v){
    std::ostringstream out;out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"{\"input_valid\":"<<(v.input_valid?"true":"false")<<",\"accepted\":"<<(v.accepted?"true":"false")
       <<",\"direction\":\""<<pi_pair_direction_name(v.direction)<<"\",\"catalogue_prior\":\""<<ligand_pi_prior_name(v.prior)
       <<"\",\"prior_relation\":";json_string(out,v.prior_relation);out<<",\"detail\":";json_string(out,v.detail);
    out<<",\"ranking_score\":";number(out,v.ranking_score);out<<",\"support_score\":";number(out,v.support_score);
    const bool local_channel=!v.channel.channel_id.empty();
    out<<",\"ranking_measure\":\""<<(local_channel?"minimum-endpoint-cross-Fock-magnitude":"weak-composition-heuristic")
       <<"\",\"ranking_unit\":\""<<(local_channel?"hartree":"dimensionless")
       <<"\",\"score_meaning\":\""<<(local_channel?"cross-Fock-strength-not-probability":"heuristic-support-not-probability")
       <<"\",\"splitting_hartree\":";number(out,v.splitting_hartree);
    out<<",\"composition_contrast\":";number(out,v.composition_contrast);
    out<<",\"complementary_composition\":"<<(v.complementary_composition?"true":"false")
       <<",\"opposite_overlap_signs\":"<<(v.opposite_overlap_signs?"true":"false")<<",\"weak\":"<<(v.weak?"true":"false")<<",\"lower\":";
    components(out,v.lower);out<<",\"upper\":";components(out,v.upper);
    out<<",\"channel_id\":";json_string(out,v.channel.channel_id);
    out<<",\"canonical_fingerprint\":";json_string(out,v.channel.canonical_fingerprint);
    out<<",\"channel_spin\":";json_string(out,v.channel.spin);
    out<<",\"operator_kind\":";json_string(out,v.channel.operator_kind);
    out<<",\"channel_direction\":";json_string(out,v.channel.direction);
    out<<",\"lower_character\":";json_string(out,v.channel.lower_character);
    out<<",\"upper_character\":";json_string(out,v.channel.upper_character);
    out<<",\"lower_cross_fock_max_hartree\":";number(out,v.channel.lower_cross_fock_max_hartree);
    out<<",\"upper_cross_fock_min_hartree\":";number(out,v.channel.upper_cross_fock_min_hartree);
    out<<",\"operator_error_hartree\":";number(out,v.channel.operator_error_hartree);
    out<<",\"operator_tolerance_hartree\":";number(out,v.channel.operator_tolerance_hartree);
    out<<",\"lower_members\":[";for(std::size_t i=0;i<v.channel.lower_members.size();++i){if(i)out<<',';out<<v.channel.lower_members[i];}
    out<<"],\"upper_members\":[";for(std::size_t i=0;i<v.channel.upper_members.size();++i){if(i)out<<',';out<<v.channel.upper_members[i];}
    out<<"],\"same_operator_verified\":"<<(v.channel.same_operator_verified?"true":"false")
       <<",\"occupations_verified\":"<<(v.channel.occupations_verified?"true":"false")
       <<",\"complete_membership_verified\":"<<(v.channel.complete_membership_verified?"true":"false")
       <<",\"energy_semantics\":\"canonical-group-energy-separation\"}";return out.str();
}
WeakCrystalFieldAssessment assess_weak_crystal_field(const PiPartnerComponents& a,
    const PiPartnerComponents& b,double split_threshold,double overlap_threshold) {
    WeakCrystalFieldAssessment result;result.first=a;result.second=b;
    result.split_threshold_hartree=split_threshold;result.overlap_threshold=overlap_threshold;
    const auto valid=[](const auto& c) {
        return std::isfinite(c.energy_hartree)&&std::isfinite(c.metal_d)&&c.metal_d>=0&&
            std::isfinite(c.metal_ligand_overlap);
    };
    result.input_valid=valid(a)&&valid(b)&&std::isfinite(split_threshold)&&split_threshold>=0&&
        std::isfinite(overlap_threshold)&&overlap_threshold>=0;
    if(!result.input_valid){result.detail="invalid-input";return result;}
    result.splitting_hartree=std::abs(b.energy_hartree-a.energy_hartree);
    if(a.metal_d<.60||b.metal_d<.60||result.splitting_hartree>split_threshold||
        std::abs(a.metal_ligand_overlap)>overlap_threshold||
        std::abs(b.metal_ligand_overlap)>overlap_threshold) {
        result.detail="weak-local-d-evidence-insufficient";return result;
    }
    result.support_score=std::clamp(split_threshold>0?
        1-result.splitting_hartree/split_threshold:1,0.0,1.0);
    result.accepted=result.support_score>=.20;
    result.detail=result.accepted?"weak-local-d-separation":"weak-gap-support-insufficient";
    return result;
}
std::string weak_crystal_field_assessment_json(const WeakCrystalFieldAssessment& v) {
    std::ostringstream out;out<<std::setprecision(std::numeric_limits<double>::max_digits10);
    out<<"{\"input_valid\":"<<(v.input_valid?"true":"false")<<",\"accepted\":"<<(v.accepted?"true":"false");
    out<<",\"splitting_hartree\":";number(out,v.splitting_hartree);
    out<<",\"support_score\":";number(out,v.support_score);
    out<<",\"score_meaning\":\"heuristic-support-not-probability\",\"split_threshold_hartree\":";
    number(out,v.split_threshold_hartree);out<<",\"overlap_threshold\":";number(out,v.overlap_threshold);
    out<<",\"detail\":";json_string(out,v.detail);out<<",\"first\":";components(out,v.first);
    out<<",\"second\":";components(out,v.second);out<<'}';return out.str();
}
}
