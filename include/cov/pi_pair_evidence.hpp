#pragma once
#include <limits>
#include <string>
#include <vector>

namespace cov {
enum class LigandPiPrior { Unresolved, SigmaOnly, Donor, Acceptor, Ambiguous };
enum class PiPairDirection { Unresolved, Donor, Acceptor, Coupled, WeakNearNonbonding };
struct PiPartnerComponents {
    double energy_hartree=0;
    double metal_d=0;
    double ligand_p=0;
    double pi_fraction=0;
    double metal_ligand_overlap=0;
};
// A counterpart is a relation within one independently verified local Fock
// channel. Composition contrast alone cannot establish this relation.
struct PiPartnerChannelEvidence {
    std::string channel_id, canonical_fingerprint, spin;
    std::vector<std::size_t> lower_members, upper_members;
    bool same_operator_verified=false, occupations_verified=false;
    bool complete_membership_verified=false;
    std::string operator_kind;
    std::string direction="unresolved";
    std::string lower_character, upper_character;
    double lower_cross_fock_max_hartree=0, upper_cross_fock_min_hartree=0;
    double operator_error_hartree=0;
    double operator_tolerance_hartree=0;
};
struct PiPartnerAssessment {
    bool input_valid=false;
    bool accepted=false;
    PiPairDirection direction=PiPairDirection::Unresolved;
    LigandPiPrior prior=LigandPiPrior::Unresolved;
    // Verified channels rank by minimum endpoint cross-Fock magnitude (hartree).
    // Weak fallback ranking/support are dimensionless heuristics, never probabilities.
    double ranking_score=std::numeric_limits<double>::quiet_NaN();
    double support_score=std::numeric_limits<double>::quiet_NaN();
    double splitting_hartree=std::numeric_limits<double>::quiet_NaN();
    double composition_contrast=std::numeric_limits<double>::quiet_NaN();
    bool complementary_composition=false;
    bool opposite_overlap_signs=false;
    bool weak=false;
    std::string prior_relation="undetermined";
    std::string detail;
    PiPartnerComponents lower;
    PiPartnerComponents upper;
    PiPartnerChannelEvidence channel;
};
struct WeakCrystalFieldAssessment {
    bool input_valid=false;
    bool accepted=false;
    double splitting_hartree=std::numeric_limits<double>::quiet_NaN();
    double support_score=std::numeric_limits<double>::quiet_NaN();
    double split_threshold_hartree=std::numeric_limits<double>::quiet_NaN();
    double overlap_threshold=std::numeric_limits<double>::quiet_NaN();
    PiPartnerComponents first;
    PiPartnerComponents second;
    std::string detail;
};
// Numerical weak-field screen only. Callers must independently establish the
// compatible local axes, spin, d-shell irreps and all orbital identities.
[[nodiscard]] WeakCrystalFieldAssessment assess_weak_crystal_field(
    const PiPartnerComponents& first,const PiPartnerComponents& second,
    double split_threshold_hartree=0.020,double overlap_threshold=0.025);
[[nodiscard]] std::string weak_crystal_field_assessment_json(const WeakCrystalFieldAssessment&);
[[nodiscard]] PiPartnerAssessment assess_pi_partner(
    const PiPartnerComponents& lower,const PiPartnerComponents& upper,LigandPiPrior,
    double weak_split_hartree=0.020,double weak_overlap=0.025);
[[nodiscard]] PiPartnerAssessment assess_pi_channel_partner(
    const PiPartnerComponents& lower,const PiPartnerComponents& upper,
    LigandPiPrior,const PiPartnerChannelEvidence&);
[[nodiscard]] const char* ligand_pi_prior_name(LigandPiPrior) noexcept;
[[nodiscard]] const char* pi_pair_direction_name(PiPairDirection) noexcept;
[[nodiscard]] std::string pi_partner_assessment_json(const PiPartnerAssessment&);
}
