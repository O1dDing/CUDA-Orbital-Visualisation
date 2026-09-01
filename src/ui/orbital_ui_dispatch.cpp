#include "cov/orbital_ui.hpp"
#include "cov/orbital_ui_text.hpp"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

namespace cov::ui {

void draw_orbital_browser_legacy(const Wavefunction& wavefunction,
                                 std::size_t selected_index,
                                 OrbitalUIState& state,
                                 Language language,
                                 float ui_scale,
                                 OrbitalUIActions& actions);

void draw_energy_diagram_legacy(const Wavefunction& wavefunction,
                                std::size_t selected_index,
                                OrbitalUIState& state,
                                Language language,
                                float ui_scale,
                                OrbitalUIActions& actions);

namespace {

constexpr ImU32 kSigmaColour = IM_COL32(57, 210, 232, 255);
constexpr ImU32 kPiColour = IM_COL32(226, 93, 220, 255);
constexpr ImU32 kDeltaColour = IM_COL32(244, 155, 65, 255);
constexpr ImU32 kPhiColour = IM_COL32(55, 199, 170, 255);
constexpr ImU32 kBondingColour = IM_COL32(77, 218, 145, 255);
constexpr ImU32 kAntibondingColour = IM_COL32(244, 93, 105, 255);
constexpr ImU32 kNonbondingColour = IM_COL32(235, 181, 65, 255);
constexpr ImU32 kMulticentreColour = IM_COL32(238, 194, 89, 255);
constexpr ImU32 kNumericColour = IM_COL32(93, 174, 255, 255);
constexpr ImU32 kUnavailableColour = IM_COL32(128, 149, 177, 255);

ImVec4 text_colour(const ImU32 colour) {
    return ImGui::ColorConvertU32ToFloat4(colour);
}

void inline_text(const std::string& value, const ImU32 colour, const bool first = false) {
    if (!first) ImGui::SameLine(0.0f, 0.0f);
    ImGui::TextColored(text_colour(colour), "%s", value.c_str());
}

void inline_plain(const std::string& value, const bool first = false) {
    if (!first) ImGui::SameLine(0.0f, 0.0f);
    ImGui::TextUnformatted(value.c_str());
}

ImU32 channel_colour(const OrbitalAngularFamily family) {
    switch (family) {
        case OrbitalAngularFamily::Sigma: return kSigmaColour;
        case OrbitalAngularFamily::Pi: return kPiColour;
        case OrbitalAngularFamily::Delta: return kDeltaColour;
        case OrbitalAngularFamily::Phi: return kPhiColour;
        default: return kUnavailableColour;
    }
}

ImU32 role_colour(const OrbitalBondingRole role) {
    switch (role) {
        case OrbitalBondingRole::Bonding: return kBondingColour;
        case OrbitalBondingRole::Antibonding: return kAntibondingColour;
        case OrbitalBondingRole::Nonbonding: return kNonbondingColour;
        default: return kUnavailableColour;
    }
}

std::string percent_text(const double value) {
    return std::to_string(static_cast<int>(std::lround(100.0 * value))) + "%";
}

struct ProvenanceCounts {
    std::size_t producer=0;
    std::size_t derived=0;
    std::size_t unavailable=0;
};

ProvenanceCounts count_symmetry(const Wavefunction& wf) {
    ProvenanceCounts out;
    for (const auto& mo:wf.orbitals) {
        switch (mo.symmetry_provenance) {
            case DataProvenance::Producer: ++out.producer; break;
            case DataProvenance::Derived: ++out.derived; break;
            default: ++out.unavailable; break;
        }
    }
    return out;
}

ProvenanceCounts count_occupation(const Wavefunction& wf) {
    ProvenanceCounts out;
    for (const auto& mo:wf.orbitals) {
        switch (mo.occupation_provenance) {
            case DataProvenance::Producer: ++out.producer; break;
            case DataProvenance::Derived: ++out.derived; break;
            default: ++out.unavailable; break;
        }
    }
    return out;
}

void provenance_strip(const Wavefunction& wf, const Language language) {
    const auto sym=count_symmetry(wf);
    const auto occ=count_occupation(wf);
    ImGui::TextDisabled(
        "%s | %s [%s/%s/%s] %zu/%zu/%zu | %s [%s/%s/%s] %zu/%zu/%zu",
        localised_wavefunction_source(wf.source,language),
        orbital_tr(OrbitalText::ProvenanceSymmetry,language),
        orbital_tr(OrbitalText::Producer,language),
        orbital_tr(OrbitalText::Derived,language),
        orbital_tr(OrbitalText::Unavailable,language),
        sym.producer,sym.derived,sym.unavailable,
        orbital_tr(OrbitalText::ProvenanceOccupation,language),
        orbital_tr(OrbitalText::Producer,language),
        orbital_tr(OrbitalText::Derived,language),
        orbital_tr(OrbitalText::Unavailable,language),
        occ.producer,occ.derived,occ.unavailable);

    const std::string point_group_suffix=wf.point_group_detected.empty()
        ?std::string{}:" "+wf.point_group_detected;
    ImGui::TextDisabled(
        "%s %s | %s %s | %s %s%s%s%s",
        orbital_tr(OrbitalText::Density,language),
        localised_data_provenance(wf.total_density_provenance,language),
        orbital_tr(OrbitalText::Overlap,language),
        localised_data_provenance(wf.ao_overlap_provenance,language),
        orbital_tr(OrbitalText::BondOrder,language),
        localised_data_provenance(wf.bond_order_provenance,language),
        wf.point_group_detected.empty()?"":" | ",
        wf.point_group_detected.empty()?"":
            orbital_tr(OrbitalText::PointGroup,language),
        point_group_suffix.c_str());

    if (!wf.enrichment_source.empty()) {
        ImGui::TextDisabled("%s",
            orbital_tr(OrbitalText::GaussianEnrichmentAttached,language));
    }
    ImGui::Spacing();
}

struct ChemistryText {
    const char* heading;
    const char* valence;
    const char* yes;
    const char* no;
    const char* ao;
    const char* pair;
    const char* family;
    const char* bonding;
    const char* multicentre;
    const char* delocalised;
    const char* members;
    const char* atoms;
    const char* electrons;
    const char* donor;
    const char* method;
    const char* contribution;
    const char* explanation;
    const char* outside;
};

ChemistryText chemistry_text(const Language language) {
    return {
        orbital_tr(OrbitalText::SelectedMOChemistry,language),
        orbital_tr(OrbitalText::ChemicalValenceManifold,language),
        orbital_tr(OrbitalText::Yes,language),
        orbital_tr(OrbitalText::No,language),
        orbital_tr(OrbitalText::ValenceAOComposition,language),
        orbital_tr(OrbitalText::AtomPairInteractions,language),
        orbital_tr(OrbitalText::OrbitalFamily,language),
        orbital_tr(OrbitalText::BondingRole,language),
        orbital_tr(OrbitalText::MulticentreFamily,language),
        orbital_tr(OrbitalText::DelocalisedPiFamily,language),
        orbital_tr(OrbitalText::MemberMOs,language),
        orbital_tr(OrbitalText::ParticipatingAtoms,language),
        orbital_tr(OrbitalText::ParticipatingElectrons,language),
        orbital_tr(OrbitalText::DonorAcceptorDirection,language),
        orbital_tr(OrbitalText::AnalysisMethod,language),
        orbital_tr(OrbitalText::MOContribution,language),
        orbital_tr(OrbitalText::MOContributionExplanation,language),
        orbital_tr(OrbitalText::OutsideMinimalReference,language),
    };
}

const char* family_symbol(const OrbitalAngularFamily family, const Language language) {
    switch (family) {
        case OrbitalAngularFamily::Sigma: return "σ";
        case OrbitalAngularFamily::Pi: return "π";
        case OrbitalAngularFamily::Delta: return "δ";
        case OrbitalAngularFamily::Phi: return "φ";
        case OrbitalAngularFamily::NotApplicable: return "N/A";
        default: return orbital_tr(OrbitalText::Mixed,language);
    }
}

std::string channel_text(const OrbitalChannelDistribution& value,
                         const Language language) {
    if (value.status==ChemistryStatus::NotApplicable) return "N/A";
    if (value.status==ChemistryStatus::Undetermined ||
        value.status==ChemistryStatus::Unavailable) return orbital_tr(OrbitalText::MixedUnd,language);
    std::ostringstream out;
    out<<family_symbol(value.dominant,language);
    if (value.status==ChemistryStatus::Percentages) {
        out<<" [σ "<<std::lround(100.0*value.sigma)
           <<"% · π "<<std::lround(100.0*value.pi)
           <<"% · δ "<<std::lround(100.0*value.delta)
           <<"% · φ "<<std::lround(100.0*value.phi);
        if (value.undetermined>0.005) {
            out<<"% · UND "<<std::lround(100.0*value.undetermined);
        }
        out<<"%]";
    }
    return out.str();
}

const char* bonding_word(const OrbitalBondingRole role,
                         const Language language) {
    return localised_orbital_bonding_role(role,language);
}

std::string bonding_text(const OrbitalBondingDistribution& value,
                          const Language language) {
    if (value.status==ChemistryStatus::NotApplicable) return "N/A";
    if (value.status==ChemistryStatus::Undetermined ||
        value.status==ChemistryStatus::Unavailable) return orbital_tr(OrbitalText::MixedUnd,language);
    std::ostringstream out;
    out<<bonding_word(value.dominant,language);
    if (value.status==ChemistryStatus::Percentages) {
        out<<" ["<<orbital_tr(OrbitalText::Bonding,language)<<" "<<std::lround(100.0*value.bonding)
           <<"% · "<<orbital_tr(OrbitalText::Antibonding,language)<<" "<<std::lround(100.0*value.antibonding)
           <<"% · "<<orbital_tr(OrbitalText::Nonbonding,language)<<" "<<std::lround(100.0*value.nonbonding);
        if (value.undetermined>0.005) {
            out<<"% · UND "<<std::lround(100.0*value.undetermined);
        }
        out<<"%]";
    }
    return out.str();
}

void draw_label_value(const char* label, const std::string& value, const ImU32 colour) {
    ImGui::Text("%s:", label);
    ImGui::SameLine();
    ImGui::TextColored(text_colour(colour), "%s", value.c_str());
}

void draw_channel_value(const char* label, const OrbitalChannelDistribution& value,
                        const Language language) {
    ImGui::Text("%s:", label);
    ImGui::SameLine();
    if (value.status == ChemistryStatus::NotApplicable) {
        ImGui::TextColored(text_colour(kUnavailableColour), "N/A");
        return;
    }
    if (value.status == ChemistryStatus::Undetermined || value.status == ChemistryStatus::Unavailable) {
        ImGui::TextColored(text_colour(kUnavailableColour), "%s",
                           orbital_tr(OrbitalText::MixedUnd,language));
        return;
    }

    ImGui::TextColored(text_colour(channel_colour(value.dominant)), "%s", family_symbol(value.dominant,language));
    if (value.status != ChemistryStatus::Percentages) return;
    inline_plain(" [");
    const auto component = [](const char* symbol, const ImU32 colour, const double fraction, const bool separator) {
        if (separator) inline_plain(" · ");
        inline_text(symbol, colour);
        inline_plain(" ");
        inline_text(percent_text(fraction), kNumericColour);
    };
    component("σ", kSigmaColour, value.sigma, false);
    component("π", kPiColour, value.pi, true);
    component("δ", kDeltaColour, value.delta, true);
    component("φ", kPhiColour, value.phi, true);
    if (value.undetermined > 0.005) {
        inline_plain(" · ");
        inline_text("UND", kUnavailableColour);
        inline_plain(" ");
        inline_text(percent_text(value.undetermined), kNumericColour);
    }
    inline_plain("]");
}

void draw_bonding_value(const char* label, const OrbitalBondingDistribution& value,
                        const Language language) {
    ImGui::Text("%s:", label);
    ImGui::SameLine();
    if (value.status == ChemistryStatus::NotApplicable) {
        ImGui::TextColored(text_colour(kUnavailableColour), "N/A");
        return;
    }
    if (value.status == ChemistryStatus::Undetermined || value.status == ChemistryStatus::Unavailable) {
        ImGui::TextColored(text_colour(kUnavailableColour), "%s",
                           orbital_tr(OrbitalText::MixedUnd,language));
        return;
    }

    ImGui::TextColored(text_colour(role_colour(value.dominant)), "%s", bonding_word(value.dominant,language));
    if (value.status != ChemistryStatus::Percentages) return;
    inline_plain(" [");
    const auto component = [](const char* name, const ImU32 colour, const double fraction, const bool separator) {
        if (separator) inline_plain(" · ");
        inline_text(name, colour);
        inline_plain(" ");
        inline_text(percent_text(fraction), kNumericColour);
    };
    component(orbital_tr(OrbitalText::Bonding,language), kBondingColour, value.bonding, false);
    component(orbital_tr(OrbitalText::Antibonding,language), kAntibondingColour, value.antibonding, true);
    component(orbital_tr(OrbitalText::Nonbonding,language), kNonbondingColour, value.nonbonding, true);
    if (value.undetermined > 0.005) {
        inline_plain(" · ");
        inline_text("UND", kUnavailableColour);
        inline_plain(" ");
        inline_text(percent_text(value.undetermined), kNumericColour);
    }
    inline_plain("]");
}

std::string superscript_number(std::size_t value) {
    static constexpr const char* digits[]={"⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹"};
    const std::string plain=std::to_string(value);
    std::string out;
    for (const char c:plain) out+=digits[c-'0'];
    return out;
}

std::string subscript_number(const int value) {
    static constexpr const char* digits[]={"₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"};
    const std::string plain=std::to_string(std::max(0,value));
    std::string out;
    for (const char c:plain) out+=digits[c-'0'];
    return out;
}

std::string multicentre_descriptor(const OrbitalChemistry& chemistry) {
    if (chemistry.multicentre_label.empty()) return "N/A";
    std::ostringstream out;
    out<<chemistry.multicentre_label;
    if (chemistry.multicentre_participating_atoms>0u) {
        out<<" · "<<chemistry.multicentre_participating_atoms<<"c/"
           <<std::max(0LL,std::llround(
                  chemistry.multicentre_participating_electrons))
           <<"e";
    }
    return out.str();
}

std::string delocalised_descriptor(const OrbitalChemistry& chemistry,
                                  const Language language) {
    if (chemistry.delocalised_family_id.empty()) return "N/A";
    std::string symbol="Π";
    symbol+=superscript_number(
        chemistry.delocalised_participating_atoms);
    symbol+=subscript_number(static_cast<int>(std::lround(
        chemistry.delocalised_participating_electrons)));
    return symbol+" · "+orbital_tr(OrbitalText::DelocalisedPi,language);
}

std::string delocalised_orbital_members(const OrbitalChemistry& chemistry) {
    std::ostringstream result;
    for (const auto orbital_index : chemistry.delocalised_family_orbitals) {
        if (result.tellp() > 0) result << ", ";
        result << "MO" << orbital_index + 1u;
    }
    return result.str();
}

std::string atom_members(const Wavefunction& wf,
                         const std::vector<std::uint32_t>& atom_indices) {
    std::ostringstream result;
    for (const auto atom_index : atom_indices) {
        if (result.tellp() > 0) result << ", ";
        if (atom_index < wf.atoms.size()) result << wf.atoms[atom_index].symbol;
        result << atom_index + 1u;
    }
    return result.str();
}

void draw_participation(const Wavefunction& wf,
                        const ChemistryText& text,
                        const std::size_t atom_count,
                        const double electron_count,
                        const std::vector<std::uint32_t>& atom_indices,
                        const ImU32 colour) {
    if (!atom_indices.empty()) {
        draw_label_value(text.atoms,
                         std::to_string(atom_count)+" · "+
                             atom_members(wf,atom_indices),
                         colour);
    }
    if (atom_count>0u) {
        draw_label_value(text.electrons,
                         std::to_string(std::max(
                             0LL,std::llround(electron_count))),
                         colour);
    }
}

void draw_selected_chemistry(const Wavefunction& wf,
                             const std::size_t selected_index,
                             const Language language) {
    if (selected_index>=wf.orbitals.size()) return;
    const auto& chemistry=wf.orbitals[selected_index].chemistry;
    const auto text=chemistry_text(language);

    ImGui::SeparatorText(text.heading);
    if (!chemistry.available) {
        ImGui::TextColored(text_colour(kUnavailableColour), "%s",
                           orbital_tr(OrbitalText::AOMetricUnavailable,language));
        return;
    }

    ImGui::Text("%s:", text.valence);
    ImGui::SameLine();
    ImGui::TextColored(text_colour(chemistry.valence_manifold ? kBondingColour : kUnavailableColour),
                       "%s", chemistry.valence_manifold ? text.yes : text.no);
    inline_plain(" · ");
    inline_text(percent_text(chemistry.valence_weight), kNumericColour);
    draw_channel_value(text.family, chemistry.channel,language);
    draw_bonding_value(text.bonding, chemistry.bonding,language);

    const bool has_multicentre=!chemistry.multicentre_label.empty();
    draw_label_value(text.multicentre,
                     multicentre_descriptor(chemistry),
                     has_multicentre?kMulticentreColour:kUnavailableColour);
    if (has_multicentre) {
        draw_participation(
            wf,text,chemistry.multicentre_participating_atoms,
            chemistry.multicentre_participating_electrons,
            chemistry.multicentre_participating_atom_indices,
            kMulticentreColour);
    }

    const bool has_delocalised=!chemistry.delocalised_family_id.empty();
    draw_label_value(text.delocalised,
                     delocalised_descriptor(chemistry,language),
                     has_delocalised?kPiColour:kUnavailableColour);
    if (has_delocalised) {
        if (!chemistry.delocalised_family_orbitals.empty()) {
            draw_label_value(text.members,
                             delocalised_orbital_members(chemistry),
                             kPiColour);
        }
        draw_participation(
            wf,text,chemistry.delocalised_participating_atoms,
            chemistry.delocalised_participating_electrons,
            chemistry.delocalised_participating_atom_indices,
            kPiColour);
    }

    ImGui::TextDisabled("%s",text.ao);
    const std::size_t ao_count=std::min<std::size_t>(
        8u,chemistry.ao_contributions.size());
    for (std::size_t i=0;i<ao_count;++i) {
        const auto& contribution=chemistry.ao_contributions[i];
        ImGui::Bullet();
        ImGui::SameLine();
        ImGui::Text("%s:", contribution.label.c_str());
        ImGui::SameLine();
        ImGui::TextColored(text_colour(kNumericColour), "%.1f%%", 100.0 * contribution.weight);
    }
    if (chemistry.unresolved_weight>0.005) {
        ImGui::Bullet();
        ImGui::SameLine();
        ImGui::TextColored(text_colour(kUnavailableColour), "%s:", text.outside);
        ImGui::SameLine();
        ImGui::TextColored(text_colour(kNumericColour), "%.1f%%", 100.0 * chemistry.unresolved_weight);
    }

    if (!chemistry.interactions.empty()) {
        ImGui::TextDisabled("%s",text.pair);
        std::vector<const OrbitalPairInteraction*> interactions;
        interactions.reserve(chemistry.interactions.size());
        for (const auto& interaction:chemistry.interactions) {
            interactions.push_back(&interaction);
        }
        std::sort(interactions.begin(),interactions.end(),
            [](const auto* a,const auto* b) {
                return std::abs(a->overlap_character)>
                       std::abs(b->overlap_character);
            });
        const std::size_t count=std::min<std::size_t>(6u,interactions.size());
        for (std::size_t i=0;i<count;++i) {
            const auto& interaction=*interactions[i];
            ImGui::Bullet();
            ImGui::SameLine();
            ImGui::Text("%s–%s ·", interaction.atom_a_label.c_str(), interaction.atom_b_label.c_str());
            ImGui::SameLine();
            ImGui::TextColored(text_colour(channel_colour(interaction.channel.dominant)), "%s",
                               channel_text(interaction.channel,language).c_str());
            inline_plain(" · ");
            inline_text(bonding_text(interaction.bonding,language), role_colour(interaction.bonding.dominant));
            inline_plain(std::string(" · ") + text.contribution + " ");
            std::ostringstream contribution_value;
            contribution_value.setf(std::ios::fixed);
            contribution_value.precision(4);
            contribution_value << std::showpos << interaction.occupied_overlap_contribution;
            inline_text(contribution_value.str(), kNumericColour);
            inline_plain(" · Mayer ");
            std::ostringstream mayer_value;
            mayer_value.setf(std::ios::fixed);
            mayer_value.precision(4);
            mayer_value << interaction.total_mayer_index;
            inline_text(mayer_value.str(), kMulticentreColour);
        }
        ImGui::TextDisabled("%s",text.explanation);
    }

    draw_label_value(text.donor, chemistry.donor_acceptor,
                     chemistry.donor_acceptor == "UND" ? kUnavailableColour : kMulticentreColour);
    const std::string method=localised_chemistry_method(chemistry.method,language);
    ImGui::TextDisabled("%s: %s · %s", text.method, method.c_str(),
                       tr(Text::Confidence,language));
    ImGui::SameLine();
    ImGui::TextColored(text_colour(kNumericColour), "%.0f%%", 100.0 * chemistry.confidence);
    if (!chemistry.note.empty()) {
        const std::string note=localised_chemistry_note(chemistry.note,language);
        ImGui::TextDisabled("%s",note.c_str());
    }
}

} // namespace

void draw_orbital_browser(const Wavefunction& wavefunction,
                          std::size_t selected_index,
                          OrbitalUIState& state,
                          Language language,
                          float ui_scale,
                          OrbitalUIActions& actions) {
    provenance_strip(wavefunction,language);
    draw_orbital_browser_legacy(
        wavefunction,selected_index,state,language,ui_scale,actions);
}

void draw_energy_diagram(const Wavefunction& wavefunction,
                         std::size_t selected_index,
                         OrbitalUIState& state,
                         Language language,
                         float ui_scale,
                         OrbitalUIActions& actions) {
    draw_energy_diagram_legacy(
        wavefunction,selected_index,state,language,ui_scale,actions);
    draw_selected_chemistry(wavefunction,selected_index,language);
}

} // namespace cov::ui
