#include "cov/orbital_ui.hpp"

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

void provenance_strip(const Wavefunction& wf) {
    const auto sym=count_symmetry(wf);
    const auto occ=count_occupation(wf);
    ImGui::TextDisabled(
        "%s | symmetry P/D/U %zu/%zu/%zu | occupation P/D/U %zu/%zu/%zu",
        wavefunction_source_name(wf.source),
        sym.producer,sym.derived,sym.unavailable,
        occ.producer,occ.derived,occ.unavailable);

    ImGui::TextDisabled(
        "density %s | overlap %s | bond order %s%s%s",
        data_provenance_name(wf.total_density_provenance),
        data_provenance_name(wf.ao_overlap_provenance),
        data_provenance_name(wf.bond_order_provenance),
        wf.point_group_detected.empty()?"":" | point group ",
        wf.point_group_detected.empty()?"":wf.point_group_detected.c_str());

    if (!wf.enrichment_source.empty()) {
        ImGui::TextDisabled("Gaussian LOG/OUT enrichment attached");
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
    switch (language) {
        case Language::ChineseSimplified:
            return {"所选 MO 化学性质","化学价层轨道组","是","否",
                    "价层 AO 组成","原子对相互作用","轨道类型",
                    "成键性质","多中心族","离域 π 族","成员 MO","参与原子","参与电子",
                    "给体 / 受体方向",
                    "分析方法","MO 贡献",
                    "MO 贡献来自重叠布居；Mayer 为总密度原子对指数。",
                    "UND / 最小价层参考之外"};
        case Language::Japanese:
            return {"選択 MO の化学的性質","化学原子価軌道群","はい","いいえ",
                    "原子価 AO 組成","原子対相互作用","軌道型",
                    "結合性","多中心族","非局在化 π 族","構成 MO","参加原子","参加電子",
                    "供与体 / 受容体",
                    "解析法","MO 寄与",
                    "MO 寄与は重なり密度由来、Mayer は全密度の原子対指数。",
                    "UND / 最小原子価参照外"};
        case Language::French:
            return {"Caractère chimique de l’OM sélectionnée",
                    "Espace orbitalaire de valence chimique","oui","non",
                    "Composition AO de valence","Interactions par paire atomique",
                    "Type orbitalaire","Caractère liant",
                    "Famille multicentrique","Famille π délocalisée",
                    "OM membres","Atomes participants","Électrons participants",
                    "Direction donneur / accepteur","Méthode",
                    "Contribution OM",
                    "La contribution OM est fondée sur la population de recouvrement ; Mayer est l’indice de paire de densité totale.",
                    "UND / hors référence de valence minimale"};
        default:
            return {"Selected MO chemistry","Chemical-valence manifold","yes","no",
                    "Valence AO composition","Atom-pair interactions",
                    "Orbital family","Bonding role",
                    "Multicentre family","Delocalised π family",
                    "Member MOs","Participating atoms","Participating electrons",
                    "Donor / acceptor direction","Method","MO contribution",
                    "MO contribution is overlap-population based; Mayer is the total density-level pair index.",
                    "UND / outside minimal valence reference"};
    }
}

const char* family_symbol(const OrbitalAngularFamily family) {
    switch (family) {
        case OrbitalAngularFamily::Sigma: return "σ";
        case OrbitalAngularFamily::Pi: return "π";
        case OrbitalAngularFamily::Delta: return "δ";
        case OrbitalAngularFamily::Phi: return "φ";
        case OrbitalAngularFamily::NotApplicable: return "N/A";
        default: return "mixed";
    }
}

std::string channel_text(const OrbitalChannelDistribution& value) {
    if (value.status==ChemistryStatus::NotApplicable) return "N/A";
    if (value.status==ChemistryStatus::Undetermined ||
        value.status==ChemistryStatus::Unavailable) return "mixed(UND)";
    std::ostringstream out;
    out<<family_symbol(value.dominant);
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

const char* bonding_word(const OrbitalBondingRole role) {
    switch (role) {
        case OrbitalBondingRole::Bonding: return "bonding";
        case OrbitalBondingRole::Antibonding: return "antibonding";
        case OrbitalBondingRole::Nonbonding: return "nonbonding";
        case OrbitalBondingRole::NotApplicable: return "N/A";
        default: return "mixed";
    }
}

std::string bonding_text(const OrbitalBondingDistribution& value) {
    if (value.status==ChemistryStatus::NotApplicable) return "N/A";
    if (value.status==ChemistryStatus::Undetermined ||
        value.status==ChemistryStatus::Unavailable) return "mixed(UND)";
    std::ostringstream out;
    out<<bonding_word(value.dominant);
    if (value.status==ChemistryStatus::Percentages) {
        out<<" [bonding "<<std::lround(100.0*value.bonding)
           <<"% · antibonding "<<std::lround(100.0*value.antibonding)
           <<"% · nonbonding "<<std::lround(100.0*value.nonbonding);
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

void draw_channel_value(const char* label, const OrbitalChannelDistribution& value) {
    ImGui::Text("%s:", label);
    ImGui::SameLine();
    if (value.status == ChemistryStatus::NotApplicable) {
        ImGui::TextColored(text_colour(kUnavailableColour), "N/A");
        return;
    }
    if (value.status == ChemistryStatus::Undetermined || value.status == ChemistryStatus::Unavailable) {
        ImGui::TextColored(text_colour(kUnavailableColour), "mixed(UND)");
        return;
    }

    ImGui::TextColored(text_colour(channel_colour(value.dominant)), "%s", family_symbol(value.dominant));
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

void draw_bonding_value(const char* label, const OrbitalBondingDistribution& value) {
    ImGui::Text("%s:", label);
    ImGui::SameLine();
    if (value.status == ChemistryStatus::NotApplicable) {
        ImGui::TextColored(text_colour(kUnavailableColour), "N/A");
        return;
    }
    if (value.status == ChemistryStatus::Undetermined || value.status == ChemistryStatus::Unavailable) {
        ImGui::TextColored(text_colour(kUnavailableColour), "mixed(UND)");
        return;
    }

    ImGui::TextColored(text_colour(role_colour(value.dominant)), "%s", bonding_word(value.dominant));
    if (value.status != ChemistryStatus::Percentages) return;
    inline_plain(" [");
    const auto component = [](const char* name, const ImU32 colour, const double fraction, const bool separator) {
        if (separator) inline_plain(" · ");
        inline_text(name, colour);
        inline_plain(" ");
        inline_text(percent_text(fraction), kNumericColour);
    };
    component("bonding", kBondingColour, value.bonding, false);
    component("antibonding", kAntibondingColour, value.antibonding, true);
    component("nonbonding", kNonbondingColour, value.nonbonding, true);
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

std::string delocalised_descriptor(const OrbitalChemistry& chemistry) {
    if (chemistry.delocalised_family_id.empty()) return "N/A";
    std::string symbol="Π";
    symbol+=superscript_number(
        chemistry.delocalised_participating_atoms);
    symbol+=subscript_number(static_cast<int>(std::lround(
        chemistry.delocalised_participating_electrons)));
    return symbol+" · delocalised-pi";
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
        ImGui::TextColored(text_colour(kUnavailableColour), "UND — AO metric unavailable");
        return;
    }

    ImGui::Text("%s:", text.valence);
    ImGui::SameLine();
    ImGui::TextColored(text_colour(chemistry.valence_manifold ? kBondingColour : kUnavailableColour),
                       "%s", chemistry.valence_manifold ? text.yes : text.no);
    inline_plain(" · ");
    inline_text(percent_text(chemistry.valence_weight), kNumericColour);
    draw_channel_value(text.family, chemistry.channel);
    draw_bonding_value(text.bonding, chemistry.bonding);

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
                     delocalised_descriptor(chemistry),
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
                               channel_text(interaction.channel).c_str());
            inline_plain(" · ");
            inline_text(bonding_text(interaction.bonding), role_colour(interaction.bonding.dominant));
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
    ImGui::TextDisabled("%s: %s · confidence", text.method, chemistry.method.c_str());
    ImGui::SameLine();
    ImGui::TextColored(text_colour(kNumericColour), "%.0f%%", 100.0 * chemistry.confidence);
    if (!chemistry.note.empty()) {
        ImGui::TextDisabled("%s",chemistry.note.c_str());
    }
}

} // namespace

void draw_orbital_browser(const Wavefunction& wavefunction,
                          std::size_t selected_index,
                          OrbitalUIState& state,
                          Language language,
                          float ui_scale,
                          OrbitalUIActions& actions) {
    provenance_strip(wavefunction);
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
