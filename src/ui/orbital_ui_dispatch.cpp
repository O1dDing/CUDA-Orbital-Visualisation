#include "cov/orbital_ui.hpp"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
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
                    "成键性质","多中心 / 离域族","给体 / 受体方向",
                    "分析方法","MO 贡献",
                    "MO 贡献来自重叠布居；Mayer 为总密度原子对指数。",
                    "UND / 最小价层参考之外"};
        case Language::Japanese:
            return {"選択 MO の化学的性質","化学原子価軌道群","はい","いいえ",
                    "原子価 AO 組成","原子対相互作用","軌道型",
                    "結合性","多中心 / 非局在化族","供与体 / 受容体",
                    "解析法","MO 寄与",
                    "MO 寄与は重なり密度由来、Mayer は全密度の原子対指数。",
                    "UND / 最小原子価参照外"};
        case Language::French:
            return {"Caractère chimique de l’OM sélectionnée",
                    "Espace orbitalaire de valence chimique","oui","non",
                    "Composition AO de valence","Interactions par paire atomique",
                    "Type orbitalaire","Caractère liant",
                    "Famille multicentrique / délocalisée",
                    "Direction donneur / accepteur","Méthode",
                    "Contribution OM",
                    "La contribution OM est fondée sur la population de recouvrement ; Mayer est l’indice de paire de densité totale.",
                    "UND / hors référence de valence minimale"};
        default:
            return {"Selected MO chemistry","Chemical-valence manifold","yes","no",
                    "Valence AO composition","Atom-pair interactions",
                    "Orbital family","Bonding role",
                    "Multicentre / delocalised family",
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

std::string family_descriptor(const OrbitalChemistry& chemistry) {
    if (chemistry.multicentre_label.empty()) return "UND";
    std::string symbol=family_symbol(chemistry.channel.dominant);
    symbol+=superscript_number(chemistry.participating_atoms);
    symbol+=subscript_number(static_cast<int>(
        std::lround(chemistry.participating_electrons)));
    return symbol+" · "+chemistry.multicentre_label;
}

void draw_selected_chemistry(const Wavefunction& wf,
                             const std::size_t selected_index,
                             const Language language) {
    if (selected_index>=wf.orbitals.size()) return;
    const auto& chemistry=wf.orbitals[selected_index].chemistry;
    const auto text=chemistry_text(language);

    ImGui::SeparatorText(text.heading);
    if (!chemistry.available) {
        ImGui::TextDisabled("UND — AO metric unavailable");
        return;
    }

    ImGui::Text("%s: %s · %.1f%%",
                text.valence,
                chemistry.valence_manifold?text.yes:text.no,
                100.0*chemistry.valence_weight);
    ImGui::Text("%s: %s",text.family,
                channel_text(chemistry.channel).c_str());
    ImGui::Text("%s: %s",text.bonding,
                bonding_text(chemistry.bonding).c_str());
    ImGui::Text("%s: %s",text.multicentre,
                family_descriptor(chemistry).c_str());

    ImGui::TextDisabled("%s",text.ao);
    const std::size_t ao_count=std::min<std::size_t>(
        8u,chemistry.ao_contributions.size());
    for (std::size_t i=0;i<ao_count;++i) {
        const auto& contribution=chemistry.ao_contributions[i];
        ImGui::BulletText("%s: %.1f%%",
                          contribution.label.c_str(),
                          100.0*contribution.weight);
    }
    if (chemistry.unresolved_weight>0.005) {
        ImGui::BulletText("%s: %.1f%%",
                          text.outside,
                          100.0*chemistry.unresolved_weight);
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
            ImGui::BulletText(
                "%s–%s · %s · %s · %s %+.4f · Mayer %.4f",
                interaction.atom_a_label.c_str(),
                interaction.atom_b_label.c_str(),
                channel_text(interaction.channel).c_str(),
                bonding_text(interaction.bonding).c_str(),
                text.contribution,
                interaction.occupied_overlap_contribution,
                interaction.total_mayer_index);
        }
        ImGui::TextDisabled("%s",text.explanation);
    }

    ImGui::Text("%s: %s",text.donor,
                chemistry.donor_acceptor.c_str());
    ImGui::TextDisabled("%s: %s · confidence %.0f%%",
                        text.method,chemistry.method.c_str(),
                        100.0*chemistry.confidence);
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
