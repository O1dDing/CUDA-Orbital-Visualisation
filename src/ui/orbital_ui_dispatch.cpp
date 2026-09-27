#include "cov/orbital_ui.hpp"
#include "cov/nbo_ui.hpp"
#include "cov/orbital_ui_text.hpp"
#include "cov/numerical_diagnostics.hpp"
#include "cov/overlap.hpp"
#include "cov/orbital_chemistry_summary.hpp"
#include "cov/validation.hpp"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
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

void continue_text(const std::string& value, const float spacing = 0.0f) {
    ImGui::SameLine(0.0f, spacing);
    // Wrapped cards have no room for a new colored run at the right edge.
    // Move the whole run to the next line before ImGui wraps its individual
    // words into the remaining sliver. Longer runs can then use the full width.
    if (ImGui::CalcTextSize(value.c_str()).x > ImGui::GetContentRegionAvail().x) {
        ImGui::NewLine();
    }
}

void inline_text(const std::string& value, const ImU32 colour, const bool first = false) {
    if (!first) continue_text(value);
    ImGui::TextColored(text_colour(colour), "%s", value.c_str());
}

void inline_plain(const std::string& value, const bool first = false) {
    if (!first) continue_text(value);
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

    const bool chinese=language==Language::ChineseSimplified;
    const auto density_status=[&](NumericalStatus status) {
        switch(status) {
            case NumericalStatus::Available:return chinese?"可用":"available";
            case NumericalStatus::MissingInput:return chinese?"信息不足":"missing information";
            case NumericalStatus::InvalidInput:return chinese?"输入无效":"invalid input";
            case NumericalStatus::Failed:return chinese?"计算失败":"calculation failed";
            default:return chinese?"未计算":"not computed";
        }
    };
    ImGui::TextWrapped(chinese?"总密度：%s（%s）；自旋密度：%s（%s）":
        "Total density: %s (%s); spin density: %s (%s)",
        density_status(wf.total_density_diagnostics.status),
        localised_data_provenance(wf.total_density_provenance,language),
        density_status(wf.spin_density_diagnostics.status),
        localised_data_provenance(wf.spin_density_provenance,language));
    if(wf.spin_density_diagnostics.occupation_model==OrbitalOccupationModel::SharedIntegerDeterminant &&
       wf.spin_density_diagnostics.status==NumericalStatus::Available) {
        ImGui::TextWrapped("%s",chinese?
            "自旋密度按共享整数占据解释：双占据成对，单占据归入 α 自旋。":
            "Spin density assumes shared integer occupations: double occupation is paired; single occupation is alpha.");
    } else if(wf.total_density_diagnostics.status==NumericalStatus::Available &&
              wf.spin_density_diagnostics.status==NumericalStatus::MissingInput) {
        ImGui::TextWrapped("%s",chinese?
            "输入可以确定总密度，但不能确定自旋分配；完整 Mayer 键级暂不可用。":
            "The input determines total density but not its spin partition; full Mayer bond orders are unavailable.");
    }

    if (!wf.enrichment_source.empty()) {
        ImGui::TextDisabled("%s",
            orbital_tr(OrbitalText::GaussianEnrichmentAttached,language));
    }
    if (wf.ao_metric_diagnostics.status!=NumericalStatus::NotComputed) {
        if (wf.ao_metric_diagnostics.status!=NumericalStatus::Available ||
            (!wf.orbitals.empty() && !orbital_metric_usable(wf))) {
            ImGui::TextWrapped("%s",chinese
                ?"基组与轨道数值检查未通过；成键和派生轨道标签暂不可用。原始轨道数据保留。"
                :"Basis/orbital numerical checks failed; bonding and derived orbital labels are unavailable. Original orbital data are retained.");
        }
        if (wf.producer_ao_overlap_basis_status==NumericalStatus::InvalidInput) {
            ImGui::TextWrapped("%s",chinese
                ?"输入中的重叠矩阵与基组不一致；分析使用独立基组积分，原矩阵保留供核查。"
                :"The input overlap matrix is inconsistent with the basis; analysis uses independent basis integrals and retains the original matrix for inspection.");
        }
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
    out<<orbital_channel_fraction_summary(value);
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
    out<<orbital_bonding_fraction_summary(value,
        orbital_tr(OrbitalText::Bonding,language),
        orbital_tr(OrbitalText::Antibonding,language),
        orbital_tr(OrbitalText::Nonbonding,language));
    return out.str();
}

void draw_label_value(const char* label, const std::string& value, const ImU32 colour) {
    ImGui::Text("%s:", label);
    continue_text(value,ImGui::GetStyle().ItemSpacing.x);
    ImGui::TextColored(text_colour(colour), "%s", value.c_str());
}

void draw_channel_value(const char* label, const OrbitalChannelDistribution& value,
                        const Language language) {
    ImGui::Text("%s:", label);
    continue_text(channel_text(value,language),ImGui::GetStyle().ItemSpacing.x);
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
    continue_text(bonding_text(value,language),ImGui::GetStyle().ItemSpacing.x);
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

const char* unresolved_canonical_assignment(const Language language) {
    switch(language) {
    case Language::ChineseSimplified:
        return "此方法未建立正则轨道归属；请查看上方关系证据";
    case Language::Japanese:
        return "この方法では正準軌道への帰属は未確定です。上の関係証拠を参照";
    case Language::French:
        return "Attribution à l'orbitale canonique non établie par cette méthode ; voir les relations ci-dessus";
    default:
        return "No canonical assignment established by this method; see relation evidence above";
    }
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

void draw_routed_subspaces(const RoutedAnalysis& routed,const std::size_t selected_index) {
    if(selected_index<routed.mo_relations.size()){
        const auto& relations=routed.mo_relations[selected_index];
        ImGui::SeparatorText("Related localized subspaces");
        if(relations.available()){
            ImGui::TextDisabled("%zu actual projected relations; source records retained in export",
                                relations.value->size());
            std::map<std::string,const RoutedLocalRelation*> registry;
            for(const auto& row:routed.local_relation_registry)
                registry.emplace(row.id,&row);
            std::vector<const RoutedRelationProjection*> ranked;
            for(const auto& row:*relations.value)ranked.push_back(&row);
            const auto strength=[](const RoutedRelationProjection* row){
                double best=0;for(const auto& weight:row->canonical_projection_weights)
                    if(weight)best=std::max(best,*weight);
                return best;
            };
            std::sort(ranked.begin(),ranked.end(),[&](auto a,auto b){
                return strength(a)>strength(b);
            });
            const auto draw_relation=[&](const RoutedRelationProjection& projection){
                const auto found=registry.find(projection.relation_id);
                if(found==registry.end())return;
                const auto& relation=*found->second;
                ImGui::BulletText("%s",relation.label.c_str());
                for(std::size_t k=0;k<relation.orbitals.size();++k)
                    if(k<projection.canonical_projection_weights.size() &&
                       projection.canonical_projection_weights[k])
                        ImGui::TextDisabled("%s %zu: %.2f%% canonical projection",
                            nbo_orbital_kind_name(relation.orbitals[k].kind),
                            relation.orbitals[k].index+1,
                            100* *projection.canonical_projection_weights[k]);
                if(relation.kind=="donor_acceptor")
                    ImGui::TextDisabled("Localized E(2); not an allocated canonical MO energy");
            };
            for(std::size_t k=0;k<std::min<std::size_t>(6,ranked.size());++k)
                draw_relation(*ranked[k]);
            if(ranked.size()>6 && ImGui::TreeNode("All related source records")){
                ImGui::BeginChild("##related_local_records",ImVec2(0,240),false);
                ImGuiListClipper clipper;clipper.Begin(static_cast<int>(ranked.size()));
                while(clipper.Step())for(int k=clipper.DisplayStart;k<clipper.DisplayEnd;++k)
                    draw_relation(*ranked[static_cast<std::size_t>(k)]);
                ImGui::EndChild();ImGui::TreePop();
            }
        }else ImGui::TextDisabled("Localized relations: %s (%s)",
            routed_status_name(relations.status),relations.reason.c_str());
    }
    bool mapped_pi=false,available_pi=false;
    for(const auto& coupling:routed.pi_couplings){
        if(!coupling.available()){
            ImGui::TextDisabled("Same-operator π coupling: %s (%s)",
                routed_status_name(coupling.status),coupling.reason.c_str());
            continue;
        }
        available_pi=true;
        const auto& record=*coupling.value;
        for(const auto& group:record.groups)
            if(std::find(group.members.begin(),group.members.end(),selected_index)!=group.members.end()){
                mapped_pi=true;
                ImGui::SeparatorText("Same-operator π subspace");
                ImGui::Text("%s · %s · %zu coupled dimensions",
                    group.character.c_str(),record.direction.c_str(),record.coupled_rank);
                ImGui::TextDisabled("Cross-Fock range %.5f to %.5f Ha; partition support %.1f%% / %.1f%%",
                    group.cross_fock_min_hartree,group.cross_fock_max_hartree,
                    100*group.centre_weight,100*group.ligand_weight);
                ImGui::TextWrapped("%s",record.direction_evidence.c_str());
            }
    }
    if(routed.pi_couplings.empty())
        ImGui::TextDisabled("Same-operator π coupling: not analysed");
    else if(available_pi && !mapped_pi)
        ImGui::TextDisabled("No verified π coupling subspace maps to this canonical MO");
}

void draw_selected_chemistry(const Wavefunction& wf,
                             const std::size_t selected_index,
                             const RoutedAnalysis* routed,
                             const Language language) {
    if (selected_index>=wf.orbitals.size()) return;
    const auto& chemistry=wf.orbitals[selected_index].chemistry;
    const auto text=chemistry_text(language);

    const bool show_chemistry=ImGui::CollapsingHeader(
        (std::string(text.heading)+"##diagram.chemistry").c_str());
    validation::item("diagram.chemistry");
    if(!show_chemistry)return;
    if(routed && selected_index<routed->mo_composition.size()){
        const auto& composition=routed->mo_composition[selected_index];
        ImGui::TextDisabled("Composition route: %s via %s · %s",
            routed_status_name(composition.status),
            routed_provider_name(composition.provider),composition.method.c_str());
        if(!composition.reason.empty())ImGui::TextWrapped("%s",composition.reason.c_str());
        if(!composition.fallback_reason.empty())
            ImGui::TextWrapped("NBO fallback: %s",composition.fallback_reason.c_str());
    }
    if(routed)draw_routed_subspaces(*routed,selected_index);
    if(routed && selected_index<routed->mo_composition.size()){
        const auto& composition=routed->mo_composition[selected_index];
        if(composition.available() && composition.provider==RoutedProvider::Nbo){
            const auto& value=*composition.value;
            ImGui::TextWrapped("%s: %s",text.method,"validated NAO projection");
            if(value.retained_norm)ImGui::Text("Retained NAO norm: %.6f",*value.retained_norm);
            else ImGui::TextDisabled("Retained NAO norm: unavailable");
            if(value.residual_norm)ImGui::Text("AO-metric reconstruction residual: %.3g",*value.residual_norm);
            else ImGui::TextDisabled("AO-metric reconstruction residual: unavailable");
            ImGui::TextDisabled("%s",value.complete?
                "Complete verified NAO expansion":"Reduced or incompletely verified NAO expansion");
            std::vector<const NboNaoGroupContribution*> shells;
            for(const auto& shell:value.shells)shells.push_back(&shell);
            std::sort(shells.begin(),shells.end(),[](auto a,auto b){return a->weight>b->weight;});
            for(std::size_t k=0;k<std::min<std::size_t>(8,shells.size());++k)
                ImGui::BulletText("%s: %.2f%%",shells[k]->label.c_str(),100*shells[k]->weight);
            return;
        }
    }
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
                     has_multicentre?multicentre_descriptor(chemistry):
                         unresolved_canonical_assignment(language),
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
                     has_delocalised?delocalised_descriptor(chemistry,language):
                         unresolved_canonical_assignment(language),
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
            continue_text(channel_text(interaction.channel,language),ImGui::GetStyle().ItemSpacing.x);
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

    draw_label_value(text.donor,
                     chemistry.donor_acceptor=="UND"?
                         unresolved_canonical_assignment(language):chemistry.donor_acceptor,
                     chemistry.donor_acceptor == "UND" ? kUnavailableColour : kMulticentreColour);
    const std::string method=localised_chemistry_method(chemistry.method,language);
    ImGui::TextDisabled("%s: %s · classification support (legacy AO model)", text.method, method.c_str());
    continue_text(percent_text(chemistry.confidence),ImGui::GetStyle().ItemSpacing.x);
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
    draw_selected_chemistry(wavefunction,selected_index,
        state.nbo_ui?state.nbo_ui->routed:nullptr,language);
}

} // namespace cov::ui
