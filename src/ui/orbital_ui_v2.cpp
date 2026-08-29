#include "cov/orbital_ui.hpp"

#include "cov/mo_diagram.hpp"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace cov::ui {
namespace {

constexpr ImU32 kSigmaColour = IM_COL32(57, 210, 232, 255);
constexpr ImU32 kPiColour = IM_COL32(226, 93, 220, 255);
constexpr ImU32 kDeltaColour = IM_COL32(244, 155, 65, 255);
constexpr ImU32 kPhiColour = IM_COL32(55, 199, 170, 255);
constexpr ImU32 kBondingColour = IM_COL32(77, 218, 145, 255);
constexpr ImU32 kAntibondingColour = IM_COL32(244, 93, 105, 255);
constexpr ImU32 kNonbondingColour = IM_COL32(235, 181, 65, 255);
constexpr ImU32 kSymmetryColour = IM_COL32(177, 125, 245, 255);
constexpr ImU32 kMulticentreColour = IM_COL32(238, 194, 89, 255);
constexpr ImU32 kNumericColour = IM_COL32(93, 174, 255, 255);
constexpr ImU32 kUnavailableColour = IM_COL32(128, 149, 177, 255);

ImVec4 text_colour(const ImU32 colour) {
    return ImGui::ColorConvertU32ToFloat4(colour);
}

void labelled_value(const char* label, const std::string& value, const ImU32 colour) {
    ImGui::Text("%s:", label);
    ImGui::SameLine();
    ImGui::TextColored(text_colour(colour), "%s", value.c_str());
}

void labelled_number(const char* label, const std::string& value) {
    labelled_value(label, value, kNumericColour);
}

std::string fixed_number(const double value, const int precision) {
    std::ostringstream result;
    result.setf(std::ios::fixed, std::ios::floatfield);
    result.precision(precision);
    result << value;
    return result.str();
}

ImU32 family_colour(const std::string_view family) {
    if (family == "sigma") return kSigmaColour;
    if (family == "pi") return kPiColour;
    if (family == "delta") return kDeltaColour;
    if (family == "phi") return kPhiColour;
    return kUnavailableColour;
}

ImU32 bonding_colour(const BondingClass value) {
    switch (value) {
        case BondingClass::Bonding: return kBondingColour;
        case BondingClass::Antibonding: return kAntibondingColour;
        case BondingClass::Nonbonding: return kNonbondingColour;
        default: return kUnavailableColour;
    }
}

std::string superscript_number_ui(const std::size_t value) {
    static constexpr const char* digits[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
    std::string result;
    for (const char digit : std::to_string(value)) result += digits[digit - '0'];
    return result;
}

std::string subscript_number_ui(const int value) {
    static constexpr const char* digits[] = {"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"};
    std::string result;
    for (const char digit : std::to_string(std::max(0, value))) result += digits[digit - '0'];
    return result;
}

std::string pi_descriptor_ui(const DelocalisedPiDescriptor& descriptor) {
    return std::string("Π") + superscript_number_ui(descriptor.participating_atoms) +
           subscript_number_ui(static_cast<int>(std::lround(descriptor.participating_electrons)));
}

struct DelocalisedLabels {
    const char* members;
    const char* atoms;
    const char* electrons;
};

DelocalisedLabels delocalised_labels(const Language language) {
    switch (language) {
        case Language::ChineseSimplified: return {"成员 MO", "参与原子", "参与电子"};
        case Language::Japanese: return {"構成 MO", "参加原子", "参加電子"};
        case Language::French: return {"OM membres", "Atomes participants", "Électrons participants"};
        default: return {"Member MOs", "Participating atoms", "Participating electrons"};
    }
}

const char* intermediate_toggle_label(const Language language) {
    switch (language) {
        case Language::ChineseSimplified: return "隐藏中间框架轨道";
        case Language::Japanese: return "中間骨格軌道を非表示";
        case Language::French: return "Masquer les OM intermédiaires";
        default: return "Hide intermediate framework MOs";
    }
}

const char* intermediate_toggle_tooltip(const Language language) {
    switch (language) {
        case Language::ChineseSimplified:
            return "隐藏次级配体中心、配体内部及极化型轨道；保留主要 d 能级、必要的 σ 成键/反键代表、当前轨道和 donor/acceptor 配对。";
        case Language::Japanese:
            return "副次的な配位子中心・配位子内部・分極軌道を隠し、主要 d 準位、必要な σ 結合/反結合代表、選択軌道と donor/acceptor 対を保持します。";
        case Language::French:
            return "Masque les OM secondaires centrées sur les ligands, internes ou de polarisation; conserve les niveaux d, les représentants sigma nécessaires, la sélection et les paires donor/acceptor.";
        default:
            return "Hide secondary ligand-centred, ligand-internal and polarisation MOs; keep the principal d levels, necessary sigma bonding/antibonding representatives, the selected MO and donor/acceptor pairs.";
    }
}

const char* filter_name(const OrbitalFilterMode mode, const Language language) {
    switch (mode) {
        case OrbitalFilterMode::All: return tr(Text::FilterAll, language);
        case OrbitalFilterMode::Occupied: return tr(Text::FilterOccupied, language);
        case OrbitalFilterMode::Virtual: return tr(Text::FilterVirtual, language);
        case OrbitalFilterMode::Core: return tr(Text::FilterCore, language);
        case OrbitalFilterMode::Valence: return tr(Text::FilterValence, language);
        default: return tr(Text::FilterAuto, language);
    }
}

const char* spin_name_ui(const Spin spin, const Language language) {
    return tr(spin == Spin::Beta ? Text::Beta : Text::Alpha, language);
}

const char* family_symbol_ui(const std::string& family) {
    if (family == "sigma") return "σ";
    if (family == "pi") return "π";
    if (family == "delta") return "δ";
    if (family == "phi") return "φ";
    return "N/A";
}

const char* bonding_ui(const BondingClass value) {
    switch (value) {
        case BondingClass::Bonding: return "bonding";
        case BondingClass::Nonbonding: return "nonbonding";
        case BondingClass::Antibonding: return "antibonding";
        default: return "N/A";
    }
}

std::string lower_ascii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](const unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool search_matches(const OrbitalMetadata& item, const char* query) {
    if (!query || !*query) return true;
    std::ostringstream haystack;
    haystack << item.raw_mo_number << ' ' << item.display_label << ' '
             << item.symmetry << ' '
             << (item.spin == Spin::Beta ? "beta" : "alpha") << ' ';
    switch (item.region) {
        case OrbitalRegion::Core: haystack << "core"; break;
        case OrbitalRegion::Valence: haystack << "valence"; break;
        default: haystack << "virtual"; break;
    }
    return lower_ascii(haystack.str()).find(lower_ascii(query)) != std::string::npos;
}

std::optional<std::size_t> previous_occupied(const std::vector<MolecularOrbital>& orbitals,
                                             std::size_t from,
                                             const Spin spin,
                                             const double threshold) {
    if (orbitals.empty()) return std::nullopt;
    if (from > orbitals.size()) from = orbitals.size();
    while (from > 0) {
        --from;
        if (orbitals[from].spin == spin && orbitals[from].occupation > threshold) return from;
    }
    return std::nullopt;
}

std::optional<std::size_t> next_virtual(const std::vector<MolecularOrbital>& orbitals,
                                        std::size_t from,
                                        const Spin spin,
                                        const double threshold) {
    for (std::size_t i = from + 1; i < orbitals.size(); ++i) {
        if (orbitals[i].spin == spin && orbitals[i].occupation <= threshold) return i;
    }
    return std::nullopt;
}

std::optional<std::size_t> frontier_for_selected_spin(const FrontierOrbitals& frontier,
                                                      const Spin spin,
                                                      const bool lumo) {
    if (!frontier.separate_spin_sets) return lumo ? frontier.lumo : frontier.homo;
    if (spin == Spin::Beta) return lumo ? frontier.beta_lumo : frontier.beta_homo;
    return lumo ? frontier.alpha_lumo : frontier.alpha_homo;
}

void quick_nav_button(const char* label,
                      const std::optional<std::size_t> target,
                      OrbitalUIActions& actions,
                      const float width) {
    if (!target) ImGui::BeginDisabled();
    if (ImGui::Button(label, ImVec2(width, 0.0f)) && target) actions.select_orbital = *target;
    if (!target) ImGui::EndDisabled();
}

void energy_unit_combo(OrbitalUIState& state) {
    constexpr EnergyUnit units[] = {
        EnergyUnit::Hartree,
        EnergyUnit::ElectronVolt,
        EnergyUnit::JoulePerMol,
        EnergyUnit::KilojoulePerMol,
        EnergyUnit::CaloriePerMol,
        EnergyUnit::KilocaloriePerMol,
    };
    if (ImGui::BeginCombo("##energy_unit", energy_unit_symbol(state.energy_unit))) {
        for (const EnergyUnit unit : units) {
            const bool selected = state.energy_unit == unit;
            if (ImGui::Selectable(energy_unit_symbol(unit), selected)) state.energy_unit = unit;
            if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
}

void filter_combo(OrbitalUIState& state, const Language language) {
    if (ImGui::BeginCombo("##orbital_filter", filter_name(state.filter.mode, language))) {
        constexpr OrbitalFilterMode modes[] = {
            OrbitalFilterMode::AutoReasonable, OrbitalFilterMode::All,
            OrbitalFilterMode::Occupied, OrbitalFilterMode::Virtual,
            OrbitalFilterMode::Core, OrbitalFilterMode::Valence,
        };
        for (const auto mode : modes) {
            const bool selected = state.filter.mode == mode;
            if (ImGui::Selectable(filter_name(mode, language), selected)) state.filter.mode = mode;
            if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
}

float rich_symmetry_width(const SymmetryNotation& notation) {
    if (notation.base.empty()) return ImGui::CalcTextSize("N/A").x;
    const float base = ImGui::CalcTextSize(notation.base.c_str()).x;
    const float sub = notation.subscript.empty() ? 0.0f : ImGui::CalcTextSize(notation.subscript.c_str()).x * 0.72f;
    const float sup = notation.superscript.empty() ? 0.0f : ImGui::CalcTextSize(notation.superscript.c_str()).x * 0.72f;
    return base + std::max(sub, sup);
}

void draw_rich_symmetry(const std::string& raw,
                        const ImU32 colour = IM_COL32(220, 228, 240, 255)) {
    const SymmetryNotation notation = parse_symmetry_notation(raw);
    if (notation.base.empty()) {
        ImGui::TextColored(text_colour(kUnavailableColour), "N/A");
        return;
    }
    const ImVec2 pos = ImGui::GetCursorScreenPos();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    ImFont* font = ImGui::GetFont();
    const float size = ImGui::GetFontSize();
    draw->AddText(font, size, pos, colour, notation.base.c_str());
    const float x = pos.x + ImGui::CalcTextSize(notation.base.c_str()).x;
    if (!notation.subscript.empty()) {
        draw->AddText(font, size * 0.72f, ImVec2(x, pos.y + size * 0.36f),
                      colour, notation.subscript.c_str());
    }
    if (!notation.superscript.empty()) {
        draw->AddText(font, size * 0.72f, ImVec2(x, pos.y - size * 0.08f),
                      colour, notation.superscript.c_str());
    }
    ImGui::Dummy(ImVec2(rich_symmetry_width(notation), size * 1.08f));
}

int group_member_index(const std::string& label) {
    const auto dash = label.find('-');
    if (dash == std::string::npos || dash + 1 >= label.size()) return 0;
    const char suffix = label[dash + 1];
    return suffix >= 'a' && suffix <= 'z' ? static_cast<int>(suffix - 'a') : 0;
}

std::size_t group_base_raw(const OrbitalMetadata& item) {
    if (item.degeneracy_size <= 1) return item.raw_mo_number;
    const auto member = static_cast<std::size_t>(std::max(0, group_member_index(item.display_label)));
    return item.raw_mo_number >= member ? item.raw_mo_number - member : item.raw_mo_number;
}

void tooltip_source_confidence(const AnnotationSource source,
                               const double confidence,
                               const bool heuristic,
                               const Language language) {
    ImGui::TextDisabled("%s: %s · %s",
                        tr(Text::ClassificationSource, language),
                        annotation_source_name(source),
                        tr(Text::Confidence, language));
    ImGui::SameLine();
    ImGui::TextColored(text_colour(kNumericColour), "%.0f%%", confidence * 100.0);
    ImGui::SameLine();
    ImGui::TextDisabled("· %s", heuristic ? "heuristic" : "direct/parsed");
}

std::string delocalised_member_list(const MODiagramData& data,
                                    const DelocalisedPiDescriptor& descriptor) {
    std::ostringstream result;
    for (const auto orbital_index : descriptor.orbital_indices) {
        if (result.tellp() > 0) result << ", ";
        if (orbital_index < data.metadata.size()) {
            result << "MO " << data.metadata[orbital_index].raw_mo_number;
        } else {
            result << "MO " << orbital_index + 1u;
        }
    }
    return result.str();
}

std::string delocalised_atom_list(const DelocalisedPiDescriptor& descriptor) {
    std::ostringstream result;
    for (const auto atom_index : descriptor.atom_indices) {
        if (result.tellp() > 0) result << ", ";
        result << atom_index + 1u;
    }
    return result.str();
}

ImU32 pi_interaction_colour(PiInteractionKind kind);

void draw_level_tooltip(const MODiagramData& data,
                        const MODiagramLevel& level,
                        const OrbitalUIState& state,
                        const Language language,
                        const std::size_t orbital_index) {
    const OrbitalMetadata& metadata=orbital_index<data.metadata.size()
        ?data.metadata[orbital_index]:level.metadata;
    std::string displayed_symmetry=metadata.symmetry;
    if (displayed_symmetry.empty() || displayed_symmetry=="?" ||
        displayed_symmetry=="N/A" || displayed_symmetry=="n/a") {
        displayed_symmetry=level.metadata.symmetry;
    }
    ImGui::BeginTooltip();
    ImGui::Text("MO %s", metadata.display_label.c_str());
    ImGui::Separator();
    labelled_number(tr(Text::RawMO, language), std::to_string(metadata.raw_mo_number));
    labelled_number(tr(Text::InternalIndex, language), std::to_string(metadata.orbital_index));
    labelled_number(tr(Text::ExactEnergy, language),
                    format_energy(metadata.energy_hartree, state.energy_unit, 6));
    labelled_number("Ha", fixed_number(metadata.energy_hartree, 10));
    labelled_number("eV", fixed_number(convert_hartree(metadata.energy_hartree, EnergyUnit::ElectronVolt), 8));
    labelled_number("J/mol", fixed_number(convert_hartree(metadata.energy_hartree, EnergyUnit::JoulePerMol), 2));
    labelled_number("kJ/mol", fixed_number(convert_hartree(metadata.energy_hartree, EnergyUnit::KilojoulePerMol), 5));
    labelled_number("cal/mol", fixed_number(convert_hartree(metadata.energy_hartree, EnergyUnit::CaloriePerMol), 3));
    labelled_number("kcal/mol", fixed_number(convert_hartree(metadata.energy_hartree, EnergyUnit::KilocaloriePerMol), 5));
    labelled_number(tr(Text::Occupation, language), fixed_number(metadata.occupation, 3));
    ImGui::Text("%s: %s", tr(Text::Spin, language), spin_name_ui(metadata.spin, language));

    ImGui::TextUnformatted(tr(Text::Symmetry, language));
    ImGui::SameLine();
    draw_rich_symmetry(displayed_symmetry, kSymmetryColour);
    if (!data.ligand_field_point_group.empty()) {
        labelled_value("local ligand field",
                       data.ligand_field_point_group+" first shell",
                       kSymmetryColour);
    }

    labelled_number(tr(Text::DegenerateSet, language), std::to_string(level.metadata.degeneracy_size));
    if (level.metadata.degeneracy_size > 1) {
        std::ostringstream members;
        const std::size_t member_count=level.member_indices.empty()
            ?level.metadata.degeneracy_size:level.member_indices.size();
        const std::size_t fallback_base=group_base_raw(level.metadata);
        for (std::size_t i = 0; i < member_count; ++i) {
            if (i) members << ", ";
            const std::size_t member_index=level.member_indices.empty()
                ?fallback_base+i-1u:level.member_indices[i];
            if (member_index<data.metadata.size()) {
                members << data.metadata[member_index].display_label;
            } else {
                members << member_index+1u;
            }
        }
        labelled_value(tr(Text::DegenerateMembers, language), members.str(), kNumericColour);
    }

    labelled_number("group occupation",fixed_number(level.total_occupation,3));
    labelled_number("metal s / p / d",
                    fixed_number(100.0*level.metal_s_weight,1)+"% / "+
                    fixed_number(100.0*level.metal_p_weight,1)+"% / "+
                    fixed_number(100.0*level.metal_d_weight,1)+"%");
    labelled_number("ligand p",fixed_number(100.0*level.ligand_p_weight,1)+"%");
    labelled_number("sigma / pi channel",
                    fixed_number(100.0*level.sigma_fraction,1)+"% / "+
                    fixed_number(100.0*level.pi_fraction,1)+"%");
    labelled_number("M-L overlap",fixed_number(level.metal_ligand_overlap,6));
    if (level.raw_data_fallback) {
        labelled_value("selection","recovered from complete raw MO block",kNumericColour);
    }
    if (level.approximate_nonbonding) {
        labelled_value("weak-field treatment","approximately nonbonding",kUnavailableColour);
    }

    const auto level_index=static_cast<std::size_t>(&level-data.levels.data());
    for (const auto& interaction:data.pi_interactions) {
        if (interaction.lower_level!=level_index &&
            interaction.upper_level!=level_index &&
            interaction.retained_level!=level_index) continue;
        labelled_value("pi interaction",
                       pi_interaction_kind_name(interaction.kind),
                       pi_interaction_colour(interaction.kind));
        labelled_number("pi splitting",
                        format_energy(interaction.splitting_hartree,
                                      state.energy_unit,6));
        labelled_number("pair confidence",fixed_number(interaction.confidence,3));
    }

    ImGui::Separator();
    labelled_value(tr(Text::OrbitalFamily, language), family_symbol_ui(level.annotation.family),
                   family_colour(level.annotation.family));
    tooltip_source_confidence(level.annotation.family_source, level.annotation.family_confidence,
                              level.annotation.family_source == AnnotationSource::Heuristic, language);
    labelled_value(tr(Text::BondingClassLabel, language), bonding_ui(level.annotation.bonding_class),
                   bonding_colour(level.annotation.bonding_class));
    tooltip_source_confidence(level.annotation.bonding_source, level.annotation.bonding_confidence,
                              level.annotation.bonding_source == AnnotationSource::Heuristic, language);

    std::string multicentre_value = level.annotation.multicentre.available
        ? level.annotation.multicentre.label : "N/A";
    ImU32 multicentre_tone = level.annotation.multicentre.available
        ? kMulticentreColour : kUnavailableColour;
    if (level.annotation.delocalised_pi.available) {
        multicentre_value = pi_descriptor_ui(level.annotation.delocalised_pi) + " · " + multicentre_value;
        multicentre_tone = kPiColour;
    }
    labelled_value(tr(Text::MulticentreBond, language), multicentre_value, multicentre_tone);
    if (level.annotation.multicentre.available) {
        tooltip_source_confidence(level.annotation.multicentre.source,
                                  level.annotation.multicentre.confidence,
                                  level.annotation.multicentre.heuristic, language);
    }
    if (level.annotation.delocalised_pi.available) {
        const auto labels = delocalised_labels(language);
        const auto& descriptor = level.annotation.delocalised_pi;
        labelled_value(tr(Text::DelocalisedPiSystem, language), pi_descriptor_ui(descriptor), kPiColour);
        labelled_value(labels.members, delocalised_member_list(data, descriptor), kPiColour);
        std::string atom_value = std::to_string(descriptor.participating_atoms);
        if (!descriptor.atom_indices.empty()) atom_value += " · " + delocalised_atom_list(descriptor);
        labelled_number(labels.atoms, atom_value);
        labelled_number(labels.electrons,
                        std::to_string(static_cast<int>(std::lround(descriptor.participating_electrons))));
        tooltip_source_confidence(level.annotation.delocalised_pi.source,
                                  level.annotation.delocalised_pi.confidence,
                                  level.annotation.delocalised_pi.heuristic, language);
    } else {
        labelled_value(tr(Text::DelocalisedPiSystem, language), "N/A", kUnavailableColour);
    }
    ImGui::EndTooltip();
}

float map_energy_y(const double energy,
                   const EnergyTransform& transform,
                   const float top,
                   const float bottom) {
    if (transform.knots.empty()) return 0.5f * (top + bottom);
    const double c0 = transform.knots.front().coordinate;
    const double c1 = transform.knots.back().coordinate;
    const double c = energy_display_coordinate(energy, transform);
    const double t = (c - c0) / std::max(1.0e-12, c1 - c0);
    return bottom - static_cast<float>(t) * (bottom - top);
}

double energy_at_fraction(const double fraction, const EnergyTransform& transform) {
    if (transform.knots.empty()) return 0.0;
    const double c0 = transform.knots.front().coordinate;
    const double c1 = transform.knots.back().coordinate;
    return energy_from_display_coordinate(c0 + fraction * (c1 - c0), transform);
}

struct DiagramPoint {
    const MODiagramLevel* level = nullptr;
    ImVec2 left{};
    ImVec2 right{};
    float y = 0.0f;
};

struct DiagramMemberPoint {
    const MODiagramLevel* level = nullptr;
    std::size_t orbital_index = 0;
    ImVec2 left{};
    ImVec2 right{};
    float y = 0.0f;
};

struct PiDiagramGroup {
    const DelocalisedPiDescriptor* descriptor = nullptr;
    std::vector<const DiagramPoint*> points;
};

std::string pi_group_key(const DelocalisedPiDescriptor& descriptor) {
    if (!descriptor.family_id.empty()) return descriptor.family_id;
    std::ostringstream key;
    key << descriptor.label << ':' << descriptor.participating_atoms << ':'
        << std::lround(descriptor.participating_electrons);
    for (const auto atom_index : descriptor.atom_indices) key << ":a" << atom_index;
    for (const auto orbital_index : descriptor.orbital_indices) key << ':' << orbital_index;
    return key.str();
}

void draw_pi_groups(ImDrawList* draw,
                    const std::vector<DiagramPoint>& points,
                    const ImVec2 canvas_max,
                    const float ui_scale) {
    std::map<std::string, PiDiagramGroup> groups;
    for (const auto& point : points) {
        if (!point.level || !point.level->annotation.delocalised_pi.available) continue;
        const auto& descriptor = point.level->annotation.delocalised_pi;
        if (descriptor.participating_atoms < 3u) continue;
        auto& group = groups[pi_group_key(descriptor)];
        if (!group.descriptor) group.descriptor = &descriptor;
        group.points.push_back(&point);
    }

    float bracket_x = canvas_max.x - 13.0f * ui_scale;
    for (const auto& [key, group] : groups) {
        (void)key;
        if (!group.descriptor || group.points.size() < 2u) continue;
        if (!group.descriptor->orbital_indices.empty() &&
            group.points.size() != group.descriptor->orbital_indices.size()) {
            continue;
        }
        float upper = group.points.front()->y;
        float lower = group.points.front()->y;
        for (const auto* point : group.points) {
            upper = std::min(upper, point->y);
            lower = std::max(lower, point->y);
        }
        upper -= 5.0f * ui_scale;
        lower += 5.0f * ui_scale;
        if (lower - upper < 14.0f * ui_scale) {
            const float middle = 0.5f * (upper + lower);
            upper = middle - 7.0f * ui_scale;
            lower = middle + 7.0f * ui_scale;
        }

        const float cap = 7.0f * ui_scale;
        draw->AddLine(ImVec2(bracket_x, upper), ImVec2(bracket_x, lower), kPiColour, 2.0f * ui_scale);
        draw->AddLine(ImVec2(bracket_x - cap, upper), ImVec2(bracket_x, upper), kPiColour, 2.0f * ui_scale);
        draw->AddLine(ImVec2(bracket_x - cap, lower), ImVec2(bracket_x, lower), kPiColour, 2.0f * ui_scale);

        const std::string label = pi_descriptor_ui(*group.descriptor);
        const ImVec2 label_size = ImGui::CalcTextSize(label.c_str());
        draw->AddText(ImVec2(bracket_x - cap - label_size.x - 5.0f * ui_scale,
                             0.5f * (upper + lower) - 0.5f * label_size.y),
                      kPiColour, label.c_str());
        bracket_x -= std::max(20.0f * ui_scale, label_size.x + 12.0f * ui_scale);
    }
}

ImU32 pi_interaction_colour(const PiInteractionKind kind) {
    switch (kind) {
        case PiInteractionKind::Donor: return IM_COL32(242,163,64,255);
        case PiInteractionKind::Acceptor: return IM_COL32(70,206,218,255);
        case PiInteractionKind::WeakNearNonbonding:
            return IM_COL32(174,154,190,255);
        default: return kPiColour;
    }
}

void draw_ligand_field_pi_interactions(
    ImDrawList* draw,
    const MODiagramData& data,
    const std::vector<DiagramPoint>& points,
    const ImVec2 canvas_max,
    const EnergyUnit unit,
    const float ui_scale) {
    float bracket_x=canvas_max.x-13.0f*ui_scale;
    for (const auto& interaction:data.pi_interactions) {
        if (interaction.retained_level>=points.size()) continue;
        const ImU32 colour=pi_interaction_colour(interaction.kind);
        const char* symbol="Δπ";
        const ImVec2 symbol_size=ImGui::CalcTextSize(symbol);
        const float cap=7.0f*ui_scale;
        float hit_top=0.0f;
        float hit_bottom=0.0f;

        if (interaction.kind==PiInteractionKind::WeakNearNonbonding ||
            !interaction.lower_visible || !interaction.upper_visible ||
            interaction.lower_level>=points.size() ||
            interaction.upper_level>=points.size()) {
            const float y=points[interaction.retained_level].y;
            hit_top=y-7.0f*ui_scale;
            hit_bottom=y+7.0f*ui_scale;
            draw->AddLine(ImVec2(bracket_x,hit_top),
                          ImVec2(bracket_x,hit_bottom),colour,2.0f*ui_scale);
            draw->AddLine(ImVec2(bracket_x-cap,hit_top),
                          ImVec2(bracket_x,hit_top),colour,2.0f*ui_scale);
            draw->AddLine(ImVec2(bracket_x-cap,hit_bottom),
                          ImVec2(bracket_x,hit_bottom),
                          colour,2.0f*ui_scale);
        } else {
            float upper=std::min(points[interaction.lower_level].y,
                                 points[interaction.upper_level].y)-5.0f*ui_scale;
            float lower=std::max(points[interaction.lower_level].y,
                                 points[interaction.upper_level].y)+5.0f*ui_scale;
            if (lower-upper<14.0f*ui_scale) {
                const float middle=0.5f*(upper+lower);
                upper=middle-7.0f*ui_scale;
                lower=middle+7.0f*ui_scale;
            }
            hit_top=upper;
            hit_bottom=lower;
            draw->AddLine(ImVec2(bracket_x,upper),ImVec2(bracket_x,lower),
                          colour,2.0f*ui_scale);
            draw->AddLine(ImVec2(bracket_x-cap,upper),ImVec2(bracket_x,upper),
                          colour,2.0f*ui_scale);
            draw->AddLine(ImVec2(bracket_x-cap,lower),ImVec2(bracket_x,lower),
                          colour,2.0f*ui_scale);
        }
        const float symbol_x=bracket_x-cap-symbol_size.x-5.0f*ui_scale;
        const float symbol_y=0.5f*(hit_top+hit_bottom)-0.5f*symbol_size.y;
        draw->AddText(ImVec2(symbol_x,symbol_y),colour,symbol);

        if (ImGui::IsItemHovered()) {
            const ImVec2 mouse=ImGui::GetIO().MousePos;
            const float pad=6.0f*ui_scale;
            if (mouse.x>=symbol_x-pad && mouse.x<=bracket_x+pad &&
                mouse.y>=hit_top-pad && mouse.y<=hit_bottom+pad) {
                ImGui::BeginTooltip();
                if (!interaction.symmetry.empty()) {
                    labelled_value("symmetry",interaction.symmetry,kSymmetryColour);
                }
                labelled_value("interaction",
                               pi_interaction_kind_name(interaction.kind),colour);
                labelled_number("splitting",
                    format_energy(interaction.splitting_hartree,unit,6));
                labelled_number("splitting (Ha)",
                    fixed_number(interaction.splitting_hartree,10));
                labelled_number("confidence",fixed_number(interaction.confidence,3));
                ImGui::EndTooltip();
            }
        }
        bracket_x-=32.0f*ui_scale;
    }
}

void draw_arrow(ImDrawList* draw, const ImVec2 start, const bool up, const ImU32 colour) {
    const float dy = up ? -13.0f : 13.0f;
    const ImVec2 tip(start.x, start.y + dy);
    draw->AddLine(start, tip, colour, 1.6f);
    const float head = up ? 4.0f : -4.0f;
    draw->AddTriangleFilled(tip, ImVec2(tip.x - 3.5f, tip.y + head),
                            ImVec2(tip.x + 3.5f, tip.y + head), colour);
}

std::string compact_metadata(const OrbitalMetadata& item, const EnergyUnit unit) {
    std::ostringstream out;
    out << "label=" << item.display_label << "; raw_mo=" << item.raw_mo_number
        << "; internal_index=" << item.orbital_index
        << "; energy_hartree=" << item.energy_hartree
        << "; energy_display=" << convert_hartree(item.energy_hartree, unit)
        << ' ' << energy_unit_symbol(unit)
        << "; occupation=" << item.occupation << "; symmetry=" << item.symmetry;
    return out.str();
}

bool same_degeneracy_settings(const DegeneracySettings& left,
                              const DegeneracySettings& right) noexcept {
    return left.tolerance_hartree==right.tolerance_hartree &&
        left.require_same_spin==right.require_same_spin &&
        left.require_compatible_symmetry==right.require_compatible_symmetry;
}

bool same_filter_settings(const OrbitalFilterSettings& left,
                          const OrbitalFilterSettings& right) noexcept {
    return left.mode==right.mode &&
        left.occupation_threshold==right.occupation_threshold &&
        left.virtual_window_hartree==right.virtual_window_hartree &&
        left.core_energy_cutoff_hartree==right.core_energy_cutoff_hartree;
}

bool same_diagram_options(const MODiagramOptions& left,
                          const MODiagramOptions& right) noexcept {
    return left.mode==right.mode &&
        left.energy_unit==right.energy_unit &&
        left.energy_axis_mode==right.energy_axis_mode &&
        same_degeneracy_settings(left.degeneracy,right.degeneracy) &&
        same_filter_settings(left.filter,right.filter) &&
        left.selected_index==right.selected_index &&
        left.neighbourhood==right.neighbourhood &&
        left.max_levels==right.max_levels &&
        left.max_virtual_levels==right.max_virtual_levels &&
        left.hide_ligand_centred_intermediates==
            right.hide_ligand_centred_intermediates &&
        left.nonlinear_minimum_gap_weight==
            right.nonlinear_minimum_gap_weight &&
        left.weak_pi_split_hartree==right.weak_pi_split_hartree &&
        left.weak_metal_ligand_overlap==
            right.weak_metal_ligand_overlap &&
        left.width==right.width && left.height==right.height &&
        left.include_hidden_in_metadata==
            right.include_hidden_in_metadata;
}

bool diagram_cache_matches(const OrbitalUIDiagramCache& cache,
                           const Wavefunction& wavefunction,
                           const MODiagramOptions& options) noexcept {
    return cache.data.has_value() && cache.options.has_value() &&
        cache.wavefunction==&wavefunction &&
        cache.orbital_data==wavefunction.orbitals.data() &&
        cache.atom_count==wavefunction.atoms.size() &&
        cache.orbital_count==wavefunction.orbitals.size() &&
        same_diagram_options(*cache.options,options);
}

bool browser_cache_matches(const OrbitalUIBrowserCache& cache,
                           const Wavefunction& wavefunction,
                           const DegeneracySettings& degeneracy,
                           const OrbitalFilterSettings& filter) noexcept {
    return cache.frontier.has_value() &&
        cache.wavefunction==&wavefunction &&
        cache.orbital_data==wavefunction.orbitals.data() &&
        cache.atom_count==wavefunction.atoms.size() &&
        cache.orbital_count==wavefunction.orbitals.size() &&
        cache.metadata.size()==wavefunction.orbitals.size() &&
        cache.degeneracy.has_value() && cache.filter.has_value() &&
        same_degeneracy_settings(*cache.degeneracy,degeneracy) &&
        same_filter_settings(*cache.filter,filter);
}

} // namespace

void draw_orbital_browser(const Wavefunction& wavefunction,
                          const std::size_t selected_index,
                          OrbitalUIState& state,
                          const Language language,
                          const float ui_scale,
                          OrbitalUIActions& actions) {
    if (wavefunction.orbitals.empty()) {
        ImGui::TextDisabled("%s", tr(Text::NoOrbitals, language));
        return;
    }

    state.degeneracy.tolerance_hartree = std::clamp(state.degeneracy.tolerance_hartree, 1.0e-9, 1.0e-2);
    state.filter.virtual_window_hartree = std::clamp(state.filter.virtual_window_hartree, 0.01, 20.0);
    if (!browser_cache_matches(state.browser_cache,wavefunction,
                               state.degeneracy,state.filter)) {
        state.browser_cache.wavefunction=&wavefunction;
        state.browser_cache.orbital_data=wavefunction.orbitals.data();
        state.browser_cache.atom_count=wavefunction.atoms.size();
        state.browser_cache.orbital_count=wavefunction.orbitals.size();
        state.browser_cache.degeneracy=state.degeneracy;
        state.browser_cache.filter=state.filter;
        state.browser_cache.frontier=find_frontier_orbitals(
            wavefunction.orbitals,state.filter.occupation_threshold);
        // Selection is a transient 3-D inspection state.  None of the table
        // metadata used below depends on its selected flag, so a sentinel
        // keeps the expensive ligand-field classification cacheable while
        // row highlighting continues to use selected_index directly.
        state.browser_cache.metadata=build_orbital_metadata(
            wavefunction,wavefunction.orbitals.size(),
            state.degeneracy,state.filter);
    }
    const FrontierOrbitals& frontier=*state.browser_cache.frontier;
    const auto& metadata=state.browser_cache.metadata;
    const Spin selected_spin = selected_index < wavefunction.orbitals.size()
        ? wavefunction.orbitals[selected_index].spin : Spin::Alpha;
    const auto homo = frontier_for_selected_spin(frontier, selected_spin, false);
    const auto lumo = frontier_for_selected_spin(frontier, selected_spin, true);
    const auto hm1 = homo ? previous_occupied(wavefunction.orbitals, *homo, selected_spin,
                                              state.filter.occupation_threshold) : std::nullopt;
    const auto lp1 = lumo ? next_virtual(wavefunction.orbitals, *lumo, selected_spin,
                                         state.filter.occupation_threshold) : std::nullopt;

    const float gap = ImGui::GetStyle().ItemSpacing.x;
    const float bw = std::max(52.0f * ui_scale, (ImGui::GetContentRegionAvail().x - 3.0f * gap) * 0.25f);
    quick_nav_button(tr(Text::HOMOMinus1, language), hm1, actions, bw); ImGui::SameLine();
    quick_nav_button(tr(Text::HOMO, language), homo, actions, bw); ImGui::SameLine();
    quick_nav_button(tr(Text::LUMO, language), lumo, actions, bw); ImGui::SameLine();
    quick_nav_button(tr(Text::LUMOPlus1, language), lp1, actions, bw);

    ImGui::TextDisabled("%s", tr(Text::Search, language));
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputTextWithHint("##orbital_search", tr(Text::Search, language),
                             state.search.data(), state.search.size());

    if (ImGui::BeginTable("##filter_controls", 2, ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_NoSavedSettings)) {
        ImGui::TableNextColumn(); ImGui::TextDisabled("%s", tr(Text::Filter, language));
        ImGui::SetNextItemWidth(-1.0f); filter_combo(state, language);
        ImGui::TableNextColumn(); ImGui::TextDisabled("%s", tr(Text::EnergyUnit, language));
        ImGui::SetNextItemWidth(-1.0f); energy_unit_combo(state);
        ImGui::EndTable();
    }

    if (state.filter.mode == OrbitalFilterMode::AutoReasonable) {
        float window = static_cast<float>(state.filter.virtual_window_hartree);
        ImGui::TextDisabled("%s (Ha)", tr(Text::HighVirtualWindow, language));
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::SliderFloat("##virtual_window", &window, 0.05f, 5.0f, "%.2f")) {
            state.filter.virtual_window_hartree = window;
        }
    }

    double tolerance = state.degeneracy.tolerance_hartree;
    ImGui::TextDisabled("%s (Ha)", tr(Text::DegeneracyTolerance, language));
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::InputDouble("##degeneracy_tol", &tolerance, 1.0e-6, 1.0e-5, "%.2e")) {
        state.degeneracy.tolerance_hartree = std::clamp(tolerance, 1.0e-9, 1.0e-2);
    }
    ImGui::Checkbox(tr(Text::GroupedLabels, language), &state.grouped_labels);

    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < metadata.size(); ++i) {
        if (metadata[i].visible && search_matches(metadata[i], state.search.data())) candidates.push_back(i);
    }
    ImGui::SameLine();
    ImGui::TextDisabled("%s: %zu / %zu", tr(Text::VisibleOrbitals, language), candidates.size(), metadata.size());

    if (ImGui::BeginTable("##orbital_table", 4,
                          ImGuiTableFlags_RowBg | ImGuiTableFlags_BordersInnerH |
                          ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchProp |
                          ImGuiTableFlags_NoSavedSettings, ImVec2(0.0f, 270.0f * ui_scale))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("MO", ImGuiTableColumnFlags_WidthFixed, 76.0f * ui_scale);
        ImGui::TableSetupColumn(tr(Text::Energy, language), ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn(tr(Text::Occupation, language), ImGuiTableColumnFlags_WidthFixed, 64.0f * ui_scale);
        ImGui::TableSetupColumn(tr(Text::Symmetry, language), ImGuiTableColumnFlags_WidthFixed, 80.0f * ui_scale);
        ImGui::TableHeadersRow();
        ImGuiListClipper clipper;
        clipper.Begin(static_cast<int>(candidates.size()));
        while (clipper.Step()) {
            for (int row = clipper.DisplayStart; row < clipper.DisplayEnd; ++row) {
                const std::size_t index = candidates[static_cast<std::size_t>(row)];
                const auto& item = metadata[index];
                const std::string label = state.grouped_labels ? item.display_label : std::to_string(item.raw_mo_number);
                ImGui::PushID(static_cast<int>(index));
                ImGui::TableNextRow(); ImGui::TableNextColumn();
                if (ImGui::Selectable(label.c_str(), index == selected_index,
                                      ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap)) {
                    actions.select_orbital = index;
                }
                ImGui::TableNextColumn();
                ImGui::TextColored(text_colour(kNumericColour), "%s",
                                   format_energy(item.energy_hartree, state.energy_unit, 5).c_str());
                ImGui::TableNextColumn();
                ImGui::TextColored(text_colour(kNumericColour), "%.2f", item.occupation);
                ImGui::TableNextColumn(); draw_rich_symmetry(item.symmetry, kSymmetryColour);
                ImGui::PopID();
            }
        }
        ImGui::EndTable();
    }

    if (selected_index < metadata.size() && ImGui::Button(tr(Text::CopyMetadata, language))) {
        const std::string text = compact_metadata(metadata[selected_index], state.energy_unit);
        ImGui::SetClipboardText(text.c_str());
    }
}

void draw_energy_diagram(const Wavefunction& wavefunction,
                         const std::size_t selected_index,
                         OrbitalUIState& state,
                         const Language language,
                         const float ui_scale,
                         OrbitalUIActions& actions) {
    if (wavefunction.orbitals.empty()) {
        ImGui::TextDisabled("%s", tr(Text::NoOrbitals, language));
        return;
    }

    MODiagramOptions options;
    options.energy_unit = state.energy_unit;
    options.energy_axis_mode = state.energy_axis_mode;
    options.degeneracy = state.degeneracy;
    options.filter = state.filter;
    // A compact ligand-field diagram is a canonical chemical summary.  The
    // MO selected for 3-D inspection must not add, remove or reprioritise its
    // rows.  The live selected_index is still used below for member
    // highlighting and CUDA dispatch; an out-of-range sentinel keeps the
    // expensive structural build independent of inspection state.
    options.selected_index = state.hide_ligand_centred_intermediates
        ?wavefunction.orbitals.size():selected_index;
    options.neighbourhood = static_cast<std::size_t>(std::max(3, state.diagram_neighbourhood));
    options.hide_ligand_centred_intermediates=
        state.hide_ligand_centred_intermediates;
    options.nonlinear_minimum_gap_weight = 0.055;
    if (!diagram_cache_matches(state.diagram_cache,wavefunction,options)) {
        state.diagram_cache.wavefunction=&wavefunction;
        state.diagram_cache.orbital_data=wavefunction.orbitals.data();
        state.diagram_cache.atom_count=wavefunction.atoms.size();
        state.diagram_cache.orbital_count=wavefunction.orbitals.size();
        state.diagram_cache.options=options;
        state.diagram_cache.data=build_mo_diagram_data(wavefunction,options);
    }
    const MODiagramData& data=*state.diagram_cache.data;

    ImGui::TextDisabled("%s", tr(Text::EnergyDiagram, language));
    ImGui::TextDisabled("%s", data.selection.summary.c_str());
    ImGui::TextDisabled("%s", tr(Text::EnergyScale, language));
    ImGui::SameLine();
    if (ImGui::RadioButton(tr(Text::LinearEnergyScale, language), state.energy_axis_mode == EnergyAxisMode::Linear)) {
        state.energy_axis_mode = EnergyAxisMode::Linear;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton(tr(Text::NonlinearFocus, language), state.energy_axis_mode == EnergyAxisMode::NonlinearFocus)) {
        state.energy_axis_mode = EnergyAxisMode::NonlinearFocus;
    }
    ImGui::SetNextItemWidth(155.0f * ui_scale);
    ImGui::SliderInt("##diagram_neighbourhood", &state.diagram_neighbourhood, 3, 32, "%d", ImGuiSliderFlags_AlwaysClamp);
    ImGui::SameLine(); ImGui::TextDisabled("%s", tr(Text::AroundSelected, language));
    if (ImGui::Checkbox(intermediate_toggle_label(language),
                        &state.hide_ligand_centred_intermediates)) {
        // The diagram is rebuilt on the next immediate-mode frame.
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("%s",intermediate_toggle_tooltip(language));
    }

    const float height = 430.0f * ui_scale;
    ImGui::InvisibleButton("##energy_diagram_canvas", ImVec2(-1.0f, height), ImGuiButtonFlags_MouseButtonLeft);
    const ImVec2 p0 = ImGui::GetItemRectMin();
    const ImVec2 p1 = ImGui::GetItemRectMax();
    ImDrawList* draw = ImGui::GetWindowDrawList();
    draw->AddRectFilled(p0, p1, IM_COL32(12, 18, 27, 235), 7.0f * ui_scale);
    draw->AddRect(p0, p1, IM_COL32(43, 58, 77, 220), 7.0f * ui_scale);
    if (data.levels.empty()) return;

    const float top = p0.y + 34.0f * ui_scale;
    const float bottom = p1.y - 24.0f * ui_scale;
    const float left = p0.x + 62.0f * ui_scale;
    const float right = p1.x - 30.0f * ui_scale;
    const float lane_span = right - left;
    const float axis_x = left - 28.0f * ui_scale;
    draw->AddLine(ImVec2(axis_x, bottom), ImVec2(axis_x, top), IM_COL32(129,148,171,255), 1.4f * ui_scale);
    draw->AddTriangleFilled(ImVec2(axis_x, top - 5.0f * ui_scale),
                            ImVec2(axis_x - 4.0f * ui_scale, top + 3.0f * ui_scale),
                            ImVec2(axis_x + 4.0f * ui_scale, top + 3.0f * ui_scale), IM_COL32(129,148,171,255));
    draw->AddText(ImVec2(p0.x + 7.0f * ui_scale, p0.y + 7.0f * ui_scale),
                  IM_COL32(129,148,171,255),
                  state.energy_axis_mode == EnergyAxisMode::Linear
                    ? tr(Text::LinearEnergyScale, language)
                    : tr(Text::NonlinearEnergyScale, language));

    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const double energy = energy_at_fraction(fraction, data.energy_transform);
        const float y = map_energy_y(energy, data.energy_transform, top, bottom);
        draw->AddLine(ImVec2(axis_x - 4.0f * ui_scale, y), ImVec2(axis_x + 4.0f * ui_scale, y),
                      IM_COL32(100,118,140,255), 1.0f);
        const std::string tick = format_energy(energy, state.energy_unit, 3);
        draw->AddText(ImVec2(axis_x + 6.0f * ui_scale, y - ImGui::GetTextLineHeight() * 0.5f),
                      kNumericColour, tick.c_str());
    }

    std::map<long long, std::vector<std::size_t>> coincident;
    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        coincident[static_cast<long long>(std::llround(data.levels[i].layout_energy_hartree * 1.0e12))].push_back(i);
    }
    std::vector<float> x(data.levels.size(), left + 0.5f * lane_span);
    for (const auto& [_, ids] : coincident) {
        const float spacing = std::min(105.0f * ui_scale, lane_span / std::max(2.0f, static_cast<float>(ids.size())));
        for (std::size_t j = 0; j < ids.size(); ++j) {
            x[ids[j]] += (static_cast<float>(j) - (static_cast<float>(ids.size()) - 1.0f) * 0.5f) * spacing;
        }
    }

    std::vector<DiagramPoint> points;
    points.reserve(data.levels.size());
    std::vector<DiagramMemberPoint> member_points;
    member_points.reserve(data.metadata.size());
    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        const auto& level = data.levels[i];
        const float y = map_energy_y(level.layout_energy_hartree, data.energy_transform, top, bottom);
        const std::size_t members=level.member_indices.empty()
            ?std::max<std::size_t>(1u,level.metadata.degeneracy_size)
            :level.member_indices.size();
        const float member_half=12.0f*ui_scale;
        const float member_spacing=31.0f*ui_scale;
        const float group_half=member_half+
            0.5f*member_spacing*static_cast<float>(members-1u);
        ImU32 colour = IM_COL32(196,210,227,255);
        if (level.metadata.region == OrbitalRegion::Virtual) colour = IM_COL32(126,143,165,255);
        if (level.homo) colour = IM_COL32(79,210,157,255);
        if (level.lumo) colour = IM_COL32(228,176,82,255);
        if (level.annotation.family != "unavailable") colour = family_colour(level.annotation.family);
        for (std::size_t member=0;member<members;++member) {
            const float member_x=x[i]+(static_cast<float>(member)-
                0.5f*static_cast<float>(members-1u))*member_spacing;
            const std::size_t member_orbital=member<level.member_indices.size()
                ?level.member_indices[member]
                :level.metadata.orbital_index+member;
            const std::size_t spin_counterpart=
                member<level.member_spin_counterparts.size()
                    ?level.member_spin_counterparts[member]
                    :wavefunction.orbitals.size();
            const bool counterpart_selected=spin_counterpart==selected_index;
            const bool member_selected=member_orbital==selected_index ||
                                       counterpart_selected;
            if (member_selected) {
                draw->AddLine(ImVec2(member_x-member_half-1.0f*ui_scale,y),
                              ImVec2(member_x+member_half+1.0f*ui_scale,y),
                              kNumericColour,5.0f*ui_scale);
            }
            draw->AddLine(ImVec2(member_x-member_half,y),
                          ImVec2(member_x+member_half,y),colour,
                          member_selected?2.8f:1.8f);
            const ElectronGlyphs electrons=member<level.member_electrons.size()
                ?level.member_electrons[member]:level.electrons;
            if (electrons.alpha>0) {
                draw_arrow(draw,ImVec2(member_x-5.0f*ui_scale,y-2.0f),
                           true,IM_COL32(234,242,252,255));
            }
            if (electrons.beta>0) {
                draw_arrow(draw,ImVec2(member_x+5.0f*ui_scale,y-2.0f),
                           false,IM_COL32(234,242,252,255));
            }
            member_points.push_back({&level,
                counterpart_selected?spin_counterpart:member_orbital,
                ImVec2(member_x-member_half,y),
                ImVec2(member_x+member_half,y),y});
        }
        points.push_back({&level,ImVec2(x[i]-group_half,y),
                          ImVec2(x[i]+group_half,y),y});
    }

    draw_pi_groups(draw, points, p1, ui_scale);
    draw_ligand_field_pi_interactions(
        draw,data,points,p1,state.energy_unit,ui_scale);

    if (ImGui::IsItemHovered()) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        const DiagramMemberPoint* nearest = nullptr;
        float best = 1.0e9f;
        for (const auto& point : member_points) {
            if (mouse.x < point.left.x - 10.0f * ui_scale || mouse.x > point.right.x + 10.0f * ui_scale) continue;
            const float distance = std::abs(mouse.y - point.y);
            if (distance < best && distance <= 8.0f * ui_scale) { best = distance; nearest = &point; }
        }
        if (nearest && nearest->level) {
            draw_level_tooltip(data,*nearest->level,state,language,
                               nearest->orbital_index);
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                actions.select_orbital=nearest->orbital_index;
            }
        }
    }

    if (ImGui::Button(tr(Text::ExportBundle, language), ImVec2(-1.0f, 0.0f))) actions.export_diagram = true;
}

} // namespace cov::ui
