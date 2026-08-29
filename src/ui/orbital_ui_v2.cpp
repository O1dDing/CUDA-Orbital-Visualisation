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

void draw_level_tooltip(const MODiagramData& data,
                        const MODiagramLevel& level,
                        const OrbitalUIState& state,
                        const Language language) {
    ImGui::BeginTooltip();
    ImGui::Text("MO %s", level.metadata.display_label.c_str());
    ImGui::Separator();
    labelled_number(tr(Text::RawMO, language), std::to_string(level.metadata.raw_mo_number));
    labelled_number(tr(Text::InternalIndex, language), std::to_string(level.metadata.orbital_index));
    labelled_number(tr(Text::ExactEnergy, language),
                    format_energy(level.metadata.energy_hartree, state.energy_unit, 6));
    labelled_number("Ha", fixed_number(level.metadata.energy_hartree, 10));
    labelled_number("eV", fixed_number(convert_hartree(level.metadata.energy_hartree, EnergyUnit::ElectronVolt), 8));
    labelled_number("J/mol", fixed_number(convert_hartree(level.metadata.energy_hartree, EnergyUnit::JoulePerMol), 2));
    labelled_number("kJ/mol", fixed_number(convert_hartree(level.metadata.energy_hartree, EnergyUnit::KilojoulePerMol), 5));
    labelled_number("cal/mol", fixed_number(convert_hartree(level.metadata.energy_hartree, EnergyUnit::CaloriePerMol), 3));
    labelled_number("kcal/mol", fixed_number(convert_hartree(level.metadata.energy_hartree, EnergyUnit::KilocaloriePerMol), 5));
    labelled_number(tr(Text::Occupation, language), fixed_number(level.metadata.occupation, 3));
    ImGui::Text("%s: %s", tr(Text::Spin, language), spin_name_ui(level.metadata.spin, language));

    ImGui::TextUnformatted(tr(Text::Symmetry, language));
    ImGui::SameLine();
    draw_rich_symmetry(level.metadata.symmetry, kSymmetryColour);

    labelled_number(tr(Text::DegenerateSet, language), std::to_string(level.metadata.degeneracy_size));
    if (level.metadata.degeneracy_size > 1) {
        std::ostringstream members;
        const std::size_t base = group_base_raw(level.metadata);
        for (std::size_t i = 0; i < level.metadata.degeneracy_size; ++i) {
            if (i) members << ", ";
            const std::size_t raw = base + i;
            if (raw > 0 && raw <= data.metadata.size()) members << data.metadata[raw - 1].display_label;
            else members << raw;
        }
        labelled_value(tr(Text::DegenerateMembers, language), members.str(), kNumericColour);
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
    const FrontierOrbitals frontier = find_frontier_orbitals(wavefunction.orbitals, state.filter.occupation_threshold);
    const auto metadata = build_orbital_metadata(wavefunction, selected_index, state.degeneracy, state.filter);
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
    options.selected_index = selected_index;
    options.neighbourhood = static_cast<std::size_t>(std::max(3, state.diagram_neighbourhood));
    options.nonlinear_minimum_gap_weight = 0.055;
    const MODiagramData data = build_mo_diagram_data(wavefunction, options);

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
    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        const auto& level = data.levels[i];
        const float y = map_energy_y(level.layout_energy_hartree, data.energy_transform, top, bottom);
        const float half = std::min(35.0f * ui_scale, lane_span * 0.18f);
        ImU32 colour = IM_COL32(196,210,227,255);
        if (level.metadata.region == OrbitalRegion::Virtual) colour = IM_COL32(126,143,165,255);
        if (level.homo) colour = IM_COL32(79,210,157,255);
        if (level.lumo) colour = IM_COL32(228,176,82,255);
        if (level.annotation.family != "unavailable") colour = family_colour(level.annotation.family);
        if (level.metadata.selected) {
            draw->AddLine(ImVec2(x[i] - half - 1.0f * ui_scale, y),
                          ImVec2(x[i] + half + 1.0f * ui_scale, y),
                          kNumericColour, 5.0f * ui_scale);
        }
        draw->AddLine(ImVec2(x[i] - half, y), ImVec2(x[i] + half, y), colour,
                      level.metadata.selected ? 2.8f : 1.8f);
        if (level.electrons.alpha > 0) draw_arrow(draw, ImVec2(x[i] - 7.0f * ui_scale, y - 2.0f), true, IM_COL32(234,242,252,255));
        if (level.electrons.beta > 0) draw_arrow(draw, ImVec2(x[i] + 7.0f * ui_scale, y - 2.0f), false, IM_COL32(234,242,252,255));
        if (level.metadata.selected) {
            draw->AddText(ImVec2(x[i] + half + 4.0f * ui_scale, y - ImGui::GetTextLineHeight() * 0.5f),
                          colour, level.metadata.display_label.c_str());
        }
        points.push_back({&level, ImVec2(x[i] - half, y), ImVec2(x[i] + half, y), y});
    }

    draw_pi_groups(draw, points, p1, ui_scale);

    if (ImGui::IsItemHovered()) {
        const ImVec2 mouse = ImGui::GetIO().MousePos;
        const DiagramPoint* nearest = nullptr;
        float best = 1.0e9f;
        for (const auto& point : points) {
            if (mouse.x < point.left.x - 10.0f * ui_scale || mouse.x > point.right.x + 10.0f * ui_scale) continue;
            const float distance = std::abs(mouse.y - point.y);
            if (distance < best && distance <= 8.0f * ui_scale) { best = distance; nearest = &point; }
        }
        if (nearest && nearest->level) {
            draw_level_tooltip(data, *nearest->level, state, language);
            if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) actions.select_orbital = nearest->level->metadata.orbital_index;
        }
    }

    if (ImGui::Button(tr(Text::ExportBundle, language), ImVec2(-1.0f, 0.0f))) actions.export_diagram = true;
}

} // namespace cov::ui
