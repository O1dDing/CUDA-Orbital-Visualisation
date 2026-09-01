#include "cov/orbital_ui_text.hpp"

#include <array>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>

namespace {

using cov::AnnotationSource;
using cov::BondingClass;
using cov::DataProvenance;
using cov::MODiagramData;
using cov::MODiagramLevel;
using cov::MODiagramMode;
using cov::OrbitalBondingRole;
using cov::PiInteractionDescriptor;
using cov::PiInteractionKind;
using cov::WavefunctionSource;
using cov::ui::Language;
using cov::ui::OrbitalText;

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "ui_orbital_localisation_smoke: " << message << '\n';
    std::exit(1);
}

void require(const bool condition, const std::string& message) {
    if (!condition) fail(message);
}

void require_contains(const std::string_view value,
                      const std::string_view needle,
                      const std::string& context) {
    require(value.find(needle) != std::string_view::npos,
            context + " missing '" + std::string(needle) + "'");
}

} // namespace

int main() {
    constexpr std::array languages{
        Language::English,
        Language::ChineseSimplified,
        Language::Japanese,
        Language::French,
    };

    for (const auto language : languages) {
        for (std::size_t key = 0;
             key < static_cast<std::size_t>(OrbitalText::Count); ++key) {
            const char* text = cov::ui::orbital_tr(
                static_cast<OrbitalText>(key), language);
            require(text != nullptr && *text != '\0',
                    "empty orbital localisation entry");
        }
        require(std::string_view(cov::ui::orbital_ui_glyph_seed(language)).size() > 100u,
                "glyph seed is unexpectedly short");
    }

    require(std::string_view(cov::ui::orbital_tr(
                OrbitalText::CoordinationGeometry,
                Language::ChineseSimplified)) == "配位几何",
            "Simplified Chinese coordination label mismatch");
    require(std::string_view(cov::ui::orbital_tr(
                OrbitalText::PiWeakNearNonbonding,
                Language::Japanese)).find("弱場") != std::string_view::npos,
            "Japanese weak-field label is not localised");
    require(std::string_view(cov::ui::orbital_tr(
                OrbitalText::GaussianEnrichmentAttached,
                Language::French)).find("Enrichissement") != std::string_view::npos,
            "French enrichment label is not localised");

    require(std::string_view(cov::ui::localised_wavefunction_source(
                WavefunctionSource::Fchk, Language::ChineseSimplified)) == "FCHK",
            "FCHK scientific identifier changed");
    require(std::string_view(cov::ui::localised_data_provenance(
                DataProvenance::Derived, Language::ChineseSimplified)) == "推导数据",
            "provenance value is not localised");
    require(std::string_view(cov::ui::localised_annotation_source(
                AnnotationSource::ParsedLabel, Language::Japanese)) == "ラベル解析",
            "annotation source is not localised");
    require(std::string_view(cov::ui::localised_bonding_class(
                BondingClass::Antibonding, Language::French)) == "antiliante",
            "bonding class is not localised");
    require(std::string_view(cov::ui::localised_orbital_bonding_role(
                OrbitalBondingRole::Nonbonding,
                Language::ChineseSimplified)) == "非键",
            "orbital bonding role is not localised");
    require(std::string_view(cov::ui::localised_pi_interaction_kind(
                PiInteractionKind::Acceptor,
                Language::Japanese)).find("アクセプター") != std::string_view::npos,
            "pi interaction is not localised");

    require(cov::ui::localised_geometry_name(
                "OC-6", "Octahedral", Language::ChineseSimplified) == "八面体形",
            "OC-6 geometry is not localised");
    require(cov::ui::localised_geometry_name(
                "SPC-10", "Sphenocorona", Language::French) == "sphénocouronne",
            "SPC-10 geometry is not localised");
    require(cov::ui::localised_geometry_name(
                "future-id", "Future geometry", Language::Japanese) == "Future geometry",
            "unknown geometry fallback changed");

    const std::string method =
        "COV FCHK S-metric minimal atomic-reference projection";
    require_contains(cov::ui::localised_chemistry_method(
                         method, Language::ChineseSimplified),
                     "S 度量", "Chinese chemistry method");
    require_contains(cov::ui::localised_chemistry_note(
                         "No stable atom-pair interaction frame; chemistry remains UND",
                         Language::French),
                     "UND", "French chemistry note");
    require(cov::ui::localised_chemistry_note(
                "future backend note", Language::ChineseSimplified) ==
                "future backend note",
            "unknown backend note fallback changed");

    MODiagramData data;
    data.mode = MODiagramMode::ValenceCentral;
    data.metadata.resize(12u);
    data.levels.resize(5u);
    data.levels[2].raw_data_fallback = true;
    data.selection.valence_occupied_count = 3u;
    data.selection.frontier_virtual_count = 2u;
    data.selection.hidden_count = 7u;
    data.selection.protected_overflow_count = 1u;
    data.ligand_field_point_group = "Oh";
    data.ligand_field_geometry_id = "OC-6";
    data.ligand_field_coordination_number = 6u;
    data.pi_interactions.resize(2u);
    data.spin_counterpart_pair_count = 4u;
    data.spin_counterpart_unmatched_visible = 1u;

    const std::string zh_summary =
        cov::ui::localised_diagram_selection_summary(
            data, Language::ChineseSimplified);
    const std::string ja_summary =
        cov::ui::localised_diagram_selection_summary(data, Language::Japanese);
    const std::string fr_summary =
        cov::ui::localised_diagram_selection_summary(data, Language::French);
    require_contains(zh_summary, "MO 图摘要", "Chinese diagram summary");
    require_contains(zh_summary, "Oh", "Chinese diagram summary");
    require_contains(ja_summary, "表示準位", "Japanese diagram summary");
    require_contains(fr_summary, "niveaux visibles", "French diagram summary");

    require_contains(cov::ui::orbital_ui_glyph_seed(
                         Language::ChineseSimplified),
                     "斯芬诺冠形", "Chinese glyph seed");
    require_contains(cov::ui::orbital_ui_glyph_seed(Language::Japanese),
                     "アクセプター", "Japanese glyph seed");
    require_contains(cov::ui::orbital_ui_glyph_seed(Language::French),
                     "tétradécaédrique", "French glyph seed");

    std::cout << "ui_orbital_localisation_smoke: four-language orbital UI catalogue verified\n";
    return 0;
}
