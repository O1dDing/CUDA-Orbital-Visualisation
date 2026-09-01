#include "cov/ui.hpp"

#include <imgui.h>

#include <cstring>
#include <iostream>

namespace {

bool font_contains(const ImFont* font, const ImWchar codepoint) {
    return font != nullptr && font->FindGlyphNoFallback(codepoint) != nullptr;
}

bool expect_glyphs(const ImFont* font,
                   const ImWchar* codepoints,
                   const std::size_t count,
                   const char* description) {
    for (std::size_t i = 0; i < count; ++i) {
        if (!font_contains(font, codepoints[i])) {
            std::cerr << "font atlas missing " << description
                      << " codepoint U+" << std::hex
                      << static_cast<unsigned int>(codepoints[i]) << std::dec << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    cov::ui::apply_theme();
    if (!cov::ui::configure_fonts(15.0f)) {
        std::cerr << "font atlas build failed\n";
        ImGui::DestroyContext();
        return 1;
    }

    constexpr cov::ui::Language languages[] = {
        cov::ui::Language::English,
        cov::ui::Language::ChineseSimplified,
        cov::ui::Language::Japanese,
        cov::ui::Language::French,
    };

    for (const auto language : languages) {
        if (!cov::ui::language_name(language) || !*cov::ui::language_name(language)) {
            std::cerr << "missing language name\n";
            ImGui::DestroyContext();
            return 2;
        }
        for (int i = 0; i < static_cast<int>(cov::ui::Text::Count); ++i) {
            const char* text = cov::ui::tr(static_cast<cov::ui::Text>(i), language);
            if (!text || !*text) {
                std::cerr << "missing localisation string at key=" << i
                          << " language=" << static_cast<int>(language) << '\n';
                ImGui::DestroyContext();
                return 3;
            }
        }
    }

    if (std::strcmp(cov::ui::tr(cov::ui::Text::Load,
                                cov::ui::Language::ChineseSimplified),
                    "加载") != 0) {
        std::cerr << "UTF-8 Chinese localisation mismatch\n";
        ImGui::DestroyContext();
        return 4;
    }
    if (std::strcmp(cov::ui::tr(cov::ui::Text::OpenFile,
                                cov::ui::Language::Japanese),
                    "ファイルを開く…") != 0) {
        std::cerr << "UTF-8 Japanese localisation mismatch\n";
        ImGui::DestroyContext();
        return 5;
    }
    if (std::strstr(cov::ui::tr(cov::ui::Text::DegeneracyTolerance,
                                cov::ui::Language::French),
                    "dégénérescence") == nullptr) {
        std::cerr << "UTF-8 French localisation mismatch\n";
        ImGui::DestroyContext();
        return 6;
    }
    if (std::strcmp(cov::ui::tr(cov::ui::Text::MediumBallStick,
                                cov::ui::Language::English),
                    "Enhanced ball-and-stick") != 0 ||
        std::strcmp(cov::ui::tr(cov::ui::Text::AroundSelected,
                                cov::ui::Language::ChineseSimplified),
                    "价电子 MO 图范围") != 0 ||
        std::strcmp(cov::ui::tr(cov::ui::Text::SimpleDiagram,
                                cov::ui::Language::French),
                    "Diagramme MO de valence") != 0) {
        std::cerr << "valence diagram localisation mismatch\n";
        ImGui::DestroyContext();
        return 7;
    }

    const char* zh_seed = cov::ui::supplemental_glyph_seed(
        cov::ui::Language::ChineseSimplified);
    const char* ja_seed = cov::ui::supplemental_glyph_seed(
        cov::ui::Language::Japanese);
    const char* scientific_seed = cov::ui::scientific_glyph_seed();
    if (std::strstr(zh_seed, "轨道材质") == nullptr ||
        std::strstr(zh_seed, "标准") == nullptr ||
        std::strstr(zh_seed, "柔和自动打光") == nullptr ||
        std::strstr(zh_seed, "CUDA 设备") == nullptr ||
        std::strstr(ja_seed, "軌道マテリアル") == nullptr ||
        std::strstr(ja_seed, "標準") == nullptr ||
        std::strstr(ja_seed, "ソフト自動照明") == nullptr ||
        std::strstr(scientific_seed, "Π⁵₆") == nullptr) {
        std::cerr << "supplemental font seed is incomplete\n";
        ImGui::DestroyContext();
        return 8;
    }

    ImFont* primary = ImGui::GetIO().Fonts->Fonts.empty()
                          ? nullptr
                          : ImGui::GetIO().Fonts->Fonts.front();
    constexpr ImWchar scientific_glyphs[] = {
        0x03A0, // Π
        0x2075, // ⁵
        0x2086, // ₆
        0x03C3, // σ
        0x03C0, // π
        0x03B4, // δ
        0x03C6, // φ
    };
    if (!expect_glyphs(primary, scientific_glyphs,
                       sizeof(scientific_glyphs) / sizeof(scientific_glyphs[0]),
                       "scientific")) {
        ImGui::DestroyContext();
        return 9;
    }

    // A runner may legitimately have no CJK font installed.  When a matching
    // OS font was found and merged, however, verify the exact glyphs that used
    // to render as question marks rather than merely checking the source seed.
    constexpr ImWchar chinese_glyphs[] = {
        0x8F68, 0x9053, 0x6750, 0x8D28, // 轨道材质
        0x6807, 0x51C6,                 // 标准
        0x67D4, 0x548C, 0x81EA, 0x52A8, 0x6253, 0x5149, // 柔和自动打光
        0x8BBE, 0x5907,                 // 设备
    };
    if (std::strstr(cov::ui::font_status(), "ZH fallback missing") == nullptr &&
        !expect_glyphs(primary, chinese_glyphs,
                       sizeof(chinese_glyphs) / sizeof(chinese_glyphs[0]),
                       "Chinese UI")) {
        ImGui::DestroyContext();
        return 10;
    }

    constexpr ImWchar japanese_glyphs[] = {
        0x8ECC, 0x9053, 0x30DE, 0x30C6, 0x30EA, 0x30A2, 0x30EB, // 軌道マテリアル
        0x6A19, 0x6E96,                                             // 標準
        0x30BD, 0x30D5, 0x30C8, 0x81EA, 0x52D5, 0x7167, 0x660E,   // ソフト自動照明
    };
    if (std::strstr(cov::ui::font_status(), "JA fallback missing") == nullptr &&
        !expect_glyphs(primary, japanese_glyphs,
                       sizeof(japanese_glyphs) / sizeof(japanese_glyphs[0]),
                       "Japanese UI")) {
        ImGui::DestroyContext();
        return 11;
    }

    ImGui::DestroyContext();
    std::cout << "ui_smoke ok\n";
    return 0;
}
