#include "cov/ui.hpp"

#include <imgui.h>

#include <cstring>
#include <iostream>

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

    ImGui::DestroyContext();
    std::cout << "ui_smoke ok\n";
    return 0;
}
