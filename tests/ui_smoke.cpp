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
        const char* title = cov::ui::tr(cov::ui::Text::AppTitle, language);
        const char* rendering = cov::ui::tr(cov::ui::Text::RenderingSection, language);
        if (!title || !*title || !rendering || !*rendering ||
            !cov::ui::language_name(language) || !*cov::ui::language_name(language)) {
            std::cerr << "missing localisation string\n";
            ImGui::DestroyContext();
            return 2;
        }
    }

    if (std::strcmp(cov::ui::tr(cov::ui::Text::Load,
                                cov::ui::Language::ChineseSimplified),
                    "加载") != 0) {
        std::cerr << "UTF-8 localisation mismatch\n";
        ImGui::DestroyContext();
        return 3;
    }

    ImGui::DestroyContext();
    std::cout << "ui_smoke ok\n";
    return 0;
}
