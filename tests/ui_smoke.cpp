#include "cov/ui.hpp"
#include "cov/nbo_aomo_ui.hpp"

#include <imgui.h>

#include <cstring>
#include <initializer_list>
#include <iostream>
#include <string_view>

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

bool expect_utf8_text(const ImFont* font,const char* seed,const char* text,
                      const char* description,bool check_atlas) {
    const std::string_view available(seed);
    const auto* p=reinterpret_cast<const unsigned char*>(text);
    while(*p) {
        const auto* start=p;
        unsigned int codepoint=0;
        std::size_t bytes=0;
        if(*p<0x80){codepoint=*p;bytes=1;}
        else if((*p&0xe0)==0xc0){codepoint=*p&0x1f;bytes=2;}
        else if((*p&0xf0)==0xe0){codepoint=*p&0x0f;bytes=3;}
        else if((*p&0xf8)==0xf0){codepoint=*p&0x07;bytes=4;}
        else return false;
        for(std::size_t i=1;i<bytes;++i) {
            if((p[i]&0xc0)!=0x80)return false;
            codepoint=(codepoint<<6)|(p[i]&0x3f);
        }
        p+=bytes;
        if(available.find(std::string_view(reinterpret_cast<const char*>(start),bytes))==
           std::string_view::npos ||
           (check_atlas && (codepoint>0xffff || !font_contains(font,
               static_cast<ImWchar>(codepoint))))) {
            std::cerr<<"font seed or atlas missing "<<description<<" codepoint U+"
                     <<std::hex<<codepoint<<std::dec<<'\n';
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
        0x03B1, // α: scoped MO member spin
        0x03B2, // β: scoped MO member spin
        0x03A0, // Π
        0x2075, // ⁵
        0x2086, // ₆
        0x03C3, // σ
        0x03C0, // π
        0x03B4, // δ
        0x03C6, // φ
        0x2212, // −: signed phase legend
        0x2192, // →: NPA colour legend
        0x2013, // –: AO/NAO–MO title
        0x00B2, // ²: numerical component weight
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
        0x5173, 0x95ED, 0x8BE6, 0x60C5, // 关闭详情
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
        0x8A73, 0x7D30, 0x9589,                                  // 詳細・閉じる
        0x8A2D,                                                   // 設: MO 図の設定
    };
    if (std::strstr(cov::ui::font_status(), "JA fallback missing") == nullptr &&
        !expect_glyphs(primary, japanese_glyphs,
                       sizeof(japanese_glyphs) / sizeof(japanese_glyphs[0]),
                       "Japanese UI")) {
        ImGui::DestroyContext();
        return 11;
    }

    // These strings are drawn by the integration graph and scene toolbar,
    // outside kStrings. Verify both exact seed membership and real atlas glyphs.
    const bool have_zh=std::strstr(cov::ui::font_status(),"ZH fallback missing")==nullptr;
    for(const char* phrase:{"相位：红 + / 蓝 −","完整带符号系数",
                            "轨道不透明度（降低可看清原子与键）",
                            "灰色点划线：连接证据存在，整数键型未确定；点选查看轨道详情和键级指数。",
                            "整体 AO/NAO–MO 图","全部带符号 NAO–NHO 系数",
                            "MO 图设置","正则 MO 能量参考图",
                            "波函数文件或计算目录（FCHK 优先；自动关联 NBO）"})
        if(!expect_utf8_text(primary,zh_seed,phrase,"Chinese integration",have_zh)){
            ImGui::DestroyContext();return 12;
        }
    const bool have_ja=std::strstr(cov::ui::font_status(),"JA fallback missing")==nullptr;
    for(const char* phrase:{"位相：赤 + / 青 −","軌道の不透明度（下げると原子と結合が見えます）",
                            "全体 AO/NAO–MO 図","すべての符号付き NAO–NHO 係数",
                            "MO 図の設定","正準 MO エネルギー参照図",
                            "波動関数ファイルまたは計算フォルダー（FCHK 優先・NBO 自動関連付け）"})
        if(!expect_utf8_text(primary,ja_seed,phrase,"Japanese integration",have_ja)){
            ImGui::DestroyContext();return 13;
        }
    if((have_zh && !expect_utf8_text(primary,zh_seed,zh_seed,"Chinese atlas seed",true)) ||
       (have_ja && !expect_utf8_text(primary,ja_seed,ja_seed,"Japanese atlas seed",true))) {
        ImGui::DestroyContext();return 14;
    }

    // The unified diagram owns its evolving labels. Exercise actual merged
    // font glyphs so a new preset or scientific symbol cannot become a '?' UI.
    for(const auto language:languages) {
        const auto seed=cov::ui::nbo_aomo_glyph_seed(language);
        const bool available=language==cov::ui::Language::ChineseSimplified?have_zh:
            language==cov::ui::Language::Japanese?have_ja:true;
        if(available&&!expect_utf8_text(primary,seed.c_str(),seed.c_str(),
                                       "unified orbital diagram",true)) {
            ImGui::DestroyContext();return 15;
        }
    }

    ImGui::DestroyContext();
    std::cout << "ui_smoke ok\n";
    return 0;
}
