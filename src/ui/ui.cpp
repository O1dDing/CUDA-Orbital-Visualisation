#include "cov/ui.hpp"

#include <imgui.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <initializer_list>
#include <string>

namespace cov::ui {
namespace {

constexpr std::size_t kTextCount = static_cast<std::size_t>(Text::Count);

struct LocalisedString {
    const char* en;
    const char* zh;
    const char* ja;
    const char* fr;
};

constexpr std::array<LocalisedString, kTextCount> kStrings{{
    {"CUDA Orbital Visualisation", "CUDA Orbital Visualisation", "CUDA Orbital Visualisation", "CUDA Orbital Visualisation"},
    {"GPU-first molecular orbital viewer", "GPU 优先的分子轨道可视化", "GPUファースト分子軌道ビューア", "Visualiseur d’orbitales moléculaires orienté GPU"},
    {"Language", "语言", "言語", "Langue"},
    {"File", "文件", "ファイル", "Fichier"},
    {"Molden file", "Molden 文件", "Molden ファイル", "Fichier Molden"},
    {"Load", "加载", "読み込む", "Charger"},
    {"Drop a .molden file anywhere on the viewport, or enter a path.", "可将 .molden 文件拖到视口任意位置，或直接输入路径。", ".molden ファイルをビューポートへドロップするか、パスを入力してください。", "Déposez un fichier .molden dans la vue ou saisissez son chemin."},
    {"Wavefunction", "波函数", "波動関数", "Fonction d’onde"},
    {"Atoms", "原子数", "原子数", "Atomes"},
    {"Shells", "壳层数", "シェル数", "Couches"},
    {"Basis functions", "基函数数", "基底関数数", "Fonctions de base"},
    {"Orbitals", "轨道数", "軌道数", "Orbitales"},
    {"Shell convention", "壳层约定", "シェル規約", "Convention des couches"},
    {"Orbital", "轨道", "軌道", "Orbitale"},
    {"Molden MO (1-based)", "Molden MO（从 1 开始）", "Molden MO（1 始まり）", "MO Molden (base 1)"},
    {"Internal index", "内部索引", "内部インデックス", "Indice interne"},
    {"Energy", "能量", "エネルギー", "Énergie"},
    {"Occupation", "占据数", "占有数", "Occupation"},
    {"Spin", "自旋", "スピン", "Spin"},
    {"Symmetry", "对称性", "対称性", "Symétrie"},
    {"Rendering", "渲染", "レンダリング", "Rendu"},
    {"Isovalue", "等值面", "等値面", "Isovaleur"},
    {"Grid", "网格", "グリッド", "Grille"},
    {"Recompute grid", "重新计算网格", "グリッドを再計算", "Recalculer la grille"},
    {"Reset camera", "重置相机", "カメラをリセット", "Réinitialiser la caméra"},
    {"Performance", "性能", "パフォーマンス", "Performances"},
    {"CUDA device", "CUDA 设备", "CUDA デバイス", "Périphérique CUDA"},
    {"Last CUDA kernel", "最近 CUDA 内核", "直近の CUDA カーネル", "Dernier noyau CUDA"},
    {"GPU resident", "GPU 常驻", "GPU 常駐", "Résident GPU"},
    {"Left-drag to orbit · mouse wheel to zoom", "按住鼠标左键旋转 · 滚轮缩放", "左ドラッグで回転 · ホイールでズーム", "Glisser gauche : rotation · molette : zoom"},
    {"Changing isovalue is instant and does not recompute the CUDA grid.", "调整等值面不会重新计算 CUDA 网格，可即时更新。", "等値面の変更では CUDA グリッドを再計算せず、即時更新します。", "Changer l’isovaleur est instantané et ne recalcule pas la grille CUDA."},
    {"Experimental MVP · scientific validation in progress", "实验性 MVP · 科学数值验证仍在进行", "実験的 MVP · 科学的検証を継続中", "MVP expérimental · validation scientifique en cours"},
    {"Ready", "就绪", "準備完了", "Prêt"},
    {"Parsing", "正在解析", "解析中", "Analyse"},
    {"Loaded", "已加载", "読み込み完了", "Chargé"},
    {"Grid updated", "网格已更新", "グリッド更新完了", "Grille mise à jour"},
    {"Error", "错误", "エラー", "Erreur"},
    {"(none)", "（无）", "（なし）", "(aucun)"},
    {"Alpha", "Alpha", "Alpha", "Alpha"},
    {"Beta", "Beta", "Beta", "Beta"},
    {"UI fonts", "界面字体", "UI フォント", "Polices UI"},

    {"Open File…", "打开文件…", "ファイルを開く…", "Ouvrir un fichier…"},
    {"Recent files", "最近文件", "最近のファイル", "Fichiers récents"},
    {"Current file", "当前文件", "現在のファイル", "Fichier actuel"},
    {"Orbital browser", "轨道浏览器", "軌道ブラウザ", "Explorateur d’orbitales"},
    {"Search", "搜索", "検索", "Rechercher"},
    {"Filter", "筛选", "フィルター", "Filtre"},
    {"Auto · reasonable", "自动 · 合理范围", "自動 · 妥当範囲", "Auto · plage raisonnable"},
    {"All", "全部", "すべて", "Toutes"},
    {"Occupied", "已占据", "占有", "Occupées"},
    {"Virtual", "虚轨道", "仮想", "Virtuelles"},
    {"Core", "内层", "内殻", "Cœur"},
    {"Valence", "价层", "価電子", "Valence"},
    {"Virtual window", "虚轨道窗口", "仮想軌道ウィンドウ", "Fenêtre virtuelle"},
    {"Degeneracy tolerance", "简并阈值", "縮退判定しきい値", "Tolérance de dégénérescence"},
    {"Grouped labels", "分组标签", "グループ表示", "Étiquettes groupées"},
    {"Raw numbering", "原始编号", "元の番号", "Numérotation brute"},
    {"Degenerate set", "简并组", "縮退組", "Groupe dégénéré"},
    {"Energy unit", "能量单位", "エネルギー単位", "Unité d’énergie"},
    {"HOMO", "HOMO", "HOMO", "HOMO"},
    {"LUMO", "LUMO", "LUMO", "LUMO"},
    {"HOMO-1", "HOMO-1", "HOMO-1", "HOMO-1"},
    {"LUMO+1", "LUMO+1", "LUMO+1", "LUMO+1"},
    {"Energy-level diagram", "能级图", "エネルギー準位図", "Diagramme des niveaux d’énergie"},
    {"Around selected", "围绕当前轨道", "選択軌道の周辺", "Autour de la sélection"},
    {"Generate MO diagram", "生成 MO 图", "MO 図を生成", "Générer le diagramme MO"},
    {"Export diagram + metadata", "导出图与元数据", "図とメタデータを書き出す", "Exporter diagramme + métadonnées"},
    {"Exported", "已导出", "書き出し完了", "Exporté"},
    {"Export failed", "导出失败", "書き出し失敗", "Échec de l’export"},
    {"Molecule style", "分子样式", "分子表示", "Style moléculaire"},
    {"Medium ball + thin stick", "中号球 + 细棍", "中サイズ球 + 細い結合", "Boules moyennes + bâtons fins"},
    {"Stick + delocalisation", "棍 + 离域虚线", "結合 + 非局在化破線", "Bâtons + délocalisation"},
    {"Atom size", "原子大小", "原子サイズ", "Taille des atomes"},
    {"Bond size", "键粗细", "結合の太さ", "Épaisseur des liaisons"},
    {"Molecule opacity", "分子透明度", "分子の不透明度", "Opacité de la molécule"},
    {"Orbital opacity", "轨道透明度", "軌道の不透明度", "Opacité de l’orbitale"},
    {"Show hydrogens", "显示氢原子", "水素を表示", "Afficher les hydrogènes"},
    {"Dashed bonds use a conservative delocalisation heuristic.", "虚线键使用保守的离域启发式判断。", "破線結合は保守的な非局在化ヒューリスティックです。", "Les liaisons en pointillés utilisent une heuristique prudente de délocalisation."},
    {"SALC not confidently available", "无法可靠获得 SALC", "SALC を信頼して決定できません", "SALC non disponible avec confiance"},
    {"Symmetry-grouped", "按对称性分组", "対称性でグループ化", "Groupé par symétrie"},
    {"Simple ordering", "简单能级排序", "単純な準位順序", "Ordre simple"},
    {"Machine metadata", "机读元数据", "機械可読メタデータ", "Métadonnées machine"},
    {"Visible orbitals", "可见轨道", "表示軌道", "Orbitales visibles"},
    {"Raw MO", "原始 MO", "元の MO", "MO brute"},
    {"Region", "区域", "領域", "Région"},
    {"Core", "内层", "内殻", "Cœur"},
    {"Valence", "价层", "価電子", "Valence"},
    {"Virtual", "虚轨道", "仮想", "Virtuelle"},
    {"Reasonable energy window", "合理能量窗口", "妥当なエネルギー範囲", "Fenêtre d’énergie raisonnable"},
    {"Native Open File is unavailable on this platform.", "当前平台不支持原生“打开文件”。", "このプラットフォームではネイティブのファイル選択を利用できません。", "La boîte de dialogue native n’est pas disponible sur cette plateforme."},
    {"Copy metadata", "复制元数据", "メタデータをコピー", "Copier les métadonnées"},
    {"No orbitals", "无轨道", "軌道がありません", "Aucune orbitale"},
}};

std::string g_font_status = "Dear ImGui default";

std::string first_existing(std::initializer_list<const char*> candidates) {
    std::error_code ec;
    for (const char* candidate : candidates) {
        if (candidate && std::filesystem::exists(candidate, ec) && !ec) {
            return candidate;
        }
        ec.clear();
    }
    return {};
}

const char* file_name_or_default(const std::string& path, const char* fallback) {
    if (path.empty()) return fallback;
    static thread_local std::string name;
    name = std::filesystem::path(path).filename().string();
    return name.c_str();
}

void merge_font_if_available(const std::string& path,
                             const float pixel_size,
                             const ImWchar* ranges,
                             bool& loaded) {
    if (path.empty()) return;
    ImFontConfig cfg{};
    cfg.MergeMode = true;
    cfg.PixelSnapH = true;
    cfg.OversampleH = 1;
    cfg.OversampleV = 1;
    if (ImGui::GetIO().Fonts->AddFontFromFileTTF(path.c_str(), pixel_size, &cfg, ranges)) {
        loaded = true;
    }
}

const char* localised(const LocalisedString& value, const Language language) noexcept {
    switch (language) {
        case Language::ChineseSimplified: return value.zh;
        case Language::Japanese: return value.ja;
        case Language::French: return value.fr;
        default: return value.en;
    }
}

} // namespace

const char* tr(const Text key, const Language language) noexcept {
    const auto k = static_cast<std::size_t>(key);
    if (k >= kTextCount) return "";
    return localised(kStrings[k], language);
}

const char* language_name(const Language language) noexcept {
    switch (language) {
        case Language::ChineseSimplified: return "简体中文";
        case Language::Japanese: return "日本語";
        case Language::French: return "Français";
        default: return "English";
    }
}

void apply_theme(const float scale) {
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();

    style.WindowPadding = ImVec2(14.0f, 14.0f);
    style.FramePadding = ImVec2(10.0f, 6.0f);
    style.CellPadding = ImVec2(8.0f, 5.0f);
    style.ItemSpacing = ImVec2(9.0f, 8.0f);
    style.ItemInnerSpacing = ImVec2(7.0f, 5.0f);
    style.IndentSpacing = 18.0f;
    style.ScrollbarSize = 11.0f;
    style.GrabMinSize = 10.0f;

    style.WindowRounding = 10.0f;
    style.ChildRounding = 8.0f;
    style.FrameRounding = 7.0f;
    style.PopupRounding = 8.0f;
    style.ScrollbarRounding = 10.0f;
    style.GrabRounding = 6.0f;
    style.TabRounding = 6.0f;

    style.WindowBorderSize = 1.0f;
    style.ChildBorderSize = 1.0f;
    style.FrameBorderSize = 0.0f;
    style.PopupBorderSize = 1.0f;

    ImVec4* c = style.Colors;
    c[ImGuiCol_Text] = ImVec4(0.90f, 0.93f, 0.97f, 1.00f);
    c[ImGuiCol_TextDisabled] = ImVec4(0.53f, 0.59f, 0.67f, 1.00f);
    c[ImGuiCol_WindowBg] = ImVec4(0.055f, 0.073f, 0.102f, 0.965f);
    c[ImGuiCol_ChildBg] = ImVec4(0.082f, 0.108f, 0.145f, 0.96f);
    c[ImGuiCol_PopupBg] = ImVec4(0.070f, 0.091f, 0.122f, 0.985f);
    c[ImGuiCol_Border] = ImVec4(0.18f, 0.23f, 0.30f, 0.95f);
    c[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_FrameBg] = ImVec4(0.112f, 0.145f, 0.195f, 1.00f);
    c[ImGuiCol_FrameBgHovered] = ImVec4(0.145f, 0.195f, 0.270f, 1.00f);
    c[ImGuiCol_FrameBgActive] = ImVec4(0.175f, 0.230f, 0.315f, 1.00f);
    c[ImGuiCol_TitleBg] = c[ImGuiCol_WindowBg];
    c[ImGuiCol_TitleBgActive] = c[ImGuiCol_WindowBg];
    c[ImGuiCol_TitleBgCollapsed] = c[ImGuiCol_WindowBg];
    c[ImGuiCol_MenuBarBg] = ImVec4(0.070f, 0.091f, 0.122f, 1.00f);
    c[ImGuiCol_ScrollbarBg] = ImVec4(0.045f, 0.060f, 0.083f, 0.80f);
    c[ImGuiCol_ScrollbarGrab] = ImVec4(0.21f, 0.27f, 0.35f, 0.95f);
    c[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.30f, 0.39f, 0.50f, 1.00f);
    c[ImGuiCol_ScrollbarGrabActive] = ImVec4(0.36f, 0.46f, 0.59f, 1.00f);
    c[ImGuiCol_CheckMark] = ImVec4(0.36f, 0.56f, 0.98f, 1.00f);
    c[ImGuiCol_SliderGrab] = ImVec4(0.36f, 0.56f, 0.98f, 0.88f);
    c[ImGuiCol_SliderGrabActive] = ImVec4(0.48f, 0.67f, 1.00f, 1.00f);
    c[ImGuiCol_Button] = ImVec4(0.18f, 0.32f, 0.58f, 0.82f);
    c[ImGuiCol_ButtonHovered] = ImVec4(0.26f, 0.43f, 0.78f, 1.00f);
    c[ImGuiCol_ButtonActive] = ImVec4(0.22f, 0.37f, 0.68f, 1.00f);
    c[ImGuiCol_Header] = ImVec4(0.18f, 0.32f, 0.58f, 0.62f);
    c[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.43f, 0.78f, 0.86f);
    c[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.43f, 0.78f, 1.00f);
    c[ImGuiCol_Separator] = ImVec4(0.18f, 0.23f, 0.30f, 0.80f);
    c[ImGuiCol_SeparatorHovered] = ImVec4(0.36f, 0.56f, 0.98f, 0.75f);
    c[ImGuiCol_SeparatorActive] = ImVec4(0.36f, 0.56f, 0.98f, 1.00f);
    c[ImGuiCol_ResizeGrip] = ImVec4(0.36f, 0.56f, 0.98f, 0.20f);
    c[ImGuiCol_ResizeGripHovered] = ImVec4(0.36f, 0.56f, 0.98f, 0.65f);
    c[ImGuiCol_ResizeGripActive] = ImVec4(0.36f, 0.56f, 0.98f, 0.90f);
    c[ImGuiCol_Tab] = ImVec4(0.10f, 0.14f, 0.20f, 1.00f);
    c[ImGuiCol_TabHovered] = ImVec4(0.26f, 0.43f, 0.78f, 0.85f);
    c[ImGuiCol_TabActive] = ImVec4(0.18f, 0.32f, 0.58f, 0.95f);
    c[ImGuiCol_TabUnfocused] = ImVec4(0.08f, 0.11f, 0.15f, 1.00f);
    c[ImGuiCol_TabUnfocusedActive] = ImVec4(0.13f, 0.20f, 0.31f, 1.00f);
    c[ImGuiCol_TableHeaderBg] = ImVec4(0.10f, 0.14f, 0.19f, 1.00f);
    c[ImGuiCol_TableBorderStrong] = ImVec4(0.18f, 0.23f, 0.30f, 0.85f);
    c[ImGuiCol_TableBorderLight] = ImVec4(0.15f, 0.19f, 0.25f, 0.65f);
    c[ImGuiCol_TextSelectedBg] = ImVec4(0.36f, 0.56f, 0.98f, 0.32f);
    c[ImGuiCol_DragDropTarget] = ImVec4(0.31f, 0.82f, 0.65f, 0.95f);
    c[ImGuiCol_NavHighlight] = ImVec4(0.36f, 0.56f, 0.98f, 0.90f);
    c[ImGuiCol_ModalWindowDimBg] = ImVec4(0.02f, 0.03f, 0.05f, 0.70f);

    style.ScaleAllSizes(std::clamp(scale, 0.85f, 2.0f));
}

bool configure_fonts(const float pixel_size) {
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();
    io.Fonts->Flags |= ImFontAtlasFlags_NoPowerOfTwoHeight;

#ifdef _WIN32
    const std::string base = first_existing({
        "C:/Windows/Fonts/segoeui.ttf",
        "C:/Windows/Fonts/arial.ttf",
    });
    const std::string chinese = first_existing({
        "C:/Windows/Fonts/msyh.ttc",
        "C:/Windows/Fonts/msyh.ttf",
        "C:/Windows/Fonts/simhei.ttf",
    });
    const std::string japanese = first_existing({
        "C:/Windows/Fonts/YuGothR.ttc",
        "C:/Windows/Fonts/meiryo.ttc",
        "C:/Windows/Fonts/msgothic.ttc",
    });
#elif defined(__APPLE__)
    const std::string base = first_existing({
        "/System/Library/Fonts/SFNS.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
    });
    const std::string chinese = first_existing({
        "/System/Library/Fonts/PingFang.ttc",
        "/Library/Fonts/NotoSansCJK-Regular.ttc",
    });
    const std::string japanese = first_existing({
        "/Library/Fonts/NotoSansCJK-Regular.ttc",
        "/System/Library/Fonts/AppleGothic.ttf",
    });
#else
    const std::string base = first_existing({
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
    });
    const std::string chinese = first_existing({
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/opentype/noto/NotoSansCJKsc-Regular.otf",
    });
    const std::string japanese = first_existing({
        "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
        "/usr/share/fonts/opentype/noto/NotoSansCJKjp-Regular.otf",
    });
#endif

    ImFontGlyphRangesBuilder latin_builder;
    latin_builder.AddRanges(io.Fonts->GetGlyphRangesDefault());
    for (const auto& row : kStrings) {
        latin_builder.AddText(row.en);
        latin_builder.AddText(row.fr);
    }
    latin_builder.AddText(language_name(Language::French));
    latin_builder.AddText("● · ³ ↑ ↓");
    ImVector<ImWchar> latin_ranges;
    latin_builder.BuildRanges(&latin_ranges);

    ImFont* primary = nullptr;
    if (!base.empty()) {
        ImFontConfig cfg{};
        cfg.OversampleH = 2;
        cfg.OversampleV = 1;
        primary = io.Fonts->AddFontFromFileTTF(base.c_str(), pixel_size, &cfg,
                                               latin_ranges.Data);
    }
    if (!primary) primary = io.Fonts->AddFontDefault();

    ImFontGlyphRangesBuilder zh_builder;
    ImFontGlyphRangesBuilder ja_builder;
    for (const auto& row : kStrings) {
        zh_builder.AddText(row.zh);
        ja_builder.AddText(row.ja);
    }
    zh_builder.AddText(language_name(Language::ChineseSimplified));
    ja_builder.AddText(language_name(Language::Japanese));

    ImVector<ImWchar> zh_ranges;
    ImVector<ImWchar> ja_ranges;
    zh_builder.BuildRanges(&zh_ranges);
    ja_builder.BuildRanges(&ja_ranges);

    bool zh_loaded = false;
    bool ja_loaded = false;
    merge_font_if_available(chinese, pixel_size, zh_ranges.Data, zh_loaded);
    merge_font_if_available(japanese, pixel_size, ja_ranges.Data, ja_loaded);

    const bool built = io.Fonts->Build();
    g_font_status = file_name_or_default(base, "ImGui default");
    g_font_status += " + ";
    g_font_status += zh_loaded ? file_name_or_default(chinese, "CJK") : "ZH fallback missing";
    g_font_status += " + ";
    g_font_status += ja_loaded ? file_name_or_default(japanese, "CJK") : "JA fallback missing";
    return built;
}

const char* font_status() noexcept {
    return g_font_status.c_str();
}

void section_title(const char* label) {
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.58f, 0.68f, 0.82f, 1.0f));
    ImGui::TextUnformatted(label);
    ImGui::PopStyleColor();
}

void begin_card(const char* id, const float height) {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(12.0f, 11.0f));
    ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.082f, 0.108f, 0.145f, 0.94f));
    ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.18f, 0.23f, 0.30f, 0.90f));
    ImGui::BeginChild(id, ImVec2(0.0f, height), true, ImGuiWindowFlags_None);
}

void end_card() {
    ImGui::EndChild();
    ImGui::PopStyleColor(2);
    ImGui::PopStyleVar();
}

void status_badge(const char* label, const Tone tone) {
    ImVec4 colour{};
    switch (tone) {
        case Tone::Success: colour = ImVec4(0.31f, 0.82f, 0.65f, 1.0f); break;
        case Tone::Danger: colour = ImVec4(0.94f, 0.42f, 0.47f, 1.0f); break;
        case Tone::Accent: colour = ImVec4(0.42f, 0.64f, 1.00f, 1.0f); break;
        default: colour = ImVec4(0.58f, 0.68f, 0.82f, 1.0f); break;
    }
    ImGui::TextColored(colour, "●  %s", label);
}

} // namespace cov::ui
