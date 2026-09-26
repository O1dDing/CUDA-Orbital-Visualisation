#include "cov/ui.hpp"

namespace cov::ui {

const char* tr_legacy(Text key, Language language) noexcept;

const char* tr(Text key, Language language) noexcept {
    switch (key) {
        case Text::MoldenPath:
            switch (language) {
                case Language::ChineseSimplified: return "波函数文件或计算目录（FCHK 优先；自动关联 NBO）";
                case Language::Japanese: return "波動関数ファイルまたは計算フォルダー（FCHK 優先・NBO 自動関連付け）";
                case Language::French: return "Fichier de fonction d’onde ou dossier de calcul (FCHK prioritaire ; association NBO automatique)";
                default: return "Wavefunction file or calculation folder (FCHK preferred; automatic NBO association)";
            }
        case Text::IdleHint:
            switch (language) {
                case Language::ChineseSimplified: return "可同时拖入波函数与 NBO 文件，或拖入计算目录；也可直接输入文件或目录路径。兼容 .fchk/.fch/.chk 和 .molden。";
                case Language::Japanese: return "波動関数と NBO ファイルをまとめて、または計算フォルダーをドロップできます。ファイルやフォルダーのパス入力も可能です。.fchk/.fch/.chk・.molden に対応。";
                case Language::French: return "Déposez ensemble les fichiers de fonction d’onde et NBO, ou un dossier de calcul ; vous pouvez aussi saisir leur chemin. Formats .fchk/.fch/.chk et .molden compatibles.";
                default: return "Drop wavefunction and NBO files together, or a calculation folder; you can also enter a file or folder path. Supported: .fchk/.fch/.chk and .molden.";
            }
        case Text::MoldenMO:
            switch (language) {
                case Language::ChineseSimplified: return "源文件 MO（从 1 开始）";
                case Language::Japanese: return "入力 MO（1 始まり）";
                case Language::French: return "MO source (base 1)";
                default: return "Source MO (1-based)";
            }
        case Text::RawMO:
            switch (language) {
                case Language::ChineseSimplified: return "源 MO";
                case Language::Japanese: return "入力 MO";
                case Language::French: return "MO source";
                default: return "Source MO";
            }
        default:
            return tr_legacy(key,language);
    }
}

} // namespace cov::ui
