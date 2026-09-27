#include "cov/ui.hpp"

namespace cov::ui {

const char* tr_legacy(Text key, Language language) noexcept;

const char* tr(Text key, Language language) noexcept {
    switch (key) {
        case Text::MoldenPath:
            switch (language) {
                case Language::ChineseSimplified: return "波函数文件（FCHK 优先；Molden 兼容）";
                case Language::Japanese: return "波動関数ファイル（FCHK 優先・Molden 互換）";
                case Language::French: return "Fichier de fonction d’onde (FCHK prioritaire ; Molden compatible)";
                default: return "Wavefunction file (FCHK preferred; Molden compatible)";
            }
        case Text::IdleHint:
            switch (language) {
                case Language::ChineseSimplified: return "可拖入 .fchk/.fch/.chk 或兼容的 .molden 文件，也可直接输入路径。";
                case Language::Japanese: return ".fchk/.fch/.chk または互換 .molden ファイルをドロップするか、パスを入力してください。";
                case Language::French: return "Déposez un fichier .fchk/.fch/.chk ou .molden compatible, ou saisissez son chemin.";
                default: return "Drop a .fchk/.fch/.chk or compatible .molden file, or enter a path.";
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
