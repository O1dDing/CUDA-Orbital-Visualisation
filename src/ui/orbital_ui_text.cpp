#include "cov/orbital_ui_text.hpp"

#include <array>
#include <sstream>
#include <string>
#include <string_view>

namespace cov::ui {
namespace {

struct LocalisedString {
    const char* en;
    const char* zh;
    const char* ja;
    const char* fr;
};

constexpr auto kStrings = std::to_array<LocalisedString>({
    {"symmetry", "对称性", "対称性", "symétrie"},
    {"occupation", "占据数", "占有数", "occupation"},
    {"density", "密度", "密度", "densité"},
    {"overlap", "重叠", "重なり", "recouvrement"},
    {"bond order", "键级", "結合次数", "ordre de liaison"},
    {"point group", "点群", "点群", "groupe ponctuel"},
    {"Gaussian LOG/OUT enrichment attached", "已附加 Gaussian LOG/OUT 增补数据", "Gaussian LOG/OUT 補足データを適用済み", "Enrichissement Gaussian LOG/OUT associé"},
    {"producer", "来源数据", "生成元データ", "producteur"},
    {"derived", "推导数据", "導出データ", "dérivé"},
    {"unavailable", "不可用", "利用不可", "indisponible"},
    {"unknown", "未知", "不明", "inconnue"},

    {"Selected MO chemistry", "所选 MO 化学性质", "選択 MO の化学的性質", "Caractère chimique de l’OM sélectionnée"},
    {"Chemical-valence manifold", "化学价层轨道组", "化学原子価軌道空間", "Espace orbitalaire de valence chimique"},
    {"yes", "是", "はい", "oui"},
    {"no", "否", "いいえ", "non"},
    {"Valence AO composition", "价层 AO 组成", "原子価 AO 組成", "Composition AO de valence"},
    {"Atom-pair interactions", "原子对相互作用", "原子対相互作用", "Interactions par paire atomique"},
    {"Orbital family", "轨道类型", "軌道型", "Famille orbitale"},
    {"Bonding role", "成键性质", "結合性", "Caractère liant"},
    {"Multicentre family", "多中心族", "多中心族", "Famille multicentrique"},
    {"Delocalised pi family", "离域 π 族", "非局在化 π 族", "Famille π délocalisée"},
    {"Member MOs", "成员 MO", "構成 MO", "OM membres"},
    {"Participating atoms", "参与原子", "参加原子", "Atomes participants"},
    {"Participating electrons", "参与电子", "参加電子", "Électrons participants"},
    {"Donor / acceptor direction", "给体 / 受体方向", "供与体 / 受容体方向", "Direction donneur / accepteur"},
    {"Analysis method", "分析方法", "解析法", "Méthode d’analyse"},
    {"MO contribution", "MO 贡献", "MO 寄与", "Contribution OM"},
    {"MO contribution is overlap-population based; Mayer is the total density-level pair index.", "MO 贡献来自重叠布居；Mayer 是总密度层级的原子对指数。", "MO 寄与は重なり密度に基づき、Mayer は全密度レベルの原子対指数です。", "La contribution OM repose sur la population de recouvrement ; Mayer est l’indice de paire au niveau de la densité totale."},
    {"UND / outside minimal valence reference", "UND / 最小价层参考之外", "UND / 最小原子価参照外", "UND / hors référence de valence minimale"},
    {"UND — AO metric unavailable", "UND — AO 度量不可用", "UND — AO 計量を利用できません", "UND — métrique AO indisponible"},
    {"delocalised-pi", "离域 π", "非局在化 π", "π délocalisée"},
    {"mixed", "混合", "混合", "mixte"},
    {"mixed (UND)", "混合（UND）", "混合（UND）", "mixte (UND)"},
    {"bonding", "成键", "結合性", "liante"},
    {"antibonding", "反键", "反結合性", "antiliante"},
    {"nonbonding", "非键", "非結合性", "non liante"},

    {"Orientation channels", "取向通道", "配向チャネル", "Canaux d’orientation"},
    {"Topology", "拓扑", "トポロジー", "Topologie"},
    {"path-delocalised pi", "链状 π 离域", "鎖状 π 非局在化", "délocalisation π en chaîne"},
    {"cyclic-delocalised pi", "环状 π 离域", "環状 π 非局在化", "délocalisation π cyclique"},
    {"branched-resonance pi network", "支化共振 π 网络", "分岐共鳴 π ネットワーク", "réseau π à résonance ramifiée"},
    {"orthogonal spiro pi network", "螺环正交 π 网络", "スピロ直交 π ネットワーク", "réseau π spiro orthogonal"},
    {"haptic metal-pi", "多触点金属–π", "多点金属–π", "contact haptique métal–π"},
    {"symmetry-equivalent pi direct sum", "对称等价 π 子体系直和", "対称等価 π 部分系の直和", "somme directe de sous-systèmes π équivalents"},
    {"atoms", "原子", "原子", "atomes"},
    {"coherence", "相干度", "コヒーレンス", "cohérence"},
    {"Hide intermediate framework MOs", "隐藏中间框架轨道", "中間骨格軌道を非表示", "Masquer les OM intermédiaires"},
    {"Hide secondary ligand-centred, ligand-internal and polarisation MOs; keep the principal d levels, necessary sigma bonding/antibonding representatives, the selected MO and donor/acceptor pairs.", "隐藏次级配体中心、配体内部及极化型轨道；保留主要 d 能级、必要的 σ 成键/反键代表、当前轨道和给体/受体配对。", "副次的な配位子中心・配位子内部・分極軌道を隠し、主要 d 準位、必要な σ 結合/反結合代表、選択軌道と供与体/受容体対を保持します。", "Masque les OM secondaires centrées sur les ligands, internes aux ligands ou de polarisation ; conserve les niveaux d principaux, les représentants σ liants/antiliants nécessaires, l’OM sélectionnée et les paires donneur/accepteur."},

    {"local ligand field", "局部配体场", "局所配位子場", "champ de ligands local"},
    {"first shell", "第一配位层", "第一配位圏", "première sphère"},
    {"coordination geometry", "配位几何", "配位幾何", "géométrie de coordination"},
    {"coordination number", "配位数", "配位数", "nombre de coordination"},
    {"geometry confidence", "几何置信度", "幾何信頼度", "confiance géométrique"},
    {"angular RMS (rad)", "角度 RMS（rad）", "角度 RMS（rad）", "RMS angulaire (rad)"},
    {"directional shape score", "方向形状评分", "方向形状スコア", "score de forme directionnelle"},
    {"radial variation", "径向变异", "動径変動", "variation radiale"},
    {"local molecular geometries", "局部分子几何数", "局所分子幾何数", "géométries moléculaires locales"},
    {"atom", "原子", "原子", "atome"},
    {"local geometry", "局部几何", "局所幾何", "géométrie locale"},
    {"not centred on this MO", "该 MO 未中心化于这些原子", "この MO は該当中心上にありません", "non centrée sur cette OM"},
    {"group occupation", "轨道组占据数", "軌道群占有数", "occupation du groupe"},
    {"metal s / p / d", "金属 s / p / d", "金属 s / p / d", "métal s / p / d"},
    {"ligand p", "配体 p", "配位子 p", "ligand p"},
    {"sigma / pi channel", "σ / π 通道", "σ / π チャネル", "canal σ / π"},
    {"M–L overlap", "M–L 重叠", "M–L 重なり", "recouvrement M–L"},
    {"selection", "选取方式", "選択", "sélection"},
    {"recovered from complete raw MO block", "从完整原始 MO 数据块恢复", "完全な生 MO ブロックから復元", "récupérée depuis le bloc complet des OM brutes"},
    {"weak-field treatment", "弱场处理", "弱場処理", "traitement de champ faible"},
    {"approximately nonbonding", "近似非键", "近似的に非結合性", "approximativement non liante"},
    {"pi interaction", "π 相互作用", "π 相互作用", "interaction π"},
    {"pi splitting", "π 劈裂", "π 分裂", "séparation π"},
    {"pair confidence", "配对置信度", "対の信頼度", "confiance de la paire"},
    {"interaction", "相互作用", "相互作用", "interaction"},
    {"splitting", "劈裂", "分裂", "séparation"},
    {"splitting (Ha)", "劈裂（Ha）", "分裂（Ha）", "séparation (Ha)"},

    {"direct", "直接数据", "直接データ", "directe"},
    {"parsed label", "标签解析", "ラベル解析", "étiquette interprétée"},
    {"derived", "推导", "導出", "dérivée"},
    {"heuristic", "启发式", "ヒューリスティック", "heuristique"},
    {"pi-donor splitting", "π 给体劈裂", "π ドナー分裂", "séparation π-donneuse"},
    {"pi-acceptor splitting", "π 受体劈裂", "π アクセプター分裂", "séparation π-acceptrice"},
    {"pi-coupled splitting", "π 耦合劈裂", "π 相互作用分裂", "séparation π couplée"},
    {"weak-field split; approximately nonbonding", "弱场劈裂；近似非键", "弱場分裂；近似的に非結合性", "séparation de champ faible ; approximativement non liante"},

    {"MO diagram summary", "MO 图摘要", "MO 図概要", "Résumé du diagramme OM"},
    {"ligand-field valence groups", "配体场价层组", "配位子場原子価群", "groupes de valence du champ de ligands"},
    {"delocalised pi family groups", "离域 π 族组", "非局在化 π 族群", "groupes de la famille π délocalisée"},
    {"multicentre active-space groups", "多中心活性空间组", "多中心活性空間群", "groupes de l’espace actif multicentrique"},
    {"visible levels", "可见能级", "表示準位", "niveaux visibles"},
    {"occupied", "占据", "占有", "occupés"},
    {"virtual", "虚轨道", "仮想", "virtuels"},
    {"hidden", "隐藏", "非表示", "masqués"},
    {"local field", "局部场", "局所場", "champ local"},
    {"geometry", "几何", "幾何", "géométrie"},
    {"pi pairs", "π 配对", "π 対", "paires π"},
    {"protected overflow", "受保护行溢出", "保護行の超過", "dépassement protégé"},
    {"spin pairs", "自旋配对", "スピン対", "paires de spin"},
    {"unmatched visible", "可见未配对", "表示中の未対応", "visibles non appariées"},
    {"raw-MO recovered groups", "原始 MO 恢复组", "生 MO 復元群", "groupes récupérés des OM brutes"},
    {"compact active space unavailable", "紧凑活性空间不可用", "コンパクト活性空間を利用できません", "espace actif compact indisponible"},
});

static_assert(kStrings.size() == static_cast<std::size_t>(OrbitalText::Count));

constexpr auto kGeometryNames = std::to_array<std::pair<std::string_view, LocalisedString>>({
    {"L-2", {"Linear", "线形", "直線形", "linéaire"}},
    {"A-2", {"Angular", "折线形", "折れ線形", "coudée"}},
    {"TP-3", {"Trigonal planar", "平面三角形", "三角平面形", "trigonale plane"}},
    {"TPY-3", {"Trigonal pyramidal", "三角锥形", "三角錐形", "pyramidale trigonale"}},
    {"TS-3", {"T-shaped", "T 形", "T 字形", "en T"}},
    {"T-4", {"Tetrahedral", "四面体形", "四面体形", "tétraédrique"}},
    {"SP-4", {"Square planar", "平面正方形", "正方形平面形", "plan carré"}},
    {"SS-4", {"Seesaw", "跷跷板形", "シーソー形", "en bascule"}},
    {"vTBPY-4", {"Trigonal-pyramidal (axially vacant trigonal bipyramid)", "三角锥形（轴向空位三角双锥）", "三角錐形（軸位空孔三方両錐）", "pyramidale trigonale (bipyramide trigonale à lacune axiale)"}},
    {"TBPY-5", {"Trigonal bipyramidal", "三角双锥形", "三方両錐形", "bipyramidale trigonale"}},
    {"SPY-5", {"Square pyramidal", "四方锥形", "四角錐形", "pyramidale carrée"}},
    {"OC-6", {"Octahedral", "八面体形", "八面体形", "octaédrique"}},
    {"TPR-6", {"Trigonal prismatic", "三棱柱形", "三角柱形", "prismatique trigonale"}},
    {"PBPY-7", {"Pentagonal bipyramidal", "五角双锥形", "五方両錐形", "bipyramidale pentagonale"}},
    {"COC-7", {"Capped octahedral", "加帽八面体形", "一冠八面体形", "octaèdre coiffé"}},
    {"CTPR-7", {"Capped trigonal prismatic", "加帽三棱柱形", "一冠三角柱形", "prisme trigonal coiffé"}},
    {"SAPR-8", {"Square antiprismatic", "四方反棱柱形", "四角反柱形", "antiprismatique carrée"}},
    {"TDD-8", {"Triangular dodecahedral", "三角十二面体形", "三角十二面体形", "dodécaèdre triangulaire"}},
    {"BTPR-8", {"Bicapped trigonal prismatic", "双帽三棱柱形", "二冠三角柱形", "prisme trigonal bicoiffé"}},
    {"CSAPR-9", {"Capped square antiprismatic", "加帽四方反棱柱形", "一冠四角反柱形", "antiprisme carré coiffé"}},
    {"TCTPR-9", {"Tricapped trigonal prismatic", "三帽三棱柱形", "三冠三角柱形", "prisme trigonal tricoiffé"}},
    {"PPR-10", {"Pentagonal prism", "五棱柱形", "五角柱形", "prisme pentagonal"}},
    {"PAPR-10", {"Pentagonal antiprism", "五方反棱柱形", "五角反柱形", "antiprisme pentagonal"}},
    {"BCSAPR-10", {"Bicapped square antiprismatic", "双帽四方反棱柱形", "二冠四角反柱形", "antiprisme carré bicoiffé"}},
    {"SPC-10", {"Sphenocorona", "斯芬诺冠形", "スフェノコロナ形", "sphénocouronne"}},
    {"TD-10", {"Tetradecahedral (2:6:2)", "十四面体形（2:6:2）", "十四面体形（2:6:2）", "tétradécaédrique (2:6:2)"}},
});

struct MachineTranslation {
    std::string_view machine;
    LocalisedString text;
};

constexpr auto kChemistryMethods = std::to_array<MachineTranslation>({
    {"COV FCHK S-metric minimal atomic-reference projection",
     {"COV FCHK S-metric minimal atomic-reference projection",
      "COV FCHK S 度量最小原子参考投影",
      "COV FCHK S 計量による最小原子参照射影",
      "Projection COV FCHK sur la référence atomique minimale selon la métrique S"}},
});

constexpr auto kChemistryNotes = std::to_array<MachineTranslation>({
    {"AO overlap metric unavailable",
     {"AO overlap metric unavailable", "AO 重叠度量不可用", "AO 重なり計量を利用できません", "Métrique de recouvrement des AO indisponible"}},
    {"COV minimal atomic reference could not be built",
     {"COV minimal atomic reference could not be built", "无法构建 COV 最小原子参考", "COV の最小原子参照を構築できませんでした", "La référence atomique minimale COV n’a pas pu être construite"}},
    {"AO-to-atom map or chemical-valence rank unavailable",
     {"AO-to-atom map or chemical-valence rank unavailable", "AO–原子映射或化学价层秩不可用", "AO–原子対応または化学原子価ランクを利用できません", "Correspondance AO–atome ou rang de valence chimique indisponible"}},
    {"No stable atom-pair interaction frame; chemistry remains UND",
     {"No stable atom-pair interaction frame; chemistry remains UND", "未获得稳定的原子对相互作用框架；化学判定保持为 UND", "安定な原子対相互作用フレームが得られないため、化学判定は UND のままです", "Aucun repère stable d’interaction par paire atomique ; l’attribution chimique reste UND"}},
    {"Outside the selected minimal chemical-valence canonical manifold",
     {"Outside the selected minimal chemical-valence canonical manifold", "位于所选最小化学价层正则轨道空间之外", "選択された最小化学原子価正準軌道空間の外側です", "Hors de l’espace canonique minimal de valence chimique sélectionné"}},
});

const char* localised(const LocalisedString& value, const Language language) noexcept {
    switch (language) {
        case Language::ChineseSimplified: return value.zh;
        case Language::Japanese: return value.ja;
        case Language::French: return value.fr;
        default: return value.en;
    }
}

template <std::size_t N>
std::string translate_machine_value(
    const std::string_view machine,
    const std::array<MachineTranslation, N>& values,
    const Language language) {
    for (const auto& value : values) {
        if (value.machine == machine) return localised(value.text, language);
    }
    return std::string(machine);
}

const char* mode_text(const MODiagramMode mode, const Language language) noexcept {
    switch (mode) {
        case MODiagramMode::DelocalisedPiFamilyOnly:
            return orbital_tr(OrbitalText::ModeDelocalisedPi, language);
        case MODiagramMode::MulticentreActiveSpaceOnly:
            return orbital_tr(OrbitalText::ModeMulticentre, language);
        default:
            return orbital_tr(OrbitalText::ModeValenceCentral, language);
    }
}

} // namespace

const char* orbital_tr(const OrbitalText key, const Language language) noexcept {
    const auto index = static_cast<std::size_t>(key);
    if (index >= kStrings.size()) return "";
    return localised(kStrings[index], language);
}

const char* localised_wavefunction_source(
    const WavefunctionSource source, const Language language) noexcept {
    switch (source) {
        case WavefunctionSource::Fchk: return "FCHK";
        case WavefunctionSource::Molden: return "Molden";
        default: return orbital_tr(OrbitalText::UnknownSource, language);
    }
}

const char* localised_data_provenance(
    const DataProvenance provenance, const Language language) noexcept {
    switch (provenance) {
        case DataProvenance::Producer:
            return orbital_tr(OrbitalText::Producer, language);
        case DataProvenance::Derived:
            return orbital_tr(OrbitalText::Derived, language);
        default:
            return orbital_tr(OrbitalText::Unavailable, language);
    }
}

const char* localised_annotation_source(
    const AnnotationSource source, const Language language) noexcept {
    switch (source) {
        case AnnotationSource::Direct:
            return orbital_tr(OrbitalText::AnnotationDirect, language);
        case AnnotationSource::ParsedLabel:
            return orbital_tr(OrbitalText::AnnotationParsedLabel, language);
        case AnnotationSource::Derived:
            return orbital_tr(OrbitalText::AnnotationDerived, language);
        case AnnotationSource::Heuristic:
            return orbital_tr(OrbitalText::AnnotationHeuristic, language);
        default:
            return orbital_tr(OrbitalText::Unavailable, language);
    }
}

const char* localised_bonding_class(
    const BondingClass value, const Language language) noexcept {
    switch (value) {
        case BondingClass::Bonding:
            return orbital_tr(OrbitalText::Bonding, language);
        case BondingClass::Nonbonding:
            return orbital_tr(OrbitalText::Nonbonding, language);
        case BondingClass::Antibonding:
            return orbital_tr(OrbitalText::Antibonding, language);
        default:
            return "N/A";
    }
}

const char* localised_orbital_bonding_role(
    const OrbitalBondingRole value, const Language language) noexcept {
    switch (value) {
        case OrbitalBondingRole::Bonding:
            return orbital_tr(OrbitalText::Bonding, language);
        case OrbitalBondingRole::Antibonding:
            return orbital_tr(OrbitalText::Antibonding, language);
        case OrbitalBondingRole::Nonbonding:
            return orbital_tr(OrbitalText::Nonbonding, language);
        case OrbitalBondingRole::NotApplicable:
            return "N/A";
        default:
            return orbital_tr(OrbitalText::Mixed, language);
    }
}

const char* localised_pi_interaction_kind(
    const PiInteractionKind kind, const Language language) noexcept {
    switch (kind) {
        case PiInteractionKind::Donor:
            return orbital_tr(OrbitalText::PiDonorSplitting, language);
        case PiInteractionKind::Acceptor:
            return orbital_tr(OrbitalText::PiAcceptorSplitting, language);
        case PiInteractionKind::WeakNearNonbonding:
            return orbital_tr(OrbitalText::PiWeakNearNonbonding, language);
        default:
            return orbital_tr(OrbitalText::PiCoupledSplitting, language);
    }
}

std::string localised_geometry_name(
    const std::string_view geometry_id,
    const std::string_view fallback_name,
    const Language language) {
    for (const auto& [id, name] : kGeometryNames) {
        if (id == geometry_id) return localised(name, language);
    }
    return std::string(fallback_name);
}

std::string localised_chemistry_method(
    const std::string_view method, const Language language) {
    return translate_machine_value(method, kChemistryMethods, language);
}

std::string localised_chemistry_note(
    const std::string_view note, const Language language) {
    return translate_machine_value(note, kChemistryNotes, language);
}

std::string localised_diagram_selection_summary(
    const MODiagramData& data, const Language language) {
    std::ostringstream out;
    if (data.selection.summary.starts_with("compact active space unavailable")) {
        out << orbital_tr(OrbitalText::CompactActiveSpaceUnavailable, language)
            << " · ";
    }
    out << orbital_tr(OrbitalText::DiagramSummary, language)
        << " · " << mode_text(data.mode, language)
        << ": " << orbital_tr(OrbitalText::VisibleLevels, language)
        << '=' << data.levels.size() << '/' << data.metadata.size()
        << "; " << orbital_tr(OrbitalText::Occupied, language)
        << '=' << data.selection.valence_occupied_count
        << "; " << orbital_tr(OrbitalText::Virtual, language)
        << '=' << data.selection.frontier_virtual_count
        << "; " << orbital_tr(OrbitalText::Hidden, language)
        << '=' << data.selection.hidden_count;

    if (!data.ligand_field_point_group.empty()) {
        out << "; " << orbital_tr(OrbitalText::LocalField, language)
            << '=' << data.ligand_field_point_group;
    }
    if (!data.ligand_field_geometry_id.empty()) {
        out << "; " << orbital_tr(OrbitalText::Geometry, language)
            << '=' << data.ligand_field_geometry_id
            << "; CN=" << data.ligand_field_coordination_number;
    }
    out << "; " << orbital_tr(OrbitalText::PiPairs, language)
        << '=' << data.pi_interactions.size()
        << "; " << orbital_tr(OrbitalText::ProtectedOverflow, language)
        << '=' << data.selection.protected_overflow_count;
    if (data.spin_counterpart_pair_count > 0u) {
        out << "; " << orbital_tr(OrbitalText::SpinPairs, language)
            << '=' << data.spin_counterpart_pair_count
            << "; " << orbital_tr(OrbitalText::UnmatchedVisible, language)
            << '=' << data.spin_counterpart_unmatched_visible;
    }
    std::size_t raw_recovered = 0u;
    for (const auto& level : data.levels) {
        if (level.raw_data_fallback) ++raw_recovered;
    }
    if (raw_recovered > 0u) {
        out << "; " << orbital_tr(OrbitalText::RawRecoveredGroups, language)
            << '=' << raw_recovered;
    }
    return out.str();
}

const char* orbital_ui_glyph_seed(const Language language) noexcept {
    static const std::array<std::string, static_cast<std::size_t>(Language::Count)> seeds = [] {
        std::array<std::string, static_cast<std::size_t>(Language::Count)> result;
        for (std::size_t language_index = 0; language_index < result.size(); ++language_index) {
            const auto current = static_cast<Language>(language_index);
            auto& seed = result[language_index];
            const auto append = [&seed](const std::string_view text) {
                if (!seed.empty()) seed.push_back(' ');
                seed.append(text);
            };
            for (const auto& value : kStrings) append(localised(value, current));
            for (const auto& [id, value] : kGeometryNames) {
                append(id);
                append(localised(value, current));
            }
            for (const auto& value : kChemistryMethods) append(localised(value.text, current));
            for (const auto& value : kChemistryNotes) append(localised(value.text, current));
        }
        return result;
    }();
    const auto index = static_cast<std::size_t>(language);
    return index < seeds.size() ? seeds[index].c_str() : "";
}

} // namespace cov::ui
