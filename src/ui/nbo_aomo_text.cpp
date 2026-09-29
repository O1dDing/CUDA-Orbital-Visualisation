#include "cov/nbo_aomo_text.hpp"
#include <cstring>
namespace cov::ui {
namespace {
const char* const rows[][4]={
    {"Filter and grouping settings","筛选与分组设置","絞り込みとグループ設定","Réglages du filtre et des groupes"},
    {"spaces","组","群","groupes"},
    {"AO/NAO decomposition unavailable; showing the canonical MO diagram.","AO/NAO 分解暂不可用，显示正则 MO 图。","AO/NAO 分解を利用できないため、正準 MO 図を表示します。","Décomposition AO/NAO indisponible ; affichage des OM canoniques."},
    {"Orbital interaction diagram","轨道相互作用图","軌道相互作用図","Diagramme d'interaction orbitale"},
    {"Illustrative layout (not energy)","侧栏示意（非能量）","側方の模式図（非エネルギー）","Schéma latéral (sans énergie)"},
    {"Reset view","重置视图","表示をリセット","Réinitialiser la vue"},
    {"Show core NAOs","显示内层 NAO","内殻 NAO を表示","Afficher les NAO de coeur"},
    {"Show Rydberg NAOs","显示里德堡 NAO","Rydberg NAO を表示","Afficher les NAO de Rydberg"},
    {"Hide H orbitals","隐藏 H 轨道","H 軌道を非表示","Masquer les orbitales H"},
    {"User atom sets (advanced)","用户原子集合（高级）","ユーザー原子集合（詳細）","Ensembles d'atomes utilisateur (avancé)"},
    {"Save atom set","保存原子集合","原子集合を保存","Enregistrer l'ensemble"},
    {"Clear draft","清空草稿","下書きを消去","Effacer le brouillon"},
    {"Orbital combinations (advanced)","轨道组合（高级）","軌道の組み合わせ（詳細）","Combinaisons orbitales (avancé)"},
    {"Add every local term","加入全部局域项","すべての項を追加","Ajouter tous les termes"},
    {"Clear terms","清空项","項を消去","Effacer les termes"},
    {"Show partial sum in 3D","三维查看部分和","部分和を3D表示","Afficher la somme partielle en 3D"},
    {"Overlay terms in 3D","三维逐项叠加","各項を3D重ね表示","Superposer les termes en 3D"},
    {"Show full canonical MO","查看完整正则 MO","正準 MO 全体を表示","Afficher l'OM canonique entière"},
    {"Export whole diagram","导出整体图","全体図を書き出す","Exporter le diagramme complet"},
    {"Light paper export","浅色论文导出","明色の論文用出力","Export clair pour article"},
    {"Overview","简要分析","概要","Vue d’ensemble"},
    {"Research analysis","研究分析","研究用解析","Analyse détaillée"},
    {"Full basis","全部轨道","全軌道","Toutes les orbitales"},
    {"View","视图","表示","Vue"},
    {"Other MO links","其他 MO 连线","他の MO の連線","Liens des autres OM"},
    {"Fit graph","适配图形","図を合わせる","Ajuster le graphe"},
    {"Find selected MO","定位所选 MO","選択 MO を表示","Centrer l’OM sélectionnée"},
    {"nonlinear spacing (not energy gaps)","非线性间距（不代表能隙）","非線形間隔（エネルギー差ではない）","espacement non linéaire (sans valeur d’écart)"},
    {"amber ticks: adaptive spacing, not energy gaps","琥珀色刻度：自适应间距，不代表能隙","琥珀色の目盛：適応的な間隔、エネルギー差ではない","graduations ambrées : espacement adaptatif, sans valeur d’écart"},
    {"linear energy scale","线性能量刻度","線形エネルギー目盛","échelle d’énergie linéaire"},
    {"Group means / shell stacks; raw energies retained.","简并组取均值，原子壳层紧凑叠放；原始能量保留。","縮退群は平均値、原子殻は近接配置。元のエネルギーは保持。","Moyennes des groupes et couches rapprochées ; énergies originales conservées."},
    {"Non-quantitative","非定量","定量外","Non quantitatif"},
    {"Set","组","群","Groupe"},
    {"Lines show orbital contributions; click a level for its own 3D orbital.","连线表示轨道组成；点击能级单独查看该轨道。","線は軌道の構成。準位をクリックすると、その軌道だけを3D表示。","Les lignes indiquent la composition ; cliquez sur un niveau pour voir son orbitale seule."},
    {"Some core or Rydberg orbitals are hidden.","部分内层或里德堡轨道已隐藏。","一部の内殻・Rydberg軌道を非表示。","Certaines orbitales de cœur ou de Rydberg sont masquées."},
    {"H orbitals are hidden from this view.","此视图已隐藏 H 轨道。","この表示では H 軌道を非表示にしています。","Les orbitales H sont masquées dans cette vue."},
    {"Protected orbital partners need more width; scroll the diagram horizontally.","简并伙伴需要更多宽度，可横向滚动查看。","縮退成分には横幅が必要です。横にスクロールできます。","Faites défiler horizontalement pour voir les partenaires dégénérés."},
    {"Scroll to explore; Ctrl+wheel zooms. Set buttons fold/expand complete canonical groups.","滚动浏览；Ctrl＋滚轮缩放；点击组标签展开或折叠。","スクロールで移動、Ctrl＋ホイールで拡縮。群ラベルで展開・折り畳み。","Défilez pour explorer ; Ctrl+molette pour zoomer ; cliquez sur un groupe pour le déplier."},
    {"Scroll to explore; Ctrl+wheel zooms. Click a real orbital to inspect it.","滚动浏览；Ctrl＋滚轮缩放；点击真实轨道查看详情。","スクロールで移動、Ctrl＋ホイールで拡縮。実際の軌道をクリックして確認。","Faites défiler ; Ctrl+molette pour zoomer ; cliquez sur une orbitale réelle pour l’examiner."},
    {"Choose the actual canonical MO","选择要查看的 MO","表示する MO を選択","Choisir l’OM à afficher"},
    {"Choose a numerical component","选择要查看的轨道分量","表示する軌道成分を選択","Choisir la composante orbitale"},
    {"An atom set organizes folded channels. It is not a fixed SALC or a single orbital. Its 3D partial sum is explicitly MO-dependent.","原子集合用于整理轨道；其三维部分和随所选 MO 变化，不是独立 SALC。","原子集合は軌道の整理用です。3D部分和は選択 MO に依存し、独立した SALC ではありません。","Un ensemble organise les orbitales. Sa somme partielle dépend de l’OM choisie ; ce n’est pas une SALC fixe."},
    {"Showing partial orbital composition.","当前显示部分轨道组成。","軌道成分の一部を表示しています。","Une partie de la composition orbitale est affichée."},
    {"Edit","编辑","編集","Modifier"},
    {"Delete","删除","削除","Supprimer"},
    {"MO partial sum in 3D","查看 MO 部分和","MO 部分和を3D表示","Somme partielle de l’OM en 3D"},
    {"Selected components: %zu / %zu","已选分量：%zu / %zu","選択した成分：%zu / %zu","Composantes choisies : %zu / %zu"},
    {"Orbital contribution; not a bond-order label","轨道组成分量，不代表键级","軌道の構成成分。結合次数ではありません","Contribution orbitale, sans attribution d’ordre de liaison"},
    {"Click to show this weighted component","点击单独查看此加权分量","クリックして、この重み付き成分だけを表示","Cliquer pour afficher cette composante pondérée seule"},
};
std::size_t column(Language language){const auto i=static_cast<std::size_t>(language);return i<4?i:0;}
}
const char* aomo_text(Language language,const char* english) {
    for(const auto& row:rows)if(std::strcmp(row[0],english)==0)return row[column(language)];
    return english;
}
std::string nbo_aomo_text_glyph_seed(Language language) {
    std::string seed="AO NAO MO SALC Gaussian αβ ↑↓ Σσπδφγ − 总全体轨道组成軌道成分Composantes orbitales de";
    seed+=" 轨道范围 价层轨道（σ 与 π） 仅离域 π 子集 仅多中心活性空间";
    seed+=" 軌道の範囲 価電子軌道（σ と π） 非局在 π 部分集合のみ 多中心活性空間のみ";
    seed+=" Ensemble orbital Orbitales de valence (σ et π) Sous-ensemble π délocalisé Espace actif multicentrique";
    for(const auto& row:rows){seed+=row[column(language)];seed+=' ';}
    return seed;
}
}
