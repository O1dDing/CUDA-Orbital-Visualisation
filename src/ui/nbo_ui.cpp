#include "cov/nbo_ui.hpp"
#include "cov/validation.hpp"
#include <imgui.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <set>
#include <map>
#include <stdexcept>
#include <vector>

namespace cov::ui {
namespace {
enum Label { Title, Path, Archive, Aonbo, Nbomo, Naomo, Aonao, Naonbo, FocusNote, Segment, Attach, Canonical,
             NboSet, Association, Renderable, ReportsOnly, Populations, Naos,
             Orbitals, Wiberg, E2, Matrices, EnergyNote, ExportPath, Export,
             Source, Threshold, Units, NoData, Count };
constexpr const char* labels[4][Count] = {
    {"NBO analysis", "NBO output path", "ARCHIVE .47 (optional)", "AONBO .37 (optional)", "NBOMO .49 (optional)", "NAOMO .51: NAO / canonical MO (optional)", "AONAO .52: AO / NAO (optional)", "NAONBO .53: NAO / NBO (optional)", "Focused NAO composition appears in the central MO diagram above.", "Analysis segment (-1: require unique)", "Attach NBO", "Canonical MOs", "NBO orbitals", "Source association", "NBO coefficients verified for rendering", "Report only: coefficients or source association unavailable", "NPA populations", "NAOs", "NBO orbital occupations and components", "Wiberg indices", "Second-order interactions", "Matrix evidence", "Diagonal Fock values are not canonical MO energies.", "Export base path", "Export NBO bundle", "Source", "Print threshold", "Units", "No data reported"},
    {"NBO 分析", "NBO 输出路径", "ARCHIVE .47（可选）", "AONBO .37（可选）", "NBOMO .49（可选）", "NAOMO .51：NAO / 正则 MO（可选）", "AONAO .52：AO / NAO（可选）", "NAONBO .53：NAO / NBO（可选）", "聚焦的 NAO 成分显示在上方中央 MO 图中。", "分析段（-1：要求唯一）", "关联 NBO", "正则 MO", "NBO 轨道", "同源关联", "NBO 系数已验证，可渲染", "仅报告：缺少系数或同源关联", "NPA 布居", "NAO", "NBO 轨道占据与成分", "Wiberg 指数", "二阶相互作用", "矩阵证据", "Fock 对角值不是正则 MO 能量。", "导出路径前缀", "导出 NBO 数据", "来源", "打印阈值", "单位", "未报告数据"},
    {"NBO 解析", "NBO 出力パス", "ARCHIVE .47（任意）", "AONBO .37（任意）", "NBOMO .49（任意）", "NAOMO .51：NAO / 正準 MO（任意）", "AONAO .52：AO / NAO（任意）", "NAONBO .53：NAO / NBO（任意）", "NAO 成分は上の中央 MO 図に表示します。", "解析区間（-1：一意のみ）", "NBO を関連付け", "正準 MO", "NBO 軌道", "同一源の照合", "描画用 NBO 係数を検証済み", "報告のみ：係数または照合が不足", "NPA 原子分布", "NAO", "NBO 軌道の占有数と成分", "Wiberg 指数", "二次相互作用", "行列の証拠", "Fock 対角値は正準 MO エネルギーではありません。", "出力先の基底名", "NBO 一式を書き出す", "出典", "印字しきい値", "単位", "報告データなし"},
    {"Analyse NBO", "Chemin de sortie NBO", "ARCHIVE .47 (facultatif)", "AONBO .37 (facultatif)", "NBOMO .49 (facultatif)", "NAOMO .51 : NAO / OM canoniques (facultatif)", "AONAO .52 : AO / NAO (facultatif)", "NAONBO .53 : NAO / NBO (facultatif)", "La composition NAO ciblée figure dans le diagramme central ci-dessus.", "Segment d’analyse (-1 : unique)", "Associer NBO", "OM canoniques", "Orbitales NBO", "Association des sources", "Coefficients NBO vérifiés pour l’affichage", "Rapport seul : coefficients ou association absents", "Populations NPA", "NAO", "Occupation et composants NBO", "Indices de Wiberg", "Interactions du second ordre", "Preuves matricielles", "Les valeurs diagonales de Fock ne sont pas des énergies OM canoniques.", "Chemin de base de l’export", "Exporter NBO", "Source", "Seuil d’impression", "Unités", "Aucune donnée rapportée"}
};
const char* trn(Language language, Label key) {
    const int row=std::clamp(static_cast<int>(language),0,3);
    return labels[row][key];
}
const char* nbo_local(Language language,const char* en,const char* zh,
                      const char* ja,const char* fr) {
    switch(language){case Language::ChineseSimplified:return zh;
        case Language::Japanese:return ja;case Language::French:return fr;
        default:return en;}
}
struct ReportText {
    const char *atom,*charge,*core,*valence,*rydberg,*total,*occupation,*fock,
               *component,*coefficient,*source_evidence,*contributions,
               *selected_canonical,*threshold,*shown,*remainder,*no_matrix,
               *mapping_error,*energy_gap,*matrix_weight,*effective_core,*explicit_pop,
               *beta_archive_notice,*canonical_evidence,*direct_fchk,*density_verified,
               *spin_density;
};
ReportText report_text(Language language) {
    switch(language){
        case Language::ChineseSimplified:return {"原子","电荷","内层","价层","里德堡","打印总布居","占据","Fock 对角 (Ha)","成分","系数","来源证据","NBOMO 贡献","选中正则 MO","显示阈值 |T|²","已显示权重","剩余权重","无完整 NBOMO 矩阵","行/列或自旋不匹配","能隙 (Ha)","原始权重 |T|²","ECP 核心电子","显式电子布居","FCHK 未提供 beta MO；beta 轨道来源为已关联的 archive 与 NBO 文件，不能当作 FCHK beta MO。","正则系数证据","直接 FCHK 系数","密度已核验","自旋密度"};
        case Language::Japanese:return {"原子","電荷","内殻","価電子","リュードベリ","印字総分布","占有数","Fock 対角 (Ha)","成分","係数","出典の証拠","NBOMO 寄与","選択した正準 MO","表示しきい値 |T|²","表示した重み","残りの重み","完全な NBOMO 行列なし","行・列またはスピンが不一致","エネルギー差 (Ha)","生の重み |T|²","ECP 内殻電子","明示的電子分布","FCHK に beta MO がありません。beta 軌道は照合済み archive と NBO ファイルに由来し、FCHK の beta MO とみなしません。","正準係数の証拠","FCHK 係数を直接確認","密度を検証済み","スピン密度"};
        case Language::French:return {"Atome","Charge","Cœur","Valence","Rydberg","Population totale imprimée","Occupation","Diagonale de Fock (Ha)","Composant","Coefficient","Preuves de provenance","Contributions NBOMO","OM canonique choisi","Seuil affiché |T|²","Poids affiché","Poids restant","Aucune matrice NBOMO complète","Lignes, colonnes ou spin incompatibles","Écart d’énergie (Ha)","Poids brut |T|²","Électrons de cœur ECP","Population explicite","Le FCHK ne contient pas d’OM bêta ; les orbitales bêta proviennent de l’archive associée et des fichiers NBO, sans identité OM bêta FCHK.","Preuve des coefficients canoniques","Coefficients FCHK directs","Densité vérifiée","Densité de spin"};
        default:return {"Atom","Charge","Core","Valence","Rydberg","Printed total population","Occupation","Fock diagonal (Ha)","Component","Coefficient","Source evidence","NBOMO contributions","Selected canonical MO","Display threshold |T|²","Shown weight","Remaining weight","No complete NBOMO matrix","Row/column or spin mismatch","Energy gap (Ha)","Raw weight |T|²","ECP core electrons","Explicit population","FCHK has no beta MO; beta orbitals come from the associated archive and NBO files, not FCHK beta MO identity.","Canonical coefficient evidence","Direct FCHK coefficients","Density verified","Spin density"};
    }
}
std::string fmt(double value) {
    std::ostringstream out; out<<std::setprecision(9)<<value;return out.str();
}
std::string source_label(const NboSource& s) {
    std::ostringstream out;
    out<<s.path;
    if(s.line_begin)out<<':'<<s.line_begin;
    if(!s.block.empty())out<<" ["<<s.block<<']';
    out<<" segment="<<s.analysis_segment;
    return out.str();
}
void input_path(const char* label, const char* id, std::array<char,2048>& value) {
    ImGui::TextWrapped("%s",label);
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputText(id,value.data(),value.size());
    validation::item(id+2);
}
void note(const std::string& s) { ImGui::TextWrapped("%s",s.c_str()); }
std::string csv(const std::string& value) {
    std::string result="\"";
    for(char c:value){if(c=='\"')result+="\"\"";else result+=c;}
    result+='\"';return result;
}
std::string json(const std::string& value) {
    std::string out="\"";
    for(unsigned char c:value){
        switch(c){case '\"':out+="\\\"";break;case '\\':out+="\\\\";break;
            case '\n':out+="\\n";break;case '\r':out+="\\r";break;
            case '\t':out+="\\t";break;
            default:if(c<0x20){const char hex[]="0123456789abcdef";out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}else out+=static_cast<char>(c);}
    }
    out+='\"';return out;
}
std::ofstream output_csv(const std::filesystem::path& path) {
    std::ofstream out(path,std::ios::binary);
    if(!out)throw std::runtime_error("Cannot write "+path.string());
    out<<std::setprecision(17);
    return out;
}
std::string opt(const std::optional<double>& value) { return value?fmt(*value):""; }
std::string shown_opt(const std::optional<double>& value) { return value?fmt(*value):"null"; }
std::string escape_xml(const std::string& value) {
    std::string out;
    for(char c:value){switch(c){case '&':out+="&amp;";break;case '<':out+="&lt;";break;case '>':out+="&gt;";break;case '"':out+="&quot;";break;default:out+=c;}}
    return out;
}
std::array<std::uint8_t,7> glyph(char c) {
    if(c>='a'&&c<='z')c=static_cast<char>(c-'a'+'A');
    switch(c){
        case '0':return {14,17,19,21,25,17,14};case '1':return {4,12,4,4,4,4,14};
        case '2':return {14,17,1,2,4,8,31};case '3':return {30,1,1,14,1,1,30};
        case '4':return {2,6,10,18,31,2,2};case '5':return {31,16,16,30,1,1,30};
        case '6':return {14,16,16,30,17,17,14};case '7':return {31,1,2,4,8,8,8};
        case '8':return {14,17,17,14,17,17,14};case '9':return {14,17,17,15,1,1,14};
        case 'A':return {14,17,17,31,17,17,17};case 'B':return {30,17,17,30,17,17,30};
        case 'C':return {14,17,16,16,16,17,14};case 'D':return {30,17,17,17,17,17,30};
        case 'E':return {31,16,16,30,16,16,31};case 'F':return {31,16,16,30,16,16,16};
        case 'G':return {14,17,16,23,17,17,15};case 'H':return {17,17,17,31,17,17,17};
        case 'I':return {14,4,4,4,4,4,14};case 'J':return {7,2,2,2,2,18,12};
        case 'K':return {17,18,20,24,20,18,17};case 'L':return {16,16,16,16,16,16,31};
        case 'M':return {17,27,21,21,17,17,17};case 'N':return {17,25,21,19,17,17,17};
        case 'O':return {14,17,17,17,17,17,14};case 'P':return {30,17,17,30,16,16,16};
        case 'Q':return {14,17,17,17,21,18,13};case 'R':return {30,17,17,30,20,18,17};
        case 'S':return {15,16,16,14,1,1,30};case 'T':return {31,4,4,4,4,4,4};
        case 'U':return {17,17,17,17,17,17,14};case 'V':return {17,17,17,17,17,10,4};
        case 'W':return {17,17,17,21,21,21,10};case 'X':return {17,17,10,4,10,17,17};
        case 'Y':return {17,17,10,4,4,4,4};case 'Z':return {31,1,2,4,8,16,31};
        case '-':return {0,0,0,31,0,0,0};case '.':return {0,0,0,0,0,12,12};
        case '/':return {1,2,2,4,8,8,16};case '+':return {0,4,4,31,4,4,0};
        case ':':return {0,12,12,0,12,12,0};case '_':return {0,0,0,0,0,0,31};
        case '?':return {14,17,1,2,4,0,4};case '=':return {0,31,0,31,0,0,0};
        case '(':return {2,4,8,8,8,4,2};case ')':return {8,4,2,2,2,4,8};
        default:return {0,0,0,0,0,0,0};
    }
}
void draw_text(std::vector<unsigned char>& rgba,int width,int height,int x,int y,
               const std::string& value,int scale=2) {
    for(unsigned char c:value){
        if(c==' '){x+=4*scale;continue;}
        const auto rows=glyph(static_cast<char>(c));
        for(int row=0;row<7;++row)for(int col=0;col<5;++col)
            if(rows[row]&(1u<<(4-col)))for(int dy=0;dy<scale;++dy)for(int dx=0;dx<scale;++dx){
                const int px=x+col*scale+dx,py=y+row*scale+dy;
                if(px<0||py<0||px>=width||py>=height)continue;
                const auto p=(static_cast<std::size_t>(py)*width+px)*4;
                rgba[p]=28;rgba[p+1]=39;rgba[p+2]=55;
            }
        x+=6*scale;
    }
}
std::uint32_t crc32(const unsigned char* data,std::size_t size) {
    std::uint32_t crc=0xffffffffu;
    for(std::size_t i=0;i<size;++i){crc^=data[i];for(int j=0;j<8;++j)crc=(crc>>1)^((crc&1)?0xedb88320u:0u);}
    return ~crc;
}
void put32(std::vector<unsigned char>& out,std::uint32_t v) {
    out.push_back(static_cast<unsigned char>(v>>24));out.push_back(static_cast<unsigned char>(v>>16));
    out.push_back(static_cast<unsigned char>(v>>8));out.push_back(static_cast<unsigned char>(v));
}
void chunk(std::vector<unsigned char>& out,const char type[4],const std::vector<unsigned char>& data) {
    put32(out,static_cast<std::uint32_t>(data.size()));
    const auto begin=out.size();out.insert(out.end(),type,type+4);out.insert(out.end(),data.begin(),data.end());
    put32(out,crc32(out.data()+begin,out.size()-begin));
}
void write_png(const std::filesystem::path& path,const std::vector<unsigned char>& rgba,int width,int height) {
    std::vector<unsigned char> raw;raw.reserve(static_cast<std::size_t>(height)*(1+width*4));
    for(int y=0;y<height;++y){raw.push_back(0);const auto start=static_cast<std::size_t>(y)*width*4;raw.insert(raw.end(),rgba.begin()+start,rgba.begin()+start+width*4);}
    std::vector<unsigned char> z={0x78,0x01};
    std::size_t pos=0;
    while(pos<raw.size()){
        const std::size_t n=std::min<std::size_t>(65535,raw.size()-pos);const bool last=pos+n==raw.size();
        z.push_back(last?1:0);z.push_back(static_cast<unsigned char>(n));z.push_back(static_cast<unsigned char>(n>>8));
        z.push_back(static_cast<unsigned char>(~n));z.push_back(static_cast<unsigned char>((~n)>>8));
        z.insert(z.end(),raw.begin()+pos,raw.begin()+pos+n);pos+=n;
    }
    std::uint32_t a=1,b=0;for(auto v:raw){a=(a+v)%65521u;b=(b+a)%65521u;}put32(z,(b<<16)|a);
    std::vector<unsigned char> png={137,80,78,71,13,10,26,10};
    std::vector<unsigned char> ihdr;put32(ihdr,width);put32(ihdr,height);ihdr.insert(ihdr.end(),{8,6,0,0,0});
    chunk(png,"IHDR",ihdr);chunk(png,"IDAT",z);chunk(png,"IEND",{});
    std::ofstream out(path,std::ios::binary);if(!out)throw std::runtime_error("Cannot write PNG");
    out.write(reinterpret_cast<const char*>(png.data()),static_cast<std::streamsize>(png.size()));
    if(!out)throw std::runtime_error("PNG write failed");
}

struct FocusShell {
    std::size_t atom=0;
    std::string symbol, key, label, type;
    int n=0, l=0;
    double weight=0;
    std::optional<double> electrons;
};
struct FocusView {
    std::string snapshot_id, status, reason;
    std::vector<std::size_t> eligible;
    std::optional<std::size_t> index;
    const NboMoDecomposition* decomposition=nullptr;
    std::vector<FocusShell> shells;
};
bool type_starts(const std::string& type,const char* prefix) {
    for(std::size_t i=0;prefix[i];++i)
        if(i>=type.size() || std::tolower(static_cast<unsigned char>(type[i]))!=prefix[i])return false;
    return true;
}
bool valence_nao(const NboNaoContribution& row) {
    return type_starts(row.type,"val");
}
const char* shell_letter(int l) {
    switch(l){case 0:return "s";case 1:return "p";case 2:return "d";
              case 3:return "f";case 4:return "g";default:return "?";}
}
std::string focus_shell_key(std::size_t atom,int n,int l,const std::string& type) {
    return std::to_string(atom)+":"+std::to_string(n)+":"+std::to_string(l)+":"+type;
}
bool shell_enabled(const FocusShell& shell,const NboFocusUIState& state) {
    const bool by_class=type_starts(shell.type,"val") ||
        (type_starts(shell.type,"cor") && state.show_core) ||
        (type_starts(shell.type,"ryd") && state.show_rydberg);
    return (by_class || state.explicitly_included_shells.contains(shell.key)) &&
           !state.hidden_shells.contains(shell.key);
}
FocusView make_focus_view(const NboDataset& dataset,
                          const MODiagramViewSnapshot* diagram,
                          const NboFocusUIState* state) {
    FocusView v;
    if(!diagram || !diagram->data.view) {
        v.status="no_central_snapshot";
        v.reason="The canonical MO diagram has no current view snapshot.";
        return v;
    }
    v.snapshot_id=diagram->data.view->id;
    std::set<std::size_t> indices;
    for(const auto& level:diagram->data.levels) {
        if(level.member_indices.empty())indices.insert(level.metadata.orbital_index);
        else indices.insert(level.member_indices.begin(),level.member_indices.end());
        indices.insert(level.member_spin_counterparts.begin(),
                       level.member_spin_counterparts.end());
    }
    v.eligible.assign(indices.begin(),indices.end());
    if(v.eligible.empty()) {
        v.status="no_central_mo";
        v.reason="The current central diagram has no MO rows.";
        return v;
    }
    if(!state || !state->canonical_index || !indices.contains(*state->canonical_index)) {
        v.status="choose_central_mo";
        v.reason="Choose an identified MO from the central diagram.";
        return v;
    }
    v.index=*state->canonical_index;
    v.decomposition=nbo_mo_decomposition(dataset,*v.index);
    if(!v.decomposition || !v.decomposition->available) {
        v.status=v.decomposition?v.decomposition->status:"decomposition_unavailable";
        v.reason=v.decomposition?v.decomposition->detail:
            "No validated NAO decomposition exists for this canonical MO.";
        return v;
    }
    std::map<std::string,FocusShell> grouped;
    for(const auto& row:v.decomposition->rows) {
        // Unknown n/l cannot form a shell node. Its full numeric row remains below.
        if(!row.principal_n || !row.angular_l)continue;
        const auto key=focus_shell_key(row.atom,*row.principal_n,*row.angular_l,row.type);
        auto& shell=grouped[key];
        shell.atom=row.atom;shell.symbol=row.symbol;shell.key=key;
        shell.type=row.type;
        shell.n=*row.principal_n;shell.l=*row.angular_l;
        shell.label=std::to_string(shell.n)+shell_letter(shell.l);
        shell.weight+=row.weight;
        if(row.electron_contribution)
            shell.electrons=shell.electrons.value_or(0)+*row.electron_contribution;
    }
    for(auto& [key,shell]:grouped)v.shells.push_back(std::move(shell));
    std::stable_sort(v.shells.begin(),v.shells.end(),[](const auto& a,const auto& b){
        if(a.atom!=b.atom)return a.atom<b.atom;
        if(a.n!=b.n)return a.n<b.n;
        return a.l<b.l;
    });
    if(v.shells.empty()) {
        v.status="no_known_nl_shells";
        v.reason="Validated decomposition has no NAO row with a known n/l shell.";
    } else {v.status="available";v.reason="Choose atoms and shells from the validated NAO decomposition.";}
    return v;
}
std::string focus_mo_label(const MODiagramViewSnapshot& diagram,std::size_t index) {
    for(const auto& level:diagram.data.levels)if(mo_diagram_level_covers_orbital(level,index)) {
        std::string label="MO "+std::to_string(index+1)+" "+level.annotation.family;
        if(level.annotation.bonding_class==BondingClass::Nonbonding)label+=" / nonbonding";
        if(level.annotation.multicentre.available)label+=" / multicentre";
        return label;
    }
    return "MO "+std::to_string(index+1);
}
std::vector<const FocusShell*> visible_focus_shells(const FocusView& v,
                                                     const NboFocusUIState& state) {
    std::vector<const FocusShell*> result;
    for(const auto& shell:v.shells)
        if(state.visible_atoms.contains(shell.atom) && shell_enabled(shell,state))
            result.push_back(&shell);
    return result;
}
std::set<std::size_t> grouped_atoms(const NboFocusUIState& state) {
    std::set<std::size_t> result;
    for(const auto& group:state.ligand_groups)
        result.insert(group.atoms.begin(),group.atoms.end());
    return result;
}
struct FocusAngularWeights {double full=0, selected=0;};
std::map<int,FocusAngularWeights> angular_weights(const FocusView& view,
        const NboFocusUIState& state,const NboFocusAtomGroup& group) {
    std::map<int,FocusAngularWeights> result;
    if(view.decomposition)for(const auto& row:view.decomposition->rows)
        if(group.atoms.contains(row.atom) && row.angular_l)
            result[*row.angular_l].full+=row.weight;
    for(const auto* shell:visible_focus_shells(view,state))
        if(group.atoms.contains(shell->atom))result[shell->l].selected+=shell->weight;
    return result;
}
} // namespace

std::string nbo_glyph_seed(Language language) {
    std::string out;
    for(const char* value:labels[std::clamp(static_cast<int>(language),0,3)]) {
        out+=value;out+=' ';
    }
    const auto r=report_text(language);
    for(const char* value:{r.atom,r.charge,r.core,r.valence,r.rydberg,r.total,
                           r.occupation,r.fock,r.component,r.coefficient,
                           r.source_evidence,r.contributions,r.selected_canonical,
                           r.threshold,r.shown,r.remainder,r.no_matrix,
                           r.mapping_error,r.energy_gap,r.matrix_weight,
                           r.effective_core,r.explicit_pop,r.beta_archive_notice,
                           r.canonical_evidence,r.direct_fchk,r.density_verified,
                           r.spin_density}){
        out+=value;out+=' ';
    }
    // Dynamic labels used by the central focus view must be present when the
    // language atlas is built, including the saved ligand-group controls.
    out += "中央 MO 的 NAO 成分 选择关注原子 图中仅显示这些原子 自定义配体原子组 明确勾选 "
           "先勾选原子 再保存为 已分组原子不能重复归属 保存新配体组 清空待保存选择 清空全部配体组 删除 "
           "恢复全部壳层 纳入 按原子 按壳层 配体按 分组 当前检查的 不在中央图集合中 中央图不会为检查而加行 "
           "权重为正交 基中的 系数 不是原始 AO 系数平方 成分图所选 与当前 检查不同 在 3D 查看此 MO "
           "选择原子后绘图 完整数值仍可展开查看 取消壳层勾选可隐藏图中的成分 不会删除数值 显示全部 数值与来源 "
           "中央 MO の NAO 成分 注目原子を選択 配位子原子グループを指定 原子を選び として保存します "
           "原子の重複所属はできません 新しい配位子群を保存 未保存の選択を消去 全群を消去 削除 "
           "全殻を復元 を含む 原子別 殻別 配位子を 別 検査中の OM は中央図の集合外です 中央図には行を追加しません "
           "重みは直交 基底での 係数 であり 元の AO 係数の二乗ではありません この OM を 3D で表示 "
           "原子を選ぶと描画します 全数値は詳細で確認できます 殻の選択を外すと図で非表示になり 数値は残ります "
           "Composants NAO des OM centraux Choisir les atomes à afficher Groupe d'atomes ligand explicite "
           "Enregistrer un groupe ligand Effacer la sélection Effacer tous les groupes Supprimer Rétablir toutes les couches "
           "Inclure Core Inclure Rydberg Par atome Par couche Ligands par s/p/d Voir cette OM en 3D";
    return out;
}

void draw_nbo_focus_view(NboFocusUIState& focus,const NboDataset& dataset,
                         const Wavefunction& canonical,
                         const MODiagramViewSnapshot& diagram,
                         Language language,float scale) {
    const bool zh=language==Language::ChineseSimplified;
    const bool ja=language==Language::Japanese;
    const bool fr=language==Language::French;
    const char* title=zh?"中央 MO 的 NAO 成分":ja?"中央 MO の NAO 成分":fr?"Composants NAO des OM centraux":"NAO components of central MOs";
    const char* atom_title=zh?"选择关注原子（图中仅显示这些原子）":ja?"注目原子を選択":fr?"Choisir les atomes à afficher":"Choose atoms to show in the graph";
    const char* ligand_title=zh?"自定义配体原子组（明确勾选）":ja?"配位子原子グループを指定":fr?"Groupe d'atomes ligand explicite":"Define ligand atom group explicitly";
    ImGui::SeparatorText(title);
    validation::anchor("nbo.focus");
    const auto current=diagram.data.view?diagram.data.view->inspected_orbital_index:std::nullopt;
    auto view=make_focus_view(dataset,&diagram,&focus);
    const auto eligible=[&](std::size_t index){return std::find(view.eligible.begin(),view.eligible.end(),index)!=view.eligible.end();};
    if(focus.diagram_id!=view.snapshot_id) {
        if(current!=focus.last_inspected && current && eligible(*current))focus.canonical_index=*current;
        else if(!focus.canonical_index || !eligible(*focus.canonical_index))
            focus.canonical_index=view.eligible.empty()?std::nullopt:
                std::optional<std::size_t>(view.eligible.front());
        focus.diagram_id=view.snapshot_id;
        focus.last_inspected=current;
        view=make_focus_view(dataset,&diagram,&focus);
    }
    ImGui::TextDisabled("Central snapshot: %s",view.snapshot_id.c_str());
    validation::field("nbo.focus.snapshot",view.snapshot_id);
    if(current && !eligible(*current)) {
        ImGui::TextWrapped("%s",zh?"当前检查的 MO 不在中央图集合中；中央图不会为检查而加行。":
            ja?"検査中の OM は中央図の集合外です。中央図には行を追加しません。":
            fr?"L'OM inspectée est hors de l'ensemble ciblé du diagramme central ; aucune ligne n'est ajoutée.":
            "The inspected MO is outside the central diagram set; inspection does not add a central row.");
    }
    if(!view.eligible.empty()) {
        const std::string preview=focus.canonical_index?focus_mo_label(diagram,*focus.canonical_index):"Choose MO";
        if(ImGui::BeginCombo("##nbo.focus.mo",preview.c_str())) {
            for(const auto index:view.eligible) {
                const auto label=focus_mo_label(diagram,index);
                if(ImGui::Selectable(label.c_str(),focus.canonical_index==index))focus.canonical_index=index;
            }
            ImGui::EndCombo();
        }
        validation::item("nbo.focus.mo");
        view=make_focus_view(dataset,&diagram,&focus);
    }
    if(view.status!="available") {
        ImGui::TextWrapped("%s: %s",view.status.c_str(),view.reason.c_str());
        validation::field("nbo.focus.status",view.status+": "+view.reason);
        return;
    }
    if(!view.index || *view.index>=canonical.orbitals.size()) {
        ImGui::TextWrapped("Canonical MO identity is unavailable in the current wavefunction.");
        validation::field("nbo.focus.status","canonical_mo_identity_unavailable");
        return;
    }
    if(current!=view.index) {
        ImGui::TextWrapped("%s",zh?"成分图所选 MO 与当前 3D 检查 MO 不同。":
            ja?"成分図の OM と 3D 検査中の OM は異なります。":
            fr?"L'OM ciblée ici diffère de l'OM inspectée en 3D.":
            "The focused MO differs from the MO currently inspected in 3D.");
        if(ImGui::Button(zh?"在 3D 查看此 MO":ja?"この OM を 3D で表示":
                         fr?"Voir cette OM en 3D":"Inspect this MO in 3D"))
            focus.pending_canonical_selection=view.index;
        validation::item("nbo.focus.inspect_3d");
    }
    ImGui::TextWrapped("%s",zh?"权重为正交 NAO 基中的 |系数|²；不是原始 AO 系数平方。":
        ja?"重みは直交 NAO 基底での |係数|² であり、元の AO 係数の二乗ではありません。":
        fr?"Les poids sont |coefficient|² dans la base NAO orthogonale, pas les carrés des coefficients AO bruts.":
        "Weights are |coefficient|² in the orthogonal NAO basis, not squared raw AO coefficients.");
    ImGui::TextWrapped("%s",atom_title);
    std::map<std::size_t,std::string> atoms;
    for(const auto& shell:view.shells)atoms.emplace(shell.atom,shell.symbol);
    if(ImGui::Button(zh?"显示全部当前原子":ja?"全原子を表示":fr?"Afficher tous les atomes":"Show all current atoms"))
        for(const auto& [atom,_]:atoms)focus.visible_atoms.insert(atom);
    validation::item("nbo.focus.show_all_atoms");
    ImGui::SameLine();
    if(ImGui::Button(zh?"隐藏全部":ja?"全て隠す":fr?"Tout masquer":"Hide all"))focus.visible_atoms.clear();
    validation::item("nbo.focus.hide_all_atoms");
    std::size_t atom_slot=0;
    for(const auto& [atom,symbol]:atoms) {
        const std::string label=symbol+" "+std::to_string(atom)+"##focus.atom."+std::to_string(atom);
        bool checked=focus.visible_atoms.contains(atom);
        if(ImGui::Checkbox(label.c_str(),&checked)) {
            if(checked)focus.visible_atoms.insert(atom);else focus.visible_atoms.erase(atom);
        }
        validation::item("nbo.focus.atom."+std::to_string(atom));
        if(++atom_slot%4)ImGui::SameLine();
    }
    ImGui::NewLine();
    ImGui::TextWrapped("%s",ligand_title);
    ImGui::TextDisabled("%s",zh?"先勾选原子，再保存为 L1、L2…；已分组原子不能重复归属。":
        ja?"原子を選び L1、L2… として保存します。原子の重複所属はできません。":
        fr?"Choisissez les atomes, puis enregistrez L1, L2… ; chaque atome appartient à un seul groupe.":
        "Choose atoms, then save L1, L2…; an atom belongs to at most one group.");
    const auto assigned=grouped_atoms(focus);
    atom_slot=0;
    for(const auto& [atom,symbol]:atoms) {
        const std::string label=symbol+" "+std::to_string(atom)+"##focus.ligand."+std::to_string(atom);
        bool checked=focus.ligand_atoms.contains(atom);
        ImGui::BeginDisabled(assigned.contains(atom));
        if(ImGui::Checkbox(label.c_str(),&checked)) {
            if(checked){focus.ligand_atoms.insert(atom);focus.visible_atoms.insert(atom);}
            else focus.ligand_atoms.erase(atom);
        }
        ImGui::EndDisabled();
        validation::item("nbo.focus.ligand."+std::to_string(atom));
        if(++atom_slot%4)ImGui::SameLine();
    }
    ImGui::NewLine();
    ImGui::BeginDisabled(focus.ligand_atoms.empty());
    if(ImGui::Button(zh?"保存新配体组":ja?"新しい配位子群を保存":fr?"Enregistrer un groupe ligand":"Save ligand group")) {
        NboFocusAtomGroup group;group.id=focus.next_group_id++;group.atoms=focus.ligand_atoms;
        focus.ligand_groups.push_back(std::move(group));focus.ligand_atoms.clear();
    }
    ImGui::EndDisabled();validation::item("nbo.focus.group.add");ImGui::SameLine();
    if(ImGui::Button(zh?"清空待保存选择":ja?"未保存の選択を消去":fr?"Effacer la sélection":"Clear draft"))focus.ligand_atoms.clear();
    validation::item("nbo.focus.group.clear_draft");ImGui::SameLine();
    if(ImGui::Button(zh?"清空全部配体组":ja?"全群を消去":fr?"Effacer tous les groupes":"Clear all groups"))
        focus.ligand_groups.clear();
    validation::item("nbo.focus.group.clear");
    for(std::size_t i=0;i<focus.ligand_groups.size();) {
        const auto& group=focus.ligand_groups[i];
        const std::string name="L"+std::to_string(group.id);
        std::string members;
        for(auto atom:group.atoms){if(!members.empty())members+=", ";
            members+=atoms.contains(atom)?atoms.at(atom)+std::to_string(atom):std::to_string(atom);}
        ImGui::Text("%s: %s",name.c_str(),members.c_str());ImGui::SameLine();
        const std::string button=(zh?"删除##":ja?"削除##":fr?"Supprimer##":"Delete##")+name;
        if(ImGui::Button(button.c_str())) {
            const auto removed_id=group.id;
            focus.ligand_groups.erase(focus.ligand_groups.begin()+static_cast<std::ptrdiff_t>(i));
            validation::item("nbo.focus.group.delete."+std::to_string(removed_id));
            continue;
        }
        validation::item("nbo.focus.group.delete."+std::to_string(group.id));
        const auto weights=angular_weights(view,focus,group);
        for(const auto& [l,w]:weights)
            ImGui::Text("  %s: %s full=%.9g; selected-shell subtotal=%.9g",
                        name.c_str(),shell_letter(l),w.full,w.selected);
        ++i;
    }
    if(ImGui::Button(zh?"恢复全部壳层":ja?"全殻を復元":fr?"Rétablir toutes les couches":"Restore all shells"))
        focus.hidden_shells.clear();
    validation::item("nbo.focus.restore_shells");
    bool core=focus.show_core;
    if(ImGui::Checkbox(zh?"纳入 Core":ja?"Core を含む":fr?"Inclure Core":"Include Core",&core))focus.show_core=core;
    validation::item("nbo.focus.include_core");ImGui::SameLine();
    bool rydberg=focus.show_rydberg;
    if(ImGui::Checkbox(zh?"纳入 Rydberg":ja?"Rydberg を含む":fr?"Inclure Rydberg":"Include Rydberg",&rydberg))focus.show_rydberg=rydberg;
    validation::item("nbo.focus.include_rydberg");
    ImGui::SameLine();
    if(ImGui::RadioButton(zh?"按原子":ja?"原子別":fr?"Par atome":"By atom",focus.group_by_atom))focus.group_by_atom=true;
    validation::item("nbo.focus.group.atom");
    ImGui::SameLine();
    if(ImGui::RadioButton(zh?"按壳层":ja?"殻別":fr?"Par couche":"By shell",!focus.group_by_atom))focus.group_by_atom=false;
    validation::item("nbo.focus.group.shell");
    ImGui::SameLine();
    if(ImGui::Checkbox(zh?"配体按 s/p/d 分组":ja?"配位子を s/p/d 別":fr?"Ligands par s/p/d":"Group ligands by s/p/d",
                       &focus.group_ligands_by_l)){}
    validation::item("nbo.focus.group.mode.ligand_l");
    if(focus.visible_atoms.empty()) {
        ImGui::TextDisabled("%s",zh?"选择原子后绘图；完整数值仍可展开查看。":
            ja?"原子を選ぶと描画します。全数値は詳細で確認できます。":
            fr?"Choisissez des atomes pour tracer ; les valeurs complètes restent disponibles.":
            "Select atoms to draw; complete numeric data remains available below.");
    } else {
        ImGui::TextWrapped("%s",zh?"取消壳层勾选可隐藏图中的成分；不会删除数值。":
            ja?"殻の選択を外すと図で非表示になり、数値は残ります。":
            fr?"Décochez une couche pour la masquer du graphe ; les valeurs numériques restent.":
            "Uncheck a shell to hide it from the graph; its numeric value remains.");
        std::map<std::size_t,double> atom_weights;
        for(const auto& shell:view.shells)if(focus.visible_atoms.contains(shell.atom) &&
            shell_enabled(shell,focus))atom_weights[shell.atom]+=shell.weight;
        for(const auto& shell:view.shells) {
            if(!focus.visible_atoms.contains(shell.atom))continue;
            bool shown=shell_enabled(shell,focus);
            const std::string label=shell.symbol+std::to_string(shell.atom)+" "+shell.label+" "+shell.type+
                "##focus.shell."+shell.key;
            if(ImGui::Checkbox(label.c_str(),&shown)) {
                if(shown){focus.hidden_shells.erase(shell.key);focus.explicitly_included_shells.insert(shell.key);}
                else{focus.hidden_shells.insert(shell.key);focus.explicitly_included_shells.erase(shell.key);}
            }
            validation::item("nbo.focus.shell."+shell.key);
            if(shown && !focus.group_by_atom) {
                ImGui::SameLine();
                ImGui::SetNextItemWidth(std::max(90.0f*scale,ImGui::GetContentRegionAvail().x-105.0f*scale));
                const auto fraction=static_cast<float>(std::clamp(shell.weight,0.0,1.0));
                ImGui::ProgressBar(fraction,ImVec2(-1.0f,0.0f),fmt(shell.weight).c_str());
            }
        }
        if(focus.group_by_atom)for(const auto& [atom,weight]:atom_weights) {
            double full=0;
            for(const auto& group:view.decomposition->atoms)if(group.atom==atom)full+=group.weight;
            ImGui::Text("%s%zu selected-shell subtotal %.9g / full atom %.9g",
                        atoms[atom].c_str(),atom,weight,full);ImGui::SameLine();
            const auto fraction=static_cast<float>(std::clamp(weight,0.0,1.0));
            ImGui::ProgressBar(fraction,ImVec2(-1.0f,0.0f),fmt(weight).c_str());
        }
        const auto graph_shells=visible_focus_shells(view,focus);
        struct LiveNode {std::string label;double selected=0;std::optional<double> full;};
        std::vector<LiveNode> graph_nodes;
        if(focus.group_ligands_by_l) {
            const auto members=grouped_atoms(focus);
            for(const auto* shell:graph_shells)if(!members.contains(shell->atom))
                graph_nodes.push_back({shell->symbol+std::to_string(shell->atom)+" "+shell->label+" "+shell->type,shell->weight,{}});
            for(const auto& group:focus.ligand_groups) {
                const auto weights=angular_weights(view,focus,group);
                std::set<int> selected_l;
                for(const auto* shell:graph_shells)if(group.atoms.contains(shell->atom))selected_l.insert(shell->l);
                for(int l:selected_l)graph_nodes.push_back({"L"+std::to_string(group.id)+" "+shell_letter(l),
                    weights.at(l).selected,weights.at(l).full});
            }
        } else if(focus.group_by_atom)for(const auto& [atom,weight]:atom_weights)
            graph_nodes.push_back({atoms[atom]+std::to_string(atom)+" selected shells",weight,{}});
        else for(const auto* shell:graph_shells)
            graph_nodes.push_back({shell->symbol+std::to_string(shell->atom)+" "+shell->label+" "+shell->type,shell->weight,{}});
        if(!graph_nodes.empty()) {
            ImGui::SeparatorText(zh?"MO → NAO 成分图":ja?"MO → NAO 成分図":
                                 fr?"Graphe OM → NAO":"MO → NAO composition graph");
            validation::anchor("nbo.focus.graph");
            const float width=std::max(340.0f*scale,ImGui::GetContentRegionAvail().x);
            const float row=60.0f*scale;
            const float height=std::max(150.0f*scale,55.0f*scale+row*graph_nodes.size());
            ImGui::InvisibleButton("##nbo.focus.graph",ImVec2(width,height));
            validation::item("nbo.focus.graph.canvas");
            auto* draw=ImGui::GetWindowDrawList();
            const auto p=ImGui::GetItemRectMin();
            const float mo_x=p.x+80.0f*scale,mo_y=p.y+height*0.5f;
            draw->AddCircleFilled(ImVec2(mo_x,mo_y),24.0f*scale,IM_COL32(64,98,160,255));
            draw->AddText(ImVec2(mo_x-22.0f*scale,mo_y-7.0f*scale),IM_COL32(255,255,255,255),
                          ("MO "+std::to_string(*view.index+1)).c_str());
            const auto& mo=canonical.orbitals[*view.index];
            const std::string mo_identity="MO "+std::to_string(*view.index+1)+" ["+
                std::string(mo.spin==Spin::Beta?"beta":"alpha")+"]  E="+fmt(mo.energy_hartree)+" Ha";
            draw->AddText(ImVec2(p.x+8.0f*scale,p.y+5.0f*scale),IM_COL32(168,188,218,255),mo_identity.c_str());
            for(std::size_t i=0;i<graph_nodes.size();++i) {
                const float y=p.y+50.0f*scale+row*i;
                const float x=p.x+std::max(170.0f*scale,width*0.38f);
                draw->AddLine(ImVec2(mo_x+24.0f*scale,mo_y),ImVec2(x-8.0f*scale,y),
                              IM_COL32(74,151,206,215),1.5f*scale);
                draw->AddCircleFilled(ImVec2(x,y),6.0f*scale,IM_COL32(74,151,206,255));
                const float tx=x+14.0f*scale;
                const float wrap=std::max(80.0f,p.x+width-tx-10.0f*scale);
                const auto& node=graph_nodes[i];
                draw->AddText(ImGui::GetFont(),ImGui::GetFontSize(),
                              ImVec2(tx,y-24.0f*scale),IM_COL32(220,231,245,255),
                              node.label.c_str(),nullptr,wrap);
                const std::string selected_text="selected="+fmt(node.selected);
                draw->AddText(ImVec2(tx,y-4.0f*scale),IM_COL32(150,215,250,255),selected_text.c_str());
                if(node.full){const auto full_text="full l="+fmt(*node.full);
                    draw->AddText(ImVec2(tx,y+16.0f*scale),IM_COL32(170,188,210,255),full_text.c_str());}
            }
        }
    }
    bool details=focus.show_full_details;
    if(ImGui::Checkbox(zh?"显示全部 NAO 数值与来源":ja?"全 NAO 数値と出典":
                       fr?"Afficher toutes les valeurs NAO et sources":"Show all NAO values and sources",&details))
        focus.show_full_details=details;
    validation::item("nbo.focus.details");
    if(focus.show_full_details && view.decomposition) {
        ImGui::Text("MO %zu | %s | %s",*view.index+1,view.decomposition->status.c_str(),
                    view.decomposition->detail.c_str());
        ImGui::Text("NAO weight sum = %s; normalization error = %s",
                    shown_opt(view.decomposition->weight_sum).c_str(),
                    shown_opt(view.decomposition->normalization_error).c_str());
        for(const auto& row:view.decomposition->rows) {
            ImGui::TextWrapped("NAO %zu %s%zu %s %s n=%s l=%s c=%.9g |c|²=%.9g e=%s | %s",
                row.nao_id,row.symbol.c_str(),row.atom,row.type.c_str(),row.angular.c_str(),
                row.principal_n?std::to_string(*row.principal_n).c_str():"?",
                row.angular_l?std::to_string(*row.angular_l).c_str():"?",
                row.coefficient,row.weight,shown_opt(row.electron_contribution).c_str(),
                source_label(row.source).c_str());
        }
    }
    validation::field("nbo.focus.status",view.status+";mo="+std::to_string(*view.index)+
        ";shells="+std::to_string(view.shells.size()));
}

void draw_selected_nbo_context(NboUIState& s,const NboIntegration& integration,
                               const Wavefunction* canonical,
                               const MODiagramViewSnapshot* diagram,Language language) {
    if(s.selected_atoms.empty()&&!s.selected_structure)return;
    ImGui::SeparatorText(nbo_local(language,"Selected in 3D","三维中选中",
        "3D で選択中","Sélection 3D"));
    validation::anchor("nbo.context");
    std::string atoms;
    for(auto atom:s.selected_atoms) {
        if(!atoms.empty())atoms+=", ";
        if(canonical&&atom<canonical->atoms.size())
            atoms+=canonical->atoms[atom].symbol+std::to_string(atom+1);
        else atoms+="#"+std::to_string(atom+1);
    }
    if(!atoms.empty())note(std::string(nbo_local(language,"Atoms: ","原子：",
        "原子：","Atomes : "))+atoms);
    if(s.selected_structure&&*s.selected_structure<integration.structure.size()) {
        const auto& evidence=integration.structure[*s.selected_structure];
        note(evidence.label+" ["+evidence.kind+"]");
        note(evidence.detail);
        if(evidence.wiberg)note("Wiberg="+fmt(*evidence.wiberg));
        if(evidence.value)note("Value="+fmt(*evidence.value)+" "+evidence.units);
        note(source_label(evidence.source));
        for(const auto& ref:evidence.orbitals) {
            const auto* orbital=nbo_orbital(integration,ref);
            if(!orbital||orbital->coefficients.empty())continue;
            const auto label=orbital->label+"##nbo.context.orbital."+orbital->id;
            if(ImGui::Button(label.c_str()))
                s.aomo.pending_selection=nbo_single_selection(integration,ref);
            validation::item("nbo.context.orbital."+orbital->id);
            ImGui::SameLine();
        }
        if(!evidence.orbitals.empty())ImGui::NewLine();
        if(evidence.orbitals.size()>=2) {
            if(ImGui::Button((std::string(nbo_local(language,"Overlay both orbitals",
                    "叠加两个真实轨道","2つの実軌道を重ねる","Superposer les deux orbitales"))+
                    "##nbo.context.overlay").c_str())) {
                NboOrbitalSelection selection;selection.dataset_id=integration.id;
                selection.label=evidence.label;selection.mode=NboSelectionMode::Overlay;
                for(const auto& ref:evidence.orbitals)selection.terms.push_back({ref,1});
                s.aomo.pending_selection=std::move(selection);
            }
            validation::item("nbo.context.overlay");
        }
    }
    const bool show=ImGui::CollapsingHeader((std::string(nbo_local(language,
        "Related real orbitals and central MOs","相关真实轨道与中央 MO",
        "関連する実軌道と中央 MO","Orbitales réelles et OM centrales liées"))+
        "##nbo.context.related").c_str(),ImGuiTreeNodeFlags_DefaultOpen);
    validation::item("nbo.context.related");
    if(!show)return;
    ImGui::BeginChild("##nbo.context.list",ImVec2(0,235),ImGuiChildFlags_Border);
    std::size_t related=0;
    for(const auto& orbital:integration.orbitals) {
        if(orbital.ref.kind==NboOrbitalKind::Canonical ||
           orbital.ref.kind==NboOrbitalKind::GaussianAO ||
           orbital.coefficients.empty())continue;
        const bool touches=std::any_of(orbital.atoms.begin(),orbital.atoms.end(),
            [&](std::size_t atom){return s.selected_atoms.contains(atom);});
        if(!touches)continue;
        ++related;
        const auto label=orbital.label+" ["+nbo_orbital_kind_name(orbital.ref.kind)+
            "]##nbo.context.typed."+orbital.id;
        if(ImGui::Selectable(label.c_str())) {
            s.aomo.pending_selection=nbo_single_selection(integration,orbital.ref);
            if(orbital.ref.kind==NboOrbitalKind::NHO)s.inspected_nho=orbital.ref;
            if(orbital.ref.kind==NboOrbitalKind::NLMO)s.inspected_nlmo=orbital.ref;
        }
        validation::item("nbo.context.typed."+orbital.id);
    }
    if(!related)note(nbo_local(language,"No verified local orbital coefficients cover these atoms.",
        "这些原子没有可用的已验证局域轨道系数。",
        "これらの原子を含む検証済み局在軌道係数がありません。",
        "Aucun coefficient orbital local vérifié ne couvre ces atomes."));
    if(diagram) {
        std::set<std::size_t> central;
        for(const auto& level:diagram->data.levels) {
            if(level.member_indices.empty())central.insert(level.metadata.orbital_index);
            else central.insert(level.member_indices.begin(),level.member_indices.end());
            central.insert(level.member_spin_counterparts.begin(),level.member_spin_counterparts.end());
        }
        if(!central.empty())ImGui::SeparatorText("Central MO links through selected NAOs");
        for(auto index:central) {
            double weight=0;bool has=false;
            for(const auto& link:integration.links) {
                if(link.canonical_index!=index||link.orbital.kind!=NboOrbitalKind::NAO||
                   !link.weight)continue;
                const auto* descriptor=nbo_orbital(integration,link.orbital);
                if(!descriptor)continue;
                if(std::any_of(descriptor->atoms.begin(),descriptor->atoms.end(),
                    [&](std::size_t atom){return s.selected_atoms.contains(atom);})) {
                    has=true;weight+=*link.weight;
                }
            }
            if(!has)continue;
            const auto descriptor=std::find_if(integration.orbitals.begin(),integration.orbitals.end(),
                [&](const auto& orbital){return orbital.ref.kind==NboOrbitalKind::Canonical &&
                    orbital.ref.index==index;});
            const auto label="MO "+std::to_string(index+1)+"  selected-atom NAO |c|²="+
                fmt(weight)+"##nbo.context.mo."+std::to_string(index);
            if(ImGui::Selectable(label.c_str())&&descriptor!=integration.orbitals.end())
                s.aomo.pending_selection=nbo_single_selection(integration,descriptor->ref);
            validation::item("nbo.context.mo."+std::to_string(index));
        }
    }
    ImGui::EndChild();
}

NboUIActions draw_nbo_panel(NboUIState& s, Language language,
                            bool canonical_loaded, bool renderable,
                            bool active, float scale,
                            const Wavefunction* canonical, std::size_t canonical_index,
                            const MODiagramViewSnapshot* diagram) {
    (void)diagram;
    NboUIActions actions;
    const auto words=report_text(language);
    ImGui::Spacing();
    section_title(trn(language,Title));
    if(s.input_discovery) {
        note(s.input_status.empty()?"Detected NBO inputs":s.input_status);
        validation::field("nbo.input.status",s.input_status);
        const auto& discovery=*s.input_discovery;
        if(discovery.selection_required || discovery.candidates.size()>1) {
            note(language==Language::ChineseSimplified?
                "检测到多个计算或分析段；请选择要关联的候选。":
                language==Language::Japanese?
                "複数の計算または解析区間が見つかりました。候補を選択してください。":
                language==Language::French?
                "Plusieurs calculs ou segments ont été trouvés. Choisissez un candidat.":
                "Multiple calculations or analysis segments were found. Choose a candidate.");
            for(std::size_t i=0;i<discovery.candidates.size();++i) {
                const auto& c=discovery.candidates[i];
                std::string label=c.label.empty()?c.id:c.label;
                if(c.analysis_segment)label+=" [step "+std::to_string(*c.analysis_segment)+"]";
                label+="##nbo.candidate."+std::to_string(i);
                if(ImGui::Selectable(label.c_str()))s.pending_candidate=i;
                validation::item("nbo.candidate."+std::to_string(i));
                if(ImGui::IsItemHovered()) {
                    ImGui::BeginTooltip();
                    ImGui::Text("Canonical: %s",c.canonical.string().c_str());
                    ImGui::Text("NBO report: %s",c.report.string().c_str());
                    for(const auto& reason:c.diagnostics)ImGui::TextWrapped("%s",reason.c_str());
                    ImGui::EndTooltip();
                }
            }
        }
        for(const auto& diagnostic:discovery.diagnostics)note(diagnostic);
        const char* advanced=language==Language::ChineseSimplified?"高级：手动指定输入":
            language==Language::Japanese?"詳細：入力を手動指定":
            language==Language::French?"Avancé : saisir les fichiers":"Advanced: enter input paths manually";
        if(ImGui::Checkbox(advanced,&s.show_advanced_inputs)){}
        validation::item("nbo.inputs.advanced");
    }
    if(!s.input_discovery || s.show_advanced_inputs) {
    input_path(trn(language,Path),"##nbo.path",s.path);
    input_path(trn(language,Archive),"##nbo.archive47",s.archive47);
    input_path(trn(language,Aonbo),"##nbo.aonbo",s.aonbo);
    input_path(trn(language,Nbomo),"##nbo.nbomo",s.nbomo);
    input_path(trn(language,Naomo),"##nbo.naomo",s.naomo);
    input_path(trn(language,Aonao),"##nbo.aonao",s.aonao);
    input_path(trn(language,Naonbo),"##nbo.naonbo",s.naonbo);
    ImGui::TextWrapped("%s",trn(language,Segment));
    ImGui::SetNextItemWidth(100.0f*scale);
    ImGui::InputInt("##nbo.segment",&s.analysis_segment);
    validation::item("nbo.segment");
    ImGui::BeginDisabled(!canonical_loaded || active);
    actions.attach=ImGui::Button(trn(language,Attach),ImVec2(-1.0f,0));
    ImGui::EndDisabled();validation::item("nbo.attach");
    }
    if(!s.error.empty())note(s.error);
    if(!s.dataset)return actions;
    const auto& d=*s.dataset;
    if(s.integration)draw_selected_nbo_context(s,*s.integration,canonical,diagram,language);
    validation::field("nbo.association",d.association.status+": "+d.association.detail);
    const bool show_source_details=!s.integration || ImGui::CollapsingHeader(
        (std::string(nbo_local(language,"Source validation details","来源核验详情",
            "出典検証の詳細","Détails de validation de la source"))+
            "##nbo.source.validation").c_str());
    validation::item("nbo.source.validation");
    if(show_source_details) {
    note(std::string(trn(language,Association))+": "+d.association.status+" — "+d.association.detail);
    for(const auto& evidence:d.association.canonical_evidence){
        note(std::string(words.canonical_evidence)+" ["+nbo_spin_name(evidence.spin)+"] "+
             evidence.coefficient_source+"; "+words.direct_fchk+"="+
             (evidence.direct_fchk_coefficients?"1":"0")+"; "+words.density_verified+"="+
             (evidence.density_verified?"1":"0"));
    }
    }
    if(s.integration) {
        const auto& integrated=*s.integration;
        if(s.inspected_dataset_id!=integrated.id) {
            s.inspected_dataset_id=integrated.id;
            s.inspected_nho.reset();s.inspected_nlmo.reset();
            s.nho_sum_owner.reset();s.nho_sum_nao_indices.clear();
        }
        ImGui::TextWrapped("%s",nbo_local(language,
            "Verified capabilities and 3D controls are available below.",
            "下方可使用已验证的能力与三维控件。",
            "検証済み機能と3D操作は以下にあります。",
            "Les fonctions vérifiées et commandes 3D figurent ci-dessous."));
        const bool show_capabilities=ImGui::CollapsingHeader(
            (std::string(nbo_local(language,"Data capability details","数据能力详情",
                "データ機能の詳細","Détails des capacités des données"))+
                "##nbo.capabilities").c_str());
        validation::item("nbo.capabilities");
        for(const char* key:{"source_association","report","aomo","aomo_full","nao","pnao","nho","nbo","nlmo",
                              "charges","wiberg","interactions","structure"}) {
            if(const auto* capability=nbo_capability(integrated,key)) {
                if(show_capabilities)note(std::string(key)+" — "+nbo_capability_state_name(capability->state)+
                    (capability->detail.empty()?"":": "+capability->detail));
                validation::field(std::string("nbo.capability.")+key,
                    std::string(nbo_capability_state_name(capability->state))+": "+capability->detail);
            }
        }
        ImGui::SeparatorText(nbo_local(language,"3D structure and orbital evidence",
            "三维结构与轨道证据","3D 構造と軌道の証拠","Preuves 3D de structure et d'orbitales"));
        const auto* charges_cap=nbo_capability(integrated,"charges");
        const auto* wiberg_cap=nbo_capability(integrated,"wiberg");
        const auto* e2_cap=nbo_capability(integrated,"interactions");
        const bool charges_available=charges_cap&&charges_cap->available();
        const bool spin_available=std::any_of(d.populations.begin(),d.populations.end(),
            [](const auto& row){return row.spin_density.has_value();});
        if((s.atom_colour_mode==1&&!charges_available) ||
           (s.atom_colour_mode==2&&!spin_available))s.atom_colour_mode=0;
        ImGui::RadioButton((std::string(nbo_local(language,"Element","元素","元素","Élément"))+
            "##nbo.atom.color.element").c_str(),&s.atom_colour_mode,0);
        validation::item("nbo.atom.color.element");ImGui::SameLine();
        ImGui::BeginDisabled(!charges_available);
        ImGui::RadioButton((std::string(nbo_local(language,"NPA charge","NPA 电荷","NPA 電荷","Charge NPA"))+
            "##nbo.atom.color.charge").c_str(),&s.atom_colour_mode,1);
        validation::item("nbo.atom.color.charge");
        if(!charges_available&&ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
            ImGui::SetTooltip("%s",charges_cap?charges_cap->detail.c_str():"NPA charges unavailable");
        ImGui::EndDisabled();ImGui::SameLine();
        ImGui::BeginDisabled(!spin_available);
        ImGui::RadioButton((std::string(nbo_local(language,"Spin","自旋","スピン","Spin"))+
            "##nbo.atom.color.spin").c_str(),&s.atom_colour_mode,2);
        validation::item("nbo.atom.color.spin");
        if(!spin_available&&ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
            ImGui::SetTooltip("No verified atom spin density is available");
        ImGui::EndDisabled();
        ImGui::BeginDisabled(!wiberg_cap||!wiberg_cap->available());
        ImGui::Checkbox((std::string(nbo_local(language,"Bond indices","键级指数","結合指数","Indices de liaison"))+
            "##nbo.overlay.bond").c_str(),&s.show_bond_indices);
        validation::item("nbo.overlay.wiberg");
        if(ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled) &&
           (!wiberg_cap||!wiberg_cap->available()))
            ImGui::SetTooltip("%s",wiberg_cap?wiberg_cap->detail.c_str():"Wiberg indices unavailable");
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(!e2_cap||!e2_cap->available());
        ImGui::Checkbox((std::string(nbo_local(language,"E(2) interactions","E(2) 相互作用",
            "E(2) 相互作用","Interactions E(2)"))+"##nbo.overlay.e2").c_str(),&s.show_e2);
        validation::item("nbo.overlay.e2");
        if(ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled) &&
           (!e2_cap||!e2_cap->available()))
            ImGui::SetTooltip("%s",e2_cap?e2_cap->detail.c_str():"E(2) interactions unavailable");
        ImGui::EndDisabled();
        const bool structure_open=ImGui::CollapsingHeader((std::string(nbo_local(language,"Structure evidence",
            "结构证据","構造の証拠","Preuves structurales"))+"##nbo.structure").c_str());
        validation::item("nbo.structure");
        if(structure_open) {
            for(std::size_t i=0;i<integrated.structure.size();++i) {
                const auto& evidence=integrated.structure[i];
                const auto label=evidence.label+" ["+evidence.kind+"]##nbo.structure."+std::to_string(i);
                if(ImGui::Selectable(label.c_str(),s.selected_structure==i)) {
                    s.selected_structure=i;
                    s.selected_atoms.clear();
                    s.selected_atoms.insert(evidence.atoms.begin(),evidence.atoms.end());
                    std::string kind=evidence.kind;
                    std::transform(kind.begin(),kind.end(),kind.begin(),[](unsigned char c){return static_cast<char>(std::tolower(c));});
                    if((kind=="e2" || kind=="interaction" || kind=="donor_acceptor") &&
                       evidence.orbitals.size()>=2) {
                        NboOrbitalSelection selection;
                        selection.dataset_id=integrated.id;
                        selection.label=evidence.label;
                        selection.mode=NboSelectionMode::Overlay;
                        for(const auto& ref:evidence.orbitals)selection.terms.push_back({ref,1});
                        s.aomo.pending_selection=std::move(selection);
                    } else if(!evidence.orbitals.empty() &&
                       (evidence.orbitals[0].kind==NboOrbitalKind::NBO ||
                        evidence.orbitals[0].kind==NboOrbitalKind::NHO ||
                        evidence.orbitals[0].kind==NboOrbitalKind::NLMO))
                        s.aomo.pending_selection=nbo_single_selection(integrated,evidence.orbitals[0]);
                }
                validation::item("nbo.structure."+std::to_string(i));
                if(ImGui::IsItemHovered())ImGui::SetTooltip("%s\n%s",evidence.detail.c_str(),
                    source_label(evidence.source).c_str());
            }
        }
        for(const auto kind:{NboOrbitalKind::NAO,NboOrbitalKind::PNAO,
                             NboOrbitalKind::NHO,NboOrbitalKind::NBO,NboOrbitalKind::NLMO}) {
            const std::string key=nbo_orbital_kind_name(kind);
            std::string capability_key=key;
            std::transform(capability_key.begin(),capability_key.end(),capability_key.begin(),
                [](unsigned char c){return static_cast<char>(std::tolower(c));});
            const auto* family_cap=nbo_capability(integrated,capability_key);
            const bool family_available=family_cap&&family_cap->available();
            ImGui::BeginDisabled(!family_available);
            const bool family_open=ImGui::CollapsingHeader((key+" "+nbo_local(language,"actual orbitals",
                "真实轨道","実軌道","orbitales réelles")+"##nbo.family."+key).c_str());
            validation::item("nbo.family."+key);
            if(!family_available&&ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                ImGui::SetTooltip("%s",family_cap?family_cap->detail.c_str():
                    "No verified coefficients for this orbital family");
            ImGui::EndDisabled();
            if(!family_available) {
                ImGui::TextDisabled("%s",family_cap?family_cap->detail.c_str():
                    "No verified coefficients for this orbital family");
                continue;
            }
            if(!family_open)continue;
            std::size_t shown=0;
            for(const auto& orbital:integrated.orbitals) {
                if(orbital.ref.kind!=kind || orbital.coefficients.empty())continue;
                ++shown;
                const std::string label=orbital.label+" ["+nbo_spin_name(orbital.ref.spin)+"]"+
                    (orbital.occupation?" occ="+fmt(*orbital.occupation):"")+
                    "##nbo.typed."+orbital.id;
                const bool selected=s.aomo.selection && s.aomo.selection->terms.size()==1 &&
                    s.aomo.selection->terms.front().orbital==orbital.ref;
                if(ImGui::Selectable(label.c_str(),selected)) {
                    s.aomo.pending_selection=nbo_single_selection(integrated,orbital.ref);
                    s.selected_atoms.clear();
                    s.selected_atoms.insert(orbital.atoms.begin(),orbital.atoms.end());
                    s.selected_structure.reset();
                    if(kind==NboOrbitalKind::NLMO)s.inspected_nlmo=orbital.ref;
                    if(kind==NboOrbitalKind::NHO)s.inspected_nho=orbital.ref;
                }
                validation::item("nbo.typed."+orbital.id);
                if(ImGui::IsItemHovered())ImGui::SetTooltip("%s\n%s\n%s",
                    orbital.detail.c_str(),orbital.energy_semantics.c_str(),
                    source_label(orbital.source).c_str());
            }
            if(!shown)note(nbo_local(language,
                "No verified real-space coefficients for this family.",
                "该轨道类没有已验证的真实空间系数。",
                "この軌道系列には検証済み実空間係数がありません。",
                "Aucun coefficient spatial vérifié pour cette famille."));
        }
        if(s.inspected_nlmo && nbo_orbital(integrated,*s.inspected_nlmo)) {
            const auto ref=*s.inspected_nlmo;
            const auto components=nbo_nlmo_components(integrated,ref);
            const auto parent=nbo_nlmo_parent(integrated,ref);
            const bool nlmo_details_open=ImGui::CollapsingHeader((std::string(nbo_local(language,
                    "NLMO main part and delocalization tail",
                    "NLMO 主成分与离域尾部","NLMO 主成分と非局在化テール",
                    "Composante principale et queue délocalisée NLMO"))+
                    "##nbo.nlmo.components").c_str(),ImGuiTreeNodeFlags_DefaultOpen);
            validation::item("nbo.nlmo.components");
            if(nlmo_details_open) {
                ImGui::Text("NLMO %zu [%s]",ref.index+1,nbo_spin_name(ref.spin));
                if(!parent)note(nbo_local(language,
                    "No printed and verified parent NBO: main/tail assignment unavailable.",
                    "没有已打印并核验的母 NBO：主成分/尾部划分不可用。",
                    "印字・検証済みの親 NBO がなく、主成分とテールは区別できません。",
                    "Aucun NBO parent imprimé et vérifié : décomposition principale/queue indisponible."));
                else {
                    std::vector<NboOrbitalTerm> tail;
                    for(const auto& component:components)if(!(component.orbital==*parent))
                        tail.push_back(component);
                    const auto main=std::find_if(components.begin(),components.end(),
                        [&](const auto& component){return component.orbital==*parent;});
                    if(main!=components.end()) {
                        const auto* parent_orbital=nbo_orbital(integrated,*parent);
                        const std::string label=(parent_orbital?parent_orbital->label:"Parent NBO")+
                            "  c="+fmt(main->coefficient)+"##nbo.nlmo.main";
                        if(ImGui::Button(label.c_str())) {
                            NboOrbitalSelection selection;selection.dataset_id=integrated.id;
                            selection.label="Verified NLMO main component";
                            selection.mode=NboSelectionMode::WeightedComponent;
                            selection.terms={*main};
                            s.aomo.pending_selection=std::move(selection);
                        }
                        validation::item("nbo.nlmo.main");
                    }
                    if(!tail.empty()) {
                        const auto label=std::string(nbo_local(language,
                            "Show actual tail partial sum","查看真实尾部部分和",
                            "実際のテール部分和を表示","Afficher la somme partielle de la queue"))+
                            "##nbo.nlmo.tail";
                        if(ImGui::Button(label.c_str())) {
                            NboOrbitalSelection selection;selection.dataset_id=integrated.id;
                            selection.label="NLMO delocalization tail (validated NBO terms)";
                            selection.mode=NboSelectionMode::PartialSum;
                            selection.terms=tail;
                            s.aomo.pending_selection=std::move(selection);
                        }
                        validation::item("nbo.nlmo.tail");
                    }
                }
                note(nbo_local(language,"All signed NBONLMO coefficients:",
                    "全部带符号 NBONLMO 系数：","すべての符号付き NBONLMO 係数：",
                    "Tous les coefficients NBONLMO signés :"));
                for(const auto& term:components) {
                    const auto* orbital=nbo_orbital(integrated,term.orbital);
                    note((orbital?orbital->label:std::string("NBO ")+std::to_string(term.orbital.index+1))+
                        "  c="+fmt(term.coefficient)+
                        (parent && term.orbital==*parent?" [verified parent]":""));
                }
            }
        }
        if(s.inspected_nho && nbo_orbital(integrated,*s.inspected_nho)) {
            const auto ref=*s.inspected_nho;
            if(!s.nho_sum_owner || !(*s.nho_sum_owner==ref)) {
                s.nho_sum_owner=ref;
                s.nho_sum_nao_indices.clear();
            }
            const bool nho_details_open=ImGui::CollapsingHeader((std::string(nbo_local(language,
                    "NHO hybrid components and 3D lobe",
                    "NHO 杂化成分与三维瓣","NHO 混成成分と 3D ローブ",
                    "Composantes hybrides NHO et lobe 3D"))+
                    "##nbo.nho.components").c_str(),ImGuiTreeNodeFlags_DefaultOpen);
            validation::item("nbo.nho.components");
            if(nho_details_open) {
                const auto terms=nbo_nho_components(integrated,ref);
                if(terms.empty())note(nbo_local(language,
                    "Validated NAO–NHO components are unavailable.",
                    "已验证的 NAO–NHO 成分不可用。",
                    "検証済み NAO–NHO 成分を利用できません。",
                    "Composantes NAO–NHO vérifiées indisponibles."));
                else {
                    std::map<char,double> angular_weights;
                    std::map<char,std::vector<NboOrbitalTerm>> angular_terms;
                    double unknown=0,total=0;
                    for(const auto& term:terms) {
                        const double weight=term.coefficient*term.coefficient;
                        total+=weight;
                        const auto row=std::find_if(d.naos.begin(),d.naos.end(),
                            [&](const NboNao& nao){return nao.id==term.orbital.index+1 &&
                                nao.spin==term.orbital.spin;});
                        char angular='?';
                        if(row!=d.naos.end())for(unsigned char c:row->angular)
                            if(std::isalpha(c)){angular=static_cast<char>(std::tolower(c));break;}
                        if(angular=='s'||angular=='p'||angular=='d'||angular=='f'||angular=='g')
                            {angular_weights[angular]+=weight;angular_terms[angular].push_back(term);}
                        else unknown+=weight;
                    }
                    ImGui::Text("NHO %zu [%s]  NAO sum |c|²=%s",
                        ref.index+1,nbo_spin_name(ref.spin),fmt(total).c_str());
                    for(const auto& [angular,weight]:angular_weights)
                        ImGui::Text("%c  |c|²=%s  (%s%% of this verified NHO)",
                            angular,fmt(weight).c_str(),fmt(total>0?100*weight/total:0).c_str());
                    if(unknown>0)ImGui::Text("?  |c|²=%s (NAO angular label unavailable)",fmt(unknown).c_str());
                    ImGui::TextWrapped("%s",nbo_local(language,
                        "Angular selections and accumulated terms are unnormalized signed NAO partial sums; percentages above describe the complete verified NHO.",
                        "角动量分组及累加项均为保留原始符号且未归一化的 NAO 部分和；上方百分比描述完整已验证 NHO。",
                        "角運動量群と累積項は符号を保ち、正規化しない NAO 部分和です。上の割合は検証済み NHO 全体を示します。",
                        "Les groupes angulaires et termes cumulés sont des sommes partielles NAO signées, sans normalisation ; les pourcentages décrivent le NHO vérifié complet."));
                    for(const auto angular:{'s','p','d','f','g'}) {
                        const auto group=angular_terms.find(angular);
                        if(group==angular_terms.end() || angular_weights[angular]<=0)continue;
                        const std::string label=std::string(1,angular)+" "+nbo_local(language,
                            "partial","部分和","部分和","somme partielle")+"  |c|²="+
                            fmt(angular_weights[angular])+"##nbo.nho.angular."+angular;
                        if(ImGui::SmallButton(label.c_str())) {
                            NboOrbitalSelection selection;selection.dataset_id=integrated.id;
                            selection.label="NHO "+std::to_string(ref.index+1)+" signed "+angular+
                                " NAO partial sum";
                            selection.mode=NboSelectionMode::PartialSum;
                            selection.terms=group->second;
                            s.aomo.pending_selection=std::move(selection);
                        }
                        validation::item(std::string("nbo.nho.angular.")+angular);
                    }
                    ImGui::Text("%s: %zu / %zu",nbo_local(language,
                        "Selected signed NAO terms","已选带符号 NAO 项",
                        "選択した符号付き NAO 項","Termes NAO signés sélectionnés"),
                        s.nho_sum_nao_indices.size(),terms.size());
                    if(ImGui::SmallButton((std::string(nbo_local(language,
                        "Select all","全选","すべて選択","Tout sélectionner"))+"##nbo.nho.sum.all").c_str()))
                        for(const auto& term:terms)s.nho_sum_nao_indices.insert(term.orbital.index);
                    validation::item("nbo.nho.sum.all");ImGui::SameLine();
                    if(ImGui::SmallButton((std::string(nbo_local(language,
                        "Clear","清空","クリア","Effacer"))+"##nbo.nho.sum.clear").c_str()))
                        s.nho_sum_nao_indices.clear();
                    validation::item("nbo.nho.sum.clear");
                    std::vector<NboOrbitalTerm> selected_terms;
                    for(const auto& term:terms)
                        if(s.nho_sum_nao_indices.contains(term.orbital.index))selected_terms.push_back(term);
                    ImGui::BeginDisabled(selected_terms.empty());
                    if(ImGui::Button((std::string(nbo_local(language,
                        "Show signed partial sum","查看带符号部分和",
                        "符号付き部分和を表示","Afficher la somme partielle signée"))+
                        "##nbo.nho.sum.partial").c_str())) {
                        NboOrbitalSelection selection;selection.dataset_id=integrated.id;
                        selection.label="NHO "+std::to_string(ref.index+1)+
                            " selected signed NAO partial sum";
                        selection.mode=NboSelectionMode::PartialSum;
                        selection.terms=selected_terms;
                        s.aomo.pending_selection=std::move(selection);
                    }
                    validation::item("nbo.nho.sum.partial");
                    if(ImGui::Button((std::string(nbo_local(language,
                        "Overlay selected real terms","分别叠加已选真实项",
                        "選択した実際の項を重ねて表示","Superposer les termes réels sélectionnés"))+
                        "##nbo.nho.sum.overlay").c_str())) {
                        NboOrbitalSelection selection;selection.dataset_id=integrated.id;
                        selection.label="NHO "+std::to_string(ref.index+1)+
                            " separate signed NAO components";
                        selection.mode=NboSelectionMode::Overlay;
                        selection.terms=selected_terms;
                        s.aomo.pending_selection=std::move(selection);
                    }
                    validation::item("nbo.nho.sum.overlay");
                    ImGui::EndDisabled();
                    if(ImGui::Button((std::string(nbo_local(language,
                        "Full actual NHO","完整真实 NHO","実際の NHO 全体","NHO réel complet"))+
                        "##nbo.nho.sum.full").c_str()))
                        s.aomo.pending_selection=nbo_single_selection(integrated,ref);
                    validation::item("nbo.nho.sum.full");
                    const bool all_nho_terms=ImGui::TreeNode((std::string(nbo_local(language,
                        "Every signed NAO–NHO coefficient","全部带符号 NAO–NHO 系数",
                        "すべての符号付き NAO–NHO 係数","Tous les coefficients NAO–NHO signés"))+
                        "##nbo.nho.all").c_str());
                    validation::item("nbo.nho.all");
                    if(all_nho_terms) {
                        for(const auto& term:terms) {
                            const auto* orbital=nbo_orbital(integrated,term.orbital);
                            bool checked=s.nho_sum_nao_indices.contains(term.orbital.index);
                            if(ImGui::Checkbox(("##nbo.nho.sum.item."+
                                std::to_string(term.orbital.index)).c_str(),&checked)) {
                                if(checked)s.nho_sum_nao_indices.insert(term.orbital.index);
                                else s.nho_sum_nao_indices.erase(term.orbital.index);
                            }
                            validation::item("nbo.nho.sum.item."+
                                std::to_string(term.orbital.index));
                            ImGui::SameLine();
                            const std::string label=(orbital?orbital->label:std::string("NAO ")+
                                std::to_string(term.orbital.index+1))+"  c="+fmt(term.coefficient)+
                                "##nbo.nho.nao."+std::to_string(term.orbital.index);
                            if(ImGui::Selectable(label.c_str())) {
                                NboOrbitalSelection selection;selection.dataset_id=integrated.id;
                                selection.label="NAO component of verified NHO";
                                selection.mode=NboSelectionMode::WeightedComponent;
                                selection.terms={term};
                                s.aomo.pending_selection=std::move(selection);
                            }
                            validation::item("nbo.nho.nao."+std::to_string(term.orbital.index));
                        }
                        ImGui::TreePop();
                    }
                }
            }
        }
    } else note(renderable?trn(language,Renderable):trn(language,ReportsOnly));
    note(std::string(trn(language,Source))+": "+source_label(d.source));
    const bool fchk_has_beta=canonical && std::any_of(canonical->orbitals.begin(),canonical->orbitals.end(),
        [](const MolecularOrbital& mo){return mo.spin==Spin::Beta;});
    if(d.association.compatible && d.archive && d.archive->open_shell && !fchk_has_beta)
        note(words.beta_archive_notice);
    ImGui::BeginDisabled(!active);
    actions.canonical_set=ImGui::Button(trn(language,Canonical),ImVec2((ImGui::GetContentRegionAvail().x-ImGui::GetStyle().ItemSpacing.x)/2,0));
    ImGui::EndDisabled();validation::item("nbo.set.canonical");
    ImGui::SameLine();ImGui::BeginDisabled(!renderable || !canonical_loaded || active);
    actions.nbo_set=ImGui::Button(trn(language,NboSet),ImVec2(-1,0));
    ImGui::EndDisabled();validation::item("nbo.set.nbo");
    note(trn(language,EnergyNote));
    note(trn(language,FocusNote));
    if(ImGui::CollapsingHeader(words.contributions)){
        note("Complete NBO-to-MO matrix report; canonical MO/NAO focus uses the central diagram above.");
        if(!canonical || canonical_index>=canonical->orbitals.size() || !d.association.compatible){
            note(words.mapping_error);
        } else {
            const auto& mo=canonical->orbitals[canonical_index];
            const auto target_spin=mo.spin==Spin::Beta?NboSpin::Beta:NboSpin::Alpha;
            const auto evidence=std::find_if(d.association.canonical_evidence.begin(),
                d.association.canonical_evidence.end(),[&](const NboCanonicalEvidence& x){return x.spin==target_spin || (target_spin==NboSpin::Alpha && x.spin==NboSpin::Total);});
            const bool direct=evidence==d.association.canonical_evidence.end() || evidence->direct_fchk_coefficients;
            const auto source_index=mo.source_orbital_index;
            note(std::string(words.selected_canonical)+": "+
                 (source_index==std::numeric_limits<std::size_t>::max()?"?":std::to_string(source_index+1))+" ["+
                 (target_spin==NboSpin::Beta?"beta":"alpha")+"]");
            const NboMatrix* matrix=nullptr;
            auto accept=[&](const NboMatrix& m){
                std::string kind=m.kind;
                std::transform(kind.begin(),kind.end(),kind.begin(),[](unsigned char c){return static_cast<char>(std::toupper(c));});
                return kind=="NBOMO";
            };
            for(const auto& m:d.matrices)if(accept(m) && m.spin==target_spin){matrix=&m;break;}
            if(!matrix && target_spin==NboSpin::Alpha)
                for(const auto& m:d.matrices)if(accept(m) && m.spin==NboSpin::Total){matrix=&m;break;}
            if(!matrix && d.archive){
                for(const auto& m:d.archive->matrices)if(accept(m) && m.spin==target_spin){matrix=&m;break;}
                if(!matrix && target_spin==NboSpin::Alpha)
                    for(const auto& m:d.archive->matrices)if(accept(m) && m.spin==NboSpin::Total){matrix=&m;break;}
            }
            if(!direct){note(words.mapping_error);}else if(!matrix){note(words.no_matrix);}else{
                std::vector<std::size_t> row_to_orbital;
                for(std::size_t i=0;i<d.orbitals.size();++i)
                    if(d.orbitals[i].spin==matrix->spin)row_to_orbital.push_back(i);
                if(source_index==std::numeric_limits<std::size_t>::max() ||
                   source_index>=matrix->columns || row_to_orbital.size()!=matrix->rows ||
                   matrix->values.size()!=matrix->rows*matrix->columns){
                    note(words.mapping_error);
                }else{
                    struct Contribution {std::size_t orbital;double weight;};
                    std::vector<Contribution> entries;
                    double total=0,shown=0;
                    for(std::size_t row=0;row<matrix->rows;++row){
                        const double t=matrix->values[row*matrix->columns+source_index];
                        const double weight=t*t;
                        total+=weight;
                        entries.push_back({row_to_orbital[row],weight});
                    }
                    std::stable_sort(entries.begin(),entries.end(),[](const auto& a,const auto& b){return a.weight>b.weight;});
                    for(const auto& entry:entries){
                        const auto& orbital=d.orbitals[entry.orbital];
                        const std::string label="#"+std::to_string(orbital.id)+
                            " ["+nbo_spin_name(orbital.spin)+"] |T|²="+fmt(entry.weight);
                        ImGui::PushID(static_cast<int>(entry.orbital));
                        if(ImGui::Selectable(label.c_str(),s.selected_orbital==entry.orbital))
                            actions.selected_orbital=entry.orbital;
                        validation::item("nbo.contribution."+std::to_string(entry.orbital));
                        if(ImGui::IsItemHovered())ImGui::SetTooltip("%s",orbital.label.c_str());
                        ImGui::PopID();
                        shown+=entry.weight;
                    }
                    note(std::string(words.shown)+": "+fmt(shown)+"; "+words.remainder+": "+fmt(total-shown)+
                         "; Σ|T|²="+fmt(total));
                    note(std::string(words.source_evidence)+": "+source_label(matrix->source));
                    validation::field("nbo.nbomo", "spin="+std::string(nbo_spin_name(matrix->spin))+
                        ";canonical_source_index="+std::to_string(source_index)+
                        ";shown="+fmt(shown)+";remaining="+fmt(total-shown));
                }
            }
        }
    }
    if(ImGui::CollapsingHeader(trn(language,Populations),ImGuiTreeNodeFlags_DefaultOpen)){
        if(d.populations.empty())note(trn(language,NoData));
        for(const auto& p:d.populations){
            const auto label=std::string(words.atom)+" "+std::to_string(p.atom)+" "+p.symbol+" ["+nbo_spin_name(p.spin)+"]  "+words.charge+"="+fmt(p.charge)+" "+words.core+"="+fmt(p.core)+" "+words.valence+"="+fmt(p.valence)+" "+words.rydberg+"="+fmt(p.rydberg)+" "+words.total+"="+fmt(p.total)+" "+words.effective_core+"="+shown_opt(p.effective_core_electrons)+" "+words.explicit_pop+"="+shown_opt(p.explicit_population)+" "+words.spin_density+"="+shown_opt(p.spin_density);
            if(ImGui::Selectable((label+"##nbo.npa."+std::to_string(p.atom)+nbo_spin_name(p.spin)).c_str(),
                p.atom>0 && s.selected_atoms.contains(p.atom-1))) {
                s.selected_atoms.clear();if(p.atom>0)s.selected_atoms.insert(p.atom-1);
                s.selected_structure.reset();
            }
            validation::item("nbo.npa.atom."+std::to_string(p.atom));
        }
    }
    if(ImGui::CollapsingHeader(trn(language,Naos))){
        if(d.naos.empty())note(trn(language,NoData));
        for(const auto& n:d.naos){
            note(std::to_string(n.id)+" "+words.atom+" "+std::to_string(n.atom)+" "+n.symbol+" "+n.angular+" "+n.type+" ["+nbo_spin_name(n.spin)+"] "+words.occupation+"="+fmt(n.occupation)+" "+words.fock+"="+shown_opt(n.energy_hartree)+" "+words.spin_density+"="+shown_opt(n.spin_density));
        }
    }
    if(ImGui::CollapsingHeader(trn(language,Orbitals),ImGuiTreeNodeFlags_DefaultOpen)){
        if(d.orbitals.empty())note(trn(language,NoData));
        for(std::size_t i=0;i<d.orbitals.size();++i){
            const auto& o=d.orbitals[i];ImGui::PushID(static_cast<int>(i));
            const std::string label="#"+std::to_string(o.id)+" ["+nbo_spin_name(o.spin)+"] "+words.occupation+"="+fmt(o.occupation);
            if(ImGui::Selectable(label.c_str(),s.selected_orbital==i))actions.selected_orbital=i;
            validation::item("nbo.orbital."+std::to_string(i));
            if(ImGui::IsItemClicked() && s.integration) {
                const NboOrbitalRef ref{NboOrbitalKind::NBO,o.spin,o.id>0?o.id-1:0};
                if(o.id>0 && nbo_orbital(*s.integration,ref))
                    s.aomo.pending_selection=nbo_single_selection(*s.integration,ref);
            }
            if(ImGui::IsItemHovered())ImGui::SetTooltip("%s %s",o.kind.c_str(),o.label.c_str());
            if(s.selected_orbital==i){
                note(o.kind+" "+o.label);
                note(std::string(words.fock)+": "+(o.diagonal_fock_hartree?fmt(*o.diagonal_fock_hartree):"null"));
                for(const auto& c:o.components)note(std::string(words.component)+" "+words.atom+" "+std::to_string(c.atom)+"  "+fmt(c.percent)+"%  "+words.coefficient+"="+fmt(c.coefficient)+"  "+c.hybrid);
            }
            ImGui::PopID();
        }
    }
    if(ImGui::CollapsingHeader(trn(language,Wiberg))){
        if(d.wiberg.empty())note(trn(language,NoData));
        for(const auto& w:d.wiberg)note(std::to_string(w.atom_a)+"–"+std::to_string(w.atom_b)+" ["+nbo_spin_name(w.spin)+"] = "+fmt(w.value));
    }
    if(ImGui::CollapsingHeader(trn(language,E2))){
        for(const auto& section:d.e2_sections)note(std::string(trn(language,Threshold))+": "+shown_opt(section.printing_threshold)+" "+section.units+" ["+nbo_spin_name(section.spin)+"] "+section.missing_reason);
        if(d.e2.empty())note(trn(language,NoData));
        for(const auto& e:d.e2)note(std::to_string(e.donor)+" → "+std::to_string(e.acceptor)+" ["+nbo_spin_name(e.spin)+"] E(2)="+fmt(e.value)+" "+e.units+" Δε="+fmt(e.energy_gap_hartree)+" Ha F="+fmt(e.fock_hartree)+" Ha  "+words.threshold+"="+shown_opt(e.printing_threshold));
    }
    if(ImGui::CollapsingHeader(trn(language,Matrices))){
        for(const auto& m:d.matrices)note(m.kind+" ["+nbo_spin_name(m.spin)+"] "+std::to_string(m.rows)+"×"+std::to_string(m.columns));
        if(d.archive)for(const auto& m:d.archive->matrices)note("ARCHIVE "+m.kind+" ["+nbo_spin_name(m.spin)+"] "+std::to_string(m.rows)+"×"+std::to_string(m.columns));
        for(const auto& c:d.cmo_summaries)note("CMO summary: "+c.block);
    }
    if(ImGui::CollapsingHeader(words.source_evidence)){
        note(source_label(d.source));
        for(const auto& m:d.matrices)note(source_label(m.source));
        if(d.archive)note(source_label(d.archive->source));
        for(const auto& c:d.cmo_summaries)note(source_label(c));
        if(s.selected_orbital<d.orbitals.size())note(source_label(d.orbitals[s.selected_orbital].source));
        if(s.selected_orbital<d.orbitals.size() && d.orbitals[s.selected_orbital].energy_source)
            note(source_label(*d.orbitals[s.selected_orbital].energy_source));
        for(const auto& section:d.e2_sections)note(source_label(section.source));
    }
    input_path(trn(language,ExportPath),"##nbo.export.path",s.export_path);
    actions.export_bundle=ImGui::Button(trn(language,Export),ImVec2(-1,0));
    validation::item("nbo.export");
    if(!s.export_status.empty())note(s.export_status);
    return actions;
}

void export_nbo_focus_bundle(const NboDataset& d,const std::filesystem::path& base,
                             const MODiagramViewSnapshot* diagram,
                             const NboFocusUIState* focus) {
    const auto file=[&](const char* suffix){auto p=base;p+=suffix;return p;};
    const NboFocusUIState empty;
    const auto& state=focus?*focus:empty;
    const auto v=make_focus_view(d,diagram,focus);
    const auto shells=visible_focus_shells(v,state);
    struct Node{std::string id,label;std::size_t atom=0;std::optional<int> n,l;double weight=0;
                std::optional<double> full_atom_weight;std::string reason;
                std::vector<std::size_t> atom_ids;std::optional<double> full_angular_weight;};
    std::vector<Node> nodes;
    if(v.status=="available") {
        if(state.group_ligands_by_l) {
            const auto members=grouped_atoms(state);
            for(const auto* shell:shells)if(!members.contains(shell->atom))
                nodes.push_back({"shell:"+shell->key,
                    shell->symbol+std::to_string(shell->atom)+" "+shell->label+" "+shell->type,
                    shell->atom,shell->n,shell->l,shell->weight,std::nullopt,
                    "ungrouped_atom_selected_shell"});
            for(const auto& group:state.ligand_groups) {
                const auto weights=angular_weights(v,state,group);
                std::set<int> selected_l;
                for(const auto* shell:shells)if(group.atoms.contains(shell->atom))selected_l.insert(shell->l);
                for(int l:selected_l) {
                    Node node;node.id="ligand:L"+std::to_string(group.id)+":l:"+std::to_string(l);
                    node.label="L"+std::to_string(group.id)+" "+shell_letter(l)+" selected-shell subtotal";
                    node.l=l;node.weight=weights.at(l).selected;
                    node.full_angular_weight=weights.at(l).full;
                    node.atom_ids.assign(group.atoms.begin(),group.atoms.end());
                    node.reason="saved_ligand_group_and_enabled_shells";
                    nodes.push_back(std::move(node));
                }
            }
        } else if(state.group_by_atom) {
            std::map<std::size_t,Node> grouped;
            for(const auto* shell:shells){auto& node=grouped[shell->atom];
                node={"atom:"+std::to_string(shell->atom),shell->symbol+std::to_string(shell->atom)+" selected-shell subtotal",
                      shell->atom,std::nullopt,std::nullopt,node.weight+shell->weight,std::nullopt,
                      "selected_atom_and_enabled_shell_subtotal"};}
            if(v.decomposition)for(const auto& group:v.decomposition->atoms)
                if(auto it=grouped.find(group.atom);it!=grouped.end())
                    it->second.full_atom_weight=it->second.full_atom_weight.value_or(0)+group.weight;
            for(auto& [_,node]:grouped)nodes.push_back(std::move(node));
        } else for(const auto* shell:shells)
            nodes.push_back({"shell:"+shell->key,shell->symbol+std::to_string(shell->atom)+" "+shell->label+" "+shell->type,
                             shell->atom,shell->n,shell->l,shell->weight,std::nullopt,
                             state.explicitly_included_shells.contains(shell->key)?"explicit_shell_selection":
                             type_starts(shell->type,"val")?"selected_atom_default_valence":
                             type_starts(shell->type,"cor")?"selected_atom_core_enabled":
                             type_starts(shell->type,"ryd")?"selected_atom_rydberg_enabled":"explicit_shell_selection"});
    }
    const std::set<std::size_t> central(v.eligible.begin(),v.eligible.end());
    {auto out=output_csv(file(".focus.csv"));
        out<<"snapshot_id,canonical_index,source_orbital_index,spin,status,available,in_central_view,nao_id,atom,symbol,type,angular,principal_n,angular_l,coefficient,weight,electron_contribution,graph_visible,source_path,source_line\n";
        for(const auto& mo:d.mo_decompositions){
            auto prefix=[&](){out<<csv(v.snapshot_id)<<','<<mo.canonical_index<<','<<mo.source_orbital_index<<','
                <<csv(nbo_spin_name(mo.spin))<<','<<csv(mo.status)<<','<<(mo.available?"true":"false")<<','
                <<(central.contains(mo.canonical_index)?"true":"false")<<',';};
            if(mo.rows.empty()){prefix();out<<",,,,,,,,,,false,,\n";continue;}
            for(const auto& row:mo.rows){prefix();
                const bool shown=v.index && *v.index==mo.canonical_index &&
                    row.principal_n && row.angular_l && state.visible_atoms.contains(row.atom) &&
                    std::any_of(shells.begin(),shells.end(),[&](const FocusShell* shell){
                        return shell->key==focus_shell_key(row.atom,*row.principal_n,*row.angular_l,row.type);});
                out<<row.nao_id<<','<<row.atom<<','<<csv(row.symbol)<<','<<csv(row.type)<<','<<csv(row.angular)<<',';
                if(row.principal_n)out<<*row.principal_n;
                out<<',';if(row.angular_l)out<<*row.angular_l;
                out<<','<<row.coefficient<<','<<row.weight<<',';
                if(row.electron_contribution)out<<*row.electron_contribution;
                out<<','<<(shown?"true":"false")<<','<<csv(row.source.path)<<','<<row.source.line_begin<<'\n';
            }
        }
        if(!out)throw std::runtime_error("NBO focus CSV write failed");
    }
    {auto out=output_csv(file(".focus.groups.csv"));
        out<<"snapshot_id,focused_canonical_index,group_id,atom_ids,angular_l,angular_label,full_weight,selected_shell_subtotal\n";
        for(const auto& group:state.ligand_groups) {
            std::string atom_list;
            for(auto atom:group.atoms){if(!atom_list.empty())atom_list+=';';atom_list+=std::to_string(atom);}
            for(const auto& [l,w]:angular_weights(v,state,group))
                out<<csv(v.snapshot_id)<<','<<(v.index?std::to_string(*v.index):"")<<','<<group.id<<','
                   <<csv(atom_list)<<','<<l<<','<<csv(shell_letter(l))<<','<<w.full<<','<<w.selected<<'\n';
        }
        if(!out)throw std::runtime_error("NBO focus group CSV write failed");
    }
    {std::ofstream out(file(".focus.json"),std::ios::binary);
        if(!out)throw std::runtime_error("Cannot write NBO focus JSON");
        out<<std::setprecision(17);
        out<<"{\"schema\":\"nbo_focus_view_v1\",\"snapshot_id\":"<<json(v.snapshot_id)
           <<",\"selection_anchor\":";
        if(diagram&&diagram->data.view)out<<diagram->data.view->selection_anchor;else out<<"null";
        out<<",\"numerical_scope\":\"all_available_canonical_decompositions\",\"numerical_file\":\".focus.csv\""
           <<",\"status\":"<<json(v.status)<<",\"reason\":"<<json(v.reason)
           <<",\"source_output\":"<<json(d.source.path)<<",\"focused_mo_index\":";
        if(v.index)out<<*v.index;else out<<"null";
        out<<",\"grouping\":"<<json(state.group_ligands_by_l?"ligand_angular":state.group_by_atom?"atom":"shell")
           <<",\"show_core\":"<<(state.show_core?"true":"false")
           <<",\"show_rydberg\":"<<(state.show_rydberg?"true":"false")
           <<",\"central_mo_indices\":[";
        for(std::size_t i=0;i<v.eligible.size();++i){if(i)out<<',';out<<v.eligible[i];}
        out<<"],\"visible_atoms\":[";bool first=true;
        for(auto atom:state.visible_atoms){if(!first)out<<',';first=false;out<<atom;}
        out<<"],\"draft_ligand_atoms\":[";first=true;
        for(auto atom:state.ligand_atoms){if(!first)out<<',';first=false;out<<atom;}
        out<<"],\"ligand_groups\":[";
        for(std::size_t i=0;i<state.ligand_groups.size();++i) {
            const auto& group=state.ligand_groups[i];if(i)out<<',';
            out<<"{\"id\":"<<json("L"+std::to_string(group.id))<<",\"atom_ids\":[";
            bool gf=true;for(auto atom:group.atoms){if(!gf)out<<',';gf=false;out<<atom;}
            out<<"],\"angular_weights\":[";
            bool wf=true;for(const auto& [l,w]:angular_weights(v,state,group)) {
                if(!wf)out<<',';wf=false;
                out<<"{\"angular_l\":"<<l<<",\"label\":"<<json(shell_letter(l))
                   <<",\"full_weight\":"<<w.full<<",\"selected_shell_subtotal\":"<<w.selected<<'}';
            }
            out<<"]}";
        }
        out<<"],\"group_numerical_file\":\".focus.groups.csv\",\"hidden_shell_keys\":[";first=true;
        for(const auto& key:state.hidden_shells){if(!first)out<<',';first=false;out<<json(key);}
        out<<"],\"explicit_shell_keys\":[";first=true;
        for(const auto& key:state.explicitly_included_shells){if(!first)out<<',';first=false;out<<json(key);}
        out<<"],\"mo_node\":";
        if(v.index && diagram) {
            const auto it=std::find_if(diagram->data.metadata.begin(),diagram->data.metadata.end(),
                [&](const auto& item){return item.orbital_index==*v.index;});
            out<<"{\"id\":"<<json("canonical_mo:"+std::to_string(*v.index))
               <<",\"canonical_index\":"<<*v.index<<",\"energy_hartree\":";
            if(it!=diagram->data.metadata.end())out<<it->energy_hartree;else out<<"null";
            out<<",\"spin\":"<<(it!=diagram->data.metadata.end()?json(it->spin==Spin::Beta?"beta":"alpha"):"null")<<'}';
        } else out<<"null";
        out<<",\"graph_nodes\":[";
        for(std::size_t i=0;i<nodes.size();++i){const auto& node=nodes[i];if(i)out<<',';
            out<<"{\"id\":"<<json(node.id)<<",\"type\":\"composition\",\"label\":"<<json(node.label)
               <<",\"atom_ids\":[";
            if(node.atom_ids.empty())out<<node.atom;
            else for(std::size_t j=0;j<node.atom_ids.size();++j){if(j)out<<',';out<<node.atom_ids[j];}
            out<<"],\"principal_n\":";
            if(node.n)out<<*node.n;else out<<"null";
            out<<",\"angular_l\":";if(node.l)out<<*node.l;else out<<"null";
            out<<",\"raw_weight\":"<<node.weight
               <<",\"full_atom_weight\":";
            if(node.full_atom_weight)out<<*node.full_atom_weight;else out<<"null";
            out<<",\"full_angular_weight\":";
            if(node.full_angular_weight)out<<*node.full_angular_weight;else out<<"null";
            out<<",\"weight_scope\":"<<json(state.group_ligands_by_l?"selected_shell_subtotal_by_ligand_l":
                state.group_by_atom?"selected_shell_subtotal":"selected_shell")
               <<",\"visibility_reason\":"<<json(node.reason)<<'}';
        }
        out<<"],\"graph_edges\":[";
        for(std::size_t i=0;i<nodes.size();++i){if(i)out<<',';
            out<<"{\"id\":"<<json("composition:mo:"+std::to_string(v.index.value_or(0))+":"+nodes[i].id)
               <<",\"type\":\"composition\",\"source\":"
               <<json("canonical_mo:"+std::to_string(v.index.value_or(0)))
               <<",\"target\":"<<json(nodes[i].id)
               <<",\"coefficient\":null,\"raw_weight\":"<<nodes[i].weight
               <<",\"evidence\":"<<json(v.decomposition?source_label(v.decomposition->matrix_source):"")<<'}';
        }
        out<<"],\"mo_rows\":[";
        for(std::size_t i=0;i<d.mo_decompositions.size();++i){const auto& mo=d.mo_decompositions[i];if(i)out<<',';
            out<<"{\"canonical_index\":"<<mo.canonical_index<<",\"source_orbital_index\":"<<mo.source_orbital_index
               <<",\"spin\":"<<json(nbo_spin_name(mo.spin))<<",\"status\":"<<json(mo.status)
               <<",\"detail\":"<<json(mo.detail)<<",\"available\":"<<(mo.available?"true":"false")
               <<",\"in_central_view\":"<<(central.contains(mo.canonical_index)?"true":"false")
               <<",\"matrix_source\":"<<json(source_label(mo.matrix_source))<<",\"rows\":[";
            for(std::size_t j=0;j<mo.rows.size();++j){const auto& row=mo.rows[j];if(j)out<<',';
                out<<"{\"nao_id\":"<<row.nao_id<<",\"atom\":"<<row.atom
                   <<",\"symbol\":"<<json(row.symbol)<<",\"type\":"<<json(row.type)
                   <<",\"angular\":"<<json(row.angular)<<",\"principal_n\":";
                if(row.principal_n)out<<*row.principal_n;else out<<"null";
                out<<",\"angular_l\":";if(row.angular_l)out<<*row.angular_l;else out<<"null";
                out<<",\"coefficient\":"<<row.coefficient<<",\"weight\":"<<row.weight
                   <<",\"electron_contribution\":";
                if(row.electron_contribution)out<<*row.electron_contribution;else out<<"null";
                out<<",\"source\":"<<json(source_label(row.source))<<'}';
            }
            out<<"]}";
        }
        out<<"]}";
        if(!out)throw std::runtime_error("NBO focus JSON write failed");
    }
    constexpr int width=1100;
    const int height=std::max(210,150+static_cast<int>(nodes.size())*32);
    std::vector<unsigned char> rgba(static_cast<std::size_t>(width)*height*4,255);
    auto fill=[&](int x0,int y0,int x1,int y1,unsigned char r,unsigned char g,unsigned char b){
        for(int y=std::max(0,y0);y<std::min(height,y1);++y)for(int x=std::max(0,x0);x<std::min(width,x1);++x){
            const auto p=(static_cast<std::size_t>(y)*width+x)*4;rgba[p]=r;rgba[p+1]=g;rgba[p+2]=b;}};
    auto line=[&](int x0,int y0,int x1,int y1){
        const int dx=std::abs(x1-x0),dy=-std::abs(y1-y0),sx=x0<x1?1:-1,sy=y0<y1?1:-1;
        int err=dx+dy;
        for(;;){fill(x0,y0,x0+2,y0+2,74,151,206);if(x0==x1&&y0==y1)break;
            const int e2=2*err;if(e2>=dy){err+=dy;x0+=sx;}if(e2<=dx){err+=dx;y0+=sy;}}
    };
    const std::string title="CANONICAL MO / ORTHOGONAL NAO COMPOSITION";
    std::string mo_descriptor;
    if(v.index && diagram)for(const auto& meta:diagram->data.metadata)
        if(meta.orbital_index==*v.index){
            mo_descriptor=" ["+std::string(meta.spin==Spin::Beta?"beta":"alpha")+"] E="+
                fmt(meta.energy_hartree)+" Ha";break;
        }
    const std::string subtitle="SNAPSHOT "+v.snapshot_id+" | MO "+
        (v.index?std::to_string(*v.index+1):"NONE")+mo_descriptor+" | "+v.status;
    std::ostringstream svg;svg<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<width
        <<"\" height=\""<<height<<"\" viewBox=\"0 0 "<<width<<' '<<height
        <<"\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/>"
        <<"<text x=\"24\" y=\"30\" font-family=\"monospace\" font-size=\"18\">"<<escape_xml(title)<<"</text>"
        <<"<text x=\"24\" y=\"54\" font-family=\"monospace\" font-size=\"12\">"<<escape_xml(subtitle)<<"</text>";
    draw_text(rgba,width,height,24,15,title,2);draw_text(rgba,width,height,24,45,subtitle,1);
    if(nodes.empty()){
        const std::string message=v.status=="available"?"NO ATOMS OR SHELLS SELECTED":v.reason;
        draw_text(rgba,width,height,24,90,message.substr(0,100),1);
        svg<<"<text x=\"24\" y=\"105\" font-family=\"monospace\" font-size=\"12\">"
           <<escape_xml(message)<<"</text>";
    } else {
        const int mo_y=75+static_cast<int>(nodes.size())*16;
        const std::string mo_label="MO "+std::to_string(v.index.value_or(0)+1);
        fill(75,mo_y-22,150,mo_y+22,64,98,160);
        draw_text(rgba,width,height,88,mo_y-4,mo_label,1);
        svg<<"<rect x=\"75\" y=\""<<mo_y-22<<"\" width=\"75\" height=\"44\" rx=\"16\" fill=\"#4062a0\"/>"
           <<"<text x=\"88\" y=\""<<mo_y+5<<"\" font-family=\"monospace\" font-size=\"12\" fill=\"white\">"
           <<escape_xml(mo_label)<<"</text>";
        for(std::size_t i=0;i<nodes.size();++i){const auto& node=nodes[i];const int y=90+static_cast<int>(i)*32;
            line(150,mo_y,350,y);fill(345,y-5,355,y+5,57,133,197);
            const std::string value=node.label+"  WEIGHT="+fmt(node.weight)+
                (node.full_atom_weight?" / FULL ATOM="+fmt(*node.full_atom_weight):"")+
                (node.full_angular_weight?" / FULL L="+fmt(*node.full_angular_weight):"");
            draw_text(rgba,width,height,370,y-4,value.substr(0,110),1);
            svg<<"<line x1=\"150\" y1=\""<<mo_y<<"\" x2=\"350\" y2=\""<<y
               <<"\" stroke=\"#4a97ce\" stroke-width=\"2\"/><circle cx=\"350\" cy=\""<<y
               <<"\" r=\"5\" fill=\"#3985c5\"/><text x=\"370\" y=\""<<y+5
               <<"\" font-family=\"monospace\" font-size=\"12\">"<<escape_xml(value)<<"</text>";
        }
    }
    const std::string footer="NAO WEIGHTS ARE NOT RAW AO COEFFICIENT SQUARES; E2 IS SEPARATE";
    draw_text(rgba,width,height,24,height-22,footer,1);
    svg<<"<text x=\"24\" y=\""<<height-10<<"\" font-family=\"monospace\" font-size=\"10\">"
       <<footer<<"</text></svg>";
    {std::ofstream out(file(".focus.svg"),std::ios::binary);if(!out)throw std::runtime_error("Cannot write NBO focus SVG");
        out<<svg.str();if(!out)throw std::runtime_error("NBO focus SVG write failed");}
    write_png(file(".focus.png"),rgba,width,height);
}

void export_nbo_bundle(const NboDataset& d,std::size_t selected,const std::filesystem::path& base,
                       const Wavefunction* canonical,std::size_t canonical_index,
                       double contribution_threshold,bool nbo_active,
                       std::size_t rendered_orbital_index,
                       const MODiagramViewSnapshot* diagram,
                       const NboFocusUIState* focus){
    (void)contribution_threshold;
    auto path=base;path.replace_extension();
    if(!path.parent_path().empty())std::filesystem::create_directories(path.parent_path());
    const auto file=[&](const char* suffix){auto p=path;p+=suffix;return p;};
    {std::ofstream out(file(".nbo.json"),std::ios::binary);if(!out)throw std::runtime_error("Cannot write NBO JSON");out<<serialize_nbo_json(d);if(!out)throw std::runtime_error("NBO JSON write failed");}
    {auto out=output_csv(file(".npa.csv"));out<<"atom,symbol,spin,charge,core,valence,rydberg,printed_total_population,effective_core_electrons,explicit_population,spin_density,source_path,source_line,analysis_segment\n";
     for(const auto& p:d.populations)out<<p.atom<<','<<csv(p.symbol)<<','<<csv(nbo_spin_name(p.spin))<<','<<p.charge<<','<<p.core<<','<<p.valence<<','<<p.rydberg<<','<<p.total<<','<<opt(p.effective_core_electrons)<<','<<opt(p.explicit_population)<<','<<opt(p.spin_density)<<','<<csv(p.source.path)<<','<<p.source.line_begin<<','<<p.source.analysis_segment<<'\n';}
    {auto out=output_csv(file(".nao.csv"));out<<"id,atom,symbol,spin,angular,type,occupation,diagonal_fock_hartree,spin_density,source_path,source_line,analysis_segment\n";
     for(const auto& n:d.naos)out<<n.id<<','<<n.atom<<','<<csv(n.symbol)<<','<<csv(nbo_spin_name(n.spin))<<','<<csv(n.angular)<<','<<csv(n.type)<<','<<n.occupation<<','<<opt(n.energy_hartree)<<','<<opt(n.spin_density)<<','<<csv(n.source.path)<<','<<n.source.line_begin<<','<<n.source.analysis_segment<<'\n';}
    {auto out=output_csv(file(".nbo.csv"));out<<"index,id,ordinal,spin,kind,label,occupation,diagonal_fock_hartree,component_atom,component_percent,component_coefficient,component_hybrid,energy_source_path,energy_source_line,source_path,source_line,analysis_segment\n";
     for(std::size_t i=0;i<d.orbitals.size();++i){const auto& o=d.orbitals[i];const auto write=[&](const NboLocalComponent* c){out<<i<<','<<o.id<<','<<o.ordinal<<','<<csv(nbo_spin_name(o.spin))<<','<<csv(o.kind)<<','<<csv(o.label)<<','<<o.occupation<<','<<opt(o.diagonal_fock_hartree)<<',';if(c)out<<c->atom<<','<<c->percent<<','<<c->coefficient<<','<<csv(c->hybrid);else out<<",,,";out<<','<<(o.energy_source?csv(o.energy_source->path):"")<<','<<(o.energy_source?std::to_string(o.energy_source->line_begin):"")<<','<<csv(o.source.path)<<','<<o.source.line_begin<<','<<o.source.analysis_segment<<'\n';};if(o.components.empty())write(nullptr);else for(const auto& c:o.components)write(&c);}}
    {auto out=output_csv(file(".wiberg.csv"));out<<"atom_a,atom_b,spin,value,source_path,source_line,analysis_segment\n";for(const auto& w:d.wiberg)out<<w.atom_a<<','<<w.atom_b<<','<<csv(nbo_spin_name(w.spin))<<','<<w.value<<','<<csv(w.source.path)<<','<<w.source.line_begin<<','<<w.source.analysis_segment<<'\n';}
    {auto out=output_csv(file(".e2.csv"));out<<"donor,acceptor,spin,value,units,printing_threshold,energy_gap_hartree,fock_hartree,source_path,source_line,analysis_segment\n";for(const auto& e:d.e2)out<<e.donor<<','<<e.acceptor<<','<<csv(nbo_spin_name(e.spin))<<','<<e.value<<','<<csv(e.units)<<','<<opt(e.printing_threshold)<<','<<e.energy_gap_hartree<<','<<e.fock_hartree<<','<<csv(e.source.path)<<','<<e.source.line_begin<<','<<e.source.analysis_segment<<'\n';}
    {auto out=output_csv(file(".e2-sections.csv"));out<<"spin,printing_threshold,units,missing_reason,source_path,source_line,analysis_segment\n";for(const auto& e:d.e2_sections)out<<csv(nbo_spin_name(e.spin))<<','<<opt(e.printing_threshold)<<','<<csv(e.units)<<','<<csv(e.missing_reason)<<','<<csv(e.source.path)<<','<<e.source.line_begin<<','<<e.source.analysis_segment<<'\n';}
    {
        std::ofstream out(file(".view.json"),std::ios::binary);
        if(!out)throw std::runtime_error("Cannot write NBO view snapshot");
        out<<std::setprecision(17);
        out<<"{\"schema\":1,\"export_dataset\":\"nbo\",\"rendered_set\":"<<json(nbo_active?"nbo":"canonical")
           <<",\"rendered_orbital_index\":"<<rendered_orbital_index
           <<",\"dataset_source\":"<<json(d.source.path)
           <<",\"association\":"<<json(d.association.status)
           <<",\"selected_nbo_index\":";
        if(selected<d.orbitals.size())out<<selected;else out<<"null";
        out<<",\"selected_nbo_id\":";
        if(selected<d.orbitals.size())out<<d.orbitals[selected].id;else out<<"null";
        out<<",\"selected_nbo_spin\":"<<(selected<d.orbitals.size()?json(nbo_spin_name(d.orbitals[selected].spin)):"null")
           <<",\"selected_nbo_occupation_electrons\":";
        if(selected<d.orbitals.size())out<<d.orbitals[selected].occupation;else out<<"null";
        out<<",\"selected_nbo_diagonal_fock_hartree\":";
        if(selected<d.orbitals.size() && d.orbitals[selected].diagonal_fock_hartree)
            out<<*d.orbitals[selected].diagonal_fock_hartree;else out<<"null";
        out<<",\"selected_nbo_fock_source\":";
        if(selected<d.orbitals.size() && d.orbitals[selected].energy_source)
            out<<json(d.orbitals[selected].energy_source->path);else out<<"null";
        out<<",\"occupation_units\":\"electrons\",\"energy_semantics\":\"NBO diagonal Fock is not canonical MO energy\""
           <<",\"canonical_mo\":";
        const bool canonical_has_beta=canonical && std::any_of(canonical->orbitals.begin(),canonical->orbitals.end(),
            [](const MolecularOrbital& mo){return mo.spin==Spin::Beta;});
        const NboMatrix* matrix=nullptr;
        std::size_t source_index=std::numeric_limits<std::size_t>::max();
        NboSpin spin=NboSpin::Total;
        if(canonical && canonical_index<canonical->orbitals.size()){
            const auto& mo=canonical->orbitals[canonical_index];
            source_index=mo.source_orbital_index;spin=mo.spin==Spin::Beta?NboSpin::Beta:NboSpin::Alpha;
            out<<"{\"list_index\":"<<canonical_index<<",\"source_index\":";
            if(source_index==std::numeric_limits<std::size_t>::max())out<<"null";else out<<source_index;
            out<<",\"spin\":"<<json(nbo_spin_name(spin))<<'}';
        }else out<<"null";
        out<<",\"canonical_beta_mo_available\":"<<(canonical_has_beta?"true":"false")
           <<",\"selected_nbo_provenance\":"
           <<json(selected<d.orbitals.size() && d.orbitals[selected].spin==NboSpin::Beta && !canonical_has_beta
                   ? "associated archive and NBO sidecars; FCHK beta MO identity unavailable"
                   : "NBO producer orbital and associated canonical dataset")
           <<",\"selected_nbo_coefficient_evidence\":";
        const NboCanonicalEvidence* selected_evidence=nullptr;
        if(selected<d.orbitals.size())for(const auto& e:d.association.canonical_evidence)
            if(e.spin==d.orbitals[selected].spin){selected_evidence=&e;break;}
        if(selected_evidence){
            out<<"{\"source\":"<<json(selected_evidence->coefficient_source)
               <<",\"direct_fchk_coefficients\":"<<(selected_evidence->direct_fchk_coefficients?"true":"false")
               <<",\"density_verified\":"<<(selected_evidence->density_verified?"true":"false")
               <<",\"detail\":"<<json(selected_evidence->detail)<<'}';
        }else out<<"null";
        out
           <<",\"nbomo_threshold\":null";
        auto is_nbomo=[](const NboMatrix& m){std::string kind=m.kind;std::transform(kind.begin(),kind.end(),kind.begin(),[](unsigned char c){return static_cast<char>(std::toupper(c));});return kind=="NBOMO";};
        bool direct=true;
        for(const auto& e:d.association.canonical_evidence)
            if(e.spin==spin || (spin==NboSpin::Alpha && e.spin==NboSpin::Total))
                direct=e.direct_fchk_coefficients;
        if(d.association.compatible && direct){
            for(const auto& m:d.matrices)if(is_nbomo(m)&&m.spin==spin){matrix=&m;break;}
            if(!matrix&&spin==NboSpin::Alpha)for(const auto& m:d.matrices)if(is_nbomo(m)&&m.spin==NboSpin::Total){matrix=&m;break;}
            if(!matrix&&d.archive){for(const auto& m:d.archive->matrices)if(is_nbomo(m)&&m.spin==spin){matrix=&m;break;}
                if(!matrix&&spin==NboSpin::Alpha)for(const auto& m:d.archive->matrices)if(is_nbomo(m)&&m.spin==NboSpin::Total){matrix=&m;break;}}
        }
        std::vector<std::size_t> row_to_orbital;
        if(matrix)for(std::size_t i=0;i<d.orbitals.size();++i)if(d.orbitals[i].spin==matrix->spin)row_to_orbital.push_back(i);
        const bool mapped=matrix&&source_index<matrix->columns&&row_to_orbital.size()==matrix->rows&&matrix->values.size()==matrix->rows*matrix->columns;
        out<<",\"canonical_direct_fchk_coefficients\":"<<(direct?"true":"false")
           <<",\"nbomo_status\":"<<json(mapped?"mapped":!direct?"archive_canonical_only":d.association.compatible?"missing_or_mismatch":"association_unverified")
           <<",\"nbomo_source\":"<<(mapped?json(matrix->source.path):"null")<<",\"contributions\":[";
        if(mapped){
            struct Entry{std::size_t orbital,row;double weight;};std::vector<Entry> entries;
            double total=0,shown=0;
            for(std::size_t row=0;row<matrix->rows;++row){const double t=matrix->values[row*matrix->columns+source_index];const double w=t*t;entries.push_back({row_to_orbital[row],row,w});total+=w;}
            std::stable_sort(entries.begin(),entries.end(),[](const auto& a,const auto& b){return a.weight>b.weight;});
            for(std::size_t i=0;i<entries.size();++i){const auto& e=entries[i];const bool visible=true;if(visible)shown+=e.weight;
                if(i)out<<',';out<<"{\"nbo_index\":"<<e.orbital<<",\"nbo_id\":"<<d.orbitals[e.orbital].id<<",\"matrix_row\":"<<e.row<<",\"weight_raw\":"<<e.weight<<",\"displayed\":"<<(visible?"true":"false")<<'}';}
            out<<"],\"total_weight_raw\":"<<total<<",\"shown_weight_raw\":"<<shown<<",\"remaining_weight_raw\":"<<total-shown;
        }else out<<"],\"total_weight_raw\":null,\"shown_weight_raw\":null,\"remaining_weight_raw\":null";
        out<<'}';if(!out)throw std::runtime_error("NBO view snapshot write failed");
    }
    // A fixed occupation view is generated from the same dataset object as JSON/CSV.
    // It intentionally has no canonical MO energy axis.
    constexpr int width=960,height=420;
    std::vector<unsigned char> rgba(static_cast<std::size_t>(width)*height*4,255);
    const auto fill=[&](int x0,int y0,int x1,int y1,unsigned char r,unsigned char g,unsigned char b){for(int y=std::max(y0,0);y<std::min(y1,height);++y)for(int x=std::max(x0,0);x<std::min(x1,width);++x){const auto p=(static_cast<std::size_t>(y)*width+x)*4;rgba[p]=r;rgba[p+1]=g;rgba[p+2]=b;}};
    const std::size_t start=d.orbitals.empty()?0:std::min(selected,d.orbitals.size()-1);
    const std::size_t count=std::min<std::size_t>(12,d.orbitals.size()-start);
    const std::string title="NBO OCCUPATION (ELECTRONS)";
    const std::string scene="SCENE: "+std::string(nbo_active?"NBO":"CANONICAL")+
        " ORBITAL INDEX "+std::to_string(rendered_orbital_index);
    const std::string source="SOURCE: "+source_label(d.source);
    const std::string selected_text=d.orbitals.empty()?"SELECTED: NONE":
        "SELECTED: NBO "+std::to_string(d.orbitals[start].id)+" "+nbo_spin_name(d.orbitals[start].spin)+
        " OCC="+fmt(d.orbitals[start].occupation)+" ELECTRONS";
    draw_text(rgba,width,height,28,12,title,2);
    draw_text(rgba,width,height,28,38,source.substr(0,145),1);
    draw_text(rgba,width,height,28,52,scene,1);
    draw_text(rgba,width,height,28,66,selected_text,1);
    std::ostringstream svg;svg<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"960\" height=\"420\" viewBox=\"0 0 960 420\"><rect width=\"960\" height=\"420\" fill=\"white\"/><text x=\"28\" y=\"28\" font-family=\"monospace\" font-size=\"16\">"<<escape_xml(title)<<"</text><text x=\"28\" y=\"48\" font-family=\"monospace\" font-size=\"10\">"<<escape_xml(source.substr(0,145))<<"</text><text x=\"28\" y=\"62\" font-family=\"monospace\" font-size=\"10\">"<<escape_xml(scene)<<"</text><text x=\"28\" y=\"76\" font-family=\"monospace\" font-size=\"10\">"<<escape_xml(selected_text)<<"</text>";
    for(std::size_t j=0;j<count;++j){
        const auto& o=d.orbitals[start+j];const int y=90+static_cast<int>(j)*24;
        const int bar=std::clamp(static_cast<int>(std::round(o.occupation/2.0*420)),0,420);
        const std::string label="NBO "+std::to_string(o.id)+" "+nbo_spin_name(o.spin);
        const std::string value=fmt(o.occupation)+" E";
        fill(430,y,430+bar,y+17,45,105,175);
        draw_text(rgba,width,height,28,y+2,label,1);
        draw_text(rgba,width,height,440+bar,y+2,value,1);
        svg<<"<text x=\"28\" y=\""<<y+12<<"\" font-family=\"monospace\" font-size=\"10\">"<<escape_xml(label)<<"</text><rect x=\"430\" y=\""<<y<<"\" width=\""<<bar<<"\" height=\"17\" fill=\"#2d69af\"/><text x=\""<<440+bar<<"\" y=\""<<y+12<<"\" font-family=\"monospace\" font-size=\"10\">"<<escape_xml(value)<<"</text>";
    }
    const std::string footer="FOCK DIAGONAL IS NOT CANONICAL MO ENERGY";
    draw_text(rgba,width,height,28,401,footer,1);
    svg<<"<text x=\"28\" y=\"408\" font-family=\"monospace\" font-size=\"10\">"<<escape_xml(footer)<<"</text></svg>";
    {std::ofstream out(file(".view.svg"),std::ios::binary);if(!out)throw std::runtime_error("Cannot write NBO SVG");out<<svg.str();if(!out)throw std::runtime_error("NBO SVG write failed");}
    write_png(file(".view.png"),rgba,width,height);
    export_nbo_focus_bundle(d,path,diagram,focus);
}
} // namespace cov::ui
