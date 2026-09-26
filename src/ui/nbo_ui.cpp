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
#include <stdexcept>
#include <vector>

namespace cov::ui {
namespace {
enum Label { Title, Path, Archive, Aonbo, Nbomo, Segment, Attach, Canonical,
             NboSet, Association, Renderable, ReportsOnly, Populations, Naos,
             Orbitals, Wiberg, E2, Matrices, EnergyNote, ExportPath, Export,
             Source, Threshold, Units, NoData, Count };
constexpr const char* labels[4][Count] = {
    {"NBO analysis", "NBO output path", "ARCHIVE .47 (optional)", "AONBO .37 (optional)", "NBOMO .49 (optional)", "Analysis segment (-1: require unique)", "Attach NBO", "Canonical MOs", "NBO orbitals", "Source association", "NBO coefficients verified for rendering", "Report only: coefficients or source association unavailable", "NPA populations", "NAOs", "NBO orbital occupations and components", "Wiberg indices", "Second-order interactions", "Matrix evidence", "Diagonal Fock values are not canonical MO energies.", "Export base path", "Export NBO bundle", "Source", "Print threshold", "Units", "No data reported"},
    {"NBO 分析", "NBO 输出路径", "ARCHIVE .47（可选）", "AONBO .37（可选）", "NBOMO .49（可选）", "分析段（-1：要求唯一）", "关联 NBO", "正则 MO", "NBO 轨道", "同源关联", "NBO 系数已验证，可渲染", "仅报告：缺少系数或同源关联", "NPA 布居", "NAO", "NBO 轨道占据与成分", "Wiberg 指数", "二阶相互作用", "矩阵证据", "Fock 对角值不是正则 MO 能量。", "导出路径前缀", "导出 NBO 数据", "来源", "打印阈值", "单位", "未报告数据"},
    {"NBO 解析", "NBO 出力パス", "ARCHIVE .47（任意）", "AONBO .37（任意）", "NBOMO .49（任意）", "解析区間（-1：一意のみ）", "NBO を関連付け", "正準 MO", "NBO 軌道", "同一源の照合", "描画用 NBO 係数を検証済み", "報告のみ：係数または照合が不足", "NPA 原子分布", "NAO", "NBO 軌道の占有数と成分", "Wiberg 指数", "二次相互作用", "行列の証拠", "Fock 対角値は正準 MO エネルギーではありません。", "出力先の基底名", "NBO 一式を書き出す", "出典", "印字しきい値", "単位", "報告データなし"},
    {"Analyse NBO", "Chemin de sortie NBO", "ARCHIVE .47 (facultatif)", "AONBO .37 (facultatif)", "NBOMO .49 (facultatif)", "Segment d’analyse (-1 : unique)", "Associer NBO", "OM canoniques", "Orbitales NBO", "Association des sources", "Coefficients NBO vérifiés pour l’affichage", "Rapport seul : coefficients ou association absents", "Populations NPA", "NAO", "Occupation et composants NBO", "Indices de Wiberg", "Interactions du second ordre", "Preuves matricielles", "Les valeurs diagonales de Fock ne sont pas des énergies OM canoniques.", "Chemin de base de l’export", "Exporter NBO", "Source", "Seuil d’impression", "Unités", "Aucune donnée rapportée"}
};
const char* trn(Language language, Label key) {
    const int row=std::clamp(static_cast<int>(language),0,3);
    return labels[row][key];
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
    return out;
}

NboUIActions draw_nbo_panel(NboUIState& s, Language language,
                            bool canonical_loaded, bool renderable,
                            bool active, float scale,
                            const Wavefunction* canonical, std::size_t canonical_index) {
    NboUIActions actions;
    const auto words=report_text(language);
    ImGui::Spacing();
    section_title(trn(language,Title));
    input_path(trn(language,Path),"##nbo.path",s.path);
    input_path(trn(language,Archive),"##nbo.archive47",s.archive47);
    input_path(trn(language,Aonbo),"##nbo.aonbo",s.aonbo);
    input_path(trn(language,Nbomo),"##nbo.nbomo",s.nbomo);
    ImGui::TextWrapped("%s",trn(language,Segment));
    ImGui::SetNextItemWidth(100.0f*scale);
    ImGui::InputInt("##nbo.segment",&s.analysis_segment);
    validation::item("nbo.segment");
    ImGui::BeginDisabled(!canonical_loaded || active);
    actions.attach=ImGui::Button(trn(language,Attach),ImVec2(-1.0f,0));
    ImGui::EndDisabled();validation::item("nbo.attach");
    if(!s.error.empty())note(s.error);
    if(!s.dataset)return actions;
    const auto& d=*s.dataset;
    note(std::string(trn(language,Association))+": "+d.association.status+" — "+d.association.detail);
    validation::field("nbo.association",d.association.status+": "+d.association.detail);
    for(const auto& evidence:d.association.canonical_evidence){
        note(std::string(words.canonical_evidence)+" ["+nbo_spin_name(evidence.spin)+"] "+
             evidence.coefficient_source+"; "+words.direct_fchk+"="+
             (evidence.direct_fchk_coefficients?"1":"0")+"; "+words.density_verified+"="+
             (evidence.density_verified?"1":"0"));
    }
    note(renderable?trn(language,Renderable):trn(language,ReportsOnly));
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
    if(ImGui::CollapsingHeader(words.contributions,ImGuiTreeNodeFlags_DefaultOpen)){
        ImGui::TextWrapped("%s",words.threshold);
        ImGui::SetNextItemWidth(130.0f*scale);
        ImGui::InputDouble("##nbo.nbomo.threshold",&s.contribution_threshold,0,0,"%.4f");
        if(!std::isfinite(s.contribution_threshold))s.contribution_threshold=0.01;
        s.contribution_threshold=std::clamp(s.contribution_threshold,0.0,1.0);
        validation::item("nbo.nbomo.threshold");
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
                    std::size_t displayed=0;
                    for(const auto& entry:entries){
                        if(entry.weight<s.contribution_threshold || displayed>=16)break;
                        const auto& orbital=d.orbitals[entry.orbital];
                        const std::string label="#"+std::to_string(orbital.id)+
                            " ["+nbo_spin_name(orbital.spin)+"] |T|²="+fmt(entry.weight);
                        ImGui::PushID(static_cast<int>(entry.orbital));
                        if(ImGui::Selectable(label.c_str(),s.selected_orbital==entry.orbital))
                            actions.selected_orbital=entry.orbital;
                        validation::item("nbo.contribution."+std::to_string(entry.orbital));
                        if(ImGui::IsItemHovered())ImGui::SetTooltip("%s",orbital.label.c_str());
                        ImGui::PopID();
                        shown+=entry.weight;++displayed;
                    }
                    note(std::string(words.shown)+": "+fmt(shown)+"; "+words.remainder+": "+fmt(total-shown)+
                         "; Σ|T|²="+fmt(total));
                    note(std::string(words.source_evidence)+": "+source_label(matrix->source));
                    validation::field("nbo.nbomo", "spin="+std::string(nbo_spin_name(matrix->spin))+
                        ";canonical_source_index="+std::to_string(source_index)+
                        ";threshold="+fmt(s.contribution_threshold)+
                        ";shown="+fmt(shown)+";remaining="+fmt(total-shown));
                }
            }
        }
    }
    if(ImGui::CollapsingHeader(trn(language,Populations),ImGuiTreeNodeFlags_DefaultOpen)){
        if(d.populations.empty())note(trn(language,NoData));
        for(const auto& p:d.populations){
            note(std::string(words.atom)+" "+std::to_string(p.atom)+" "+p.symbol+" ["+nbo_spin_name(p.spin)+"]  "+words.charge+"="+fmt(p.charge)+" "+words.core+"="+fmt(p.core)+" "+words.valence+"="+fmt(p.valence)+" "+words.rydberg+"="+fmt(p.rydberg)+" "+words.total+"="+fmt(p.total)+" "+words.effective_core+"="+shown_opt(p.effective_core_electrons)+" "+words.explicit_pop+"="+shown_opt(p.explicit_population)+" "+words.spin_density+"="+shown_opt(p.spin_density));
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

void export_nbo_bundle(const NboDataset& d,std::size_t selected,const std::filesystem::path& base,
                       const Wavefunction* canonical,std::size_t canonical_index,
                       double contribution_threshold,bool nbo_active,
                       std::size_t rendered_orbital_index){
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
           <<",\"nbomo_threshold\":"<<std::clamp(contribution_threshold,0.0,1.0);
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
            std::size_t displayed=0;
            for(std::size_t i=0;i<entries.size();++i){const auto& e=entries[i];const bool visible=e.weight>=contribution_threshold&&displayed<16;if(visible){shown+=e.weight;++displayed;}
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
}
} // namespace cov::ui
