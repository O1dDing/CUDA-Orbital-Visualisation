#include "cov/nbo_aomo_ui.hpp"
#include "cov/molecule_style.hpp"
#include "cov/validation.hpp"
#include <imgui.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace cov::ui {
namespace {
std::string quote(const std::string& value) {
    std::string out="\"";
    for(unsigned char c:value) {
        switch(c) {case '"':out+="\\\"";break;case '\\':out+="\\\\";break;
            case '\n':out+="\\n";break;case '\r':out+="\\r";break;case '\t':out+="\\t";break;
            default:if(c<32){const char* hex="0123456789abcdef";out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}
                else out+=static_cast<char>(c);}
    }
    return out+'"';
}
std::string csv(const std::string& value) {
    std::string out="\"";for(char c:value){if(c=='"')out+='"';out+=c;}return out+'"';
}
std::string xml(const std::string& value) {
    std::string out;
    for(char c:value) switch(c) {
        case '&':out+="&amp;";break;case '<':out+="&lt;";break;
        case '>':out+="&gt;";break;case '"':out+="&quot;";break;
        default:out+=c;
    }
    return out;
}
std::array<unsigned char,7> raster_glyph(char c) {
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
        case '(':return {2,4,8,8,8,4,2};case ')':return {8,4,2,2,2,4,8};
        case '[':return {14,8,8,8,8,8,14};case ']':return {14,2,2,2,2,2,14};
        default:return {0,0,0,0,0,0,0};
    }
}
std::string number(double v) {std::ostringstream out;out<<std::setprecision(9)<<v;return out.str();}
using RefKey=std::tuple<int,int,std::size_t>;
RefKey key(const NboOrbitalRef& r){return {static_cast<int>(r.kind),static_cast<int>(r.spin),r.index};}
bool in(const std::vector<std::size_t>& values,std::size_t n) {
    return std::find(values.begin(),values.end(),n)!=values.end();
}
bool contains_nocase(std::string value,const std::string& needle) {
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char c){return static_cast<char>(std::tolower(c));});
    return value.find(needle)!=std::string::npos;
}
float node_width(const NboAomoNode& node) {
    return node.width;
}
std::vector<std::size_t> central_indices(const MODiagramViewSnapshot& diagram) {
    std::vector<std::size_t> indices;
    for(const auto& level:diagram.data.levels) {
        if(level.member_indices.empty())indices.push_back(level.metadata.orbital_index);
        else indices.insert(indices.end(),level.member_indices.begin(),level.member_indices.end());
        indices.insert(indices.end(),level.member_spin_counterparts.begin(),level.member_spin_counterparts.end());
    }
    std::sort(indices.begin(),indices.end());indices.erase(std::unique(indices.begin(),indices.end()),indices.end());
    return indices;
}
void suggest_terminal_fragment_groups(NboAomoUIState& state,
    const NboIntegration& data,const Wavefunction& canonical,
    const std::vector<std::size_t>& central) {
    state.suggested_fragments_initialized=true;
    const auto bonds=analyse_bonds(canonical);
    std::vector<std::vector<std::size_t>> neighbors(canonical.atoms.size());
    for(const auto& bond:bonds) {
        if(bond.atom_a>=neighbors.size()||bond.atom_b>=neighbors.size())continue;
        neighbors[bond.atom_a].push_back(bond.atom_b);
        neighbors[bond.atom_b].push_back(bond.atom_a);
    }
    std::map<std::pair<std::size_t,int>,std::set<std::size_t>> siblings;
    for(std::size_t atom=0;atom<neighbors.size();++atom)
        if(neighbors[atom].size()==1)
            siblings[{neighbors[atom].front(),canonical.atoms[atom].atomic_number}].insert(atom);
    for(const auto& [_,atoms]:siblings) {
        if(atoms.size()<2)continue;
        const std::vector<std::size_t> members(atoms.begin(),atoms.end());
        const bool applicable=std::any_of(central.begin(),central.end(),[&](std::size_t mo){
            return !nbo_fragment_selection(data,mo,members,NboOrbitalKind::NAO,false).terms.empty();});
        if(applicable)state.fragment_groups.push_back({state.next_fragment_id++,atoms,true});
    }
}
std::string atom_label(const Wavefunction& canonical,std::size_t atom) {
    if(atom>=canonical.atoms.size())return "Unassigned";
    if(!canonical.atoms[atom].symbol.empty())
        return canonical.atoms[atom].symbol+std::to_string(atom+1);
    static constexpr const char* symbols[]={"?","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar"};
    const auto z=canonical.atoms[atom].atomic_number;
    return std::string(z>0 && z<static_cast<int>(std::size(symbols))?symbols[z]:"Z")+
        (z>=static_cast<int>(std::size(symbols))?std::to_string(z):"")+std::to_string(atom+1);
}
std::optional<NboOrbitalRef> canonical_ref(const NboIntegration& data,std::size_t index) {
    for(const auto& orbital:data.orbitals)
        if(orbital.ref.kind==NboOrbitalKind::Canonical && orbital.ref.index==index)
            return orbital.ref;
    return std::nullopt;
}
std::string row_id(const NboOrbitalRef& ref) {
    return std::string(nbo_orbital_kind_name(ref.kind))+":"+nbo_spin_name(ref.spin)+":"+std::to_string(ref.index);
}
std::string source_name(const NboSource& s) {
    return s.path+(s.line_begin?":"+std::to_string(s.line_begin):"")+(s.block.empty()?"":" ["+s.block+"]");
}
const char* lt(Language l,const char* en,const char* zh,const char* ja,const char* fr) {
    switch(l){case Language::ChineseSimplified:return zh;case Language::Japanese:return ja;
        case Language::French:return fr;default:return en;}
}
std::vector<std::size_t> level_members(const MODiagramLevel& level) {
    std::vector<std::size_t> members=level.member_indices;
    if(members.empty())members.push_back(level.metadata.orbital_index);
    members.insert(members.end(),level.member_spin_counterparts.begin(),level.member_spin_counterparts.end());
    std::sort(members.begin(),members.end());members.erase(std::unique(members.begin(),members.end()),members.end());
    return members;
}
NboAomoViewSnapshot make_snapshot(const NboAomoUIState& state,const NboIntegration& data,
                                  const Wavefunction& canonical,const MODiagramViewSnapshot& diagram,
                                  float plot_width) {
    NboAomoViewSnapshot view;
    view.integration_id=data.id;
    view.mo_snapshot_id=diagram.data.view?diagram.data.view->id:"no-mo-snapshot";
    view.id=view.mo_snapshot_id+":aomo:"+std::to_string(state.revision);
    view.basis_kind=state.basis_kind;
    view.mo_energy_axis_mode=energy_axis_mode_name(diagram.options.energy_axis_mode);
    view.mo_energy_axis_detail="Only the canonical MO column uses this "+
        std::string(view.mo_energy_axis_mode)+" energy transform in hartree; AO/NAO rows are not placed on that energy axis";
    view.zoom=state.zoom;view.pan_x=state.pan_x;view.pan_y=state.pan_y;
    view.show_core=state.show_core;view.show_rydberg=state.show_rydberg;
    view.selection=state.selection;view.fragment_groups=state.fragment_groups;
    view.sum_component_ids.assign(state.sum_component_ids.begin(),state.sum_component_ids.end());
    if(const auto* cap=nbo_capability(data,"aomo")){view.capability_status=nbo_capability_state_name(cap->state);view.capability_detail=cap->detail;}
    view.central_mo_indices=central_indices(diagram);
    const auto inspected=diagram.data.view?diagram.data.view->inspected_orbital_index:std::nullopt;
    view.focused_canonical_index=state.focused_canonical_index && in(view.central_mo_indices,*state.focused_canonical_index)
        ?*state.focused_canonical_index:inspected && in(view.central_mo_indices,*inspected)
            ?*inspected:view.central_mo_indices.empty()?0:view.central_mo_indices.front();
    if(const auto* decomposition=nbo_mo_decomposition(data.dataset,view.focused_canonical_index)) {
        view.focused_projection_weight=decomposition->weight_sum;
        view.focused_projection_residual_norm=decomposition->projection_residual_norm;
    }
    std::map<std::size_t,std::vector<const NboOrbitalDescriptor*>> atoms;
    for(const auto& orbital:data.orbitals) {
        if(orbital.ref.kind!=state.basis_kind)continue;
        const bool related=std::any_of(data.links.begin(),data.links.end(),[&](const auto& link){
            return link.orbital==orbital.ref && in(view.central_mo_indices,link.canonical_index);});
        if(!related)continue;
        if(state.basis_kind==NboOrbitalKind::NAO &&
           ((contains_nocase(orbital.label," cor")&&!state.show_core) ||
            (contains_nocase(orbital.label," ryd")&&!state.show_rydberg))) {
            ++view.hidden_class_count;continue;
        }
        atoms[orbital.atoms.empty()?canonical.atoms.size():orbital.atoms.front()].push_back(&orbital);
    }
    std::map<RefKey,std::size_t> basis_nodes;
    const float row_step=std::max(39.0f,ImGui::GetTextLineHeight()+13.0f);
    const float box_height=row_step-5.0f;
    float y=38;
    for(auto& [atom,orbitals]:atoms) {
        std::sort(orbitals.begin(),orbitals.end(),[](const auto* a,const auto* b){return a->ref.index<b->ref.index;});
        NboAomoNode header;header.id="atom:"+std::to_string(atom);header.label=atom_label(canonical,atom);
        header.detail="Group header; expand to inspect actual orbitals";header.atoms={atom};
        header.x=12;header.y=y;header.height=box_height;header.group_header=true;
        view.nodes.push_back(std::move(header));y+=row_step;
        if(state.collapsed_atoms.contains(atom)) {view.hidden_basis_count+=orbitals.size();continue;}
        for(const auto* orbital:orbitals) {
            NboAomoNode node;node.id=orbital->id.empty()?row_id(orbital->ref):orbital->id;
            node.label=orbital->label.empty()?row_id(orbital->ref):orbital->label;
            node.detail=orbital->label+"; "+orbital->detail;
            node.energy_semantics=orbital->energy_semantics;
            if(orbital->ref.kind==NboOrbitalKind::NAO) {
                const auto it=std::find_if(data.dataset.naos.begin(),data.dataset.naos.end(),
                    [&](const auto& row){return row.id==orbital->ref.index+1 &&
                        row.spin==orbital->ref.spin;});
                if(it!=data.dataset.naos.end())
                    node.label=it->symbol+std::to_string(it->atom)+" "+
                        it->angular+" / NAO "+std::to_string(it->id);
            } else if(orbital->ref.kind==NboOrbitalKind::GaussianAO) {
                node.label=atom_label(canonical,orbital->atoms.empty()?canonical.atoms.size():
                    orbital->atoms.front())+" / AO "+std::to_string(orbital->ref.index+1);
            }
            if(orbital->ref.spin!=NboSpin::Total)
                node.label+=" ["+std::string(nbo_spin_name(orbital->ref.spin))+"]";
            node.orbital=orbital->ref;node.atoms=orbital->atoms;
            node.energy_hartree=orbital->energy_hartree;node.occupation=orbital->occupation;
            node.x=24;node.y=y;node.height=box_height;
            basis_nodes[key(orbital->ref)]=view.nodes.size();
            view.nodes.push_back(std::move(node));y+=row_step;
        }
        y+=8;
    }
    std::map<std::size_t,std::size_t> mo_nodes;
    std::map<std::size_t,std::vector<std::size_t>> level_nodes;
    y=38;
    for(std::size_t level_index=0;level_index<diagram.data.levels.size();++level_index) {
        const auto& level=diagram.data.levels[level_index];
        const auto members=level_members(level);
        const bool grouped=members.size()>1;
        if(grouped) {
            NboAomoNode header;header.id="degenerate:"+std::to_string(level_index);
            header.label="Degenerate set ("+std::to_string(members.size())+")";
            header.detail="Group, not a single physical orbital; expand or inspect all members";
            header.x=std::max(350.0f,plot_width-195.0f);header.y=y;
            header.height=box_height;header.group_header=true;
            level_nodes[level_index].push_back(view.nodes.size());
            view.nodes.push_back(std::move(header));y+=row_step;
            if(state.collapsed_levels.contains(level_index)){view.hidden_mo_count+=members.size();continue;}
        }
        for(auto index:members) {
            if(index>=canonical.orbitals.size())continue;
            const auto& mo=canonical.orbitals[index];
            NboAomoNode node;node.id="canonical_mo:"+std::to_string(index);
            node.label="MO "+std::to_string(index+1);
            node.detail="Canonical MO; energy is Gaussian canonical eigenvalue";
            node.energy_semantics="canonical MO energy";
            node.orbital=canonical_ref(data,index);
            node.available=node.orbital.has_value();
            if(node.orbital)node.detail+="; typed source spin="+std::string(nbo_spin_name(node.orbital->spin));
            else node.detail+="; typed canonical source descriptor unavailable";
            std::ostringstream energy_label;energy_label<<std::fixed<<std::setprecision(5)<<mo.energy_hartree;
            node.label+=" ["+(node.orbital?std::string(nbo_spin_name(node.orbital->spin)):
                std::string(mo.spin==Spin::Beta?"beta":"unknown"))+"] "+
                energy_label.str()+" Ha";
            node.canonical_index=index;node.energy_hartree=mo.energy_hartree;node.occupation=mo.occupation;
            node.composition_available=std::any_of(data.links.begin(),data.links.end(),[&](const auto& link){
                return link.canonical_index==index && link.orbital.kind==state.basis_kind;});
            if(!node.composition_available)node.detail+="; no validated "+
                std::string(nbo_orbital_kind_name(state.basis_kind))+" coefficients for this spin/MO";
            node.x=std::max(370.0f,plot_width-170.0f);node.y=y;node.height=box_height;
            mo_nodes[index]=view.nodes.size();
            level_nodes[level_index].push_back(view.nodes.size());
            view.nodes.push_back(std::move(node));y+=row_step;
        }
        y+=8;
    }
    std::vector<std::size_t> level_order;
    for(std::size_t i=0;i<diagram.data.levels.size();++i)level_order.push_back(i);
    const auto& transform=diagram.data.energy_transform;
    std::stable_sort(level_order.begin(),level_order.end(),[&](std::size_t a,std::size_t b){
        return energy_display_coordinate(diagram.data.levels[a].layout_energy_hartree,transform)>
               energy_display_coordinate(diagram.data.levels[b].layout_energy_hartree,transform);
    });
    double minimum=0,maximum=1;
    if(!transform.knots.empty()){minimum=transform.knots.front().coordinate;
        maximum=transform.knots.back().coordinate;}
    float previous_bottom=20;
    const float span=std::max(300.0f,row_step*static_cast<float>(level_order.size())+80);
    for(auto level_index:level_order) {
        const double c=energy_display_coordinate(
            diagram.data.levels[level_index].layout_energy_hartree,transform);
        const double t=(c-minimum)/std::max(1e-12,maximum-minimum);
        const float desired=38+static_cast<float>(1-std::clamp(t,0.0,1.0))*span;
        float row_y=std::max(desired,previous_bottom+8);
        for(auto node_index:level_nodes[level_index]){
            view.nodes[node_index].y=row_y;
            row_y+=row_step;
        }
        previous_bottom=row_y;
    }
    for(const auto& link:data.links) {
        const auto source=basis_nodes.find(key(link.orbital));
        const auto target=mo_nodes.find(link.canonical_index);
        if(source==basis_nodes.end()||target==mo_nodes.end())continue;
        NboAomoEdge edge;edge.id="component:"+row_id(link.orbital)+":mo:"+std::to_string(link.canonical_index);
        edge.source_node=source->second;edge.target_node=target->second;
        edge.coefficient=link.coefficient;edge.weight=link.weight;edge.source=link.source;
        view.edges.push_back(std::move(edge));
    }
    // The middle column contains only actual numerical partial sums for a
    // chosen canonical MO. It never asserts a symmetry-adapted SALC.
    float fy=mo_nodes.contains(view.focused_canonical_index)?
        view.nodes[mo_nodes.at(view.focused_canonical_index)].y:55.0f;
    for(const auto& group:state.fragment_groups) {
        std::vector<std::size_t> members(group.atoms.begin(),group.atoms.end());
        const auto selection=nbo_fragment_selection(data,view.focused_canonical_index,members,state.basis_kind,false);
        if(selection.terms.empty())continue;
        NboAomoNode node;node.id="fragment:F"+std::to_string(group.id)+":mo:"+
            std::to_string(view.focused_canonical_index);
        node.label="F"+std::to_string(group.id)+" partial MO "+std::to_string(view.focused_canonical_index+1);
        node.detail=group.suggested?
            "Suggested from same-element terminal atoms sharing one displayed structural neighbour; numerical partial sum, not a symmetry SALC":
            "User-defined numerical fragment partial sum; not a symmetry SALC";
        node.fragment_group_id=group.id;node.canonical_index=view.focused_canonical_index;
        if(state.basis_kind==NboOrbitalKind::NAO) {
            double norm2=0;for(const auto& term:selection.terms)norm2+=term.coefficient*term.coefficient;
            node.metric_norm2=norm2;
        } else {
            const auto partial_view=make_nbo_selection_view(data,canonical,selection);
            if(partial_view.available && !partial_view.metric_norm2.empty())
                node.metric_norm2=partial_view.metric_norm2.front();
        }
        node.detail+="; actual partial metric norm²="+
            (node.metric_norm2?number(*node.metric_norm2):"unavailable");
        node.atoms=members;node.x=std::max(205.0f,plot_width*0.38f);
        node.y=fy;node.height=box_height;
        node.label="F"+std::to_string(group.id)+" ";
        bool first_atom=true;
        for(auto atom:members){if(!first_atom)node.label+='+';first_atom=false;
            node.label+=atom_label(canonical,atom);}
        node.label+=" / MO "+std::to_string(view.focused_canonical_index+1);
        const auto fragment_index=view.nodes.size();view.nodes.push_back(std::move(node));fy+=row_step+8;
        for(const auto& term:selection.terms)if(const auto it=basis_nodes.find(key(term.orbital));it!=basis_nodes.end()) {
            NboAomoEdge edge;edge.id="fragment-component:"+row_id(term.orbital)+":F"+std::to_string(group.id);
            edge.source_node=it->second;edge.target_node=fragment_index;edge.coefficient=term.coefficient;
            view.edges.push_back(std::move(edge));
        }
        if(const auto it=mo_nodes.find(view.focused_canonical_index);it!=mo_nodes.end()) {
            NboAomoEdge edge;edge.id="fragment-to-mo:F"+std::to_string(group.id);
            edge.source_node=fragment_index;edge.target_node=it->second;edge.coefficient=1;
            view.edges.push_back(std::move(edge));
        }
    }
    float left_end=0,right_start=plot_width;
    for(auto& node:view.nodes) {
        const float base=node.fragment_group_id?145.0f:node.group_header?120.0f:145.0f;
        node.width=std::max(base,ImGui::CalcTextSize(node.label.c_str()).x+16.0f);
        if((node.orbital && node.orbital->kind==NboOrbitalKind::Canonical) ||
           node.id.rfind("degenerate:",0)==0) {
            node.x=std::max(250.0f,plot_width-node.width-10.0f);
            right_start=std::min(right_start,node.x);
        } else if(!node.fragment_group_id)left_end=std::max(left_end,node.x+node.width);
    }
    for(auto& node:view.nodes)if(node.fragment_group_id)
        node.x=std::max(left_end+12.0f,(left_end+right_start-node.width)*0.5f);
    return view;
}
struct EdgeAppearance {
    unsigned char red=75,green=180,blue=234,alpha=8;
    float width=0.65f;
};
struct GraphAppearance {
    const NboAomoViewSnapshot& view;
    std::set<RefKey> selected_basis;
    std::map<RefKey,double> basis_to_focus;
    std::map<std::size_t,double> mo_from_selection;
    double max_focus=0,max_selected=0;
    explicit GraphAppearance(const NboAomoViewSnapshot& snapshot):view(snapshot) {
        if(view.selection)for(const auto& term:view.selection->terms)
            if(term.orbital.kind==view.basis_kind)selected_basis.insert(key(term.orbital));
        for(const auto& edge:view.edges) {
            const auto& a=view.nodes[edge.source_node],&b=view.nodes[edge.target_node];
            if(!a.orbital||!b.canonical_index)continue;
            const double magnitude=std::abs(edge.coefficient);
            if(*b.canonical_index==view.focused_canonical_index) {
                basis_to_focus[key(*a.orbital)]=std::max(basis_to_focus[key(*a.orbital)],magnitude);
                max_focus=std::max(max_focus,magnitude);
            }
            if(selected_basis.contains(key(*a.orbital))) {
                mo_from_selection[*b.canonical_index]=
                    std::max(mo_from_selection[*b.canonical_index],magnitude);
                max_selected=std::max(max_selected,magnitude);
            }
        }
    }
    static float relative(double value,double maximum) {
        return maximum>0?static_cast<float>(std::sqrt(value/maximum)):0.0f;
    }
    float node_strength(const NboAomoNode& node) const {
        if(node.orbital&&view.selection)for(const auto& term:view.selection->terms)
            if(term.orbital==*node.orbital)return 1.0f;
        if(node.canonical_index) {
            if(*node.canonical_index==view.focused_canonical_index)return 1.0f;
            if(const auto it=mo_from_selection.find(*node.canonical_index);
               it!=mo_from_selection.end())return relative(it->second,max_selected);
        }
        if(node.orbital)if(const auto it=basis_to_focus.find(key(*node.orbital));
                         it!=basis_to_focus.end())return relative(it->second,max_focus);
        return 0.0f;
    }
    EdgeAppearance edge(const NboAomoEdge& link) const {
        const auto& a=view.nodes[link.source_node],&b=view.nodes[link.target_node];
        const bool active=(b.canonical_index && *b.canonical_index==view.focused_canonical_index) ||
            (a.orbital && selected_basis.contains(key(*a.orbital)));
        double maximum=0;
        if(b.canonical_index && *b.canonical_index==view.focused_canonical_index)maximum=max_focus;
        else if(a.orbital && selected_basis.contains(key(*a.orbital)))maximum=max_selected;
        float intensity=relative(std::abs(link.coefficient),maximum);
        if(a.fragment_group_id)intensity=a.metric_norm2 && max_focus>0?
            static_cast<float>(std::min(1.0,std::sqrt(std::max(0.0,*a.metric_norm2))/max_focus)):0.0f;
        EdgeAppearance appearance;
        if(link.coefficient<0){appearance.red=242;appearance.green=145;appearance.blue=143;}
        appearance.alpha=static_cast<unsigned char>(active?14+225*intensity:8+65*intensity);
        appearance.width=0.65f+2.3f*intensity;
        return appearance;
    }
};
} // namespace

bool draw_nbo_aomo_diagram(NboAomoUIState& state,const NboIntegration& data,
                           const Wavefunction& canonical,const MODiagramViewSnapshot& diagram,
                           Language language,float scale) {
    const auto* capability=nbo_capability(data,"aomo");
    if(!capability || !capability->available()) {
        state.status=capability?capability->detail:"AO–MO decomposition data missing";
        state.drawn_snapshot.reset();return false;
    }
    if(state.source_id!=data.id){state=NboAomoUIState{};state.source_id=data.id;}
    const auto indices=central_indices(diagram);
    if(indices.empty()){state.status="The central MO diagram has no real members";state.drawn_snapshot.reset();return false;}
    const auto inspected=diagram.data.view?diagram.data.view->inspected_orbital_index:std::nullopt;
    if(inspected!=state.last_inspected){
        if(inspected && in(indices,*inspected))state.focused_canonical_index=*inspected;
        state.last_inspected=inspected;++state.revision;
    }
    if(!state.focused_canonical_index || !in(indices,*state.focused_canonical_index))
        state.focused_canonical_index=indices.front();
    if(!state.suggested_fragments_initialized)
        suggest_terminal_fragment_groups(state,data,canonical,indices);
    ImGui::SeparatorText(lt(language,"Whole AO/NAO–MO diagram","整体 AO/NAO–MO 图","全体 AO/NAO–MO 図","Diagramme AO/NAO–OM complet"));
    const bool nao=state.basis_kind==NboOrbitalKind::NAO;
    if(ImGui::RadioButton("NAO##aomo.basis",nao)){state.basis_kind=NboOrbitalKind::NAO;++state.revision;}
    validation::item("aomo.basis.nao");ImGui::SameLine();
    if(ImGui::RadioButton("Gaussian AO##aomo.basis",!nao)){state.basis_kind=NboOrbitalKind::GaussianAO;++state.revision;}
    validation::item("aomo.basis.ao");
    ImGui::SameLine();
    if(ImGui::Button("−##aomo.zoom")){state.zoom=std::max(0.35f,state.zoom/1.2f);++state.revision;}
    validation::item("aomo.zoom.out");ImGui::SameLine();
    if(ImGui::Button("+##aomo.zoom")){state.zoom=std::min(3.0f,state.zoom*1.2f);++state.revision;}
    validation::item("aomo.zoom.in");ImGui::SameLine();
    if(ImGui::Button(lt(language,"Reset view","重置视图","表示をリセット","Réinitialiser la vue"))){
        state.zoom=1;state.pan_x=0;state.pan_y=0;++state.revision;}
    validation::item("aomo.view.reset");
    auto draw_graph_options=[&]() {
    if(state.basis_kind==NboOrbitalKind::NAO) {
        if(ImGui::Checkbox(lt(language,"Show core NAOs","显示内层 NAO","内殻 NAO を表示",
            "Afficher les NAO de coeur"),&state.show_core))++state.revision;
        validation::item("aomo.show_core");ImGui::SameLine();
        if(ImGui::Checkbox(lt(language,"Show Rydberg NAOs","显示里德堡 NAO",
            "Rydberg NAO を表示","Afficher les NAO de Rydberg"),&state.show_rydberg))++state.revision;
        validation::item("aomo.show_rydberg");
    }
    if(ImGui::CollapsingHeader(lt(language,"Graph explanation and fragment groups",
        "图示说明与片段组","図の説明とフラグメント群",
        "Explication du graphe et groupes de fragments"))) {
    ImGui::TextWrapped("%s",lt(language,
        "Left: real Gaussian AO or orthogonal NAO. Middle: numerical fragment partial sums only. Right: unchanged canonical MOs. Lines retain signed coefficients; no display cutoff.",
        "左：真实 Gaussian AO 或正交 NAO；中：仅数值片段部分和；右：不变的正则 MO。连线保留带符号系数，无显示截断。",
        "左：実際の Gaussian AO または直交 NAO。中央：数値的な部分和のみ。右：元の正準 MO。符号付き係数を保持します。",
        "Gauche : AO Gaussian ou NAO orthogonales réelles ; centre : sommes partielles numériques ; droite : OM canoniques inchangées."));
    ImGui::TextDisabled("%s",lt(language,
        "Drag empty canvas to pan; mouse wheel to zoom. Atom and degenerate headers expand; only real orbital rows select 3D.",
        "拖动空白处平移，滚轮缩放。原子与简并标题可展开；只有真实轨道行选取三维。",
        "空白をドラッグして移動、ホイールで拡大。原子と縮退の見出しは展開できます。",
        "Glissez le fond pour déplacer, molette pour zoomer. Les en-têtes développent les orbitales réelles."));
    ImGui::TextWrapped("%s",lt(language,"Numerical fragment groups (not asserted SALCs)",
        "数值片段组（不声称严格 SALC）","数値的フラグメント群（厳密な SALC ではありません）",
        "Groupes de fragments numériques (pas des SALC affirmées)"));
    std::set<std::size_t> available_atoms;
    for(const auto& orbital:data.orbitals)if(orbital.ref.kind==state.basis_kind)
        available_atoms.insert(orbital.atoms.begin(),orbital.atoms.end());
    std::size_t slot=0;
    for(auto atom:available_atoms) {
        bool checked=state.draft_fragment_atoms.contains(atom);
        const auto label=atom_label(canonical,atom)+"##aomo.fragment.atom."+std::to_string(atom);
        if(ImGui::Checkbox(label.c_str(),&checked)){
            if(checked)state.draft_fragment_atoms.insert(atom);else state.draft_fragment_atoms.erase(atom);
            ++state.revision;
        }
        validation::item("aomo.fragment.atom."+std::to_string(atom));
        if(++slot%5)ImGui::SameLine();
    }
    if(slot)ImGui::NewLine();
    ImGui::BeginDisabled(state.draft_fragment_atoms.empty());
    if(ImGui::Button(lt(language,"Save fragment","保存片段","フラグメントを保存","Enregistrer le fragment"))){
        if(state.editing_fragment_id) {
            const auto it=std::find_if(state.fragment_groups.begin(),state.fragment_groups.end(),
                [&](const auto& group){return group.id==*state.editing_fragment_id;});
            if(it!=state.fragment_groups.end()){it->atoms=state.draft_fragment_atoms;it->suggested=false;}
        } else state.fragment_groups.push_back({state.next_fragment_id++,state.draft_fragment_atoms,false});
        state.draft_fragment_atoms.clear();state.editing_fragment_id.reset();++state.revision;
    }
    ImGui::EndDisabled();validation::item("aomo.fragment.add");ImGui::SameLine();
    if(ImGui::Button(lt(language,"Clear draft","清空待保存选择","未保存の選択を消去","Effacer la sélection")))
        {state.draft_fragment_atoms.clear();state.editing_fragment_id.reset();}
    validation::item("aomo.fragment.clear_draft");
    for(std::size_t i=0;i<state.fragment_groups.size();) {
        const auto& group=state.fragment_groups[i];
        std::string label="F"+std::to_string(group.id)+" (";
        for(auto atom:group.atoms)label+=atom_label(canonical,atom)+" ";label+=")";
        if(group.suggested)label+=" [suggested]";
        ImGui::Text("%s",label.c_str());ImGui::SameLine();
        if(ImGui::SmallButton(("Edit##aomo.fragment."+std::to_string(group.id)).c_str())) {
            state.draft_fragment_atoms=group.atoms;state.editing_fragment_id=group.id;
        }
        validation::item("aomo.fragment.edit."+std::to_string(group.id));
        ImGui::SameLine();
        if(ImGui::SmallButton(("Delete##aomo.fragment."+std::to_string(group.id)).c_str())){
            state.fragment_groups.erase(state.fragment_groups.begin()+static_cast<std::ptrdiff_t>(i));
            ++state.revision;continue;
        }
        validation::item("aomo.fragment.delete."+std::to_string(group.id));++i;
    }
    }
    };
    const float canvas_width=std::max(540.0f*scale,ImGui::GetContentRegionAvail().x);
    const auto snapshot=std::make_shared<const NboAomoViewSnapshot>(
        make_snapshot(state,data,canonical,diagram,canvas_width));
    state.drawn_snapshot=snapshot;
    state.status=snapshot->capability_status+": "+snapshot->capability_detail;
    validation::field("aomo.snapshot",snapshot->id);
    validation::anchor("aomo.graph");
    float max_node_bottom=0;
    for(const auto& node:snapshot->nodes)max_node_bottom=std::max(max_node_bottom,node.y+node.height);
    const float canvas_height=std::clamp(max_node_bottom+25.0f,420.0f,620.0f);
    ImGui::InvisibleButton("##aomo.canvas",ImVec2(canvas_width,canvas_height));
    validation::item("aomo.graph.canvas");
    const auto origin=ImGui::GetItemRectMin();
    const auto canvas_max=ImGui::GetItemRectMax();
    const auto transform=[&](const NboAomoNode& node){return ImVec2(origin.x+state.pan_x+node.x*state.zoom,
                                                              origin.y+state.pan_y+node.y*state.zoom);};
    const bool hover=ImGui::IsItemHovered();
    auto* draw=ImGui::GetWindowDrawList();
    draw->PushClipRect(origin,ImGui::GetItemRectMax(),true);
    draw->AddRectFilled(origin,ImGui::GetItemRectMax(),IM_COL32(24,31,44,255),6*scale);
    draw->AddText(ImVec2(origin.x+12,origin.y+8),IM_COL32(144,166,189,255),
        "AO/NAO: no shared energy axis");
    const auto axis_label="Canonical E (Ha): "+snapshot->mo_energy_axis_mode;
    const float label_width=ImGui::CalcTextSize(axis_label.c_str()).x;
    draw->AddText(ImVec2(std::max(origin.x+220,canvas_max.x-label_width-12),origin.y+8),
        IM_COL32(144,166,189,255),axis_label.c_str());
    const GraphAppearance appearance(*snapshot);
    for(const auto& edge:snapshot->edges){
        const auto& a=snapshot->nodes[edge.source_node];const auto& b=snapshot->nodes[edge.target_node];
        const auto style=appearance.edge(edge);
        const auto pa=transform(a),pb=transform(b);
        draw->AddLine(ImVec2(pa.x+node_width(a)*state.zoom,pa.y+a.height*0.5f*state.zoom),
            ImVec2(pb.x,pb.y+b.height*0.5f*state.zoom),
            IM_COL32(style.red,style.green,style.blue,style.alpha),
            style.width*scale);
    }
    std::optional<std::size_t> clicked;
    float best=1e9f;
    for(std::size_t i=0;i<snapshot->nodes.size();++i){
        const auto& node=snapshot->nodes[i];const auto p=transform(node);
        const float width=node_width(node)*state.zoom;
        const float height=node.height*state.zoom;
        const ImVec2 q(p.x+width,p.y+height);
        const ImVec2 hit_min(std::max(origin.x,p.x),std::max(origin.y,p.y));
        const ImVec2 hit_max(std::min(canvas_max.x,q.x),std::min(canvas_max.y,q.y));
        if(hit_min.x<hit_max.x && hit_min.y<hit_max.y)
            validation::hit("aomo.node."+node.id,hit_min,hit_max);
        const float strength=appearance.node_strength(node);
        const bool selected_node=strength>=0.99f;
        const int green=static_cast<int>(58+65*strength);
        const int blue=static_cast<int>(82+84*strength);
        draw->AddRectFilled(p,q,node.group_header?IM_COL32(63,77,101,255):
            !node.composition_available?IM_COL32(63,64,69,255):
            IM_COL32(43,green,blue,255),4*scale);
        draw->AddRect(p,q,selected_node?IM_COL32(122,223,255,255):IM_COL32(91,115,146,255),4*scale);
        draw->PushClipRect(p,q,true);
        draw->AddText(ImVec2(p.x+5*state.zoom,p.y+(height-ImGui::GetTextLineHeight())*0.5f),
            IM_COL32(235,241,248,255),node.label.c_str());
        draw->PopClipRect();
        if(hover && ImGui::IsMouseHoveringRect(p,q) && width*height<best){clicked=i;best=width*height;
            ImGui::SetTooltip("%s\n%s\n%s%s%s",node.label.c_str(),node.detail.c_str(),
                node.energy_semantics.c_str(),node.energy_hartree?" (Ha) = " : "",
                node.energy_hartree?number(*node.energy_hartree).c_str():"");}
    }
    draw->PopClipRect();
    if(hover && ImGui::GetIO().MouseWheel!=0){
        state.zoom=std::clamp(state.zoom*(ImGui::GetIO().MouseWheel>0?1.12f:0.89f),0.35f,3.0f);++state.revision;
    }
    if(hover && !clicked && ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
        const auto delta=ImGui::GetIO().MouseDelta;state.pan_x+=delta.x;state.pan_y+=delta.y;++state.revision;
    }
    if(hover && clicked && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        const auto& node=snapshot->nodes[*clicked];
        if(node.group_header){
            if(node.id.rfind("atom:",0)==0){const auto atom=static_cast<std::size_t>(std::stoull(node.id.substr(5)));
                if(!state.collapsed_atoms.erase(atom))state.collapsed_atoms.insert(atom);++state.revision;}
            else if(node.id.rfind("degenerate:",0)==0){const auto level=static_cast<std::size_t>(std::stoull(node.id.substr(11)));
                if(ImGui::GetIO().KeyShift) {
                    NboOrbitalSelection selection;selection.dataset_id=data.id;
                    selection.label="Canonical degenerate set "+std::to_string(level+1);
                    selection.mode=NboSelectionMode::Overlay;
                    for(const auto index:level_members(diagram.data.levels[level]))
                        if(index<canonical.orbitals.size())
                            if(const auto ref=canonical_ref(data,index))
                                selection.terms.push_back({*ref,1});
                    if(!selection.terms.empty())state.pending_selection=std::move(selection);
                } else {if(!state.collapsed_levels.erase(level))state.collapsed_levels.insert(level);++state.revision;}}
        } else if(node.fragment_group_id && node.canonical_index) {
            state.pending_selection=nbo_fragment_selection(data,*node.canonical_index,node.atoms,state.basis_kind,false);
        } else if(node.orbital) {
            state.pending_selection=nbo_single_selection(data,*node.orbital);
            if(node.canonical_index)state.focused_canonical_index=*node.canonical_index;
        }
    }
    auto activate_edge=[&](const NboAomoEdge& edge) {
        const auto& a=snapshot->nodes[edge.source_node];
        const auto& b=snapshot->nodes[edge.target_node];
        if(!b.canonical_index)return;
        if(a.fragment_group_id) {
            state.pending_selection=nbo_fragment_selection(data,*b.canonical_index,
                a.atoms,state.basis_kind,false);
        } else if(a.orbital) {
            const auto it=std::find_if(data.links.begin(),data.links.end(),[&](const auto& link){
                return link.orbital==*a.orbital && link.canonical_index==*b.canonical_index;});
            if(it!=data.links.end())state.pending_selection=nbo_component_selection(data,*it);
        }
    };
    std::vector<std::pair<float,const NboAomoEdge*>> nearby_edges;
    const auto mouse=ImGui::GetIO().MousePos;
    for(const auto& edge:snapshot->edges) {
        const auto& a=snapshot->nodes[edge.source_node];
        const auto& b=snapshot->nodes[edge.target_node];
        const auto p=transform(a),q=transform(b);
        const ImVec2 start(p.x+node_width(a)*state.zoom,p.y+a.height*0.5f*state.zoom),
            end(q.x,q.y+b.height*0.5f*state.zoom);
        const ImVec2 mid((start.x+end.x)*0.5f,(start.y+end.y)*0.5f);
        const ImVec2 hit_min(std::max(origin.x,mid.x-5*scale),std::max(origin.y,mid.y-5*scale));
        const ImVec2 hit_max(std::min(canvas_max.x,mid.x+5*scale),std::min(canvas_max.y,mid.y+5*scale));
        if(hit_min.x<hit_max.x && hit_min.y<hit_max.y)
            validation::hit("aomo.edge."+edge.id,hit_min,hit_max);
        if(!hover||clicked||!b.canonical_index)continue;
        const float dx=end.x-start.x,dy=end.y-start.y;
        const float t=std::clamp(((mouse.x-start.x)*dx+(mouse.y-start.y)*dy)/
            std::max(1.0f,dx*dx+dy*dy),0.0f,1.0f);
        const float distance=std::hypot(mouse.x-(start.x+t*dx),mouse.y-(start.y+t*dy));
        if(distance<=8*scale)nearby_edges.push_back({distance,&edge});
    }
    if(!nearby_edges.empty()) {
        std::sort(nearby_edges.begin(),nearby_edges.end(),[](const auto& a,const auto& b){
            if(std::abs(a.first-b.first)>0.001f)return a.first<b.first;
            return a.second->id<b.second->id;});
        const auto& edge=*nearby_edges.front().second;
        const auto& source_node=snapshot->nodes[edge.source_node];
        if(source_node.fragment_group_id)
            ImGui::SetTooltip("External factor=%s; actual partial metric norm²=%s\n%s",
                number(edge.coefficient).c_str(),
                source_node.metric_norm2?number(*source_node.metric_norm2).c_str():"unavailable",
                source_node.detail.c_str());
        else ImGui::SetTooltip("c=%s%s%s\n%s",
            number(edge.coefficient).c_str(),edge.weight?"  |c|²=":"",
            edge.weight?number(*edge.weight).c_str():"",
            source_name(edge.source).c_str());
        if(ImGui::IsMouseClicked(ImGuiMouseButton_Left) ||
           ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
            state.ambiguous_edge_ids.clear();
            for(const auto& candidate:nearby_edges)
                if(candidate.first-nearby_edges.front().first<=1.25f*scale)
                    state.ambiguous_edge_ids.push_back(candidate.second->id);
            if(state.ambiguous_edge_ids.size()==1)activate_edge(edge);
            else ImGui::OpenPopup("##aomo.overlapping.edges");
        }
    }
    if(ImGui::BeginPopup("##aomo.overlapping.edges")) {
        ImGui::TextUnformatted("Choose a numerical component");
        for(const auto& id:state.ambiguous_edge_ids) {
            const auto it=std::find_if(snapshot->edges.begin(),snapshot->edges.end(),
                [&](const auto& edge){return edge.id==id;});
            if(it==snapshot->edges.end())continue;
            const auto& a=snapshot->nodes[it->source_node];
            const auto& b=snapshot->nodes[it->target_node];
            const auto label=a.label+" → "+b.label+"  c="+number(it->coefficient)+"##"+id;
            if(ImGui::Selectable(label.c_str())){activate_edge(*it);ImGui::CloseCurrentPopup();}
            validation::item("aomo.edge.choice."+id);
        }
        ImGui::EndPopup();
    }
    draw_graph_options();
    if(snapshot->focused_projection_weight && snapshot->focused_projection_residual_norm) {
        const auto* full_cap=nbo_capability(data,"aomo_full");
        const auto status=full_cap&&full_cap->available()?"verified full local span":
            "verified local-subspace projection; retained norm may be below 1";
        ImGui::TextWrapped("MO %zu: NAO Σ|c|²=%s; AO-metric residual norm=%s — %s",
            snapshot->focused_canonical_index+1,
            number(*snapshot->focused_projection_weight).c_str(),
            number(*snapshot->focused_projection_residual_norm).c_str(),status);
        validation::field("aomo.projection",std::to_string(snapshot->focused_canonical_index)+
            ";weight="+number(*snapshot->focused_projection_weight)+
            ";residual_norm="+number(*snapshot->focused_projection_residual_norm));
    }
    ImGui::TextWrapped("%s",lt(language,
        "Blue/red lines: positive/negative signed coefficients. Intensity follows continuous |c|; faint numerical terms remain selectable. NAO |c|² is an orthogonal-basis weight; Gaussian AO c² is not a population.",
        "蓝/红线代表正/负系数；明暗随 |c| 连续变化，微弱数值项仍可选。正交 NAO 的 |c|² 是权重；Gaussian AO 的 c² 不是布居。",
        "青/赤の線は正/負の係数。濃淡は |c| に連続的に従い、微小項も選択可能。直交 NAO の |c|² は重みですが、Gaussian AO の c² は分布ではありません。",
        "Les lignes bleues/rouges indiquent des coefficients positifs/négatifs. Leur intensité suit |c| sans seuil ; les petits termes restent sélectionnables. |c|² des NAO orthogonales est un poids, c² des AO Gaussian n'est pas une population."));
    ImGui::TextDisabled("%s",lt(language,
        "Click a real node for 3D; click a link for its signed weighted component.",
        "点击真实节点查看三维；点击连线查看带符号系数的分量。",
        "実軌道をクリックして3D表示。線をクリックして符号付き成分を確認。",
        "Cliquez sur une orbitale pour la 3D ; cliquez sur un lien pour sa composante signée."));
    ImGui::Text("%zu by atom collapse, %zu by Core/Rydberg class, %zu by degeneracy collapse",
        snapshot->hidden_basis_count,snapshot->hidden_class_count,snapshot->hidden_mo_count);
    if(ImGui::Checkbox(lt(language,"Full signed coefficients","完整带符号系数",
        "符号付き係数の全表","Tous les coefficients signés"),&state.show_full_numeric))++state.revision;
    validation::item("aomo.full_numeric");
    if(state.sum_canonical_index!=snapshot->focused_canonical_index ||
       state.sum_basis_kind!=state.basis_kind) {
        state.sum_component_ids.clear();
        state.sum_canonical_index=snapshot->focused_canonical_index;
        state.sum_basis_kind=state.basis_kind;
    }
    if(state.show_full_numeric) {
        validation::anchor("aomo.sum.controls");
        ImGui::TextWrapped("MO %zu; %s; %s",snapshot->focused_canonical_index+1,
            nbo_orbital_kind_name(state.basis_kind),
            state.basis_kind==NboOrbitalKind::NAO?"orthonormal: |c|² is a weight":
            "Gaussian AO: c² is not an orbital weight in a nonorthogonal basis");
        const auto links=nbo_links_for_mo(data,snapshot->focused_canonical_index,state.basis_kind);
        if(ImGui::Button(lt(language,"Add every local term","加入全部局域项",
            "すべての項を追加","Ajouter tous les termes"))) {
            for(const auto& link:links)state.sum_component_ids.insert(row_id(link.orbital));
        }
        validation::item("aomo.sum.all");ImGui::SameLine();
        if(ImGui::Button(lt(language,"Clear terms","清空项",
            "項を消去","Effacer les termes")))state.sum_component_ids.clear();
        validation::item("aomo.sum.clear");
        std::vector<NboOrbitalTerm> terms;
        for(const auto& link:links)if(state.sum_component_ids.contains(row_id(link.orbital)))
            terms.push_back({link.orbital,link.coefficient});
        ImGui::BeginDisabled(terms.empty());
        if(ImGui::Button(lt(language,"Show partial sum in 3D","三维查看部分和",
            "部分和を3D表示","Afficher la somme partielle en 3D"))) {
            NboOrbitalSelection selection;selection.dataset_id=data.id;
            selection.mode=NboSelectionMode::PartialSum;
            selection.target_canonical_index=snapshot->focused_canonical_index;
            selection.label="Signed partial sum of "+std::to_string(terms.size())+
                " actual "+nbo_orbital_kind_name(state.basis_kind)+" components of MO "+
                std::to_string(snapshot->focused_canonical_index+1);
            selection.terms=terms;selection.normalize=false;
            state.pending_selection=std::move(selection);
        }
        validation::item("aomo.sum.partial");ImGui::SameLine();
        if(ImGui::Button(lt(language,"Overlay terms in 3D","三维逐项叠加",
            "各項を3D重ね表示","Superposer les termes en 3D"))) {
            NboOrbitalSelection selection;selection.dataset_id=data.id;
            selection.mode=NboSelectionMode::Overlay;
            selection.target_canonical_index=snapshot->focused_canonical_index;
            selection.label="Separate signed components of MO "+
                std::to_string(snapshot->focused_canonical_index+1);
            selection.terms=terms;selection.normalize=false;
            state.pending_selection=std::move(selection);
        }
        validation::item("aomo.sum.overlay");
        ImGui::EndDisabled();
        if(const auto canonical=canonical_ref(data,snapshot->focused_canonical_index)) {
            if(ImGui::Button(lt(language,"Show full canonical MO","查看完整正则 MO",
                "正準 MO 全体を表示","Afficher l'OM canonique entière")))
                state.pending_selection=nbo_single_selection(data,*canonical);
            validation::item("aomo.sum.full_mo");
        }
        ImGui::TextDisabled("%zu / %zu terms selected; raw amplitudes, no normalization or coefficient cutoff",
            terms.size(),links.size());
        if(state.basis_kind==NboOrbitalKind::NAO)
            ImGui::TextWrapped("All NAO terms represent the validated local-subspace projection; use Full canonical MO for any residual outside that span.");
        validation::anchor("aomo.coefficients");
        ImGui::BeginChild("##aomo.coefficients",ImVec2(0,240*scale),ImGuiChildFlags_Border);
        for(std::size_t i=0;i<links.size();++i) {
            const auto& link=links[i];
            const auto* descriptor=nbo_orbital(data,link.orbital);
            const std::string label=(descriptor?descriptor->label:row_id(link.orbital))+
                "  c="+number(link.coefficient)+(link.weight?"  |c|²="+number(*link.weight):"")+
                "##aomo.coefficient."+std::to_string(i);
            bool checked=state.sum_component_ids.contains(row_id(link.orbital));
            if(ImGui::Checkbox(("##aomo.sum.item."+row_id(link.orbital)).c_str(),&checked)) {
                if(checked)state.sum_component_ids.insert(row_id(link.orbital));
                else state.sum_component_ids.erase(row_id(link.orbital));
            }
            validation::item("aomo.sum.item."+row_id(link.orbital));
            ImGui::SameLine();
            if(ImGui::Selectable(label.c_str(),false))
                state.pending_selection=nbo_component_selection(data,link);
            validation::item("aomo.component."+row_id(link.orbital));
            if(ImGui::IsItemHovered())ImGui::SetTooltip("%s",source_name(link.source).c_str());
        }
        ImGui::EndChild();
        validation::field("aomo.coefficients.count",std::to_string(links.size()));
    }
    if(ImGui::Button(lt(language,"Export whole diagram","导出整体图","全体図を書き出す","Exporter le diagramme complet"))){
        state.export_requested=true;}
    validation::item("aomo.export");
    if(!state.export_status.empty())ImGui::TextWrapped("%s",state.export_status.c_str());
    return true;
}

NboAomoExportResult export_nbo_aomo_bundle(const NboAomoViewSnapshot& view,
    const NboIntegration& data,const std::filesystem::path& base) {
    NboAomoExportResult result;
    const auto file=[&](const char* suffix){auto path=base;path+=suffix;return path;};
    result.json_path=file(".aomo.json");result.csv_path=file(".aomo.csv");
    result.svg_path=file(".aomo.svg");result.png_path=file(".aomo.png");
    if(view.integration_id!=data.id || view.capability_status!="available") {
        result.error="AO–MO snapshot does not match verified integration";return result;
    }
    try {
        std::set<std::string> visible_nodes;
        std::set<RefKey> visible_orbitals;
        for(const auto& node:view.nodes){visible_nodes.insert(node.id);
            if(node.orbital)visible_orbitals.insert(key(*node.orbital));}
        {
            std::ofstream out(result.csv_path,std::ios::binary);
            if(!out)throw std::runtime_error("Cannot write AO–MO CSV");
            out<<std::setprecision(17);
            out<<"snapshot_id,integration_id,mo_snapshot_id,in_central_view,visible_link,basis_kind,basis_index,spin,canonical_index,coefficient,orthonormal_weight,nao_projection_weight,ao_metric_residual_norm,source_path,source_line,source_block\n";
            for(const auto& link:data.links){
                const bool central=in(view.central_mo_indices,link.canonical_index);
                const auto target="canonical_mo:"+std::to_string(link.canonical_index);
                const bool shown=visible_orbitals.contains(key(link.orbital)) && visible_nodes.contains(target);
                out<<csv(view.id)<<','<<csv(data.id)<<','<<csv(view.mo_snapshot_id)<<','
                   <<(central?"true":"false")<<','<<(shown?"true":"false")<<','
                   <<csv(nbo_orbital_kind_name(link.orbital.kind))<<','<<link.orbital.index<<','
                   <<csv(nbo_spin_name(link.orbital.spin))<<','<<link.canonical_index<<','
                   <<link.coefficient<<',';
                if(link.weight)out<<*link.weight;
                const auto* decomposition=nbo_mo_decomposition(data.dataset,link.canonical_index);
                out<<',';
                if(decomposition&&decomposition->weight_sum)out<<*decomposition->weight_sum;
                out<<',';
                if(decomposition&&decomposition->projection_residual_norm)
                    out<<*decomposition->projection_residual_norm;
                out<<','<<csv(link.source.path)<<','<<link.source.line_begin<<','
                   <<csv(link.source.block)<<'\n';
            }
            if(!out)throw std::runtime_error("AO–MO CSV write failed");
        }
        result.csv=true;
        {
            std::ofstream out(result.json_path,std::ios::binary);
            if(!out)throw std::runtime_error("Cannot write AO–MO JSON");
            out<<std::setprecision(17);
            out<<"{\"schema\":\"cov_aomo_view_v1\",\"snapshot_id\":"<<quote(view.id)
               <<",\"integration_id\":"<<quote(data.id)<<",\"mo_snapshot_id\":"<<quote(view.mo_snapshot_id)
               <<",\"mo_energy_axis_mode\":"<<quote(view.mo_energy_axis_mode)
               <<",\"mo_energy_axis_detail\":"<<quote(view.mo_energy_axis_detail)
               <<",\"basis_kind\":"<<quote(nbo_orbital_kind_name(view.basis_kind))
               <<",\"capability_status\":"<<quote(view.capability_status)
               <<",\"capability_detail\":"<<quote(view.capability_detail)
               <<",\"focused_canonical_index\":"<<view.focused_canonical_index
               <<",\"focused_projection_weight\":";
            if(view.focused_projection_weight)out<<*view.focused_projection_weight;else out<<"null";
            out<<",\"focused_projection_residual_norm\":";
            if(view.focused_projection_residual_norm)out<<*view.focused_projection_residual_norm;else out<<"null";
            out
               <<",\"full_numerical_file\":\".aomo.csv\",\"full_numerical_scope\":\"all_verified_links\",\"zoom\":"<<view.zoom
               <<",\"pan\":["<<view.pan_x<<','<<view.pan_y<<"],\"hidden_basis_count\":"<<view.hidden_basis_count
               <<",\"hidden_class_count\":"<<view.hidden_class_count
               <<",\"show_core\":"<<(view.show_core?"true":"false")
               <<",\"show_rydberg\":"<<(view.show_rydberg?"true":"false")
               <<",\"hidden_mo_count\":"<<view.hidden_mo_count<<",\"central_mo_indices\":[";
            for(std::size_t i=0;i<view.central_mo_indices.size();++i) {
                if(i)out<<',';out<<view.central_mo_indices[i];
            }
            out<<"],\"sum_component_ids\":[";
            for(std::size_t i=0;i<view.sum_component_ids.size();++i){
                if(i)out<<',';out<<quote(view.sum_component_ids[i]);
            }
            out<<"],\"fragment_groups\":[";
            for(std::size_t i=0;i<view.fragment_groups.size();++i) {
                if(i)out<<',';const auto& group=view.fragment_groups[i];
                out<<"{\"id\":"<<group.id<<",\"suggested\":"
                   <<(group.suggested?"true":"false")<<",\"atoms0\":[";
                bool first=true;for(auto atom:group.atoms){if(!first)out<<',';first=false;out<<atom;}
                out<<"]}";
            }
            out<<"],\"nodes\":[";
            for(std::size_t i=0;i<view.nodes.size();++i) {
                if(i)out<<',';const auto& node=view.nodes[i];
                out<<"{\"id\":"<<quote(node.id)<<",\"label\":"<<quote(node.label)
                   <<",\"detail\":"<<quote(node.detail)
                   <<",\"energy_semantics\":"<<quote(node.energy_semantics)
                   <<",\"x\":"<<node.x<<",\"y\":"<<node.y
                   <<",\"width\":"<<node.width<<",\"height\":"<<node.height
                   <<",\"group_header\":"<<(node.group_header?"true":"false")
                   <<",\"available\":"<<(node.available?"true":"false")
                   <<",\"composition_available\":"<<(node.composition_available?"true":"false")
                   <<",\"atoms0\":[";
                for(std::size_t j=0;j<node.atoms.size();++j){if(j)out<<',';out<<node.atoms[j];}
                out<<"],\"orbital\":";
                if(node.orbital)out<<"{\"kind\":"<<quote(nbo_orbital_kind_name(node.orbital->kind))
                    <<",\"spin\":"<<quote(nbo_spin_name(node.orbital->spin))
                    <<",\"index\":"<<node.orbital->index<<"}";else out<<"null";
                out<<",\"canonical_index\":";
                if(node.canonical_index)out<<*node.canonical_index;else out<<"null";
                out<<",\"energy_hartree\":";
                if(node.energy_hartree)out<<*node.energy_hartree;else out<<"null";
                out<<",\"occupation\":";
                if(node.occupation)out<<*node.occupation;else out<<"null";
                out<<",\"metric_norm2\":";
                if(node.metric_norm2)out<<*node.metric_norm2;else out<<"null";
                out<<'}';
            }
            out<<"],\"edges\":[";
            for(std::size_t i=0;i<view.edges.size();++i) {
                if(i)out<<',';const auto& edge=view.edges[i];
                out<<"{\"id\":"<<quote(edge.id)<<",\"source\":"<<quote(view.nodes.at(edge.source_node).id)
                   <<",\"target\":"<<quote(view.nodes.at(edge.target_node).id)
                   <<",\"coefficient\":"<<edge.coefficient<<",\"weight\":";
                if(edge.weight)out<<*edge.weight;else out<<"null";
                out<<",\"source_path\":"<<quote(edge.source.path)<<",\"source_line\":"<<edge.source.line_begin
                   <<",\"external_factor_only\":"<<(view.nodes.at(edge.source_node).fragment_group_id?"true":"false")<<'}';
            }
            out<<"],\"selection\":";
            if(view.selection) {
                out<<"{\"dataset_id\":"<<quote(view.selection->dataset_id)
                   <<",\"label\":"<<quote(view.selection->label)
                   <<",\"mode\":"<<static_cast<int>(view.selection->mode)<<",\"terms\":[";
                for(std::size_t i=0;i<view.selection->terms.size();++i) {
                    if(i)out<<',';const auto& term=view.selection->terms[i];
                    out<<"{\"kind\":"<<quote(nbo_orbital_kind_name(term.orbital.kind))
                       <<",\"spin\":"<<quote(nbo_spin_name(term.orbital.spin))
                       <<",\"index\":"<<term.orbital.index
                       <<",\"coefficient\":"<<term.coefficient<<'}';
                }
                out<<"]}";
            }else out<<"null";
            out<<"}";
            if(!out)throw std::runtime_error("AO–MO JSON write failed");
        }
        result.json=true;
        float max_y=0,max_x=0;
        for(const auto& node:view.nodes){max_y=std::max(max_y,node.y+node.height+20);
            max_x=std::max(max_x,node.x+node_width(node)+30);}
        const auto width=std::max(540,static_cast<int>(std::ceil(max_x)));
        const auto height=std::max(480,static_cast<int>(std::ceil(max_y+35)));
        const GraphAppearance appearance(view);
        {
            std::ofstream out(result.svg_path,std::ios::binary);
            if(!out)throw std::runtime_error("Cannot write AO–MO SVG");
            out<<"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""<<width
               <<"\" height=\""<<height<<"\" viewBox=\"0 0 "<<width<<' '<<height<<"\">\n";
            out<<"<rect width=\"100%\" height=\"100%\" fill=\"#18202d\"/>\n";
            out<<"<text x=\"28\" y=\"32\" fill=\"#d6e6f6\" font-size=\"16\">"
               <<"AO/NAO (no shared E axis) to canonical MO | E (Ha): "
               <<xml(view.mo_energy_axis_mode)<<" | "<<xml(view.mo_snapshot_id)<<"</text>\n";
            for(const auto& edge:view.edges) {
                const auto& a=view.nodes.at(edge.source_node);const auto& b=view.nodes.at(edge.target_node);
                const auto style=appearance.edge(edge);
                const char* stroke=edge.coefficient<0?"#f2918f":"#4bb4ea";
                out<<"<line x1=\""<<a.x+node_width(a)<<"\" y1=\""<<a.y+a.height*0.5f
                   <<"\" x2=\""<<b.x<<"\" y2=\""<<b.y+b.height*0.5f
                   <<"\" stroke=\""<<stroke<<"\" stroke-opacity=\""
                   <<static_cast<double>(style.alpha)/255.0<<"\" stroke-width=\""<<style.width
                   <<"\"><title>";
                if(a.fragment_group_id)out<<"external factor="<<edge.coefficient
                    <<"; actual partial metric norm²="<<(a.metric_norm2?number(*a.metric_norm2):"unavailable");
                else out<<"c="<<edge.coefficient<<" "<<xml(source_name(edge.source));
                out<<"</title></line>\n";
            }
            for(const auto& node:view.nodes) {
                const float strength=appearance.node_strength(node);
                const int green=static_cast<int>(58+65*strength),blue=static_cast<int>(82+84*strength);
                std::ostringstream fill;fill<<"#"<<std::hex<<std::setfill('0')
                    <<std::setw(2)<<43<<std::setw(2)<<green<<std::setw(2)<<blue;
                out<<"<g><title>"<<xml(node.detail)<<"</title><rect x=\""<<node.x<<"\" y=\""<<node.y
                   <<"\" width=\""<<node.width<<"\" height=\""<<node.height<<"\" rx=\"4\" fill=\""
                   <<(node.group_header?"#3f4d65":!node.composition_available?"#3f4045":fill.str())
                   <<"\" stroke=\""<<(strength>=0.99f?"#7adfff":"#5b7392")
                   <<"\"/><text x=\""<<node.x+5<<"\" y=\""<<node.y+node.height*0.5f+4
                   <<"\" fill=\"white\" font-size=\"12\">"<<xml(node.label)<<"</text></g>\n";
            }
            out<<"</svg>\n";
            if(!out)throw std::runtime_error("AO–MO SVG write failed");
        }
        result.svg=true;
        // PNG uses the same frozen graph geometry and IDs as SVG/JSON.
        if(height>16000)throw std::runtime_error("AO–MO PNG exceeds supported raster height");
        std::vector<unsigned char> rgba(static_cast<std::size_t>(width)*height*4,255);
        for(std::size_t p=0;p<rgba.size();p+=4){rgba[p]=24;rgba[p+1]=32;rgba[p+2]=45;}
        auto pixel=[&](int x,int y,unsigned char r,unsigned char g,unsigned char b){
            if(x<0||y<0||x>=width||y>=height)return;
            const auto p=(static_cast<std::size_t>(y)*width+x)*4;rgba[p]=r;rgba[p+1]=g;rgba[p+2]=b;
        };
        auto blend=[&](int x,int y,const EdgeAppearance& style){
            if(x<0||y<0||x>=width||y>=height)return;
            const auto p=(static_cast<std::size_t>(y)*width+x)*4;
            const float t=static_cast<float>(style.alpha)/255.0f;
            rgba[p]=static_cast<unsigned char>(rgba[p]*(1-t)+style.red*t);
            rgba[p+1]=static_cast<unsigned char>(rgba[p+1]*(1-t)+style.green*t);
            rgba[p+2]=static_cast<unsigned char>(rgba[p+2]*(1-t)+style.blue*t);
        };
        auto line=[&](float ax,float ay,float bx,float by,const EdgeAppearance& style){
            const int steps=std::max(1,static_cast<int>(std::ceil(std::hypot(bx-ax,by-ay))));
            for(int i=0;i<=steps;++i){const float t=static_cast<float>(i)/steps;
                const int x=static_cast<int>(std::round(ax+(bx-ax)*t));
                const int y=static_cast<int>(std::round(ay+(by-ay)*t));
                const float radius=style.width*0.5f;
                for(int dy=-2;dy<=2;++dy)for(int dx=-2;dx<=2;++dx)
                    if(dx*dx+dy*dy<=radius*radius+0.5f)blend(x+dx,y+dy,style);
            }
        };
        for(const auto& edge:view.edges) {
            const auto& a=view.nodes.at(edge.source_node);const auto& b=view.nodes.at(edge.target_node);
            line(a.x+node_width(a),a.y+a.height*0.5f,b.x,b.y+b.height*0.5f,
                appearance.edge(edge));
        }
        for(const auto& node:view.nodes) {
            const int x=static_cast<int>(node.x),y=static_cast<int>(node.y);
            const int nw=static_cast<int>(std::ceil(node.width));
            const int nh=static_cast<int>(std::ceil(node.height));
            const float strength=appearance.node_strength(node);
            const unsigned char red=node.group_header?63:!node.composition_available?63:43;
            const unsigned char green=node.group_header?77:!node.composition_available?64:
                static_cast<unsigned char>(58+65*strength);
            const unsigned char blue=node.group_header?101:!node.composition_available?69:
                static_cast<unsigned char>(82+84*strength);
            for(int yy=y;yy<y+nh;++yy)for(int xx=x;xx<x+nw;++xx)
                pixel(xx,yy,red,green,blue);
            const unsigned char border_r=strength>=0.99f?122:91;
            const unsigned char border_g=strength>=0.99f?223:115;
            const unsigned char border_b=strength>=0.99f?255:146;
            for(int xx=x;xx<x+nw;++xx){pixel(xx,y,border_r,border_g,border_b);
                pixel(xx,y+nh-1,border_r,border_g,border_b);}
            for(int yy=y;yy<y+nh;++yy){pixel(x,yy,border_r,border_g,border_b);
                pixel(x+nw-1,yy,border_r,border_g,border_b);}
            // Full Unicode labels remain in SVG and JSON. The raster's built-in
            // ASCII font renders each concise visible label.
            const auto& label=node.label;
            int tx=x+4;
            for(char c:label){
                if(tx+6>=x+nw)break;
                const auto glyph=raster_glyph(c);
                for(int row=0;row<7;++row)for(int col=0;col<5;++col)
                    if(glyph[row]&(1u<<(4-col)))pixel(tx+col,y+8+row,228,240,249);
                tx+=6;
            }
        }
        int title_x=14;
        const std::string title="AO/NAO NO SHARED E AXIS  MO E HA "+view.mo_energy_axis_mode;
        for(char c:title) {
            if(title_x+6>=width)break;
            const auto glyph=raster_glyph(c);
            for(int row=0;row<7;++row)for(int col=0;col<5;++col)
                if(glyph[row]&(1u<<(4-col)))pixel(title_x+col,14+row,213,230,246);
            title_x+=6;
        }
        auto put32=[](std::vector<unsigned char>& out,std::uint32_t v){
            out.push_back(static_cast<unsigned char>(v>>24));out.push_back(static_cast<unsigned char>(v>>16));
            out.push_back(static_cast<unsigned char>(v>>8));out.push_back(static_cast<unsigned char>(v));
        };
        auto crc=[](const unsigned char* bytes,std::size_t size){
            std::uint32_t v=0xffffffffu;
            for(std::size_t i=0;i<size;++i){v^=bytes[i];for(int j=0;j<8;++j)v=(v>>1)^((v&1)?0xedb88320u:0u);}
            return ~v;
        };
        auto chunk=[&](std::vector<unsigned char>& out,const char type[4],const std::vector<unsigned char>& bytes){
            put32(out,static_cast<std::uint32_t>(bytes.size()));const auto start=out.size();
            out.insert(out.end(),type,type+4);out.insert(out.end(),bytes.begin(),bytes.end());
            put32(out,crc(out.data()+start,out.size()-start));
        };
        std::vector<unsigned char> raw;raw.reserve(static_cast<std::size_t>(height)*(1+width*4));
        for(int y=0;y<height;++y){raw.push_back(0);const auto p=static_cast<std::size_t>(y)*width*4;
            raw.insert(raw.end(),rgba.begin()+p,rgba.begin()+p+width*4);}
        std::vector<unsigned char> z={0x78,0x01};
        for(std::size_t pos=0;pos<raw.size();){
            const auto n=std::min<std::size_t>(65535,raw.size()-pos);const bool last=pos+n==raw.size();
            z.push_back(last?1:0);z.push_back(static_cast<unsigned char>(n));
            z.push_back(static_cast<unsigned char>(n>>8));
            z.push_back(static_cast<unsigned char>(~n));
            z.push_back(static_cast<unsigned char>((~n)>>8));
            z.insert(z.end(),raw.begin()+pos,raw.begin()+pos+n);pos+=n;
        }
        std::uint32_t a=1,b=0;for(auto byte:raw){a=(a+byte)%65521u;b=(b+a)%65521u;}
        put32(z,(b<<16)|a);
        std::vector<unsigned char> png={137,80,78,71,13,10,26,10};
        std::vector<unsigned char> ihdr;put32(ihdr,width);put32(ihdr,height);
        ihdr.insert(ihdr.end(),{8,6,0,0,0});
        chunk(png,"IHDR",ihdr);chunk(png,"IDAT",z);chunk(png,"IEND",{});
        std::ofstream out(result.png_path,std::ios::binary);
        if(!out)throw std::runtime_error("Cannot write AO–MO PNG");
        out.write(reinterpret_cast<const char*>(png.data()),static_cast<std::streamsize>(png.size()));
        if(!out)throw std::runtime_error("AO–MO PNG write failed");
        result.png=true;
    }catch(const std::exception& error){result.error=error.what();}
    return result;
}

std::optional<std::size_t> suggest_initial_nbo_index(const NboIntegration& data,
    const MODiagramViewSnapshot* diagram) {
    if(!diagram || !diagram->data.view)return std::nullopt;
    const auto target=diagram->data.view->inspected_orbital_index;
    if(!target)return std::nullopt;
    const auto links=nbo_links_for_mo(data,*target,NboOrbitalKind::NBO);
    const NboMoLink* best=nullptr;
    for(const auto& link:links) {
        if(!link.weight)continue;
        const auto it=std::find_if(data.dataset.orbitals.begin(),data.dataset.orbitals.end(),
            [&](const auto& orbital){return orbital.id==link.orbital.index+1 &&
                orbital.spin==link.orbital.spin && orbital.occupation>0 &&
                orbital.kind!="CR" && orbital.kind!="CR*";});
        if(it==data.dataset.orbitals.end())continue;
        if(!best || *link.weight>*best->weight)best=&link;
    }
    if(!best)return std::nullopt;
    for(std::size_t i=0;i<data.dataset.orbitals.size();++i)
        if(data.dataset.orbitals[i].id==best->orbital.index+1 &&
           data.dataset.orbitals[i].spin==best->orbital.spin)return i;
    return std::nullopt;
}

} // namespace cov::ui
