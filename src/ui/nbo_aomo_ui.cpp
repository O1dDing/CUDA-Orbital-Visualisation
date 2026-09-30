#include "cov/nbo_aomo_ui.hpp"
#include "cov/nbo_aomo_hover.hpp"
#include "cov/nbo_aomo_text.hpp"
#include "cov/point_group_catalog.hpp"
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
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace cov::ui {
namespace {
void same_line_for(const char* next_label) {
    ImGui::SameLine();
    const float required=ImGui::CalcTextSize(next_label,nullptr,true).x+
        ImGui::GetFrameHeight()+ImGui::GetStyle().ItemInnerSpacing.x;
    if(ImGui::GetContentRegionAvail().x<required)ImGui::NewLine();
}
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
std::string number(double v) {std::ostringstream out;out<<std::setprecision(9)<<v;return out.str();}
std::string hex_rgb(const std::array<std::uint8_t,3>& rgb) {
    std::ostringstream out;out<<'#'<<std::hex<<std::setfill('0');
    for(const auto value:rgb)out<<std::setw(2)<<static_cast<unsigned>(value);
    return out.str();
}
std::string percent(double value) {
    std::ostringstream out;out<<std::fixed<<std::setprecision(3)<<value*100.0;
    return out.str();
}
std::array<std::string,3> coverage_lines(const NboAomoViewSnapshot& view) {
    if(!view.focused_display_weight)return {"","",""};
    return {
        aomo_text(view.language,"Lines show orbital contributions."),
        (!view.show_core || !view.show_rydberg)?aomo_text(view.language,
            "Some core or Rydberg orbitals are hidden."):std::string{},
        view.hide_h_orbitals?aomo_text(view.language,"H orbitals are hidden from this view."):
            std::string{}
    };
}
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
// A principal shell is accepted only when the producer actually names it.
std::string principal_shell(const std::string& type) {
    for(std::size_t i=0;i<type.size();++i) {
        if(!std::isdigit(static_cast<unsigned char>(type[i])))continue;
        const auto begin=i;
        while(i<type.size() && std::isdigit(static_cast<unsigned char>(type[i])))++i;
        if(i<type.size() && std::string("spdfgh").find(type[i])!=std::string::npos)
            return type.substr(begin,i-begin+1);
    }
    return {};
}
std::size_t verified_irrep_dimension(const std::string& group,const std::string& irrep) {
    const auto* definition=find_point_group(group);
    if(!definition)return 0;
    const auto normalized=[](const std::string& value) {
        std::string key;
        for(unsigned char c:value)if(!std::isspace(c) && c!='_')
            key+=static_cast<char>(std::tolower(c));
        return key;
    };
    const auto key=normalized(irrep);
    for(const auto& row:definition->irreps)
        if(normalized(std::string(row.label))==key)return row.dimension;
    return 0;
}
// The producer's Val(nd) and Val((n+1)s) rows identify a local metal
// interaction frame. A Ryd((n+1)p) row on that same atom and spin belongs to
// the frame even though its NBO class remains Rydberg.
using AtomSpin=std::pair<std::size_t,NboSpin>;
std::set<std::pair<AtomSpin,std::string>> framework_p_shells(const NboIntegration& data) {
    std::map<AtomSpin,std::set<std::string>> valence;
    for(const auto& nao:data.dataset.naos)if(contains_nocase(nao.type,"val")) {
        const auto shell=principal_shell(nao.type);
        if(!shell.empty())valence[{nao.atom,nao.spin}].insert(shell);
    }
    std::set<std::pair<AtomSpin,std::string>> result;
    for(const auto& [atom_spin,shells]:valence)
        for(const auto& shell:shells)if(shell.size()>1 && shell.back()=='d') {
            const auto n=std::stoi(shell.substr(0,shell.size()-1));
            if(shells.contains(std::to_string(n+1)+"s"))
                result.insert({atom_spin,std::to_string(n+1)+"p"});
        }
    return result;
}
template<class Paint>
void electron_strokes(const NboAomoNode& node,Paint paint) {
    if(!node.occupation_on_bar || node.occupation_label.empty())return;
    const bool paired=node.occupation_label=="↑↓";
    const float cy=node.occupation_y+node.occupation_height*0.5f;
    const float half=node.occupation_height*0.5f;
    const auto arrow=[&](float x,bool down) {
        const float tip=cy+(down?half:-half),tail=cy+(down?-half:half);
        paint(x,tail,x,tip);paint(x,tip,x-3,tip+(down?-4:4));
        paint(x,tip,x+3,tip+(down?-4:4));
    };
    const float cx=node.occupation_x+node.occupation_width*0.5f;
    arrow(cx-(paired?5.0f:0.0f),node.occupation_label=="↓");
    if(paired)arrow(cx+5.0f,true);
}
template<class Paint>
void dashed_segments(float ax,float ay,float bx,float by,Paint paint) {
    const float length=std::hypot(bx-ax,by-ay);
    if(length<0.01f)return;
    for(float d=0;d<length;d+=10.0f) {
        const float a=d/length,b=std::min(d+5.0f,length)/length;
        paint(ax+(bx-ax)*a,ay+(by-ay)*a,ax+(bx-ax)*b,ay+(by-ay)*b);
    }
}
template<class Paint>
void for_each_node_text_background(const NboAomoViewSnapshot& view,Paint paint) {
    for(const auto& node:view.nodes) {
        if(node.group_header)continue; // disclosure controls already have a fill
        const auto rectangle=[&](const char* role,float x,float y,float width,float height) {
            if(width<=0 || height<=0)return;
            const bool qualitative=node.lane!=NboAomoLane::Centre &&
                y>=view.qualitative_band_y;
            paint(node.id,role,x-2.0f,y-1.0f,width+4.0f,height+2.0f,qualitative);
        };
        rectangle("label",node.label_x,node.label_y,node.label_width,node.label_height);
        if(!node.occupation_label.empty())
            rectangle("occupation",node.occupation_x,node.occupation_y,
                node.occupation_width,node.occupation_height);
    }
}
std::vector<std::size_t> level_members(const MODiagramLevel& level,
    std::size_t canonical_count) {
    std::vector<std::size_t> members;
    const auto append=[&](std::size_t index) {
        // UDFT spatial rows use SIZE_MAX for an unmatched opposite-spin slot.
        // It is not a canonical member; likewise reject any out-of-range index.
        if(index<canonical_count)members.push_back(index);
    };
    if(level.member_indices.empty())append(level.metadata.orbital_index);
    else for(const auto index:level.member_indices)append(index);
    for(const auto index:level.member_spin_counterparts)append(index);
    std::sort(members.begin(),members.end());
    members.erase(std::unique(members.begin(),members.end()),members.end());
    return members;
}
std::vector<std::size_t> central_indices(const MODiagramViewSnapshot& diagram,
    std::size_t canonical_count) {
    std::vector<std::size_t> indices;
    for(const auto& level:diagram.data.levels) {
        const auto members=level_members(level,canonical_count);
        indices.insert(indices.end(),members.begin(),members.end());
    }
    std::sort(indices.begin(),indices.end());indices.erase(std::unique(indices.begin(),indices.end()),indices.end());
    return indices;
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

// Layout is derived from immutable orbital identities and the diagram's energy
// transform. In particular, changing the focused MO never redefines a side
// orbital or moves a quantitative energy level.
NboAomoViewSnapshot make_unified_snapshot(const NboAomoUIState& state,
    const NboIntegration& data,const Wavefunction& canonical,
    const MODiagramViewSnapshot& diagram,float available_width) {
    NboAomoViewSnapshot view;
    view.integration_id=data.id;
    view.language=state.language;
    view.mo_snapshot_id=diagram.data.view?diagram.data.view->id:"no-mo-snapshot";
    view.id=view.mo_snapshot_id+":unified:"+std::to_string(state.revision)+
        ":lang:"+std::to_string(static_cast<int>(state.language));
    view.basis_kind=state.basis_kind;
    view.preset=state.preset;view.overview=state.overview;
    view.illustrative_side_layout=state.preset==NboAomoPreset::Teaching &&
        state.illustrative_side_layout;
    view.paper_export=state.paper_export;
    view.salc_model=state.salc_model;
    view.mo_energy_axis_mode=energy_axis_mode_name(diagram.options.energy_axis_mode);
    view.display_energy_unit=energy_unit_symbol(diagram.options.energy_unit);
    view.energy_unit=diagram.options.energy_unit;
    view.mo_energy_axis_detail=view.illustrative_side_layout?
        "Canonical eigenvalues use the numeric transform; all side positions are explicitly illustrative, including sides with verified raw expectation energies":
        "Canonical eigenvalues and verified same-operator side expectations use the declared transform; other side levels are in a non-quantitative band";
    view.zoom=state.zoom;view.pan_x=state.pan_x;view.pan_y=state.pan_y;
    view.show_core=state.show_core;view.show_rydberg=state.show_rydberg;
    view.hide_h_orbitals=state.hide_h_orbitals;
    view.selection=state.selection;view.fragment_groups=state.fragment_groups;
    view.selected_side_node_id=state.selected_side_node_id;
    view.sum_component_ids.assign(state.sum_component_ids.begin(),state.sum_component_ids.end());
    if(const auto* cap=nbo_capability(data,"aomo")) {
        view.capability_status=nbo_capability_state_name(cap->state);
        view.capability_detail=cap->detail;
    }
    view.central_mo_indices=central_indices(diagram,canonical.orbitals.size());
    const auto inspected=diagram.data.view?diagram.data.view->inspected_orbital_index:std::nullopt;
    view.focused_canonical_index=state.focused_canonical_index && in(view.central_mo_indices,*state.focused_canonical_index)
        ?*state.focused_canonical_index:inspected && in(view.central_mo_indices,*inspected)
            ?*inspected:view.central_mo_indices.empty()?0:view.central_mo_indices.front();
    // Focusing a virtual MO never changes an explicit basis-class filter.
    const bool allow_rydberg=state.show_rydberg;
    const auto framework_p=framework_p_shells(data);
    const auto is_h_only=[&](const std::vector<std::size_t>& atoms) {
        return state.hide_h_orbitals && !atoms.empty() &&
            std::all_of(atoms.begin(),atoms.end(),[&](auto atom) {
                return atom<canonical.atoms.size() && canonical.atoms[atom].atomic_number==1;
            });
    };
    const auto keep_framework_p=[&](std::size_t zero_atom,NboSpin spin,
                                    const std::string& type) {
        return contains_nocase(type,"ryd") &&
            framework_p.contains({{zero_atom+1,spin},principal_shell(type)});
    };
    if(const auto* decomposition=nbo_mo_decomposition(data.dataset,view.focused_canonical_index)) {
        view.focused_projection_weight=decomposition->weight_sum;
        view.focused_projection_residual_norm=decomposition->projection_residual_norm;
    }
    std::size_t widest_level=1;
    for(const auto& level:diagram.data.levels)
        widest_level=std::max(widest_level,level_members(level,canonical.orbitals.size()).size());
    const float required_centre=62.0f*static_cast<float>(widest_level)+
        14.0f*static_cast<float>(widest_level-1);
    // Reserve the entire central manifold. Larger spin/degenerate groups make
    // the canvas horizontally scrollable instead of invading a side column.
    const float width=std::max({540.0f,available_width,
        widest_level>1?required_centre+600.0f:540.0f});
    view.canvas_width=width;
    // One canvas type size for labels, axis ticks, titles and supporting text.
    // DPI/view magnification is applied uniformly by the viewport transform.
    view.label_font_size=18.0f;
    const float row_gap=view.label_font_size+36.0f;
    const float numeric_top=78.0f;
    float numeric_span=state.preset==NboAomoPreset::Teaching?240.0f:
        state.preset==NboAomoPreset::Research?360.0f:480.0f;
    const float left_x=width<800?165.0f:195.0f;
    const float centre_x=width*0.5f-31.0f;
    const float right_x=width<800?width-165.0f:width-225.0f;
    std::vector<double> axis_energies;
    const auto include_energy=[&](double energy) {
        if(!std::isfinite(energy))return;
        axis_energies.push_back(energy);
    };
    for(const auto& level:diagram.data.levels)
        for(const auto index:level_members(level,canonical.orbitals.size()))
            if(index<canonical.orbitals.size())include_energy(canonical.orbitals[index].energy_hartree);
    if(!view.illustrative_side_layout && state.basis_kind==NboOrbitalKind::NAO &&
        state.salc_model && state.salc_model->available)
        for(const auto& orbital:state.salc_model->orbitals) {
            if(is_h_only(orbital.atoms) ||
                (contains_nocase(orbital.type,"cor")&&!state.show_core) ||
                (contains_nocase(orbital.type,"ryd")&&!allow_rydberg &&
                 !(orbital.atoms.size()==1 && keep_framework_p(
                    orbital.atoms.front(),orbital.spin,orbital.type))))continue;
            if(orbital.energy_hartree)include_energy(*orbital.energy_hartree);
        }
    // Build one transform from BOTH orbital definitions. Reusing a central-only
    // transform and extending its range to side expectations crushes the MO
    // manifold. Very close energies share a knot neighbourhood, not a false
    // energy value: interpolation still uses every node's unmodified energy.
    std::sort(axis_energies.begin(),axis_energies.end());
    std::vector<double> knot_energies;
    for(const auto energy:axis_energies)
        if(knot_energies.empty() || energy-knot_energies.back()>1e-7)
            knot_energies.push_back(energy);
    // A near-degenerate-only set still needs a reversible axis for its exact
    // distinct endpoints. A singleton knot is reserved for exactly one value.
    if(!axis_energies.empty() && knot_energies.size()==1 &&
        axis_energies.back()>axis_energies.front())
        knot_energies.push_back(axis_energies.back());
    view.energy_transform=build_energy_transform(knot_energies,
        diagram.options.energy_axis_mode,0.0);
    if(diagram.options.energy_axis_mode==EnergyAxisMode::NonlinearFocus) {
        view.energy_tick_screen_rgb={225,173,89};
        view.energy_tick_export_rgb=view.paper_export?
            std::array<std::uint8_t,3>{140,91,19}:view.energy_tick_screen_rgb;
        view.energy_tick_semantics="adaptive-nonlinear-amber";
    } else {
        view.energy_tick_export_rgb=view.paper_export?
            std::array<std::uint8_t,3>{81,101,121}:view.energy_tick_screen_rgb;
    }
    if(view.energy_transform.knots.size()==1 && !axis_energies.empty() &&
        axis_energies.back()>axis_energies.front()) {
        // The shared interpolation helper treats <=1e-12 Ha as coincident.
        // Bracket that extreme case instead of silently using a constant map.
        const double pad=axis_energies.back()-axis_energies.front()<=1e-12?1e-7:0;
        const double a=axis_energies.front()-pad,b=axis_energies.back()+pad;
        const bool linear=diagram.options.energy_axis_mode==EnergyAxisMode::Linear;
        view.energy_transform.knots={{a,linear?a:0.0},{b,linear?b:1.0}};
    }
    if(diagram.options.energy_axis_mode==EnergyAxisMode::NonlinearFocus &&
        view.energy_transform.knots.size()>1) {
        auto& knots=view.energy_transform.knots;
        struct AxisConstraint {std::size_t lower;double gap;};
        std::vector<std::vector<AxisConstraint>> incoming(knots.size());
        // Complete canonical groups are successive reading rows. Cross-lane
        // energy neighbours impose no artificial separation. Close side
        // members use horizontal slots, keeping one common numeric mapping.
        std::vector<std::pair<double,double>> central_bands;
        for(const auto& level:diagram.data.levels) {
            double first=std::numeric_limits<double>::infinity();
            double last=-first;
            for(const auto index:level_members(level,canonical.orbitals.size()))if(index<canonical.orbitals.size()) {
                first=std::min(first,canonical.orbitals[index].energy_hartree);
                last=std::max(last,canonical.orbitals[index].energy_hartree);
            }
            if(std::isfinite(first))central_bands.emplace_back(first,last);
        }
        std::sort(central_bands.begin(),central_bands.end());
        const auto knot_index=[&](double energy) {
            const auto it=std::upper_bound(knots.begin(),knots.end(),energy+1e-10,
                [](double e,const auto& knot){return e<knot.energy_hartree;});
            return it==knots.begin()?std::size_t{0}:
                static_cast<std::size_t>(std::distance(knots.begin(),it)-1);
        };
        for(std::size_t i=1;i<central_bands.size();++i) {
            const auto a=knot_index(central_bands[i-1].second);
            const auto b=knot_index(central_bands[i].first);
            if(b>a)incoming[b].push_back({a,row_gap});
        }
        const auto original=knots;
        knots.front().coordinate=0;
        for(std::size_t i=1;i<knots.size();++i) {
            // Tiny cross-column gaps remain visible but do not each consume a
            // text row. The DAG constraints reserve rows only where needed.
            knots[i].coordinate=knots[i-1].coordinate+std::max(1.5,
                (original[i].coordinate-original[i-1].coordinate)*180.0);
            for(const auto& constraint:incoming[i])
                knots[i].coordinate=std::max(knots[i].coordinate,
                    knots[constraint.lower].coordinate+constraint.gap);
        }
        numeric_span=std::max(numeric_span,static_cast<float>(knots.back().coordinate));
    }
    const auto& transform=view.energy_transform;
    double low=axis_energies.empty()?0:energy_display_coordinate(axis_energies.front(),transform);
    double high=axis_energies.empty()?1:energy_display_coordinate(axis_energies.back(),transform);
    if(high-low<1e-12){low-=0.5;high+=0.5;}
    if(diagram.options.energy_axis_mode==EnergyAxisMode::NonlinearFocus) {
        view.mo_energy_axis_mode="shared piecewise nonlinear (spacing is not an energy gap)";
    }
    view.mo_energy_axis_detail+=
        "; shared piecewise-linear interpolation of exported energy_transform knots; "
        "knot neighbourhood tolerance 1e-7 Ha, raw node energies unchanged; "
        "same-spin canonical groups use mean display energy; atomic shells use mean display energy plus a symmetric compact stack; "
        "side values are molecular operator expectations, central values are canonical eigenvalues";
    view.axis_coordinate_min=low;view.axis_coordinate_max=high;
    view.numeric_top=numeric_top;view.numeric_span=numeric_span;
    view.qualitative_band_y=numeric_top+numeric_span+66.0f;
    const int tick_count=std::max(5,static_cast<int>(numeric_span/110.0f));
    for(int i=0;i<=tick_count;++i) {
        const double coordinate=low+(high-low)*static_cast<double>(i)/tick_count;
        view.energy_ticks.push_back({energy_from_display_coordinate(coordinate,transform),
            numeric_top+static_cast<float>(1.0-static_cast<double>(i)/tick_count)*numeric_span});
    }
    const auto energy_y=[&](double energy) {
        const double coordinate=energy_display_coordinate(energy,transform);
        const double t=(coordinate-low)/(high-low);
        return numeric_top+static_cast<float>(1-t)*numeric_span;
    };
    const auto canonical_label=[&](std::size_t index) {
        const auto* name=state.names&&index<state.names->canonical.size()?&state.names->canonical[index]:nullptr;
        return canonical_mo_display_label(canonical,index,name);
    };
    std::map<std::size_t,std::size_t> mo_nodes;
    for(std::size_t li=0;li<diagram.data.levels.size();++li) {
        const auto& level=diagram.data.levels[li];
        const auto members=level_members(level,canonical.orbitals.size());
        if(members.empty())continue;
        std::map<std::size_t,std::pair<std::size_t,bool>> spatial_pairs;
        for(std::size_t j=0;j<level.member_indices.size() &&
            j<level.member_spin_counterparts.size();++j) {
            const auto primary=level.member_indices[j];
            const auto counterpart=level.member_spin_counterparts[j];
            if(primary>=canonical.orbitals.size() ||
                counterpart>=canonical.orbitals.size() || primary==counterpart ||
                canonical.orbitals[primary].spin==canonical.orbitals[counterpart].spin)
                continue;
            spatial_pairs[primary]={counterpart,true};
            spatial_pairs[counterpart]={primary,false};
        }
        std::map<Spin,std::pair<double,std::size_t>> spin_means;
        for(const auto index:members) {
            auto& mean=spin_means[canonical.orbitals[index].spin];
            mean.first+=canonical.orbitals[index].energy_hartree;++mean.second;
        }
        // Every real member stays visible and independently selectable. The
        // frozen layout later centres this complete row around the canvas.
        std::size_t mi=0;
        for(const auto index:members) {
            if(index>=canonical.orbitals.size())continue;
            const auto& mo=canonical.orbitals[index];
            NboAomoNode node;
            node.id="canonical_mo:"+std::to_string(index);
            node.label=canonical_label(index);
            node.individual_label=node.label;
            node.detail="Gaussian canonical eigenvalue; occupation="+number(mo.occupation);
            if(const auto pair=spatial_pairs.find(index);pair!=spatial_pairs.end()) {
                const auto other=pair->second.first;
                node.spatial_pair_id="spatial-pair:"+
                    std::to_string(std::min(index,other))+":"+
                    std::to_string(std::max(index,other));
                const auto& opposite=canonical.orbitals[other];
                const bool occupied=mo.occupation>0.5;
                const bool opposite_occupied=opposite.occupation>0.5;
                const std::string occupancy=occupied && opposite_occupied?" ↑↓":
                    occupied?(mo.spin==Spin::Beta?" ↓":" ↑"):
                    opposite_occupied?(opposite.spin==Spin::Beta?" ↓":" ↑"):"";
                const auto primary=pair->second.second?index:other;
                const auto partner=pair->second.second?other:index;
                const auto strip_spin=[](std::string label) {
                    for(const auto suffix:{" [alpha]"," [beta]"}) {
                        const std::string ending=suffix;
                        if(label.size()>=ending.size() &&
                            label.compare(label.size()-ending.size(),ending.size(),ending)==0) {
                            label.resize(label.size()-ending.size());break;
                        }
                    }
                    return label;
                };
                const auto primary_name=strip_spin(canonical_label(primary));
                const auto partner_name=strip_spin(canonical_label(partner));
                if(primary_name==partner_name)
                    node.spatial_pair_label=primary_name+" α/β"+occupancy;
                else node.spatial_pair_label=canonical_label(primary)+" / "+canonical_label(partner)+occupancy;
                if(pair->second.second)node.label=node.spatial_pair_label;
                else node.label.clear();
                node.detail+="; matched alpha/beta spatial display pair with "+
                    canonical_mo_source_label(canonical,other)+
                    "; each member retains its own canonical eigenvalue and selection. "
                    "Within a degenerate subspace, this display correspondence need not be unique.";
            }
            if(state.names && index<state.names->canonical.size()) {
                const auto& name=state.names->canonical[index];
                node.symmetry_irrep=name.irrep;node.symmetry_ordinal=name.ordinal;
                node.symmetry_name_verified=name.verified;node.name_detail=name.detail;
                node.detail+="; original MO "+std::to_string(index+1)+"; "+name.detail;
            }
            node.energy_semantics="canonical eigenvalue";
            node.orbital=canonical_ref(data,index);
            node.available=node.orbital.has_value();
            node.canonical_index=index;node.energy_hartree=mo.energy_hartree;
            const auto mean=spin_means.at(mo.spin);
            node.display_energy_hartree=mean.first/static_cast<double>(mean.second);
            node.display_energy_semantics=mean.second>1?"same-spin canonical group mean":"individual canonical eigenvalue";
            node.display_group_id="canonical-level:"+std::to_string(li)+
                (mo.spin==Spin::Beta?":beta":":alpha");
            node.detail+="; diagram energy="+number(*node.display_energy_hartree)+
                " Ha ("+node.display_energy_semantics+"); original eigenvalue retained";
            node.subspace_id="canonical-level:"+std::to_string(li);
            node.occupation=mo.occupation;node.quantitative_energy=true;
            node.lane=NboAomoLane::Centre;
            node.x=centre_x+static_cast<float>(mi++);
            node.y=energy_y(*node.display_energy_hartree);node.width=62;node.height=22;
            mo_nodes[index]=view.nodes.size();
            view.nodes.push_back(std::move(node));
        }
    }
    std::map<RefKey,std::size_t> basis_nodes;
    std::set<std::size_t> folded_atoms;
    std::vector<std::size_t> side_nodes;
    const auto add_side=[&](NboAomoNode node) {
        const auto index=view.nodes.size();
        side_nodes.push_back(index);view.nodes.push_back(std::move(node));
        return index;
    };
    const bool use_salc=state.basis_kind==NboOrbitalKind::NAO &&
        state.salc_model && state.salc_model->available;
    std::vector<std::size_t> salc_nodes(use_salc?state.salc_model->orbitals.size():0,
        std::numeric_limits<std::size_t>::max());
    std::vector<bool> class_filtered(salc_nodes.size(),false),group_folded(salc_nodes.size(),false);
    std::map<std::string,std::vector<std::size_t>> subspace_members;
    std::set<std::string> expanded;
    if(use_salc) {
        std::map<std::string,double> subspace_weight;
        std::size_t eligible_count=0;
        for(std::size_t i=0;i<state.salc_model->orbitals.size();++i) {
            const auto& orbital=state.salc_model->orbitals[i];
            class_filtered[i]=is_h_only(orbital.atoms) ||
                (contains_nocase(orbital.type,"cor")&&!state.show_core) ||
                (contains_nocase(orbital.type,"ryd")&&!allow_rydberg &&
                 !(orbital.atoms.size()==1 && keep_framework_p(
                    orbital.atoms.front(),orbital.spin,orbital.type)));
            if(class_filtered[i])continue;
            subspace_members[orbital.subspace_id].push_back(i);++eligible_count;
        }
        for(const auto& link:state.salc_model->links)
            if(link.canonical_index==view.focused_canonical_index &&
                link.side_index<class_filtered.size() && !class_filtered[link.side_index])
                subspace_weight[state.salc_model->orbitals[link.side_index].subspace_id]
                    +=link.weight;
        std::vector<std::string> groups;
        for(const auto& [id,_]:subspace_members)groups.push_back(id);
        std::stable_sort(groups.begin(),groups.end(),[&](const auto& a,const auto& b){
            return subspace_weight[a]>subspace_weight[b];
        });
        // The graph shows every eligible verified side orbital by default;
        // only an explicit user fold may replace a subspace with a header.
        const std::size_t node_budget=std::numeric_limits<std::size_t>::max();
        std::size_t shown=0;
        for(const auto& id:groups) {
            const bool manual_open=state.expanded_subspaces.contains(id);
            const bool manual_fold=state.collapsed_subspaces.contains(id);
            const bool open=!manual_fold && (manual_open ||
                eligible_count<=node_budget || shown==0 ||
                shown+subspace_members[id].size()<=node_budget);
            if(open){expanded.insert(id);shown+=subspace_members[id].size();}
            else for(const auto index:subspace_members[id])group_folded[index]=true;
        }
        using SalcDisplayKey=std::tuple<std::string,std::string,std::size_t,NboSpin>;
        std::map<SalcDisplayKey,std::vector<std::size_t>> named_partners;
        if(state.names)
            for(std::size_t index=0;index<state.salc_model->orbitals.size() &&
                index<state.names->salc.size();++index) {
                const auto& name=state.names->salc[index];
                const auto& member=state.salc_model->orbitals[index];
                if(class_filtered[index] || group_folded[index] || !name.verified ||
                    name.ordinal==0 || name.irrep.empty())continue;
                named_partners[{member.subspace_id,name.irrep,name.ordinal,member.spin}]
                    .push_back(index);
            }
        for(std::size_t i=0;i<state.salc_model->orbitals.size();++i) {
            const auto& orbital=state.salc_model->orbitals[i];
            if(class_filtered[i]) {
                if(is_h_only(orbital.atoms))++view.hidden_h_count;
                else ++view.hidden_class_count;
                continue;
            }
            if(group_folded[i]){++view.hidden_group_count;continue;}
            const auto fragment=std::find_if(state.salc_model->fragments.begin(),
                state.salc_model->fragments.end(),[&](const auto& f){return f.id==orbital.fragment_id;});
            NboAomoNode node;
            node.id="salc:"+orbital.id;node.label=orbital.label;
            node.detail=orbital.label+"; "+orbital.detail;
            node.energy_semantics=orbital.energy_semantics;
            if(orbital.terms.size()==1 && !orbital.atoms.empty() &&
                orbital.atoms.front()<canonical.atoms.size()) {
                std::string angular=orbital.angular;
                if(contains_nocase(angular,"px"))angular="2px";
                else if(contains_nocase(angular,"py"))angular="2py";
                else if(contains_nocase(angular,"pz"))angular="2pz";
                else if(contains_nocase(angular,"2s") ||
                    contains_nocase(orbital.label,"2s"))angular="2s";
                node.label=canonical.atoms[orbital.atoms.front()].symbol+
                    std::to_string(orbital.atoms.front()+1)+" "+angular;
            }
            else if(orbital.symmetry_adapted && orbital.atoms.size()==2 &&
                orbital.terms.size()==2) {
                const bool in_phase=orbital.terms[0].coefficient*
                    orbital.terms[1].coefficient>0;
                const bool hydrogens=std::all_of(orbital.atoms.begin(),
                    orbital.atoms.end(),[&](auto atom){return atom<canonical.atoms.size() &&
                        canonical.atoms[atom].atomic_number==1;});
                node.label=(hydrogens?"H pair":"SALC")+
                    std::string(in_phase?" +":" -");
            } else if(orbital.terms.size()>1) {
                const bool same_element=!orbital.atoms.empty() &&
                    std::all_of(orbital.atoms.begin(),orbital.atoms.end(),[&](auto atom){
                        return atom<canonical.atoms.size() &&
                            canonical.atoms[atom].atomic_number==
                                canonical.atoms[orbital.atoms.front()].atomic_number;});
                const std::string family=same_element?
                    canonical.atoms[orbital.atoms.front()].symbol+
                        std::to_string(orbital.atoms.size()):"SALC";
                const std::string shell=contains_nocase(orbital.type,"1s")?"1s":
                    contains_nocase(orbital.type,"2s")?"2s":
                    contains_nocase(orbital.type,"2p")?"2p":orbital.angular;
                const auto space=orbital.subspace_id.rfind("space");
                node.label=family+" "+shell+" S"+
                    (space==std::string::npos?std::to_string(i+1):
                        orbital.subspace_id.substr(space+5));
            }
            node.salc_index=i;node.subspace_id=orbital.subspace_id;
            node.atoms=orbital.atoms;
            node.energy_hartree=orbital.energy_hartree;node.occupation=orbital.occupation;
            node.quantitative_energy=orbital.energy_hartree.has_value() &&
                !view.illustrative_side_layout;
            node.lane=fragment!=state.salc_model->fragments.end() && fragment->side?
                NboAomoLane::Right:NboAomoLane::Left;
            node.x=node.lane==NboAomoLane::Right?right_x:left_x;
            node.width=width<800?60.0f:94.0f;node.height=22;
            if(orbital.partner_dimension>1) {
                node.x+=static_cast<float>(orbital.partner_index)*32.0f;
                node.width=25.0f;
                const bool same_element=!orbital.atoms.empty() &&
                    std::all_of(orbital.atoms.begin(),orbital.atoms.end(),[&](auto atom){
                        return atom<canonical.atoms.size() &&
                            canonical.atoms[atom].atomic_number==
                                canonical.atoms[orbital.atoms.front()].atomic_number;});
                const std::string family=same_element?
                    canonical.atoms[orbital.atoms.front()].symbol+
                        std::to_string(orbital.atoms.size()):"SALC";
                const std::string shell=contains_nocase(orbital.type,"1s")?"1s":
                    contains_nocase(orbital.type,"2s")?"2s":
                    contains_nocase(orbital.type,"2p")?"2p":orbital.angular;
                node.label=family+" "+shell+" #"+
                    std::to_string(orbital.partner_index+1);
            }
            if(orbital.spin!=NboSpin::Total)
                node.label+=" ["+std::string(nbo_spin_name(orbital.spin))+"]";
            if(orbital.atoms.size()>1 && state.names && i<state.names->salc.size()) {
                const auto& name=state.names->salc[i];
                node.label=name.verified&&name.ordinal?name.label:
                    std::string(orbital.symmetry_adapted?"SALC ":"Fragment orbital ")+
                    std::to_string(i+1)+" ["+nbo_spin_name(orbital.spin)+"]";
                node.symmetry_irrep=name.irrep;node.symmetry_ordinal=name.ordinal;
                node.symmetry_name_verified=name.verified;node.name_detail=name.detail;
                node.detail+="; "+name.detail;
            }
            if(orbital.atoms.size()==1 && orbital.atoms.front()<canonical.atoms.size()) {
                const auto shell=principal_shell(orbital.type);
                if(!shell.empty()) {
                    node.display_group_id="atomic-shell:"+std::to_string(orbital.atoms.front())+
                        ":"+nbo_spin_name(orbital.spin)+":"+shell;
                    node.shell_label=canonical.atoms[orbital.atoms.front()].symbol+
                        std::to_string(orbital.atoms.front()+1)+" "+shell;
                    if(orbital.spin!=NboSpin::Total)node.shell_label+=
                        " ["+std::string(nbo_spin_name(orbital.spin))+"]";
                }
            }
            else if(orbital.atoms.size()>1 && !orbital.subspace_id.empty()) {
                const auto subspace=std::find_if(state.salc_model->subspaces.begin(),
                    state.salc_model->subspaces.end(),[&](const auto& space) {
                        return space.id==orbital.subspace_id;
                    });
                bool shared_name=state.names && i<state.names->salc.size() &&
                    subspace!=state.salc_model->subspaces.end() &&
                    subspace->symmetry_verified;
                std::size_t dimension=0;
                if(shared_name) {
                    const auto& name=state.names->salc[i];
                    dimension=verified_irrep_dimension(state.salc_model->point_group,name.irrep);
                    const auto key=SalcDisplayKey{orbital.subspace_id,name.irrep,
                        name.ordinal,orbital.spin};
                    const auto it=named_partners.find(key);
                    shared_name=name.verified && name.ordinal>0 && dimension>0 &&
                        it!=named_partners.end() && it->second.size()==dimension;
                    if(shared_name)for(const auto member:it->second)
                        if(state.names->salc[member].label!=name.label) {
                            shared_name=false;break;
                        }
                }
                if(subspace!=state.salc_model->subspaces.end() &&
                    subspace->symmetry_verified && subspace->spin==orbital.spin &&
                    shared_name) {
                    const auto& name=state.names->salc[i];
                    node.display_group_id="verified-salc:"+orbital.subspace_id+
                        ":"+name.irrep+":"+std::to_string(name.ordinal)+
                        ":"+nbo_spin_name(orbital.spin);
                    node.shell_label=name.label;
                }
            }
            salc_nodes[i]=add_side(std::move(node));
        }
        struct FoldedFamily {std::size_t members=0,subspaces=0,first=0;double weight=0;};
        std::map<std::string,FoldedFamily> folded_families;
        std::map<std::size_t,FoldedFamily> folded_user_families;
        for(const auto& [id,members]:subspace_members)if(!expanded.contains(id) &&
            !members.empty()) {
            const auto& first=state.salc_model->orbitals[members.front()];
            const auto fragment=std::find_if(state.salc_model->fragments.begin(),
                state.salc_model->fragments.end(),[&](const auto& f){return f.id==first.fragment_id;});
            double weight=0;
            for(const auto index:members)for(const auto& link:state.salc_model->links)
                if(link.side_index==index &&
                    link.canonical_index==view.focused_canonical_index)weight+=link.weight;
            const auto user_group=std::find_if(state.fragment_groups.begin(),
                state.fragment_groups.end(),[&](const auto& group){
                    return !first.atoms.empty() &&
                        std::all_of(first.atoms.begin(),first.atoms.end(),[&](auto atom){
                            return group.atoms.contains(atom);});});
            if(user_group!=state.fragment_groups.end() &&
                !state.expanded_user_fragments.contains(user_group->id)) {
                auto& family=folded_user_families[user_group->id];
                if(family.subspaces==0)family.first=members.front();
                family.members+=members.size();++family.subspaces;
                family.weight+=weight;
                continue;
            }
            if(!state.expanded_fragments.contains(first.fragment_id)) {
                auto& family=folded_families[first.fragment_id];
                if(family.subspaces==0)family.first=members.front();
                family.members+=members.size();++family.subspaces;
                family.weight+=weight;
                continue;
            }
            NboAomoNode header;
            header.id="salc-group:"+id;header.subspace_id=id;
            header.label="+ "+first.type+" ("+std::to_string(members.size())+")";
            header.detail=(fragment!=state.salc_model->fragments.end()?
                fragment->label:std::string("fragment"))+": collapsed subspace, not an orbital; click to show actual members. "
                "Focused raw orthogonal weight="+percent(weight)+"%.";
            header.group_header=true;header.composition_available=false;
            header.lane=fragment!=state.salc_model->fragments.end() && fragment->side?
                NboAomoLane::Right:NboAomoLane::Left;
            header.x=header.lane==NboAomoLane::Right?right_x:left_x;
            header.width=width<800?100.0f:140.0f;header.height=22;
            add_side(std::move(header));
        }
        for(const auto& [id,family]:folded_families) {
            const auto& first=state.salc_model->orbitals[family.first];
            const auto fragment=std::find_if(state.salc_model->fragments.begin(),
                state.salc_model->fragments.end(),[&](const auto& f){return f.id==id;});
            const bool same_element=!first.atoms.empty() &&
                std::all_of(first.atoms.begin(),first.atoms.end(),[&](auto atom){
                    return atom<canonical.atoms.size() &&
                        canonical.atoms[atom].atomic_number==
                            canonical.atoms[first.atoms.front()].atomic_number;});
            const std::string family_label=same_element?
                canonical.atoms[first.atoms.front()].symbol+
                    std::to_string(first.atoms.size()):
                fragment!=state.salc_model->fragments.end()?
                    fragment->label:std::string("atom family");
            NboAomoNode header;
            header.id="fragment-group:"+id;header.group_header=true;
            header.label="+ "+family_label+" ("+
                std::to_string(family.subspaces)+" "+aomo_text(state.language,"spaces")+")";
            header.detail=(fragment!=state.salc_model->fragments.end()?
                fragment->label:family_label)+
                ": folded atom/equivalent-atom family, not one orbital or an asserted chemical bond. "
                "Click to list complete subspaces. Focused raw orthogonal weight="+
                percent(family.weight)+"%.";
            header.atoms=first.atoms;header.composition_available=false;
            header.lane=fragment!=state.salc_model->fragments.end() && fragment->side?
                NboAomoLane::Right:NboAomoLane::Left;
            header.x=header.lane==NboAomoLane::Right?right_x:left_x;
            header.width=width<800?100.0f:140.0f;header.height=22;
            add_side(std::move(header));
        }
        for(const auto& [id,family]:folded_user_families) {
            const auto group=std::find_if(state.fragment_groups.begin(),
                state.fragment_groups.end(),[&](const auto& g){return g.id==id;});
            if(group==state.fragment_groups.end())continue;
            const auto& first=state.salc_model->orbitals[family.first];
            NboAomoNode header;
            header.id="user-fragment:"+std::to_string(id);header.group_header=true;
            header.label="+ F"+std::to_string(id)+" ("+
                std::to_string(family.subspaces)+" "+aomo_text(state.language,"spaces")+")";
            header.detail="User-defined atom set used for folding, not a fixed orbital or SALC. "
                "Click to list constituent atom families; focused raw orthogonal weight="+
                percent(family.weight)+"%. The separate partial-sum command is MO-dependent.";
            header.atoms.assign(group->atoms.begin(),group->atoms.end());
            header.composition_available=false;
            const auto fragment=std::find_if(state.salc_model->fragments.begin(),
                state.salc_model->fragments.end(),[&](const auto& f){
                    return f.id==first.fragment_id;});
            header.lane=fragment!=state.salc_model->fragments.end() && fragment->side?
                NboAomoLane::Right:NboAomoLane::Left;
            header.x=header.lane==NboAomoLane::Right?right_x:left_x;
            header.width=width<800?100.0f:140.0f;header.height=22;
            add_side(std::move(header));
        }
    } else {
        std::set<std::size_t> related;
        for(const auto& link:data.links)if(link.orbital.kind==state.basis_kind &&
            in(view.central_mo_indices,link.canonical_index))related.insert(link.orbital.index);
        std::size_t primary_atom=canonical.atoms.size();
        for(const auto& orbital:data.orbitals)if(orbital.ref.kind==state.basis_kind &&
            related.contains(orbital.ref.index) && !orbital.atoms.empty())
            primary_atom=std::min(primary_atom,orbital.atoms.front());
        const auto atom_of=[&](const NboOrbitalDescriptor& orbital){
            return orbital.atoms.empty()?canonical.atoms.size():orbital.atoms.front();
        };
        std::map<std::size_t,std::size_t> atom_counts;
        std::map<std::size_t,double> atom_strength;
        for(const auto& orbital:data.orbitals)if(orbital.ref.kind==state.basis_kind &&
            related.contains(orbital.ref.index)) {
            const bool core=contains_nocase(orbital.label," cor");
            const bool ryd=contains_nocase(orbital.label," ryd");
            const bool frame_p=orbital.ref.kind==NboOrbitalKind::NAO &&
                orbital.atoms.size()==1 && [&] {
                    const auto row=std::find_if(data.dataset.naos.begin(),data.dataset.naos.end(),
                        [&](const auto& nao) {return nao.id==orbital.ref.index+1 &&
                            nao.spin==orbital.ref.spin;});
                    return row!=data.dataset.naos.end() && keep_framework_p(
                        orbital.atoms.front(),row->spin,row->type);
                }();
            if(is_h_only(orbital.atoms) || (core&&!state.show_core) ||
               (ryd&&!state.show_rydberg&&!frame_p))continue;
            ++atom_counts[atom_of(orbital)];
        }
        for(const auto& link:data.links)if(link.orbital.kind==state.basis_kind &&
            link.canonical_index==view.focused_canonical_index) {
            const auto* orbital=nbo_orbital(data,link.orbital);
            if(!orbital)continue;
            double strength=link.weight.value_or(0);
            if(state.basis_kind==NboOrbitalKind::GaussianAO &&
                link.canonical_index<canonical.orbitals.size()) {
                const auto& coefficients=canonical.orbitals[link.canonical_index].coefficients;
                const auto n=coefficients.size(),row=link.orbital.index;
                if(row<n && canonical.ao_overlap.size()==n*n &&
                    canonical.ao_overlap[row*n+row]>0) {
                    double projection=0;
                    for(std::size_t j=0;j<n;++j)
                        projection+=canonical.ao_overlap[row*n+j]*coefficients[j];
                    strength=projection*projection/canonical.ao_overlap[row*n+row];
                }
            }
            atom_strength[atom_of(*orbital)]+=strength;
        }
        std::vector<std::size_t> atom_order;
        for(const auto& [atom,_]:atom_counts)atom_order.push_back(atom);
        std::stable_sort(atom_order.begin(),atom_order.end(),[&](auto a,auto b){
            return atom_strength[a]>atom_strength[b];
        });
        const std::size_t node_budget=std::numeric_limits<std::size_t>::max();
        std::size_t shown=0;
        std::set<std::size_t> expanded_atoms;
        for(const auto atom:atom_order) {
            const bool open=!state.collapsed_atoms.contains(atom) &&
                (state.expanded_atoms.contains(atom)||shown==0||
                    shown+atom_counts[atom]<=node_budget);
            if(open){expanded_atoms.insert(atom);shown+=atom_counts[atom];}
            else folded_atoms.insert(atom);
        }
        for(const auto& orbital:data.orbitals) {
            if(orbital.ref.kind!=state.basis_kind || !related.contains(orbital.ref.index))continue;
            const bool core=contains_nocase(orbital.label," cor");
            const bool ryd=contains_nocase(orbital.label," ryd");
            const bool frame_p=orbital.ref.kind==NboOrbitalKind::NAO &&
                orbital.atoms.size()==1 && [&] {
                    const auto row=std::find_if(data.dataset.naos.begin(),data.dataset.naos.end(),
                        [&](const auto& nao) {return nao.id==orbital.ref.index+1 &&
                            nao.spin==orbital.ref.spin;});
                    return row!=data.dataset.naos.end() && keep_framework_p(
                        orbital.atoms.front(),row->spin,row->type);
                }();
            if(is_h_only(orbital.atoms) || (core&&!state.show_core) ||
               (ryd&&!state.show_rydberg&&!frame_p)){
                if(is_h_only(orbital.atoms))++view.hidden_h_count;
                else ++view.hidden_class_count;
                continue;
            }
            if(folded_atoms.contains(atom_of(orbital))) {
                ++view.hidden_basis_count;continue;
            }
            NboAomoNode node;
            node.id=orbital.id.empty()?row_id(orbital.ref):orbital.id;
            node.label=orbital.label.empty()?row_id(orbital.ref):orbital.label;
            // Atom identity comes from the canonical structure, while NAO
            // angular identity is matched by both producer number and spin.
            if(orbital.atoms.size()==1 && orbital.atoms.front()<canonical.atoms.size()) {
                const auto atom=orbital.atoms.front();
                node.label=canonical.atoms[atom].symbol+std::to_string(atom+1);
                if(orbital.ref.kind==NboOrbitalKind::NAO) {
                    const auto row=std::find_if(data.dataset.naos.begin(),data.dataset.naos.end(),
                        [&](const auto& nao){return nao.id==orbital.ref.index+1 && nao.spin==orbital.ref.spin;});
                    if(row!=data.dataset.naos.end() && !row->angular.empty())node.label+=" "+row->angular;
                    node.label+=" / NAO ";
                } else node.label+=" / AO ";
                node.label+=std::to_string(orbital.ref.index+1);
                const auto shell=orbital.label.find(" / shell ");
                if(orbital.ref.kind==NboOrbitalKind::GaussianAO && shell!=std::string::npos)
                    node.label+=orbital.label.substr(shell);
            }
            if(orbital.ref.spin!=NboSpin::Total)
                node.label+=" ["+std::string(nbo_spin_name(orbital.ref.spin))+"]";
            node.detail=orbital.label+"; "+orbital.detail;node.energy_semantics=orbital.energy_semantics;
            node.orbital=orbital.ref;node.atoms=orbital.atoms;
            if(orbital.atoms.size()==1 && orbital.atoms.front()<canonical.atoms.size()) {
                std::string shell,shell_key;
                if(orbital.ref.kind==NboOrbitalKind::NAO) {
                    const auto row=std::find_if(data.dataset.naos.begin(),data.dataset.naos.end(),
                        [&](const auto& nao){return nao.id==orbital.ref.index+1 && nao.spin==orbital.ref.spin;});
                    if(row!=data.dataset.naos.end()) {shell=principal_shell(row->type);shell_key=shell;}
                } else {
                    // Gaussian shell identity is a contraction index, not a
                    // hydrogenic principal quantum number: never invent n.
                    for(std::size_t si=0;si<canonical.shells.size();++si) {
                        const auto& sh=canonical.shells[si];
                        const auto count=sh.pure?2*sh.angular_momentum+1:
                            (sh.angular_momentum+1)*(sh.angular_momentum+2)/2;
                        for(std::size_t bi=sh.basis_offset;bi<sh.basis_offset+count &&
                            bi<canonical.gaussian_ao_transform.size();++bi)
                            if(canonical.gaussian_ao_transform[bi].source_index==orbital.ref.index &&
                                sh.angular_momentum<6) {
                                shell=std::string(1,"spdfgh"[sh.angular_momentum])+" (shell "+std::to_string(si+1)+")";
                                shell_key="gaussian:"+std::to_string(si);
                            }
                    }
                }
                if(!shell.empty()) {
                    node.display_group_id="atomic-shell:"+std::to_string(orbital.atoms.front())+
                        ":"+nbo_spin_name(orbital.ref.spin)+":"+shell_key;
                    node.shell_label=canonical.atoms[orbital.atoms.front()].symbol+
                        std::to_string(orbital.atoms.front()+1)+" "+shell;
                    if(orbital.ref.spin!=NboSpin::Total)node.shell_label+=
                        " ["+std::string(nbo_spin_name(orbital.ref.spin))+"]";
                }
            }
            node.energy_hartree.reset(); // unverified side Fock rows are not axis values
            node.occupation=orbital.occupation;
            node.lane=orbital.atoms.empty()||orbital.atoms.front()==primary_atom?
                NboAomoLane::Left:NboAomoLane::Right;
            node.x=node.lane==NboAomoLane::Right?right_x:left_x;
            node.width=width<800?60.0f:94.0f;node.height=22;
            basis_nodes[key(orbital.ref)]=add_side(std::move(node));
        }
        for(const auto atom:folded_atoms) {
            NboAomoNode header;
            header.id="atom:"+std::to_string(atom);header.group_header=true;
            header.label="+ "+(atom<canonical.atoms.size()?
                canonical.atoms[atom].symbol+std::to_string(atom+1):
                std::string("Other"))+" ("+std::to_string(atom_counts[atom])+")";
            header.detail="Collapsed atom group, not an orbital; click to show all actual basis orbitals. "
                "AO projection strengths are non-additive.";
            header.atoms={atom};header.composition_available=false;
            header.lane=atom==primary_atom?NboAomoLane::Left:NboAomoLane::Right;
            header.x=header.lane==NboAomoLane::Right?right_x:left_x;
            header.width=width<800?100.0f:140.0f;header.height=22;
            add_side(std::move(header));
        }
    }
    std::size_t qualitative_left=0,qualitative_right=0;
    for(const auto index:side_nodes) {
        auto& node=view.nodes[index];
        if(node.quantitative_energy && node.energy_hartree)node.y=energy_y(*node.energy_hartree);
        else {
            auto& row=node.lane==NboAomoLane::Right?qualitative_right:qualitative_left;
            node.y=view.qualitative_band_y+38.0f+static_cast<float>(row++)*(row_gap+8.0f);
        }
    }
    if(use_salc) {
        for(std::size_t i=0;i<state.salc_model->links.size();++i) {
            const auto& link=state.salc_model->links[i];
            if(link.side_index>=salc_nodes.size() ||
                salc_nodes[link.side_index]==std::numeric_limits<std::size_t>::max())continue;
            const auto target=mo_nodes.find(link.canonical_index);
            if(target==mo_nodes.end())continue;
            NboAomoEdge edge;
            edge.id="salc-component:"+state.salc_model->orbitals[link.side_index].id+
                ":mo:"+std::to_string(link.canonical_index);
            edge.source_node=salc_nodes[link.side_index];edge.target_node=target->second;
            edge.salc_link_index=i;edge.coefficient=link.coefficient;edge.weight=link.weight;
            view.edges.push_back(std::move(edge));
        }
    } else {
        for(const auto& link:data.links) {
            const auto source=basis_nodes.find(key(link.orbital));
            const auto target=mo_nodes.find(link.canonical_index);
            if(source==basis_nodes.end()||target==mo_nodes.end())continue;
            NboAomoEdge edge;
            edge.id="component:"+row_id(link.orbital)+":mo:"+std::to_string(link.canonical_index);
            edge.source_node=source->second;edge.target_node=target->second;
            edge.coefficient=link.coefficient;edge.weight=link.weight;edge.source=link.source;
            if(state.basis_kind==NboOrbitalKind::GaussianAO &&
                link.canonical_index<canonical.orbitals.size()) {
                const auto n=canonical.orbitals[link.canonical_index].coefficients.size();
                if(link.orbital.index<n && canonical.ao_overlap.size()==n*n) {
                    const auto row=link.orbital.index;
                    const double norm=canonical.ao_overlap[row*n+row];
                    if(norm>0) {
                        double projection=0;
                        for(std::size_t j=0;j<n;++j)
                            projection+=canonical.ao_overlap[row*n+j]*
                                canonical.orbitals[link.canonical_index].coefficients[j];
                        edge.projection_strength_nonadditive=projection*projection/norm;
                    }
                }
            }
            view.edges.push_back(std::move(edge));
        }
    }
    // A measured matrix/representation error informs numerical-zero handling.
    // This is a display bound, never a claim that a tiny chemical interaction
    // is symmetry-forbidden. Full raw coefficients remain in export/detail.
    double measured_error=0;
    if(state.salc_model)measured_error=std::max(
        state.salc_model->orthogonality_error,state.salc_model->representation_error);
    const double zero=std::max(1e-10,std::min(1e-8,10*measured_error));
    view.numerical_zero_bound=zero;
    view.numerical_zero_reason="max(1e-10 floating-point display floor, capped measured NAO/symmetry error x10); raw links retained";
    std::set<std::size_t> focus_members{view.focused_canonical_index};
    std::map<std::size_t,std::size_t> mo_level;
    for(std::size_t li=0;li<diagram.data.levels.size();++li) {
        const auto members=level_members(diagram.data.levels[li],canonical.orbitals.size());
        for(auto mo:members)mo_level[mo]=li;
        if(in(members,view.focused_canonical_index))
            focus_members.insert(members.begin(),members.end());
    }
    // Account for class suppression separately from reading suppression and
    // from the local-space residual. These are orthogonal-basis weights only.
    double hidden_core=0,hidden_ryd=0,hidden_valence=0,hidden_group=0,hidden_h=0;
    if(use_salc) {
        for(const auto& link:state.salc_model->links) {
            if(link.canonical_index!=view.focused_canonical_index ||
                link.side_index>=salc_nodes.size() ||
                salc_nodes[link.side_index]!=std::numeric_limits<std::size_t>::max())
                continue;
            if(group_folded[link.side_index]){hidden_group+=link.weight;continue;}
            const auto& side=state.salc_model->orbitals[link.side_index];
            if(is_h_only(side.atoms)){hidden_h+=link.weight;continue;}
            const auto& type=side.type;
            if(contains_nocase(type,"cor"))hidden_core+=link.weight;
            else if(contains_nocase(type,"ryd"))hidden_ryd+=link.weight;
            else hidden_valence+=link.weight;
        }
    } else if(state.basis_kind==NboOrbitalKind::NAO) {
        for(const auto& link:data.links)if(link.orbital.kind==NboOrbitalKind::NAO &&
            link.canonical_index==view.focused_canonical_index && link.weight &&
            !basis_nodes.contains(key(link.orbital))) {
            const auto* orbital=nbo_orbital(data,link.orbital);
            if(orbital && !orbital->atoms.empty() &&
                folded_atoms.contains(orbital->atoms.front())) {
                hidden_group+=*link.weight;continue;
            }
            if(orbital && is_h_only(orbital->atoms)){
                hidden_h+=*link.weight;continue;
            }
            if(orbital && contains_nocase(orbital->label,"cor"))hidden_core+=*link.weight;
            else if(orbital && contains_nocase(orbital->label,"ryd"))hidden_ryd+=*link.weight;
            else hidden_valence+=*link.weight;
        }
    }
    if(state.basis_kind==NboOrbitalKind::NAO) {
        view.focused_hidden_core_weight=hidden_core;
        view.focused_hidden_rydberg_weight=hidden_ryd;
        view.focused_hidden_valence_weight=hidden_valence;
        view.focused_hidden_group_weight=hidden_group;
        view.focused_hidden_h_weight=hidden_h;
    }
    struct EdgeUnit {
        std::vector<std::size_t> edges;
        double strength=0;
        bool focused=false,selected=false;
    };
    std::map<std::string,EdgeUnit> units;
    for(std::size_t i=0;i<view.edges.size();++i) {
        const auto& edge=view.edges[i];
        if(std::abs(edge.coefficient)<=zero) {
            ++view.hidden_numeric_zero_count;continue;
        }
        const auto& source=view.nodes[edge.source_node];
        const auto mo=*view.nodes[edge.target_node].canonical_index;
        const bool focused=focus_members.contains(mo);
        const bool selected=source.id==state.selected_side_node_id;
        if(!focused && !selected && !state.overview)continue;
        std::string channel=source.id;
        if(use_salc && source.salc_index)
            channel=state.salc_model->orbitals[*source.salc_index].subspace_id;
        const std::string group_key=std::to_string(mo_level[mo])+":"+channel;
        auto& unit=units[group_key];
        unit.edges.push_back(i);
        unit.strength+=edge.weight.value_or(
            edge.projection_strength_nonadditive.value_or(std::abs(edge.coefficient)));
        unit.focused|=focused;
        unit.selected|=selected;
    }
    std::vector<EdgeUnit*> ranked;
    for(auto& [_,unit]:units)ranked.push_back(&unit);
    std::stable_sort(ranked.begin(),ranked.end(),[](const auto* a,const auto* b){
        if(a->selected!=b->selected)return a->selected;
        if(a->focused!=b->focused)return a->focused;
        return a->strength>b->strength;
    });
    const std::size_t budget=state.preset==NboAomoPreset::Full?
        std::numeric_limits<std::size_t>::max():
        std::max<std::size_t>(6,static_cast<std::size_t>(view.canvas_width/
            (state.preset==NboAomoPreset::Research?34.0f:52.0f)));
    std::size_t used=0;
    double shown_weight=0,hidden_weight=0;
    for(const auto* unit:ranked) {
        // Keep complete channel and degenerate-member groups together.
        const bool show=state.preset==NboAomoPreset::Full || unit->selected ||
            used==0 || used+unit->edges.size()<=budget;
        if(show)used+=unit->edges.size();
        for(const auto index:unit->edges) {
            auto& edge=view.edges[index];edge.visible=show;
            const auto mo=view.nodes[edge.target_node].canonical_index;
            if(mo && *mo==view.focused_canonical_index && edge.weight) {
                if(show)shown_weight+=*edge.weight;
                else hidden_weight+=*edge.weight;
            }
            if(!show)++view.hidden_readability_count;
        }
    }
    if(state.basis_kind==NboOrbitalKind::NAO) {
        view.focused_display_weight=shown_weight;
        view.focused_hidden_weight=hidden_weight;
    }
    // Frozen display geometry: mean-energy MO rows and compact atomic shells.
    // Neither operation modifies raw energies, coefficients or field identities.
    const auto text_width=[&](const std::string& text) {
        return ImGui::GetFont()->CalcTextSizeA(view.label_font_size,
            10000.0f,0.0f,text.c_str()).x;
    };
    std::map<std::string,std::vector<std::size_t>> shells;
    for(std::size_t i=0;i<view.nodes.size();++i) {
        auto& node=view.nodes[i];
        if(node.individual_label.empty())node.individual_label=node.label;
        if(!node.display_energy_hartree && node.quantitative_energy)
            node.display_energy_hartree=node.energy_hartree;
        if(!node.shell_label.empty() && !node.group_header)
            shells[node.display_group_id].push_back(i);
    }
    const float shell_pitch=view.label_font_size*0.85f;
    if(diagram.options.energy_axis_mode==EnergyAxisMode::NonlinearFocus &&
       !view.illustrative_side_layout && view.energy_transform.knots.size()>1) {
        // The first pass reserves canonical rows. Once actual visible side
        // shells exist, their typography can reserve space on that SAME
        // reversible numeric axis. No raw expectation or eigenvalue moves.
        struct AxisBlock {double energy;float top,bottom;};
        std::array<std::vector<AxisBlock>,3> side_blocks;
        std::set<std::string> grouped;
        for(const auto& node:view.nodes) {
            const auto lane=static_cast<std::size_t>(node.lane);
            if(lane==1 || node.group_header || !node.quantitative_energy)continue;
            if(!node.shell_label.empty()) {
                if(!grouped.insert(node.display_group_id).second)continue;
                const auto& members=shells.at(node.display_group_id);
                if(!std::all_of(members.begin(),members.end(),[&](auto index) {
                    const auto& member=view.nodes[index];
                    return member.quantitative_energy && member.energy_hartree &&
                        std::isfinite(*member.energy_hartree);
                }))continue;
                double energy=0;
                for(const auto index:members)
                    energy+=*view.nodes[index].energy_hartree;
                energy/=static_cast<double>(members.size());
                const float half_span=static_cast<float>(members.size()-1)*shell_pitch*0.5f;
                side_blocks[lane].push_back({energy,-half_span,half_span+22.0f});
            } else if(node.display_energy_hartree || node.energy_hartree) {
                side_blocks[lane].push_back({node.display_energy_hartree?
                    *node.display_energy_hartree:*node.energy_hartree,0.0f,22.0f});
            }
        }
        struct AxisConstraint {std::size_t low,high;double gap;};
        std::vector<AxisConstraint> side_constraints;
        auto& knots=view.energy_transform.knots;
        const auto anchor=[&](double energy) {
            auto it=std::lower_bound(knots.begin(),knots.end(),energy,
                [](const auto& knot,double value){return knot.energy_hartree<value;});
            if(it!=knots.end() && std::abs(it->energy_hartree-energy)<1e-12)
                return static_cast<std::size_t>(it-knots.begin());
            if(it!=knots.begin() && std::abs((it-1)->energy_hartree-energy)<1e-12)
                return static_cast<std::size_t>((it-1)-knots.begin());
            const double coordinate=energy_display_coordinate(energy,view.energy_transform);
            return static_cast<std::size_t>(knots.insert(it,{energy,coordinate})-knots.begin());
        };
        struct AxisCluster {double energy;float top,bottom;};
        std::array<std::vector<AxisCluster>,3> clusters;
        for(const auto lane:{0,2}) {
            auto& blocks=side_blocks[lane];
            std::stable_sort(blocks.begin(),blocks.end(),[](const auto& a,const auto& b) {
                return a.energy<b.energy;
            });
            for(const auto& block:blocks) {
                auto& rows=clusters[lane];
                if(rows.empty() || block.energy-rows.back().energy>1e-7)
                    rows.push_back({block.energy,block.top,block.bottom});
                else {
                    rows.back().top=std::min(rows.back().top,block.top);
                    rows.back().bottom=std::max(rows.back().bottom,block.bottom);
                }
            }
        }
        for(const auto lane:{0,2})
            for(const auto& row:clusters[lane])anchor(row.energy);
        for(const auto lane:{0,2}) {
            const auto& rows=clusters[lane];
            for(std::size_t i=1;i<rows.size();++i) {
                const auto a=anchor(rows[i-1].energy),b=anchor(rows[i].energy);
                if(b>a)side_constraints.push_back({a,b,
                    static_cast<double>(rows[i].bottom-rows[i-1].top+3.0f)});
            }
        }
        if(!side_constraints.empty()) {
            const auto original=knots;
            std::vector<std::vector<AxisConstraint>> incoming(knots.size());
            for(const auto& constraint:side_constraints)
                incoming[constraint.high].push_back(constraint);
            const auto spaced=[&](double fraction) {
                std::vector<double> positions(knots.size());
                positions[0]=original[0].coordinate;
                for(std::size_t i=1;i<knots.size();++i) {
                    positions[i]=positions[i-1]+std::max(1e-6,
                        original[i].coordinate-original[i-1].coordinate);
                    for(const auto& constraint:incoming[i])
                        positions[i]=std::max(positions[i],
                            positions[constraint.low]+constraint.gap*fraction);
                }
                return positions;
            };
            auto positions=spaced(1.0);
            const double base_span=original.back().coordinate-original.front().coordinate;
            const double limit=std::max(960.0,base_span+480.0);
            if(positions.back()-positions.front()>limit) {
                // Exceptionally dense displays retain a bounded shared axis;
                // residual collisions use the real independent side columns.
                double lo_fraction=0,hi_fraction=1;
                for(int i=0;i<28;++i) {
                    const double mid=(lo_fraction+hi_fraction)*0.5;
                    if(const auto trial=spaced(mid);
                       trial.back()-trial.front()<=limit)lo_fraction=mid;
                    else hi_fraction=mid;
                }
                positions=spaced(lo_fraction);
            }
            for(std::size_t i=0;i<knots.size();++i)knots[i].coordinate=positions[i];
            low=energy_display_coordinate(axis_energies.front(),transform);
            high=energy_display_coordinate(axis_energies.back(),transform);
            numeric_span=std::max(numeric_span,static_cast<float>(high-low));
            view.axis_coordinate_min=low;view.axis_coordinate_max=high;
            view.numeric_span=numeric_span;
            view.qualitative_band_y=numeric_top+numeric_span+66.0f;
            view.energy_ticks.clear();
            const int count=std::max(5,static_cast<int>(numeric_span/110.0f));
            for(int i=0;i<=count;++i) {
                const double coordinate=low+(high-low)*static_cast<double>(i)/count;
                view.energy_ticks.push_back({energy_from_display_coordinate(coordinate,transform),
                    numeric_top+static_cast<float>(1.0-static_cast<double>(i)/count)*numeric_span});
            }
            for(auto& node:view.nodes)
                if(node.quantitative_energy && node.display_energy_hartree)
                    node.y=energy_y(*node.display_energy_hartree);
            view.mo_energy_axis_detail+="; side typography constraints share the reversible axis; "
                "exceptionally dense rows use bounded horizontal fallback";
        }
    }
    for(const auto& [id,members]:shells) {
        double sum=0;bool all_quantitative=true;
        float centre_y=0;
        for(const auto index:members) {
            const auto& node=view.nodes[index];centre_y+=node.y;
            all_quantitative&=node.quantitative_energy && node.energy_hartree.has_value();
            if(node.energy_hartree)sum+=*node.energy_hartree;
        }
        centre_y/=static_cast<float>(members.size());
        const double mean=sum/static_cast<double>(members.size());
        if(all_quantitative)centre_y=energy_y(mean);
        for(std::size_t j=0;j<members.size();++j) {
            auto& node=view.nodes[members[j]];
            node.shell_member_index=j;node.shell_member_count=members.size();
            node.display_offset_y=(static_cast<float>(j)-
                static_cast<float>(members.size()-1)*0.5f)*shell_pitch;
            node.display_energy_hartree=all_quantitative?std::optional<double>(mean):std::nullopt;
            node.display_energy_semantics=all_quantitative?
                "same-atom same-shell same-spin mean with symmetric display stack":
                "non-quantitative atomic shell stack; no verified energy";
            node.y=centre_y+node.display_offset_y;
            node.label=j==0?node.shell_label:"";
            node.detail+="; "+node.display_energy_semantics+
                "; individual member: "+node.individual_label;
        }
    }
    for(const auto& node:view.nodes)if(node.quantitative_energy)
        view.qualitative_band_y=std::max(view.qualitative_band_y,
            node.y+node.height+view.label_font_size+24.0f);
    for(const auto lane:{NboAomoLane::Left,NboAomoLane::Right}) {
        float cursor=view.qualitative_band_y+38.0f;
        std::set<std::string> laid_out;
        for(auto& node:view.nodes) {
            if(node.lane!=lane || node.quantitative_energy)continue;
            if(!node.shell_label.empty()) {
                if(!laid_out.insert(node.display_group_id).second)continue;
                const auto& members=shells.at(node.display_group_id);
                const float span=static_cast<float>(members.size()-1)*shell_pitch;
                for(const auto index:members)view.nodes[index].y=
                    cursor+span*0.5f+view.nodes[index].display_offset_y;
                cursor+=span+view.label_font_size+20.0f;
            } else {node.y=cursor;cursor+=view.label_font_size+36.0f;}
        }
    }
    // Each shell is one indivisible typography block. Coincident unrelated
    // blocks get columns; a shell's members always share exactly one x.
    struct DisplayBlock {std::vector<std::size_t> nodes;float top=0,bottom=0,width=0;};
    std::array<std::vector<DisplayBlock>,3> blocks;
    std::vector<std::size_t> node_slot(view.nodes.size(),0);
    std::array<float,3> slot_width{86,86,86};
    std::array<std::size_t,3> slot_count{1,1,1};
    std::set<std::string> placed_shells;
    for(std::size_t i=0;i<view.nodes.size();++i) {
        auto& node=view.nodes[i];
        node.label_width=text_width(node.label);node.label_height=view.label_font_size;
        const auto lane=static_cast<std::size_t>(node.lane);
        if(node.group_header) {
            slot_width[lane]=std::max(slot_width[lane],node.label_width+20.0f);
            continue;
        }
        const bool integer_side_occupation=node.occupation &&
            (std::abs(*node.occupation-1.0)<1e-6 || std::abs(*node.occupation-2.0)<1e-6);
        if(node.occupation && *node.occupation>0.01 &&
            (node.canonical_index || integer_side_occupation)) {
            bool beta=node.orbital && node.orbital->spin==NboSpin::Beta;
            if(node.canonical_index)beta=canonical.orbitals[*node.canonical_index].spin==Spin::Beta;
            else if(node.salc_index && view.salc_model)beta=
                view.salc_model->orbitals[*node.salc_index].spin==NboSpin::Beta;
            node.occupation_label=*node.occupation>1.5?"↑↓":beta?"↓":"↑";
            node.occupation_on_bar=true;node.occupation_width=node.occupation_label=="↑↓"?16.0f:6.0f;
            node.occupation_height=node.shell_label.empty()?22.0f:
                std::min(20.0f,view.label_font_size*0.75f);
        }
    }
    for(std::size_t i=0;i<view.nodes.size();++i) {
        const auto& node=view.nodes[i];
        if(node.group_header)continue;
        const auto lane=static_cast<std::size_t>(node.lane);
        DisplayBlock block;
        if(!node.shell_label.empty()) {
            if(!placed_shells.insert(node.display_group_id).second)continue;
            block.nodes=shells.at(node.display_group_id);
        } else block.nodes={i};
        block.top=std::numeric_limits<float>::infinity();block.bottom=-block.top;
        for(const auto index:block.nodes) {
            const auto& member=view.nodes[index];
            const bool side_label=lane!=1;
            // Side labels and occupation marks lie within the physical bar
            // envelope. Padding here would merge distinct, visible rows into
            // one overlap component and needlessly stagger their columns.
            block.top=std::min(block.top,member.y-(side_label?0.0f:3.0f));
            block.bottom=std::max(block.bottom,member.y+22.0f+
                (side_label?0.0f:view.label_font_size+8.0f));
            block.width=std::max(block.width,side_label?
                std::max(member.label_width,text_width(member.shell_label))+
                    view.orbital_bar_width+28.0f:
                std::max(view.orbital_bar_width,member.label_width)+16.0f);
        }
        slot_width[lane]=std::max(slot_width[lane],block.width);
        blocks[lane].push_back(std::move(block));
    }
    for(std::size_t lane=0;lane<3;++lane) {
        auto& order=blocks[lane];
        std::stable_sort(order.begin(),order.end(),[](const auto& a,const auto& b){return a.top<b.top;});
        std::vector<float> slot_end;
        for(const auto& block:order) {
            std::size_t slot=0;
            while(slot<slot_end.size() && slot_end[slot]>block.top)++slot;
            if(slot==slot_end.size())slot_end.push_back(0);
            slot_end[slot]=block.bottom;
            for(const auto index:block.nodes)node_slot[index]=slot;
        }
        slot_count[lane]=std::max<std::size_t>(1,slot_end.size());
    }
    std::vector<std::size_t> band_order;
    for(std::size_t i=0;i<view.nodes.size();++i)
        if(view.nodes[i].lane==NboAomoLane::Centre && !view.nodes[i].group_header)
            band_order.push_back(i);
    std::stable_sort(band_order.begin(),band_order.end(),[&](auto a,auto b) {
        return view.nodes[a].y<view.nodes[b].y;
    });
    float max_band_width=0,last_band_y=-std::numeric_limits<float>::infinity();
    std::map<std::string,float> band_units;
    const auto finish_band=[&]() {
        float total=0;
        for(const auto& [_,unit_width]:band_units)total+=unit_width;
        if(band_units.size()>1)total+=16.0f*static_cast<float>(band_units.size()-1);
        max_band_width=std::max(max_band_width,total);
        band_units.clear();
    };
    for(const auto index:band_order) {
        const auto& node=view.nodes[index];
        if(node.y-last_band_y>view.label_font_size+24.0f)finish_band();
        const auto id=node.spatial_pair_id.empty()?node.id:node.spatial_pair_id;
        band_units[id]=std::max(band_units[id],
            std::max(view.orbital_bar_width,node.label_width));
        last_band_y=node.y;
    }
    finish_band();
    // Equal outer-label reserve keeps the actual side bar envelopes mirrored
    // even when left and right text lengths differ.
    slot_width[0]=slot_width[2]=std::max(slot_width[0],slot_width[2]);
    for(std::size_t lane=0;lane<3;++lane)
        view.lane_width[lane]=slot_width[lane]*static_cast<float>(slot_count[lane]);
    view.lane_width[1]=std::max(view.lane_width[1],max_band_width);
    const float side_budget=std::max(view.lane_width[0],view.lane_width[2]);
    // Both sides reserve the same physical envelope around the MO manifold.
    view.lane_width[0]=view.lane_width[2]=side_budget;
    std::vector<float> side_slot_offset(view.nodes.size(),0);
    for(const auto lane:{0,2}) {
        // Overlap is transitive: A can touch B and B touch C even if A and C
        // do not meet. Keep that whole component on one centred column grid.
        const auto& order=blocks[lane];
        for(std::size_t begin=0;begin<order.size();) {
            std::size_t end=begin+1;
            float bottom=order[begin].bottom;
            auto first=node_slot[order[begin].nodes.front()];
            auto last=first;
            while(end<order.size() && order[end].top<bottom) {
                bottom=std::max(bottom,order[end].bottom);
                first=std::min(first,node_slot[order[end].nodes.front()]);
                last=std::max(last,node_slot[order[end].nodes.front()]);
                ++end;
            }
            const float occupied=static_cast<float>(last-first+1)*slot_width[lane];
            const float offset=(side_budget-occupied)*0.5f-
                static_cast<float>(first)*slot_width[lane];
            for(std::size_t i=begin;i<end;++i)
                for(const auto index:order[i].nodes)side_slot_offset[index]=offset;
            begin=end;
        }
    }
    const float gutter=std::clamp((std::max(720.0f,available_width)-
        (2*side_budget+view.lane_width[1]+180.0f))*0.25f,46.0f,110.0f);
    float axis_gutter=90.0f;
    for(const auto& tick:view.energy_ticks)
        axis_gutter=std::max(axis_gutter,8.0f+
            text_width(format_energy(tick.energy_hartree,view.energy_unit,3))+16.0f);
    view.lane_x[0]=axis_gutter;
    view.lane_x[1]=view.lane_x[0]+side_budget+gutter;
    view.lane_x[2]=view.lane_x[1]+view.lane_width[1]+gutter;
    view.canvas_width=view.lane_x[2]+side_budget+axis_gutter;
    for(std::size_t i=0;i<view.nodes.size();++i) {
        auto& node=view.nodes[i];const auto lane=static_cast<std::size_t>(node.lane);
        const float column=view.lane_x[lane]+side_slot_offset[i]+
            static_cast<float>(node_slot[i])*slot_width[lane];
        node.x=column;
        node.width=node.group_header?std::max(94.0f,node.label_width+12.0f):view.orbital_bar_width;
        node.height=node.group_header?view.label_font_size+8.0f:22.0f;
        if(lane!=1 && node.group_header)
            node.x=view.lane_x[lane]+(side_budget-node.width)*0.5f;
        node.label_x=node.x+(node.group_header?5.0f:0.0f);
        node.label_y=node.y+(node.group_header?4.0f:25.0f);
        if(lane==2 && !node.group_header) {
            node.label_x=node.x+node.width+12.0f;
            node.label_y=node.y-node.display_offset_y+11.0f-view.label_font_size*0.5f;
        } else if(lane==0 && !node.group_header) {
            node.x=column+slot_width[lane]-node.width;
            node.label_x=node.x-node.label_width-12.0f;
            node.label_y=node.y-node.display_offset_y+11.0f-view.label_font_size*0.5f;
        }
        node.occupation_x=node.x+(node.width-node.occupation_width)*0.5f;
        node.occupation_y=node.y+11.0f-node.occupation_height*0.5f;
    }
    // Close energy rows share a horizontal typography band. This changes only
    // x spacing; every bar stays at its own energy-axis y coordinate.
    std::vector<std::size_t> central_order;
    for(std::size_t i=0;i<view.nodes.size();++i)if(
        view.nodes[i].lane==NboAomoLane::Centre && !view.nodes[i].group_header)
        central_order.push_back(i);
    std::stable_sort(central_order.begin(),central_order.end(),[&](auto a,auto b) {
        return view.nodes[a].y<view.nodes[b].y;
    });
    struct CentralBand {float last_y=0;std::vector<std::size_t> members;};
    std::vector<CentralBand> central_bands;
    for(const auto index:central_order) {
        const float y=view.nodes[index].y;
        if(central_bands.empty() ||
            y-central_bands.back().last_y>view.label_font_size+24.0f)
            central_bands.push_back({y,{}});
        central_bands.back().last_y=y;
        central_bands.back().members.push_back(index);
    }
    for(const auto& band:central_bands) {
        const auto& members=band.members;
        std::vector<std::size_t> units;
        std::map<std::string,std::size_t> pair_unit;
        for(const auto index:members) {
            const auto& node=view.nodes[index];
            if(node.spatial_pair_id.empty()) {units.push_back(index);continue;}
            if(const auto it=pair_unit.find(node.spatial_pair_id);it!=pair_unit.end()) {
                if(!node.label.empty())units[it->second]=index;
            } else {
                pair_unit[node.spatial_pair_id]=units.size();
                units.push_back(index);
            }
        }
        float row_width=units.size()>1?16.0f*static_cast<float>(units.size()-1):0.0f;
        for(const auto index:units)
            row_width+=std::max(view.orbital_bar_width,view.nodes[index].label_width);
        float x=view.canvas_width*0.5f-row_width*0.5f;
        for(const auto index:units) {
            auto& node=view.nodes[index];
            const float unit_width=std::max(view.orbital_bar_width,node.label_width);
            node.x=x+(unit_width-node.width)*0.5f;
            node.label_x=x+(unit_width-node.label_width)*0.5f;
            node.occupation_x=node.x+(node.width-node.occupation_width)*0.5f;
            x+=unit_width+16.0f;
        }
    }
    std::map<std::string,std::vector<std::size_t>> paired_nodes;
    for(std::size_t i=0;i<view.nodes.size();++i)
        if(!view.nodes[i].spatial_pair_id.empty())
            paired_nodes[view.nodes[i].spatial_pair_id].push_back(i);
    for(const auto& [_,members]:paired_nodes)if(members.size()==2) {
        auto& first=view.nodes[members[0]];
        auto& second=view.nodes[members[1]];
        auto* labelled=first.label.empty()?&second:&first;
        const float x=labelled->x;
        first.x=second.x=x;
        // Coincident alpha/beta bars separate horizontally; their exact
        // energy-axis y coordinates never move to make room for a pair.
        if(std::abs(first.y-second.y)<first.height) {
            first.x=x-34.0f;second.x=x+34.0f;
        }
        first.occupation_x=first.x+(first.width-first.occupation_width)*0.5f;
        second.occupation_x=second.x+(second.width-second.occupation_width)*0.5f;
        first.occupation_y=first.y+11.0f-first.occupation_height*0.5f;
        second.occupation_y=second.y+11.0f-second.occupation_height*0.5f;
        labelled->label_y=std::min(first.y,second.y)-view.label_font_size-4.0f;
    }
    // A counterpart may be far from its partner on the true energy axis.
    // Pair labels are positioned after horizontal bands, so reserve their
    // actual rectangles here without moving any orbital energy line.
    std::vector<std::size_t> centre_labels;
    for(std::size_t i=0;i<view.nodes.size();++i)
        if(view.nodes[i].lane==NboAomoLane::Centre &&
           !view.nodes[i].group_header && !view.nodes[i].label.empty())
            centre_labels.push_back(i);
    std::stable_sort(centre_labels.begin(),centre_labels.end(),[&](auto a,auto b) {
        return view.nodes[a].label_y>view.nodes[b].label_y;
    });
    std::vector<std::size_t> placed_labels;
    for(const auto index:centre_labels) {
        auto& node=view.nodes[index];
        bool moved;
        do {
            moved=false;
            for(const auto earlier:placed_labels) {
                const auto& other=view.nodes[earlier];
                const bool x_overlap=node.label_x<other.label_x+other.label_width &&
                    other.label_x<node.label_x+node.label_width;
                const bool y_overlap=node.label_y<other.label_y+other.label_height &&
                    other.label_y<node.label_y+node.label_height;
                if(x_overlap && y_overlap) {
                    node.label_y=other.label_y-node.label_height-4.0f;
                    moved=true;
                }
            }
        } while(moved);
        placed_labels.push_back(index);
    }
    float top_padding=0;
    for(const auto& node:view.nodes)if(!node.group_header)
        top_padding=std::max(top_padding,78.0f-std::min(node.y,node.label_y));
    if(top_padding>0) {
        for(auto& node:view.nodes) {
            node.y+=top_padding;node.label_y+=top_padding;node.occupation_y+=top_padding;
        }
        for(auto& tick:view.energy_ticks)tick.y+=top_padding;
        view.numeric_top+=top_padding;view.qualitative_band_y+=top_padding;
    }
    view.canvas_height=view.numeric_top+numeric_span+view.label_font_size+36.0f;
    for(const auto& node:view.nodes)
        view.canvas_height=std::max(view.canvas_height,
            std::max(node.label_y+node.label_height,node.y+node.height)+22.0f);
    view.footer_y=view.canvas_height;
    const auto caption=[&](const std::string& role,const std::string& text,float x,float y) {
        view.captions.push_back({role,text,x,y,text_width(text),view.label_font_size});
    };
    const auto wrapped=[&](const std::string& role,const std::string& text,float x,float& y) {
        std::string line;
        for(std::size_t i=0;i<text.size();) {
            const auto c=static_cast<unsigned char>(text[i]);
            const std::size_t length=c<0x80?1:c<0xe0?2:c<0xf0?3:4;
            const auto next=text.substr(i,length);i+=length;
            if(!line.empty() && text_width(line+next)>view.canvas_width-x-12.0f) {
                const auto space=line.rfind(' ');
                if(space!=std::string::npos && space>line.size()/2) {
                    caption(role,line.substr(0,space),x,y);line=line.substr(space+1);
                } else {caption(role,line,x,y);line.clear();}
                y+=view.label_font_size+5.0f;
            }
            if(!line.empty() || next!=" ")line+=next;
        }
        if(!line.empty()){caption(role,line,x,y);y+=view.label_font_size+5.0f;}
    };
    caption("heading","AO / NAO",view.lane_x[0],7);
    caption("heading","MO",view.canvas_width*0.5f-text_width("MO")*0.5f,7);
    caption("heading","SALC / AO",view.lane_x[2],7);
    float header_y=30;
    wrapped("axis","E ("+view.display_energy_unit+"): "+
        (diagram.options.energy_axis_mode==EnergyAxisMode::NonlinearFocus?
            aomo_text(state.language,"Nonlinear energy axis"):
            aomo_text(state.language,"Energy axis")),12,header_y);
    wrapped("semantics",aomo_text(state.language,"Grouped levels use mean energy."),12,header_y);
    // Reserve the caption area even when a narrow chart wraps its header.
    const float caption_padding=std::max(0.0f,header_y+8.0f-view.numeric_top);
    if(caption_padding>0) {
        for(auto& node:view.nodes){node.y+=caption_padding;node.label_y+=caption_padding;node.occupation_y+=caption_padding;}
        for(auto& tick:view.energy_ticks)tick.y+=caption_padding;
        view.numeric_top+=caption_padding;view.qualitative_band_y+=caption_padding;
        view.footer_y+=caption_padding;
    }
    if(std::any_of(view.nodes.begin(),view.nodes.end(),[](const auto& node){return node.lane!=NboAomoLane::Centre && !node.quantitative_energy;}))
        for(const auto lane:{0,2})caption("nonquantitative",aomo_text(state.language,"Non-quantitative"),view.lane_x[lane],view.qualitative_band_y+5);
    float footer_y=view.footer_y;
    for(const auto& line:coverage_lines(view))if(!line.empty())wrapped("coverage",line,12,footer_y);
    view.canvas_height=footer_y+12;
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
            if(!b.canonical_index)continue;
            const double magnitude=edge.projection_strength_nonadditive.value_or(
                std::abs(edge.coefficient));
            if(*b.canonical_index==view.focused_canonical_index) {
                if(a.orbital)basis_to_focus[key(*a.orbital)]=
                    std::max(basis_to_focus[key(*a.orbital)],magnitude);
                max_focus=std::max(max_focus,magnitude);
            }
            if((a.orbital && selected_basis.contains(key(*a.orbital))) ||
                a.id==view.selected_side_node_id) {
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
        if(node.id==view.selected_side_node_id)return 1.0f;
        if(node.salc_index) {
            for(const auto& edge:view.edges)if(edge.visible && edge.source_node<view.nodes.size() &&
                &view.nodes[edge.source_node]==&node &&
                view.nodes[edge.target_node].canonical_index==view.focused_canonical_index)
                return relative(std::abs(edge.coefficient),max_focus);
        }
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
            (a.orbital && selected_basis.contains(key(*a.orbital))) ||
            a.id==view.selected_side_node_id;
        double maximum=0;
        if(b.canonical_index && *b.canonical_index==view.focused_canonical_index)maximum=max_focus;
        else if((a.orbital && selected_basis.contains(key(*a.orbital))) ||
            a.id==view.selected_side_node_id)maximum=max_selected;
        const double strength=link.projection_strength_nonadditive.value_or(
            std::abs(link.coefficient));
        float intensity=relative(strength,maximum);
        if(a.fragment_group_id)intensity=a.metric_norm2 && max_focus>0?
            static_cast<float>(std::min(1.0,std::sqrt(std::max(0.0,*a.metric_norm2))/max_focus)):0.0f;
        EdgeAppearance appearance;
        if(link.coefficient<0){appearance.red=242;appearance.green=145;appearance.blue=143;}
        appearance.alpha=static_cast<unsigned char>(active?85+65*intensity:35+25*intensity);
        appearance.width=0.7f+0.45f*intensity;
        return appearance;
    }
};
} // namespace

std::string nbo_aomo_glyph_seed(Language language) {
    return nbo_aomo_text_glyph_seed(language)+nbo_aomo_hover_glyph_seed(language);
}

bool prepare_nbo_aomo_state(NboAomoUIState& state,const NboIntegration& data,
                           const Wavefunction& canonical) {
    const auto* capability=nbo_capability(data,"aomo");
    if(!capability || !capability->available()) {
        state.status=capability?capability->detail:"AO–MO decomposition data missing";
        state.drawn_snapshot.reset();return false;
    }
    const auto attachment_id=data.id+":"+data.canonical_fingerprint;
    if(state.source_id!=attachment_id){state=NboAomoUIState{};state.source_id=attachment_id;}
    if(!state.salc_model)
        state.salc_model=std::make_shared<const NboSalcModel>(build_nbo_salc_model(canonical,data));
    if(!state.names || state.names_model!=state.salc_model.get()) {
        state.names=std::make_shared<const NboAomoNames>(
            build_nbo_aomo_names(canonical,data,state.salc_model.get()));
        state.names_model=state.salc_model.get();
    }
    return true;
}

bool draw_nbo_aomo_diagram(NboAomoUIState& state,const NboIntegration& data,
                           const Wavefunction& canonical,const MODiagramViewSnapshot& diagram,
                           Language language,float scale) {
    if(!prepare_nbo_aomo_state(state,data,canonical))return false;
    state.language=language;
    const auto indices=central_indices(diagram,canonical.orbitals.size());
    if(indices.empty()){state.status="The central MO diagram has no real members";state.drawn_snapshot.reset();return false;}
    const auto inspected=diagram.data.view?diagram.data.view->inspected_orbital_index:std::nullopt;
    if(inspected!=state.last_inspected){
        if(inspected && in(indices,*inspected)) {
            state.focused_canonical_index=*inspected;
            state.selected_side_node_id.clear();
        }
        state.last_inspected=inspected;++state.revision;
    }
    if(!state.focused_canonical_index || !in(indices,*state.focused_canonical_index))
        state.focused_canonical_index=indices.front();
    ImGui::SeparatorText(aomo_text(language,"Orbital interaction diagram"));
    // Capture the complete content row before SameLine controls change the
    // cursor; fitting must use the graph width, never a toolbar remainder.
    const float graph_available_width=std::max(1.0f,ImGui::GetContentRegionAvail().x-28*scale);
    const char* preset_names[]={aomo_text(language,"Overview"),
        aomo_text(language,"Research analysis"),aomo_text(language,"Full basis")};
    constexpr const char* preset_items[]={"aomo.preset.teaching","aomo.preset.research","aomo.preset.full"};
    const char* preset_caption=aomo_text(language,"View");
    const float preset_label_width=ImGui::CalcTextSize(preset_caption).x+ImGui::GetStyle().ItemInnerSpacing.x;
    ImGui::SetNextItemWidth(std::max(80*scale,std::min(260*scale,ImGui::GetContentRegionAvail().x-preset_label_width)));
    const bool preset_open=ImGui::BeginCombo((std::string(preset_caption)+"###aomo.preset").c_str(),
        preset_names[static_cast<int>(state.preset)]);
    validation::item("aomo.preset");
    if(preset_open) {
        for(int preset=0;preset<3;++preset) {
            const bool selected=static_cast<int>(state.preset)==preset;
            if(ImGui::Selectable(preset_names[preset],selected) && !selected) {
                state.preset=static_cast<NboAomoPreset>(preset);++state.revision;
                if(state.preset==NboAomoPreset::Full) {
                    state.show_core=true;state.show_rydberg=true;
                } else if(state.preset==NboAomoPreset::Teaching) {
                    state.show_core=false;state.show_rydberg=false;
                    state.overview=true;
                }
            }
            validation::item(preset_items[preset]);
            if(selected)ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    same_line_for(aomo_text(language,"Other MO links"));
    if(ImGui::Checkbox((std::string(aomo_text(language,"Other MO links"))+"###aomo.overview").c_str(),&state.overview))++state.revision;
    validation::item("aomo.overview");
    if(state.preset==NboAomoPreset::Teaching) {
        if(ImGui::Checkbox(aomo_text(language,"Illustrative layout (not energy)"),
            &state.illustrative_side_layout))++state.revision;
        validation::item("aomo.illustrative_side_layout");
    }
    const bool nao=state.basis_kind==NboOrbitalKind::NAO;
    if(ImGui::RadioButton("NAO##aomo.basis",nao)){state.basis_kind=NboOrbitalKind::NAO;++state.revision;}
    validation::item("aomo.basis.nao");same_line_for("Gaussian AO");
    if(ImGui::RadioButton("Gaussian AO##aomo.basis",!nao)){state.basis_kind=NboOrbitalKind::GaussianAO;++state.revision;}
    validation::item("aomo.basis.ao");
    same_line_for("−");
    if(ImGui::Button("−##aomo.zoom")){state.zoom=std::max(0.35f,state.zoom/1.2f);++state.revision;}
    validation::item("aomo.zoom.out");same_line_for("+");
    if(ImGui::Button("+##aomo.zoom")){state.zoom=std::min(3.0f,state.zoom*1.2f);++state.revision;}
    validation::item("aomo.zoom.in");same_line_for(aomo_text(language,"Reset view"));
    if(ImGui::Button(aomo_text(language,"Reset view"))){
        state.zoom=1;state.pan_x=0;state.pan_y=0;++state.revision;}
    validation::item("aomo.view.reset");
    same_line_for(aomo_text(language,"Fit graph"));
    if(ImGui::Button((std::string(aomo_text(language,"Fit graph"))+"###aomo.fit").c_str()) && state.drawn_snapshot) {
        state.zoom=std::clamp(std::min(
            graph_available_width/state.drawn_snapshot->canvas_width,
            620*scale/state.drawn_snapshot->canvas_height),0.35f,1.8f);
        state.pan_x=0;state.pan_y=0;++state.revision;
    }
    validation::item("aomo.view.fit");
    const bool scroll_to_focus=ImGui::Button((std::string(aomo_text(language,"Find selected MO"))+"###aomo.find").c_str());
    validation::item("aomo.focus.find");
    auto draw_graph_options=[&]() {
        if(state.basis_kind==NboOrbitalKind::NAO) {
            if(ImGui::Checkbox(aomo_text(language,"Show core NAOs"),&state.show_core))++state.revision;
            validation::item("aomo.show_core");same_line_for(aomo_text(language,"Show Rydberg NAOs"));
            if(ImGui::Checkbox(aomo_text(language,"Show Rydberg NAOs"),&state.show_rydberg))++state.revision;
            validation::item("aomo.show_rydberg");
        }
        if(ImGui::Checkbox(aomo_text(language,"Hide H orbitals"),&state.hide_h_orbitals))
            ++state.revision;
        validation::item("aomo.hide_h_orbitals");
        if(ImGui::CollapsingHeader(aomo_text(language,"User atom sets (advanced)"))) {
            ImGui::TextWrapped("%s",aomo_text(language,"Atom sets organize the view. Their partial sums depend on the selected MO, rather than defining a fixed SALC."));
            std::set<std::size_t> available_atoms;
            for(const auto& orbital:data.orbitals)
                if(orbital.ref.kind==state.basis_kind)
                    available_atoms.insert(orbital.atoms.begin(),orbital.atoms.end());
            std::size_t slot=0;
            for(const auto atom:available_atoms) {
                if(atom>=canonical.atoms.size())continue;
                bool checked=state.draft_fragment_atoms.contains(atom);
                const auto label=canonical.atoms[atom].symbol+std::to_string(atom+1)+
                    "##aomo.fragment.atom."+std::to_string(atom);
                if(ImGui::Checkbox(label.c_str(),&checked)) {
                    if(checked)state.draft_fragment_atoms.insert(atom);
                    else state.draft_fragment_atoms.erase(atom);
                    ++state.revision;
                }
                validation::item("aomo.fragment.atom."+std::to_string(atom));
                if(++slot%5)ImGui::SameLine();
            }
            if(slot)ImGui::NewLine();
            ImGui::BeginDisabled(state.draft_fragment_atoms.empty());
            if(ImGui::Button(aomo_text(language,"Save atom set"))) {
                if(state.editing_fragment_id) {
                    const auto it=std::find_if(state.fragment_groups.begin(),
                        state.fragment_groups.end(),[&](const auto& group){
                            return group.id==*state.editing_fragment_id;});
                    if(it!=state.fragment_groups.end()) {
                        it->atoms=state.draft_fragment_atoms;it->suggested=false;
                    }
                } else state.fragment_groups.push_back({state.next_fragment_id++,
                    state.draft_fragment_atoms,false});
                state.draft_fragment_atoms.clear();state.editing_fragment_id.reset();
                ++state.revision;
            }
            ImGui::EndDisabled();validation::item("aomo.fragment.add");ImGui::SameLine();
            if(ImGui::Button(aomo_text(language,"Clear draft"))) {
                state.draft_fragment_atoms.clear();state.editing_fragment_id.reset();
                ++state.revision;
            }
            validation::item("aomo.fragment.clear_draft");
            for(std::size_t i=0;i<state.fragment_groups.size();) {
                const auto& group=state.fragment_groups[i];
                std::string label="F"+std::to_string(group.id)+" (";
                for(auto atom:group.atoms)if(atom<canonical.atoms.size())
                    label+=canonical.atoms[atom].symbol+std::to_string(atom+1)+" ";
                label+=")";
                ImGui::TextUnformatted(label.c_str());ImGui::SameLine();
                if(ImGui::SmallButton((std::string(aomo_text(language,"Edit"))+"###aomo.fragment.edit."+
                    std::to_string(group.id)).c_str())) {
                    state.draft_fragment_atoms=group.atoms;
                    state.editing_fragment_id=group.id;++state.revision;
                }
                validation::item("aomo.fragment.edit."+std::to_string(group.id));
                ImGui::SameLine();
                if(ImGui::SmallButton((std::string(aomo_text(language,"Delete"))+"###aomo.fragment.delete."+
                    std::to_string(group.id)).c_str())) {
                    state.expanded_user_fragments.erase(group.id);
                    state.fragment_groups.erase(state.fragment_groups.begin()+
                        static_cast<std::ptrdiff_t>(i));
                    ++state.revision;continue;
                }
                validation::item("aomo.fragment.delete."+std::to_string(group.id));
                ImGui::SameLine();
                const std::vector<std::size_t> atoms(group.atoms.begin(),group.atoms.end());
                const auto partial=nbo_fragment_selection(data,
                    state.focused_canonical_index.value_or(0),atoms,state.basis_kind,false);
                ImGui::BeginDisabled(partial.terms.empty());
                if(ImGui::SmallButton((std::string(aomo_text(language,"MO partial sum in 3D"))+"###aomo.fragment.sum."+
                    std::to_string(group.id)).c_str()))state.pending_selection=partial;
                ImGui::EndDisabled();
                validation::item("aomo.fragment.sum."+std::to_string(group.id));
                ++i;
            }
        }
    };
    const float canvas_width=std::max(540.0f*scale,
        ImGui::GetContentRegionAvail().x-28.0f*scale);
    auto initial_snapshot=make_unified_snapshot(state,data,canonical,diagram,canvas_width);
    if(state.initial_fit_pending) {
        // The details panel changes width after a calculation is attached.
        // Fit until two consecutive frames agree on its available width;
        // otherwise the first, narrower layout leaves compact diagrams tiny.
        // Larger diagrams keep readable native size and remain scrollable.
        state.zoom=std::clamp(std::min(canvas_width/initial_snapshot.canvas_width,
            620*scale/initial_snapshot.canvas_height),1.0f,1.8f);
        initial_snapshot.zoom=state.zoom;
        state.initial_fit_stable_frames=
            std::abs(state.initial_fit_viewport_width-canvas_width)<0.5f
            ?state.initial_fit_stable_frames+1:0;
        state.initial_fit_viewport_width=canvas_width;
        state.initial_fit_pending=state.initial_fit_stable_frames<2;
    }
    const auto snapshot=std::make_shared<const NboAomoViewSnapshot>(std::move(initial_snapshot));
    state.drawn_snapshot=snapshot;
    state.status=snapshot->capability_status+": "+snapshot->capability_detail;
    validation::field("aomo.snapshot",snapshot->id);
    validation::anchor("aomo.graph");
    const float canvas_height=snapshot->canvas_height*state.zoom;
    ImGui::BeginChild("##aomo.graph.scroll",ImVec2(0,std::min(660.0f*scale,canvas_height+28.0f*scale)),
        ImGuiChildFlags_Border,ImGuiWindowFlags_HorizontalScrollbar);
    if(scroll_to_focus) {
        const auto focused=std::find_if(snapshot->nodes.begin(),snapshot->nodes.end(),
            [&](const auto& node){return node.canonical_index==snapshot->focused_canonical_index;});
        if(focused!=snapshot->nodes.end()) {
            ImGui::SetScrollY(std::max(0.0f,focused->y*state.zoom-180.0f*scale));
            ImGui::SetScrollX(std::max(0.0f,focused->x*state.zoom-
                ImGui::GetContentRegionAvail().x*0.5f));
        }
    }
    ImGui::InvisibleButton("##aomo.canvas",ImVec2(snapshot->canvas_width*state.zoom,canvas_height));
    validation::item("aomo.graph.canvas");
    const auto origin=ImGui::GetItemRectMin();
    const auto canvas_max=ImGui::GetItemRectMax();
    const auto transform=[&](const NboAomoNode& node){return ImVec2(origin.x+state.pan_x+node.x*state.zoom,
                                                              origin.y+state.pan_y+node.y*state.zoom);};
    const auto point=[&](float x,float y){return ImVec2(origin.x+state.pan_x+x*state.zoom,
        origin.y+state.pan_y+y*state.zoom);};
    const bool hover=ImGui::IsItemHovered();
    auto* draw=ImGui::GetWindowDrawList();
    const auto text_at=[&](float x,float y,ImU32 colour,const std::string& text,float size=14.0f){
        draw->AddText(ImGui::GetFont(),size*state.zoom,point(x,y),colour,text.c_str());
    };
    draw->PushClipRect(origin,ImGui::GetItemRectMax(),true);
    draw->AddRectFilled(origin,ImGui::GetItemRectMax(),IM_COL32(24,31,44,255),6*scale);
    const bool has_nonquant=std::any_of(snapshot->nodes.begin(),snapshot->nodes.end(),
        [](const auto& node){return node.lane!=NboAomoLane::Centre &&
            !node.quantitative_energy;});
    if(has_nonquant) {
        for(const auto lane:{0,2}) {
            draw->AddRectFilled(point(snapshot->lane_x[lane]-8,snapshot->qualitative_band_y),
                point(snapshot->lane_x[lane]+snapshot->lane_width[lane],snapshot->footer_y-8),
                IM_COL32(35,41,52,255));
        }
    }
    for(const auto& caption:snapshot->captions) {
        const auto colour=caption.role=="axis" &&
            snapshot->energy_transform.mode==EnergyAxisMode::NonlinearFocus?
            IM_COL32(snapshot->energy_tick_screen_rgb[0],snapshot->energy_tick_screen_rgb[1],
                snapshot->energy_tick_screen_rgb[2],255):IM_COL32(160,180,202,255);
        text_at(caption.x,caption.y,colour,caption.text,snapshot->label_font_size);
    }
    for(const auto& tick:snapshot->energy_ticks) {
        const float y=origin.y+state.pan_y+(tick.y+11.0f)*state.zoom;
        draw->AddLine(ImVec2(origin.x+20,y),ImVec2(canvas_max.x-20,y),
            IM_COL32(72,87,108,32),1.0f*scale);
        const auto label=format_energy(tick.energy_hartree,diagram.options.energy_unit,3);
        text_at(8,tick.y+11.0f-snapshot->label_font_size,
            IM_COL32(snapshot->energy_tick_screen_rgb[0],snapshot->energy_tick_screen_rgb[1],
                snapshot->energy_tick_screen_rgb[2],255),label,snapshot->label_font_size);
    }
    const GraphAppearance appearance(*snapshot);
    for(const auto& edge:snapshot->edges){
        if(!edge.visible)continue;
        const auto& a=snapshot->nodes[edge.source_node];const auto& b=snapshot->nodes[edge.target_node];
        const auto style=appearance.edge(edge);
        const auto pa=transform(a),pb=transform(b);
        const float ax=a.lane==NboAomoLane::Right?pa.x:pa.x+node_width(a)*state.zoom;
        const float bx=a.lane==NboAomoLane::Right?pb.x+node_width(b)*state.zoom:pb.x;
        dashed_segments(ax/state.zoom,(pa.y+a.height*0.5f*state.zoom)/state.zoom,
            bx/state.zoom,(pb.y+b.height*0.5f*state.zoom)/state.zoom,
            [&](float x1,float y1,float x2,float y2) {
                draw->AddLine(ImVec2(x1*state.zoom,y1*state.zoom),
                    ImVec2(x2*state.zoom,y2*state.zoom),
                    IM_COL32(style.red,style.green,style.blue,style.alpha),style.width*scale);
            });
    }
    // Paint every text mask after edges, before any orbital bars or labels.
    // This removes edge/glyph crossings without shifting a physical level or
    // allowing a later node's mask to erase an earlier node's orbital bar.
    for_each_node_text_background(*snapshot,[&](const std::string&,const char*,
        float x,float y,float width,float height,bool qualitative) {
        draw->AddRectFilled(point(x,y),point(x+width,y+height),qualitative?
            IM_COL32(35,41,52,255):IM_COL32(24,31,44,255));
    });
    std::optional<std::size_t> clicked;
    float best=1e9f;
    for(std::size_t i=0;i<snapshot->nodes.size();++i){
        const auto& node=snapshot->nodes[i];const auto p=transform(node);
        const float width=node_width(node)*state.zoom;
        const float height=node.height*state.zoom;
        const ImVec2 q(p.x+width,p.y+height);
        const float line_y=p.y+height*0.5f;
        const float hit_half=node.shell_member_count>1?
            snapshot->label_font_size*0.85f*0.4f*state.zoom:height*0.5f;
        const ImVec2 hit_min(std::max(origin.x,p.x),std::max(origin.y,line_y-hit_half));
        const ImVec2 hit_max(std::min(canvas_max.x,q.x),std::min(canvas_max.y,line_y+hit_half));
        if(hit_min.x<hit_max.x && hit_min.y<hit_max.y)
            validation::hit("aomo.node."+node.id,hit_min,hit_max);
        const float strength=appearance.node_strength(node);
        const bool selected_node=strength>=0.99f;
        const auto colour=selected_node?IM_COL32(122,223,255,255):
            node.group_header?IM_COL32(143,161,187,255):
            node.quantitative_energy?IM_COL32(228,239,249,255):IM_COL32(174,189,207,255);
        const bool folded_group=node.group_header;
        if(folded_group) {
            draw->AddRectFilled(p,q,IM_COL32(57,70,91,255),4*scale);
            draw->AddRect(p,q,IM_COL32(110,137,166,255),4*scale);
        } else draw->AddLine(ImVec2(p.x,line_y),ImVec2(q.x,line_y),colour,
            snapshot->orbital_bar_stroke_width*state.zoom);
        text_at(node.label_x,node.label_y,colour,node.label,snapshot->label_font_size);
        if(node.occupation_on_bar)
            electron_strokes(node,[&](float x1,float y1,float x2,float y2) {
                draw->AddLine(point(x1,y1),point(x2,y2),colour,1.5f*state.zoom);
            });
        else if(!node.occupation_label.empty())
            text_at(node.occupation_x,node.occupation_y,IM_COL32(194,209,225,255),
                node.occupation_label,snapshot->label_font_size);
        const bool bar_hover=hit_min.x<hit_max.x && hit_min.y<hit_max.y &&
            ImGui::IsMouseHoveringRect(hit_min,hit_max);
        const bool label_hover=!node.label.empty() && ImGui::IsMouseHoveringRect(
            point(node.label_x,node.label_y),point(node.label_x+node.label_width,node.label_y+node.label_height));
        const float distance=bar_hover?std::abs(ImGui::GetIO().MousePos.y-line_y):1000.0f;
        if(hover && (bar_hover || label_hover) && distance<best){clicked=i;best=distance;
            ImGui::BeginTooltip();
            ImGui::PushTextWrapPos(ImGui::GetFontSize()*38.0f);
            const auto lines=nbo_aomo_hover_lines(node,canonical,data,state.salc_model.get(),language);
            std::string hover_text;
            for(const auto& line:lines) {
                ImGui::TextUnformatted(line.c_str());
                if(!hover_text.empty())hover_text+='\n';hover_text+=line;
            }
            validation::field("aomo.hover",hover_text);
            ImGui::PopTextWrapPos();ImGui::EndTooltip();}
    }
    draw->PopClipRect();
    ImGui::EndChild();
    if(snapshot->canvas_width>ImGui::GetContentRegionAvail().x+2.0f)
        ImGui::TextDisabled("%s",aomo_text(language,"Protected orbital partners need more width; scroll the diagram horizontally."));
    ImGui::TextDisabled("%s",aomo_text(language,"Scroll to explore; Ctrl+wheel zooms. Click an orbital to inspect it."));
    if(hover && ImGui::GetIO().KeyCtrl && ImGui::GetIO().MouseWheel!=0){
        state.zoom=std::clamp(state.zoom*(ImGui::GetIO().MouseWheel>0?1.12f:0.89f),0.35f,3.0f);++state.revision;
    }
    if(hover && !clicked && ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
        const auto delta=ImGui::GetIO().MouseDelta;state.pan_x+=delta.x;state.pan_y+=delta.y;++state.revision;
    }
    if(hover && clicked && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
        const auto& node=snapshot->nodes[*clicked];
        if(!node.group_header && node.lane!=NboAomoLane::Centre &&
            !node.subspace_id.empty()) {
            state.expanded_subspaces.erase(node.subspace_id);
            state.collapsed_subspaces.insert(node.subspace_id);
            ++state.revision;
        } else if(!node.group_header && node.orbital &&
            node.lane!=NboAomoLane::Centre && !node.atoms.empty()) {
            state.expanded_atoms.erase(node.atoms.front());
            state.collapsed_atoms.insert(node.atoms.front());
            ++state.revision;
        }
    }
    bool choose_overlapping_mo=false;
    if(hover && clicked && ImGui::IsMouseClicked(ImGuiMouseButton_Left) &&
        snapshot->nodes[*clicked].canonical_index) {
        state.ambiguous_node_ids.clear();
        for(const auto& node:snapshot->nodes)if(node.canonical_index) {
            const auto p=transform(node);
            if(ImGui::IsMouseHoveringRect(p,
                ImVec2(p.x+node_width(node)*state.zoom,
                    p.y+node.height*state.zoom)))
                state.ambiguous_node_ids.push_back(node.id);
        }
        if(state.ambiguous_node_ids.size()>1) {
            choose_overlapping_mo=true;
            ImGui::OpenPopup("##aomo.overlapping.nodes");
        }
    }
    if(hover && clicked && ImGui::IsMouseClicked(ImGuiMouseButton_Left) &&
        !choose_overlapping_mo) {
        const auto& node=snapshot->nodes[*clicked];
        if(node.group_header){
            if(node.id.rfind("salc-group:",0)==0) {
                state.collapsed_subspaces.erase(node.subspace_id);
                state.expanded_subspaces.insert(node.subspace_id);
                ++state.revision;
            } else if(node.id.rfind("atom:",0)==0){const auto atom=static_cast<std::size_t>(std::stoull(node.id.substr(5)));
                state.collapsed_atoms.erase(atom);state.expanded_atoms.insert(atom);++state.revision;}
            else if(node.id.rfind("fragment-group:",0)==0){
                state.expanded_fragments.insert(node.id.substr(std::string("fragment-group:").size()));
                ++state.revision;}
            else if(node.id.rfind("user-fragment:",0)==0){
                state.expanded_user_fragments.insert(static_cast<std::size_t>(
                    std::stoull(node.id.substr(std::string("user-fragment:").size()))));
                ++state.revision;}
        } else if(node.salc_index && state.salc_model) {
            state.pending_selection=nbo_salc_selection(*state.salc_model,*node.salc_index);
            state.selected_side_node_id=node.id;++state.revision;
        } else if(node.fragment_group_id && node.canonical_index) {
            state.pending_selection=nbo_fragment_selection(data,*node.canonical_index,node.atoms,state.basis_kind,false);
        } else if(node.orbital) {
            state.pending_selection=nbo_single_selection(data,*node.orbital);
            if(node.lane!=NboAomoLane::Centre) {
                state.selected_side_node_id=node.id;++state.revision;
            }
            if(node.canonical_index) {
                state.focused_canonical_index=*node.canonical_index;
                state.selected_side_node_id.clear();++state.revision;
            }
        }
    }
    if(ImGui::BeginPopup("##aomo.overlapping.nodes")) {
        ImGui::TextUnformatted(aomo_text(language,"Choose an MO"));
        for(const auto& id:state.ambiguous_node_ids) {
            const auto it=std::find_if(snapshot->nodes.begin(),snapshot->nodes.end(),
                [&](const auto& node){return node.id==id;});
            if(it==snapshot->nodes.end() || !it->canonical_index)continue;
            const auto label=it->label+"  E="+
                (it->energy_hartree?number(*it->energy_hartree):"unknown")+
                " Ha##"+id;
            if(ImGui::Selectable(label.c_str())) {
                if(it->orbital)state.pending_selection=nbo_single_selection(data,*it->orbital);
                state.focused_canonical_index=*it->canonical_index;
                state.selected_side_node_id.clear();
                ++state.revision;
                ImGui::CloseCurrentPopup();
            }
            validation::item("aomo.node.choice."+id);
        }
        ImGui::EndPopup();
    }
    auto activate_edge=[&](const NboAomoEdge& edge) {
        const auto& a=snapshot->nodes[edge.source_node];
        const auto& b=snapshot->nodes[edge.target_node];
        if(!b.canonical_index)return;
        if(edge.salc_link_index && state.salc_model &&
            *edge.salc_link_index<state.salc_model->links.size()) {
            state.pending_selection=nbo_salc_component_selection(*state.salc_model,
                state.salc_model->links[*edge.salc_link_index]);
        } else if(a.fragment_group_id) {
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
        if(!edge.visible)continue;
        const auto& a=snapshot->nodes[edge.source_node];
        const auto& b=snapshot->nodes[edge.target_node];
        const auto p=transform(a),q=transform(b);
        const ImVec2 start(a.lane==NboAomoLane::Right?p.x:p.x+node_width(a)*state.zoom,
            p.y+a.height*0.5f*state.zoom),
            end(a.lane==NboAomoLane::Right?q.x+node_width(b)*state.zoom:q.x,
                q.y+b.height*0.5f*state.zoom);
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
        ImGui::BeginTooltip();ImGui::PushTextWrapPos(ImGui::GetFontSize()*32);
        ImGui::TextUnformatted((source_node.individual_label+" → "+
            snapshot->nodes[edge.target_node].individual_label).c_str());
        ImGui::TextUnformatted(aomo_text(language,"Orbital contribution; not a bond-order label"));
        ImGui::TextUnformatted(aomo_text(language,"Click to show this weighted component"));
        ImGui::PopTextWrapPos();ImGui::EndTooltip();
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
        ImGui::TextUnformatted(aomo_text(language,"Choose a numerical component"));
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
    if(snapshot->focused_projection_weight && snapshot->focused_projection_residual_norm)
        validation::field("aomo.projection",std::to_string(snapshot->focused_canonical_index)+
            ";weight="+number(*snapshot->focused_projection_weight)+
            ";residual_norm="+number(*snapshot->focused_projection_residual_norm));
    if(ImGui::Checkbox(aomo_text(language,"Orbital combinations (advanced)"),&state.show_full_numeric))++state.revision;
    validation::item("aomo.full_numeric");
    if(state.sum_canonical_index!=snapshot->focused_canonical_index ||
       state.sum_basis_kind!=state.basis_kind) {
        state.sum_component_ids.clear();
        state.sum_canonical_index=snapshot->focused_canonical_index;
        state.sum_basis_kind=state.basis_kind;
    }
    if(state.show_full_numeric) {
        validation::anchor("aomo.sum.controls");
        const auto links=nbo_links_for_mo(data,snapshot->focused_canonical_index,state.basis_kind);
        if(ImGui::Button(aomo_text(language,"Add every local term"))) {
            for(const auto& link:links)state.sum_component_ids.insert(row_id(link.orbital));
        }
        validation::item("aomo.sum.all");ImGui::SameLine();
        if(ImGui::Button(aomo_text(language,"Clear terms")))state.sum_component_ids.clear();
        validation::item("aomo.sum.clear");
        std::vector<NboOrbitalTerm> terms;
        for(const auto& link:links)if(state.sum_component_ids.contains(row_id(link.orbital)))
            terms.push_back({link.orbital,link.coefficient});
        ImGui::BeginDisabled(terms.empty());
        if(ImGui::Button(aomo_text(language,"Show partial sum in 3D"))) {
            NboOrbitalSelection selection;selection.dataset_id=data.id;
            selection.mode=NboSelectionMode::PartialSum;
            selection.semantic_kind="canonical_local_partial_sum";
            selection.group_id="mo:"+std::to_string(snapshot->focused_canonical_index)+
                ":"+nbo_orbital_kind_name(state.basis_kind);
            selection.target_canonical_index=snapshot->focused_canonical_index;
            selection.label=std::string(aomo_text(language,"Selected-component sum of "))+
                canonical_mo_display_label(canonical,snapshot->focused_canonical_index,
                    state.names&&snapshot->focused_canonical_index<state.names->canonical.size()
                        ?&state.names->canonical[snapshot->focused_canonical_index]:nullptr);
            selection.terms=terms;selection.normalize=false;
            state.pending_selection=std::move(selection);
        }
        validation::item("aomo.sum.partial");ImGui::SameLine();
        if(ImGui::Button(aomo_text(language,"Overlay terms in 3D"))) {
            NboOrbitalSelection selection;selection.dataset_id=data.id;
            selection.mode=NboSelectionMode::Overlay;
            selection.semantic_kind="canonical_local_component_overlay";
            selection.group_id="mo:"+std::to_string(snapshot->focused_canonical_index)+
                ":"+nbo_orbital_kind_name(state.basis_kind);
            selection.target_canonical_index=snapshot->focused_canonical_index;
            selection.label=std::string(aomo_text(language,"Selected-component overlay of "))+
                canonical_mo_display_label(canonical,snapshot->focused_canonical_index,
                    state.names&&snapshot->focused_canonical_index<state.names->canonical.size()
                        ?&state.names->canonical[snapshot->focused_canonical_index]:nullptr);
            selection.terms=terms;selection.normalize=false;
            state.pending_selection=std::move(selection);
        }
        validation::item("aomo.sum.overlay");
        ImGui::EndDisabled();
        if(const auto canonical=canonical_ref(data,snapshot->focused_canonical_index)) {
            if(ImGui::Button(aomo_text(language,"Show full canonical MO")))
                state.pending_selection=nbo_single_selection(data,*canonical);
            validation::item("aomo.sum.full_mo");
        }
        ImGui::TextDisabled(aomo_text(language,"Selected components: %zu / %zu"),
            terms.size(),links.size());
        const bool omitted_terms=std::any_of(links.begin(),links.end(),[&](const auto& link){
            return !state.sum_component_ids.contains(row_id(link.orbital))&&std::abs(link.coefficient)>1e-12;
        });
        const bool outside_local_span=state.basis_kind==NboOrbitalKind::NAO&&
            snapshot->focused_projection_residual_norm&&*snapshot->focused_projection_residual_norm>2e-5;
        if(omitted_terms||outside_local_span)
            ImGui::TextWrapped("%s",aomo_text(language,"Showing partial orbital composition."));
        validation::anchor("aomo.coefficients");
        ImGui::BeginChild("##aomo.coefficients",ImVec2(0,240*scale),ImGuiChildFlags_Border);
        for(std::size_t i=0;i<links.size();++i) {
            const auto& link=links[i];
            const auto* descriptor=nbo_orbital(data,link.orbital);
            const std::string label=(descriptor?descriptor->label:row_id(link.orbital))+
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
            if(ImGui::IsItemHovered())ImGui::SetTooltip("%s",aomo_text(language,"Click to show this weighted component"));
        }
        ImGui::EndChild();
        validation::field("aomo.coefficients.count",std::to_string(links.size()));
    }
    if(ImGui::Button(aomo_text(language,"Export whole diagram"))){
        state.export_requested=true;}
    validation::item("aomo.export");
    same_line_for(aomo_text(language,"Light paper export"));
    if(ImGui::Checkbox(aomo_text(language,"Light paper export"),&state.paper_export))++state.revision;
    validation::item("aomo.paper_export");
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
        result.error=aomo_text(view.language,"The diagram does not match the loaded data.");return result;
    }
    try {
        std::set<std::string> visible_nodes;
        std::set<RefKey> visible_orbitals;
        for(const auto& node:view.nodes){visible_nodes.insert(node.id);
            if(node.orbital)visible_orbitals.insert(key(*node.orbital));}
        {
            std::ofstream out(result.csv_path,std::ios::binary);
            if(!out)throw std::runtime_error(aomo_text(view.language,"Could not save the diagram data."));
            out<<std::setprecision(17);
            out<<"snapshot_id,integration_id,mo_snapshot_id,in_central_view,visible_link,basis_kind,basis_index,spin,canonical_index,coefficient,orthonormal_weight,nao_projection_weight,ao_metric_residual_norm,source_path,source_line,source_block,diagram_source_id,displayed_on_canvas,projection_strength_nonadditive,record_kind,source_energy_hartree,source_display_energy_hartree,source_display_offset_y,source_display_group_id,target_energy_hartree,target_display_energy_hartree,target_display_offset_y,target_display_group_id\n";
            const auto display_columns=[&](const NboAomoEdge* edge,
                std::optional<double> source_energy,std::size_t canonical_index) {
                const NboAomoNode* source=edge?&view.nodes[edge->source_node]:nullptr;
                const NboAomoNode* target=nullptr;
                for(const auto& node:view.nodes)if(node.canonical_index==canonical_index){target=&node;break;}
                out<<',';if(source_energy)out<<*source_energy;
                out<<',';if(source && source->display_energy_hartree)out<<*source->display_energy_hartree;
                out<<',';if(source)out<<source->display_offset_y;
                out<<','<<csv(source?source->display_group_id:"");
                out<<',';
                if(target && target->energy_hartree)out<<*target->energy_hartree;
                else for(const auto& orbital:data.orbitals)
                    if(orbital.ref.kind==NboOrbitalKind::Canonical && orbital.ref.index==canonical_index &&
                        orbital.energy_hartree){out<<*orbital.energy_hartree;break;}
                out<<',';if(target && target->display_energy_hartree)out<<*target->display_energy_hartree;
                out<<',';if(target)out<<target->display_offset_y;
                out<<','<<csv(target?target->display_group_id:"")<<'\n';
            };
            for(const auto& link:data.links){
                const bool central=in(view.central_mo_indices,link.canonical_index);
                const auto target="canonical_mo:"+std::to_string(link.canonical_index);
                const bool shown=visible_orbitals.contains(key(link.orbital)) && visible_nodes.contains(target);
                const auto edge_id="component:"+row_id(link.orbital)+":mo:"+
                    std::to_string(link.canonical_index);
                const auto edge=std::find_if(view.edges.begin(),view.edges.end(),
                    [&](const auto& candidate){return candidate.id==edge_id;});
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
                   <<csv(link.source.block)<<','<<csv(edge_id)<<','
                   <<(edge!=view.edges.end()&&edge->visible?"true":"false")<<',';
                if(edge!=view.edges.end()&&edge->projection_strength_nonadditive)
                    out<<*edge->projection_strength_nonadditive;
                out<<",raw_basis_link";
                const auto* source=nbo_orbital(data,link.orbital);
                display_columns(edge!=view.edges.end()?&*edge:nullptr,
                    source?source->energy_hartree:std::nullopt,link.canonical_index);
            }
            if(view.salc_model)for(const auto& link:view.salc_model->links) {
                if(link.side_index>=view.salc_model->orbitals.size())continue;
                const auto& orbital=view.salc_model->orbitals[link.side_index];
                const auto edge_id="salc-component:"+orbital.id+":mo:"+
                    std::to_string(link.canonical_index);
                const auto edge=std::find_if(view.edges.begin(),view.edges.end(),
                    [&](const auto& candidate){return candidate.id==edge_id;});
                out<<csv(view.id)<<','<<csv(data.id)<<','<<csv(view.mo_snapshot_id)<<','
                   <<(in(view.central_mo_indices,link.canonical_index)?"true":"false")<<','
                   <<(edge!=view.edges.end()?"true":"false")<<",SALC,"<<link.side_index<<','
                   <<csv(nbo_spin_name(orbital.spin))<<','<<link.canonical_index<<','
                   <<link.coefficient<<','<<link.weight<<',';
                if(const auto* decomposition=nbo_mo_decomposition(data.dataset,link.canonical_index)) {
                    if(decomposition->weight_sum)out<<*decomposition->weight_sum;
                }
                out<<',';
                if(const auto* decomposition=nbo_mo_decomposition(data.dataset,link.canonical_index)) {
                    if(decomposition->projection_residual_norm)out<<*decomposition->projection_residual_norm;
                }
                out<<','<<csv("")<<",0,"<<csv("")<<','<<csv(edge_id)<<','
                   <<(edge!=view.edges.end()&&edge->visible?"true":"false")
                   <<",,fixed_salc_link";
                display_columns(edge!=view.edges.end()?&*edge:nullptr,
                    orbital.energy_hartree,link.canonical_index);
            }
            if(!out)throw std::runtime_error(aomo_text(view.language,"Could not save the diagram data."));
        }
        result.csv=true;
        {
            std::ofstream out(result.json_path,std::ios::binary);
            if(!out)throw std::runtime_error(aomo_text(view.language,"Could not save the diagram data."));
            out<<std::setprecision(17);
            out<<"{\"schema\":\"cov_aomo_unified_view_v3\",\"snapshot_id\":"<<quote(view.id)
               <<",\"integration_id\":"<<quote(data.id)<<",\"mo_snapshot_id\":"<<quote(view.mo_snapshot_id)
                <<",\"mo_energy_axis_mode\":"<<quote(view.mo_energy_axis_mode)
                <<",\"display_energy_unit\":"<<quote(view.display_energy_unit)
                <<",\"mo_energy_axis_detail\":"<<quote(view.mo_energy_axis_detail)
                <<",\"preset\":"<<static_cast<int>(view.preset)
                <<",\"overview\":"<<(view.overview?"true":"false")
                <<",\"illustrative_side_layout\":"<<(view.illustrative_side_layout?"true":"false")
                <<",\"export_theme\":"<<quote(view.paper_export?"light_paper":"dark_screen")
                <<",\"auto_rydberg_expanded\":"<<(view.auto_rydberg_expanded?"true":"false")
                <<",\"canvas_width\":"<<view.canvas_width
                <<",\"canvas_height\":"<<view.canvas_height
                <<",\"qualitative_band_y\":"<<view.qualitative_band_y
                <<",\"label_font_size\":"<<view.label_font_size
                <<",\"orbital_bar_width\":"<<view.orbital_bar_width
                <<",\"orbital_bar_stroke_width\":"<<view.orbital_bar_stroke_width
                <<",\"footer_y\":"<<view.footer_y
                <<",\"numerical_zero_bound\":"<<view.numerical_zero_bound
                <<",\"numerical_zero_reason\":"<<quote(view.numerical_zero_reason)
               <<",\"basis_kind\":"<<quote(nbo_orbital_kind_name(view.basis_kind))
               <<",\"capability_status\":"<<quote(view.capability_status)
               <<",\"capability_detail\":"<<quote(view.capability_detail)
               <<",\"focused_canonical_index\":"<<view.focused_canonical_index
               <<",\"focused_projection_weight\":";
            if(view.focused_projection_weight)out<<*view.focused_projection_weight;else out<<"null";
            out<<",\"captions\":[";
            for(std::size_t i=0;i<view.captions.size();++i) {
                if(i)out<<',';const auto& caption=view.captions[i];
                out<<"{\"role\":"<<quote(caption.role)<<",\"text\":"<<quote(caption.text)
                   <<",\"x\":"<<caption.x<<",\"y\":"<<caption.y
                   <<",\"width\":"<<caption.width<<",\"height\":"<<caption.height<<'}';
            }
            out<<']';
            out<<",\"focused_projection_residual_norm\":";
            if(view.focused_projection_residual_norm)out<<*view.focused_projection_residual_norm;else out<<"null";
            out<<",\"focused_projection_residual_status\":"<<quote(
                view.focused_projection_residual_norm?"available":"unavailable");
            out<<",\"focused_display_weight\":";
            if(view.focused_display_weight)out<<*view.focused_display_weight;else out<<"null";
            out<<",\"focused_hidden_readability_weight\":";
            if(view.focused_hidden_weight)out<<*view.focused_hidden_weight;else out<<"null";
            out<<",\"focused_hidden_group_weight\":";
            if(view.focused_hidden_group_weight)out<<*view.focused_hidden_group_weight;
            else out<<"null";
            out<<",\"focused_hidden_class_weights\":{\"core\":";
            if(view.focused_hidden_core_weight)out<<*view.focused_hidden_core_weight;else out<<"null";
            out<<",\"rydberg\":";
            if(view.focused_hidden_rydberg_weight)out<<*view.focused_hidden_rydberg_weight;else out<<"null";
            out<<",\"valence\":";
            if(view.focused_hidden_valence_weight)out<<*view.focused_hidden_valence_weight;else out<<"null";
            out<<",\"hydrogen\":";
            if(view.focused_hidden_h_weight)out<<*view.focused_hidden_h_weight;else out<<"null";
            out<<"},\"hidden_numeric_zero_count\":"<<view.hidden_numeric_zero_count
               <<",\"hidden_readability_count\":"<<view.hidden_readability_count;
            out<<",\"hidden_group_count\":"<<view.hidden_group_count;
            out
               <<",\"full_numerical_file\":\".aomo.csv\",\"full_numerical_scope\":\"all_verified_links\",\"zoom\":"<<view.zoom
               <<",\"pan\":["<<view.pan_x<<','<<view.pan_y<<"],\"hidden_basis_count\":"<<view.hidden_basis_count
               <<",\"hidden_class_count\":"<<view.hidden_class_count
               <<",\"hidden_h_count\":"<<view.hidden_h_count
               <<",\"show_core\":"<<(view.show_core?"true":"false")
               <<",\"show_rydberg\":"<<(view.show_rydberg?"true":"false")
               <<",\"hide_h_orbitals\":"<<(view.hide_h_orbitals?"true":"false")
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
            out<<"],\"energy_tick_semantics\":"<<quote(view.energy_tick_semantics)
               <<",\"energy_tick_screen_color\":"<<quote(hex_rgb(view.energy_tick_screen_rgb))
               <<",\"energy_tick_export_color\":"<<quote(hex_rgb(view.energy_tick_export_rgb))
               <<",\"energy_transform\":{\"mode\":"
                <<quote(energy_axis_mode_name(view.energy_transform.mode))
                <<",\"interpolation\":\"piecewise_linear_in_energy\""
                <<",\"knot_neighbourhood_tolerance_hartree\":1e-7"
                <<",\"coordinate_min\":"<<view.axis_coordinate_min
                <<",\"coordinate_max\":"<<view.axis_coordinate_max
                <<",\"numeric_top\":"<<view.numeric_top
                <<",\"numeric_span\":"<<view.numeric_span
                <<",\"line_offset_y\":11,\"minimum_gap_weight\":"<<view.energy_transform.minimum_gap_weight
                <<",\"focus_hartree\":"<<view.energy_transform.focus_hartree
                <<",\"scale_hartree\":"<<view.energy_transform.scale_hartree
                <<",\"knots\":[";
            for(std::size_t i=0;i<view.energy_transform.knots.size();++i) {
                if(i)out<<',';const auto& knot=view.energy_transform.knots[i];
                out<<"{\"energy_hartree\":"<<knot.energy_hartree<<",\"coordinate\":"<<knot.coordinate<<'}';
            }
            out<<"]},\"lanes\":[";
            for(std::size_t i=0;i<3;++i) {
                if(i)out<<',';
                out<<"{\"x\":"<<view.lane_x[i]<<",\"width\":"<<view.lane_width[i]<<'}';
            }
            out<<"],\"energy_ticks\":[";
            for(std::size_t i=0;i<view.energy_ticks.size();++i) {
                if(i)out<<',';
                out<<"{\"energy_hartree\":"<<view.energy_ticks[i].energy_hartree
                   <<",\"y\":"<<view.energy_ticks[i].y<<'}';
            }
            out<<"],\"nodes\":[";
            for(std::size_t i=0;i<view.nodes.size();++i) {
                if(i)out<<',';const auto& node=view.nodes[i];
                out<<"{\"id\":"<<quote(node.id)<<",\"label\":"<<quote(node.label)
                   <<",\"detail\":"<<quote(node.detail)
                   <<",\"symmetry_irrep\":"<<quote(node.symmetry_irrep)
                   <<",\"symmetry_ordinal\":"<<node.symmetry_ordinal
                   <<",\"symmetry_name_verified\":"<<(node.symmetry_name_verified?"true":"false")
                   <<",\"name_detail\":"<<quote(node.name_detail)
                   <<",\"individual_label\":"<<quote(node.individual_label)
                   <<",\"spatial_pair_id\":"<<quote(node.spatial_pair_id)
                   <<",\"spatial_pair_label\":"<<quote(node.spatial_pair_label)
                   <<",\"display_group_id\":"<<quote(node.display_group_id)
                   <<",\"display_energy_semantics\":"<<quote(node.display_energy_semantics)
                   <<",\"display_offset_y\":"<<node.display_offset_y
                   <<",\"shell_label\":"<<quote(node.shell_label)
                   <<",\"shell_member_index\":"<<node.shell_member_index
                   <<",\"shell_member_count\":"<<node.shell_member_count
                   <<",\"subspace_id\":"<<quote(node.subspace_id)
                   <<",\"energy_semantics\":"<<quote(node.energy_semantics)
                   <<",\"x\":"<<node.x<<",\"y\":"<<node.y
                   <<",\"label_y\":"<<node.label_y
                   <<",\"label_bbox\":{\"x\":"<<node.label_x<<",\"y\":"<<node.label_y
                   <<",\"width\":"<<node.label_width<<",\"height\":"<<node.label_height<<'}'
                   <<",\"occupation_label\":"<<quote(node.occupation_label)
                   <<",\"occupation_on_bar\":"<<(node.occupation_on_bar?"true":"false")
                   <<",\"occupation_bbox\":{\"x\":"<<node.occupation_x<<",\"y\":"<<node.occupation_y
                   <<",\"width\":"<<node.occupation_width<<",\"height\":"<<node.occupation_height<<'}'
                   <<",\"renders_energy_line\":"<<(node.group_header?"false":"true")
                   <<",\"lane\":"<<static_cast<int>(node.lane)
                   <<",\"quantitative_energy\":"<<(node.quantitative_energy?"true":"false")
                   <<",\"width\":"<<node.width<<",\"height\":"<<node.height
                   <<",\"group_header\":"<<(node.group_header?"true":"false")
                   <<",\"available\":"<<(node.available?"true":"false")
                   <<",\"composition_available\":"<<(node.composition_available?"true":"false")
                   <<",\"atoms0\":[";
                for(std::size_t j=0;j<node.atoms.size();++j){if(j)out<<',';out<<node.atoms[j];}
                out<<"],\"member_canonical_indices\":[";
                for(std::size_t j=0;j<node.member_canonical_indices.size();++j){
                    if(j)out<<',';out<<node.member_canonical_indices[j];
                }
                out<<"],\"orbital\":";
                if(node.orbital)out<<"{\"kind\":"<<quote(nbo_orbital_kind_name(node.orbital->kind))
                    <<",\"spin\":"<<quote(nbo_spin_name(node.orbital->spin))
                    <<",\"index\":"<<node.orbital->index<<"}";else out<<"null";
                out<<",\"canonical_index\":";
                if(node.canonical_index)out<<*node.canonical_index;else out<<"null";
                out<<",\"salc_index\":";
                if(node.salc_index)out<<*node.salc_index;else out<<"null";
                out<<",\"energy_hartree\":";
                if(node.energy_hartree)out<<*node.energy_hartree;else out<<"null";
                out<<",\"display_energy_hartree\":";
                if(node.display_energy_hartree)out<<*node.display_energy_hartree;else out<<"null";
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
                    <<",\"coefficient\":"<<edge.coefficient
                    <<",\"visible\":"<<(edge.visible?"true":"false")
                    <<",\"weight\":";
                if(edge.weight)out<<*edge.weight;else out<<"null";
                out<<",\"projection_strength_nonadditive\":";
                if(edge.projection_strength_nonadditive)
                    out<<*edge.projection_strength_nonadditive;else out<<"null";
                out<<",\"source_path\":"<<quote(edge.source.path)<<",\"source_line\":"<<edge.source.line_begin
                   <<",\"external_factor_only\":"<<(view.nodes.at(edge.source_node).fragment_group_id?"true":"false")<<'}';
            }
            out<<"],\"salc_model\":";
            if(view.salc_model)out<<serialize_nbo_salc_json(*view.salc_model);
            else out<<"null";
            out<<",\"selection\":";
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
            if(!out)throw std::runtime_error(aomo_text(view.language,"Could not save the diagram data."));
        }
        result.json=true;
        float max_y=0,max_x=0;
        for(const auto& node:view.nodes){max_y=std::max(max_y,node.label_y+node.label_height+20);
            max_x=std::max(max_x,std::max(node.x+node.width,
                node.label_x+node.label_width)+20);}
        const auto width=std::max(static_cast<int>(std::ceil(view.canvas_width)),
            std::max(540,static_cast<int>(std::ceil(max_x))));
        const auto height=std::max(480,static_cast<int>(std::ceil(std::max(max_y,view.canvas_height))));
        const GraphAppearance appearance(view);
        const bool paper=view.paper_export;
        {
            std::ofstream out(result.svg_path,std::ios::binary);
            if(!out)throw std::runtime_error(aomo_text(view.language,"Could not save the diagram image."));
            out<<std::setprecision(9);
            out<<"<svg xmlns=\"http://www.w3.org/2000/svg\" font-family=\"Segoe UI, Arial, sans-serif\" width=\""<<width
               <<"\" height=\""<<height<<"\" viewBox=\"0 0 "<<width<<' '<<height<<"\">\n";
            out<<"<rect width=\"100%\" height=\"100%\" fill=\""
               <<(paper?"#ffffff":"#18202d")<<"\"/>\n";
            const bool has_nonquant=std::any_of(view.nodes.begin(),view.nodes.end(),
                [](const auto& node){return node.lane!=NboAomoLane::Centre &&
                    !node.quantitative_energy;});
            if(has_nonquant) {
                for(const auto lane:{0,2}) {
                    out<<"<rect x=\""<<view.lane_x[lane]-8<<"\" y=\""<<view.qualitative_band_y
                       <<"\" width=\""<<view.lane_width[lane]+8<<"\" height=\""
                       <<view.footer_y-8-view.qualitative_band_y<<"\" fill=\""
                       <<(paper?"#eef2f7":"#232934")<<"\"/>\n";
                }
            }
            for(const auto& caption:view.captions)
                out<<"<text class=\"caption\" data-role=\""<<xml(caption.role)
                   <<"\" x=\""<<caption.x<<"\" y=\""<<caption.y
                   <<"\" dominant-baseline=\"text-before-edge\" fill=\""
                   <<(caption.role=="axis" && view.energy_transform.mode==EnergyAxisMode::NonlinearFocus?
                       hex_rgb(view.energy_tick_export_rgb):(paper?"#516579":"#a0b4ca"))
                   <<"\" font-size=\""<<view.label_font_size<<"\">"<<xml(caption.text)<<"</text>\n";
            for(const auto& tick:view.energy_ticks)
                out<<"<g><line x1=\"20\" x2=\""<<width-20<<"\" y1=\""
                   <<tick.y+11<<"\" y2=\""<<tick.y+11
                   <<"\" stroke=\""<<(paper?"#aab9c8":"#48576c")
                   <<"\" stroke-opacity=\"0.12\"/>"
                   <<"<text x=\"8\" y=\""<<tick.y+11-view.label_font_size
                   <<"\" dominant-baseline=\"text-before-edge\" fill=\""
                   <<hex_rgb(view.energy_tick_export_rgb)
                   <<"\" font-size=\""<<view.label_font_size<<"\">"
                   <<xml(format_energy(tick.energy_hartree,view.energy_unit,3))
                   <<"</text></g>\n";
            for(const auto& edge:view.edges) {
                if(!edge.visible)continue;
                const auto& a=view.nodes.at(edge.source_node);const auto& b=view.nodes.at(edge.target_node);
                const auto style=appearance.edge(edge);
                const char* stroke=edge.coefficient<0?
                    (paper?"#b5444c":"#f2918f"):
                    (paper?"#177cae":"#4bb4ea");
                out<<"<line x1=\""<<(a.lane==NboAomoLane::Right?a.x:a.x+node_width(a))
                   <<"\" y1=\""<<a.y+a.height*0.5f
                   <<"\" x2=\""<<(a.lane==NboAomoLane::Right?b.x+node_width(b):b.x)
                   <<"\" y2=\""<<b.y+b.height*0.5f
                   <<"\" stroke=\""<<stroke<<"\" stroke-opacity=\""
                   <<static_cast<double>(style.alpha)/255.0<<"\" stroke-width=\""<<style.width
                   <<"\" stroke-dasharray=\"5 5\"><title>";
                out<<xml(aomo_text(view.language,"Coefficient: "))<<edge.coefficient;
                if(a.fragment_group_id && a.metric_norm2)
                    out<<"; "<<xml(aomo_text(view.language,"Component norm squared: "))<<number(*a.metric_norm2);
                else out<<" "<<xml(source_name(edge.source));
                out<<"</title></line>\n";
            }
            for_each_node_text_background(view,[&](const std::string& node_id,const char* role,
                float x,float y,float width,float height,bool qualitative) {
                out<<"<rect class=\"text-background\" data-node=\""<<xml(node_id)
                   <<"\" data-role=\""<<role<<"\" x=\""<<x<<"\" y=\""<<y
                   <<"\" width=\""<<width<<"\" height=\""<<height
                   <<"\" fill=\""<<(qualitative?(paper?"#eef2f7":"#232934"):
                       (paper?"#ffffff":"#18202d"))<<"\" fill-opacity=\"1\"/>\n";
            });
            for(const auto& node:view.nodes) {
                const float strength=appearance.node_strength(node);
                const int green=static_cast<int>(58+65*strength),blue=static_cast<int>(82+84*strength);
                const auto color=paper?
                    (strength>=0.99f?"#126b9d":node.group_header?"#4b6076":"#24384b"):
                    (strength>=0.99f?"#7adfff":node.group_header?"#8fa1bb":"#e4eff9");
                const float ly=node.y+node.height*0.5f;
                const bool folded_group=node.group_header;
                out<<"<g id=\""<<xml(node.id)<<"\" data-display-group=\""<<xml(node.display_group_id)
                   <<"\" data-individual-label=\""<<xml(node.individual_label)
                   <<"\" data-spatial-pair=\""<<xml(node.spatial_pair_id)
                   <<"\" data-energy-hartree=\""<<(node.energy_hartree?number(*node.energy_hartree):"")
                   <<"\" data-display-energy-hartree=\""<<(node.display_energy_hartree?number(*node.display_energy_hartree):"")
                   <<"\" data-display-offset-y=\""<<node.display_offset_y<<"\"><title>"<<xml(node.individual_label.empty()?node.label:node.individual_label)
                   <<"</title>";
                if(folded_group)
                    out<<"<rect x=\""<<node.x<<"\" y=\""<<node.y
                       <<"\" width=\""<<node.width<<"\" height=\""<<node.height
                       <<"\" rx=\"4\" fill=\""<<(paper?"#e0e8ef":"#39465b")
                       <<"\" stroke=\""<<(paper?"#7d91a5":"#6e89a6")<<"\"/>";
                else out<<"<line x1=\""<<node.x<<"\" y1=\""<<ly
                   <<"\" x2=\""<<node.x+node.width<<"\" y2=\""<<ly
                   <<"\" stroke=\""<<color<<"\" stroke-width=\""<<view.orbital_bar_stroke_width<<"\"/>";
                out<<"<text x=\""<<node.label_x
                   <<"\" y=\""<<node.label_y
                   <<"\" dominant-baseline=\"text-before-edge\" textLength=\""<<node.label_width
                   <<"\" lengthAdjust=\"spacingAndGlyphs\" fill=\""<<color<<"\" font-size=\""<<view.label_font_size<<"\">"
                   <<xml(node.label)<<"</text>";
                if(node.occupation_on_bar) {
                    out<<"<g class=\"electron-arrows\" data-occupation=\""<<xml(node.occupation_label)
                        <<"\" stroke=\""<<color<<"\" stroke-width=\"1.5\">";
                    electron_strokes(node,[&](float x1,float y1,float x2,float y2) {
                        out<<"<line x1=\""<<x1<<"\" y1=\""<<y1<<"\" x2=\""<<x2
                            <<"\" y2=\""<<y2<<"\"/>";
                    });
                    out<<"</g>";
                } else if(!node.occupation_label.empty())
                    out<<"<text x=\""<<node.occupation_x
                       <<"\" y=\""<<node.occupation_y
                       <<"\" dominant-baseline=\"text-before-edge\" textLength=\""<<node.occupation_width
                       <<"\" lengthAdjust=\"spacingAndGlyphs\" fill=\""<<(paper?"#21364a":"#e8ebf1")
                       <<"\" font-size=\""<<view.label_font_size<<"\">"
                       <<xml(node.occupation_label)<<"</text>";
                out<<"</g>\n";
            }
            out<<"</svg>\n";
            if(!out)throw std::runtime_error(aomo_text(view.language,"Could not save the diagram image."));
        }
        result.svg=true;
        // PNG uses the same frozen graph geometry and IDs as SVG/JSON.
        if(height>16000)throw std::runtime_error(aomo_text(view.language,"The image exceeds the supported height."));
        std::vector<unsigned char> rgba(static_cast<std::size_t>(width)*height*4,255);
        for(std::size_t p=0;p<rgba.size();p+=4){
            rgba[p]=paper?255:24;rgba[p+1]=paper?255:32;rgba[p+2]=paper?255:45;}
        auto pixel=[&](int x,int y,unsigned char r,unsigned char g,unsigned char b){
            if(x<0||y<0||x>=width||y>=height)return;
            const auto p=(static_cast<std::size_t>(y)*width+x)*4;rgba[p]=r;rgba[p+1]=g;rgba[p+2]=b;
        };
        // Raster text uses the same live ImGui font atlas as the on-screen
        // canvas, including the loaded CJK and scientific glyphs.
        unsigned char* atlas=nullptr;int atlas_w=0,atlas_h=0;
        ImGui::GetIO().Fonts->GetTexDataAsAlpha8(&atlas,&atlas_w,&atlas_h);
        const ImFont* font=ImGui::GetFont();
        std::size_t missing_glyphs=0;
        auto draw_text=[&](float x,float y,const std::string& value,
                           float size,unsigned char r,unsigned char g,unsigned char b) {
            if(!font||!atlas||atlas_w<=0||atlas_h<=0)return;
            const float scale=size/std::max(1.0f,font->FontSize);
            float cursor=x;
            for(std::size_t pos=0;pos<value.size();) {
                const auto byte=static_cast<unsigned char>(value[pos]);
                unsigned codepoint=byte;std::size_t length=1;
                if((byte&0xe0)==0xc0){codepoint=byte&0x1f;length=2;}
                else if((byte&0xf0)==0xe0){codepoint=byte&0x0f;length=3;}
                else if((byte&0xf8)==0xf0){codepoint=byte&0x07;length=4;}
                if(pos+length>value.size()){length=1;codepoint='?';}
                else for(std::size_t k=1;k<length;++k)
                    codepoint=(codepoint<<6)|(static_cast<unsigned char>(value[pos+k])&0x3f);
                pos+=length;
                const auto* glyph=font->FindGlyphNoFallback(static_cast<ImWchar>(codepoint));
                if(!glyph){++missing_glyphs;glyph=font->FindGlyph('?');}
                if(!glyph)continue;
                const int x0=static_cast<int>(std::floor(cursor+glyph->X0*scale));
                const int x1=static_cast<int>(std::ceil(cursor+glyph->X1*scale));
                const int y0=static_cast<int>(std::floor(y+glyph->Y0*scale));
                const int y1=static_cast<int>(std::ceil(y+glyph->Y1*scale));
                const int sx0=static_cast<int>(glyph->U0*atlas_w);
                const int sx1=static_cast<int>(glyph->U1*atlas_w);
                const int sy0=static_cast<int>(glyph->V0*atlas_h);
                const int sy1=static_cast<int>(glyph->V1*atlas_h);
                for(int yy=std::max(0,y0);yy<std::min(height,y1);++yy)
                    for(int xx=std::max(0,x0);xx<std::min(width,x1);++xx) {
                        const int sx=std::clamp(sx0+(xx-x0)*(sx1-sx0)/std::max(1,x1-x0),0,atlas_w-1);
                        const int sy=std::clamp(sy0+(yy-y0)*(sy1-sy0)/std::max(1,y1-y0),0,atlas_h-1);
                        const float alpha=static_cast<float>(atlas[sy*atlas_w+sx])/255.0f;
                        const auto p=(static_cast<std::size_t>(yy)*width+xx)*4;
                        rgba[p]=static_cast<unsigned char>(rgba[p]*(1-alpha)+r*alpha);
                        rgba[p+1]=static_cast<unsigned char>(rgba[p+1]*(1-alpha)+g*alpha);
                        rgba[p+2]=static_cast<unsigned char>(rgba[p+2]*(1-alpha)+b*alpha);
                    }
                cursor+=glyph->AdvanceX*scale;
            }
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
        if(std::any_of(view.nodes.begin(),view.nodes.end(),[](const auto& node){
            return node.lane!=NboAomoLane::Centre && !node.quantitative_energy;})) {
            for(const auto lane:{0,2}) {
                for(int yy=static_cast<int>(view.qualitative_band_y);yy<view.footer_y-8;++yy)
                    for(int xx=static_cast<int>(view.lane_x[lane]-8);
                        xx<view.lane_x[lane]+view.lane_width[lane];++xx)
                        if(paper)pixel(xx,yy,238,242,247);
                        else pixel(xx,yy,35,41,52);
            }
        }
        for(const auto& tick:view.energy_ticks) {
            EdgeAppearance grid;
            grid.red=paper?170:72;grid.green=paper?185:87;
            grid.blue=paper?200:108;
            grid.alpha=32;grid.width=1;
            line(20,tick.y+11,width-20,tick.y+11,grid);
        }
        for(const auto& edge:view.edges) {
            if(!edge.visible)continue;
            const auto& a=view.nodes.at(edge.source_node);const auto& b=view.nodes.at(edge.target_node);
            auto style=appearance.edge(edge);
            if(paper) {
                style.red=edge.coefficient<0?181:23;
                style.green=edge.coefficient<0?68:124;
                style.blue=edge.coefficient<0?76:174;
            }
            dashed_segments(a.lane==NboAomoLane::Right?a.x:a.x+node_width(a),
                a.y+a.height*0.5f,
                a.lane==NboAomoLane::Right?b.x+node_width(b):b.x,
                b.y+b.height*0.5f,
                [&](float x1,float y1,float x2,float y2){line(x1,y1,x2,y2,style);});
        }
        for_each_node_text_background(view,[&](const std::string&,const char*,
            float x,float y,float w,float h,bool qualitative) {
            const unsigned char red=qualitative?(paper?238:35):(paper?255:24);
            const unsigned char green=qualitative?(paper?242:41):(paper?255:32);
            const unsigned char blue=qualitative?(paper?247:52):(paper?255:45);
            for(int yy=static_cast<int>(std::floor(y));yy<static_cast<int>(std::ceil(y+h));++yy)
                for(int xx=static_cast<int>(std::floor(x));xx<static_cast<int>(std::ceil(x+w));++xx)
                    pixel(xx,yy,red,green,blue);
        });
        for(const auto& node:view.nodes) {
            const int x=static_cast<int>(node.x),y=static_cast<int>(node.y);
            const int nw=static_cast<int>(std::ceil(node.width));
            const int nh=static_cast<int>(std::ceil(node.height));
            const float strength=appearance.node_strength(node);
            const unsigned char red=paper?(strength>=0.99f?18:node.group_header?75:36):
                (strength>=0.99f?122:node.group_header?143:228);
            const unsigned char green=paper?(strength>=0.99f?107:node.group_header?96:56):
                (strength>=0.99f?223:node.group_header?161:239);
            const unsigned char blue=paper?(strength>=0.99f?157:node.group_header?118:75):
                (strength>=0.99f?255:node.group_header?187:249);
            const int line_y=y+nh/2;
            const bool folded_group=node.group_header;
            if(folded_group) {
                for(int yy=y;yy<y+nh;++yy)for(int xx=x;xx<x+nw;++xx)
                    if(paper)pixel(xx,yy,224,232,239);
                    else pixel(xx,yy,57,70,91);
            } else {
                const float centre=node.y+node.height*0.5f;
                const float lo=centre-view.orbital_bar_stroke_width*0.5f;
                const float hi=centre+view.orbital_bar_stroke_width*0.5f;
                EdgeAppearance stroke;stroke.red=red;stroke.green=green;stroke.blue=blue;
                for(int yy=static_cast<int>(std::floor(lo));yy<static_cast<int>(std::ceil(hi));++yy) {
                    const float coverage=std::max(0.0f,std::min(hi,static_cast<float>(yy+1))-
                        std::max(lo,static_cast<float>(yy)));
                    stroke.alpha=static_cast<unsigned char>(std::round(255*coverage));
                    for(int xx=x;xx<x+nw;++xx)blend(xx,yy,stroke);
                }
            }
            draw_text(node.label_x,node.label_y,
                node.label,view.label_font_size,red,green,blue);
            if(node.occupation_on_bar) {
                EdgeAppearance stroke;stroke.red=red;stroke.green=green;stroke.blue=blue;
                stroke.alpha=255;stroke.width=1.5f;
                electron_strokes(node,[&](float x1,float y1,float x2,float y2){line(x1,y1,x2,y2,stroke);});
            } else if(!node.occupation_label.empty())
                draw_text(node.occupation_x,node.occupation_y,
                    node.occupation_label,view.label_font_size,
                    paper?33:232,paper?54:235,paper?74:241);
        }
        for(const auto& caption:view.captions)
            draw_text(caption.x,caption.y,caption.text,view.label_font_size,
                caption.role=="axis" && view.energy_transform.mode==EnergyAxisMode::NonlinearFocus?
                    view.energy_tick_export_rgb[0]:(paper?81:160),
                caption.role=="axis" && view.energy_transform.mode==EnergyAxisMode::NonlinearFocus?
                    view.energy_tick_export_rgb[1]:(paper?101:180),
                caption.role=="axis" && view.energy_transform.mode==EnergyAxisMode::NonlinearFocus?
                    view.energy_tick_export_rgb[2]:(paper?121:202));
        for(const auto& tick:view.energy_ticks)
            draw_text(8,tick.y+11-view.label_font_size,
                format_energy(tick.energy_hartree,view.energy_unit,3),view.label_font_size,
                view.energy_tick_export_rgb[0],view.energy_tick_export_rgb[1],
                view.energy_tick_export_rgb[2]);
        if(missing_glyphs)
            draw_text(12,height-20,aomo_text(view.language,"Some characters could not be drawn."),
                view.label_font_size,255,181,100);
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
        if(!out)throw std::runtime_error(aomo_text(view.language,"Could not save the diagram image."));
        out.write(reinterpret_cast<const char*>(png.data()),static_cast<std::streamsize>(png.size()));
        if(!out)throw std::runtime_error(aomo_text(view.language,"Could not save the diagram image."));
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
