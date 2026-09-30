#include "cov/nbo_aomo_labels.hpp"
#include "cov/d2h_orbital_characters.hpp"
#include "cov/open_profile.hpp"
#include "cov/orbital_symmetry_scope.hpp"
#include "cov/orbital_symmetry.hpp"
#include "cov/point_group_catalog.hpp"
#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>

namespace cov::ui {
namespace {
constexpr double character_tolerance=2e-4;
constexpr double metric_tolerance=2e-5;
constexpr double energy_order_tolerance=2e-5;
thread_local std::size_t canonical_name_revision=0;
using Vec=std::array<double,3>;
using Mat=std::array<double,9>;
double dot(const Vec& a,const Vec& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double det(const Mat& a){return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);}
double trace(const Mat& a){return a[0]+a[4]+a[8];}
double matrix_error(const Mat& a,const Mat& b){double e=0;for(std::size_t i=0;i<9;++i)e=std::max(e,std::abs(a[i]-b[i]));return e;}
Mat multiply(const Mat& a,const Mat& b){Mat c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];return c;}
const Mat identity{1,0,0,0,1,0,0,0,1};
const Mat inversion{-1,0,0,0,-1,0,0,0,-1};
std::string normalized(std::string s){std::string o;for(unsigned char c:s)if(!std::isspace(c)&&c!='_')o+=char(std::tolower(c));return o;}
std::string finite_irrep(const std::string& group,const std::string& label){
    const auto* pg=find_point_group(group);if(!pg)return {};
    for(const auto& ir:pg->irreps)if(normalized(std::string(ir.label))==normalized(label))return std::string(ir.label);
    return {};
}
std::size_t dimension(const std::string& group,const std::string& label){
    const auto* pg=find_point_group(group);if(pg)for(const auto& ir:pg->irreps)if(normalized(std::string(ir.label))==normalized(label))return ir.dimension;
    return 0;
}
struct CanonicalActionCache {
    const Wavefunction* wavefunction=nullptr;
    std::vector<double> packed,metric_packed;
    std::vector<std::vector<double>> transformed,metric_transformed;
};
struct Frame {
    std::string group, detail;
    std::vector<SymmetryOperation> ops;
    std::size_t e=0, inv=0, rotation=0, mirror0=0, mirror1=0, yz=0;
    bool valid=false, mirrors_named=false, linear=false;
    Vec axis{};
    unsigned lmax=0, rotation_order=0;
    orbital_characters::D2hAxes d2h_axes;
    mutable std::shared_ptr<CanonicalActionCache> canonical_action;
};
Frame frame_for(const Wavefunction& w,const NboSalcModel* model){
    Frame f;
    if(!model||!model->group_verified)return f;
    f.group=model->point_group;f.linear=f.group=="Dinfh"||f.group=="Cinfv";
    if(model->used_group!=f.group&&(!f.linear||model->used_group!="finite sampling subgroup of "+f.group))return f;
    f.ops=model->operations;const auto n=f.ops.size();const auto* pg=find_point_group(f.group);
    // The native Dnh classifier supports general n, beyond the intentionally
    // smaller central-metal display catalogue (which omits e.g. D6h).
    if(!f.linear&&f.group.size()>2&&(f.group.front()=='C'||f.group.front()=='D')){
        std::size_t i=1;while(i+1<f.group.size()&&std::isdigit(static_cast<unsigned char>(f.group[i]))){f.rotation_order=10*f.rotation_order+unsigned(f.group[i]-'0');++i;}if(i+1!=f.group.size()||f.rotation_order<2||f.rotation_order>64)f.rotation_order=0;
    }
    const auto expected_order=pg?std::size_t(pg->order):f.rotation_order&&f.group.front()=='D'&&f.group.back()=='h'?4*f.rotation_order:f.rotation_order&&f.group.front()=='C'&&f.group.back()=='v'?2*f.rotation_order:0;
    if(!n||(!f.linear&&(!expected_order||expected_order!=n)))return f;
    if(f.linear){
        const auto geometry=analyse_molecular_symmetry(w);
        if(!geometry.linear||geometry.point_group!=f.group)return f;
        double best=0;for(const auto& a:w.atoms){Vec v{a.x-geometry.centre_bohr[0],a.y-geometry.centre_bohr[1],a.z-geometry.centre_bohr[2]};if(dot(v,v)>best){best=dot(v,v);f.axis=v;}}
        if(best<1e-20)return f;for(auto& v:f.axis)v/=std::sqrt(best);
        for(const auto& shell:w.shells)f.lmax=std::max(f.lmax,unsigned(shell.angular_momentum));
        // Sampling must resolve every supported angular frequency, without aliasing.
        if(f.lmax>4||n<(f.group=="Dinfh"?4:2)*std::max(3u,2*f.lmax+1))return f;
    }
    // SALC closure products intentionally carry only matrix/permutation data;
    // kind/order/axis metadata must not be used here.
    std::vector<std::size_t> mirrors,rotations;bool have_e=false,have_i=false;
    for(std::size_t i=0;i<n;++i){const auto& op=f.ops[i];
        if(op.atom_permutation.size()!=w.atoms.size())return f;
        std::set<std::size_t> seen;for(std::size_t a=0;a<w.atoms.size();++a){const auto b=op.atom_permutation[a];if(b>=w.atoms.size()||w.atoms[a].atomic_number!=w.atoms[b].atomic_number||!seen.insert(b).second)return f;}
        Mat trans{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)trans[3*a+b]=op.matrix[3*b+a];
        if(matrix_error(multiply(trans,op.matrix),identity)>metric_tolerance)return f;
        if(matrix_error(op.matrix,identity)<metric_tolerance){f.e=i;have_e=true;}
        else if(matrix_error(op.matrix,inversion)<metric_tolerance){f.inv=i;have_i=true;}
        else if(std::abs(det(op.matrix)+1)<metric_tolerance&&std::abs(trace(op.matrix)-1)<metric_tolerance)mirrors.push_back(i);
        else if(std::abs(det(op.matrix)-1)<metric_tolerance&&std::abs(trace(op.matrix)+1)<metric_tolerance)rotations.push_back(i);
    }
    if(!have_e)return f;
    // Recheck matrix closure, not just the advertised point-group string.
    for(const auto& a:f.ops)for(const auto& b:f.ops){const auto ab=multiply(a.matrix,b.matrix);if(std::none_of(f.ops.begin(),f.ops.end(),[&](const auto& c){return matrix_error(ab,c.matrix)<metric_tolerance;}))return f;}
    f.valid=true;
    if(f.group=="C1")f.valid=n==1;
    else if(f.group=="Ci")f.valid=n==2&&have_i;
    else if(f.group=="Cs"){f.valid=n==2&&mirrors.size()==1;if(f.valid)f.mirror0=mirrors[0];}
    else if(f.group=="C2"||f.group=="C2h"){f.valid=rotations.size()==1&&(f.group=="C2"||(have_i&&mirrors.size()==1));if(f.valid)f.rotation=rotations[0];}
    else if(f.group=="C2v"){
        f.valid=n==4&&mirrors.size()==2&&rotations.size()==1;if(!f.valid)return f;
        f.mirror0=mirrors[0];f.mirror1=mirrors[1];f.rotation=rotations[0];
        // Name the molecular plane yz by unique maximal on-plane atom support.
        // This covaries under rigid rotations and atom permutations. A tie is
        // explicitly unresolved: A1/A2 still work, but B1/B2 must not be guessed.
        std::size_t support[2]{};Vec normals[2]{};
        for(int k=0;k<2;++k){const auto& op=f.ops[mirrors[k]];
            for(std::size_t a=0;a<w.atoms.size();++a)if(op.atom_permutation[a]==a)++support[k];
            double best=0;for(int r=0;r<3;++r){Vec v{};for(int c=0;c<3;++c)v[c]=identity[3*r+c]-op.matrix[3*r+c];if(dot(v,v)>best){best=dot(v,v);normals[k]=v;}}
            if(best>0)for(auto& x:normals[k])x/=std::sqrt(best);
        }
        if(std::abs(dot(normals[0],normals[1]))>metric_tolerance){f.valid=false;return f;}
        f.mirrors_named=support[0]!=support[1];
        if(f.mirrors_named){f.yz=mirrors[support[0]>support[1]?0:1];const auto& normal=normals[support[0]>support[1]?0:1];std::ostringstream s;s<<"C2v: yz is the uniquely maximal fixed-atom mirror; x normal in input frame=("<<normal[0]<<','<<normal[1]<<','<<normal[2]<<"); z is C2 axis; B1 transforms as x, B2 as y";f.detail=s.str();}
        else f.detail="C2v mirror-axis convention unresolved (equal fixed-atom support); B1/B2 not assigned";
    }
    if(f.group=="D2h"){const auto geometry=analyse_molecular_symmetry(w);f.d2h_axes=orbital_characters::d2h_axes(w,f.ops,geometry.centre_bohr,geometry.tolerance_bohr);f.detail=f.d2h_axes.detail;}
    if(f.detail.empty())f.detail=f.linear?"Validated angular-momentum-resolving sampling of "+f.group:"Validated full finite "+f.group+" operation matrices";
    return f;
}
Frame canonical_frame(const Wavefunction& w){
    // Use the same matrix, closure, character and metric gates as an attached
    // SALC frame. A geometry label alone does not certify orbital symmetry.
    const auto geometry=analyse_molecular_symmetry(w);
    NboSalcModel frame;
    frame.point_group=frame.used_group=geometry.point_group;
    frame.group_verified=!geometry.operations.empty();
    frame.operations=geometry.operations;
    if(!geometry.linear&&frame.group_verified){
        // Geometry supplies validated generators, not necessarily every group
        // element. Complete their finite closure before the full-frame gate.
        SymmetryOperation e;e.atom_permutation.resize(w.atoms.size());
        std::iota(e.atom_permutation.begin(),e.atom_permutation.end(),0);
        frame.operations={e};
        for(std::size_t i=0;i<frame.operations.size();++i)for(const auto& generator:geometry.operations){
            if(generator.atom_permutation.size()!=w.atoms.size()){frame.group_verified=false;return frame_for(w,&frame);}
            SymmetryOperation product;product.matrix=multiply(frame.operations[i].matrix,generator.matrix);
            product.atom_permutation.resize(w.atoms.size());
            for(std::size_t a=0;a<w.atoms.size();++a){const auto b=generator.atom_permutation[a];if(b>=w.atoms.size()){frame.group_verified=false;return frame_for(w,&frame);}product.atom_permutation[a]=frame.operations[i].atom_permutation[b];}
            if(std::any_of(frame.operations.begin(),frame.operations.end(),[&](const auto& op){return op.atom_permutation==product.atom_permutation&&matrix_error(op.matrix,product.matrix)<metric_tolerance;}))continue;
            if(frame.operations.size()>=256){frame.group_verified=false;return frame_for(w,&frame);}
            double mapping_error=0;
            for(std::size_t a=0;a<w.atoms.size();++a){const auto& from=w.atoms[a];const auto& to=w.atoms[product.atom_permutation[a]];const Vec x{from.x-geometry.centre_bohr[0],from.y-geometry.centre_bohr[1],from.z-geometry.centre_bohr[2]},y{to.x-geometry.centre_bohr[0],to.y-geometry.centre_bohr[1],to.z-geometry.centre_bohr[2]};Vec difference{};for(int r=0;r<3;++r){difference[r]=-y[r];for(int c=0;c<3;++c)difference[r]+=product.matrix[3*r+c]*x[c];}mapping_error=std::max(mapping_error,std::sqrt(dot(difference,difference)));}
            if(mapping_error>geometry.tolerance_bohr){frame.group_verified=false;return frame_for(w,&frame);}
            product.max_mapping_error_bohr=mapping_error;frame.operations.push_back(std::move(product));
        }
    }
    if(geometry.linear){
        Vec axis{};double length2=0;
        for(const auto& atom:w.atoms){
            Vec v{atom.x-geometry.centre_bohr[0],atom.y-geometry.centre_bohr[1],atom.z-geometry.centre_bohr[2]};
            if(dot(v,v)>length2){length2=dot(v,v);axis=v;}
        }
        if(length2<=1e-20)return frame_for(w,&frame);
        for(auto& x:axis)x/=std::sqrt(length2);
        const Vec seed=std::abs(axis[0])<.8?Vec{1,0,0}:Vec{0,1,0};
        Vec normal{axis[1]*seed[2]-axis[2]*seed[1],axis[2]*seed[0]-axis[0]*seed[2],axis[0]*seed[1]-axis[1]*seed[0]};
        const double normal_length=std::sqrt(dot(normal,normal));
        for(auto& x:normal)x/=normal_length;
        Mat reflection=identity;
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)reflection[3*a+b]-=2*normal[a]*normal[b];
        unsigned lmax=0;for(const auto& shell:w.shells)lmax=std::max(lmax,unsigned(shell.angular_momentum));
        // Match the validated SALC sampling resolution; no angular order is
        // inferred from molecule identity or an orbital number.
        const unsigned order=std::max(3u,2*lmax+1);
        const auto inversion_op=std::find_if(geometry.operations.begin(),geometry.operations.end(),
            [](const auto& op){return op.kind==SymmetryOperationKind::Inversion;});
        const bool centrosymmetric=geometry.point_group=="Dinfh";
        if(centrosymmetric&&inversion_op==geometry.operations.end())return frame_for(w,&frame);
        frame.operations.clear();
        const Mat cross{0,-axis[2],axis[1],axis[2],0,-axis[0],-axis[1],axis[0],0};
        for(unsigned k=0;k<order;++k){
            const double angle=2*3.14159265358979323846*k/order,c=std::cos(angle),s=std::sin(angle);
            Mat rotation{};
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                rotation[3*a+b]=c*identity[3*a+b]+(1-c)*axis[a]*axis[b]+s*cross[3*a+b];
            for(int mirror=0;mirror<2;++mirror)for(int invert=0;invert<(centrosymmetric?2:1);++invert){
                SymmetryOperation op;op.matrix=mirror?multiply(rotation,reflection):rotation;
                if(invert){for(auto& x:op.matrix)x=-x;op.atom_permutation=inversion_op->atom_permutation;}
                else{op.atom_permutation.resize(w.atoms.size());std::iota(op.atom_permutation.begin(),op.atom_permutation.end(),0);}
                frame.operations.push_back(std::move(op));
            }
        }
        frame.used_group="finite sampling subgroup of "+geometry.point_group;
    }
    return frame_for(w,&frame);
}
// For a linear molecule O(2) characters are 2 cos(m theta), with an
// inversion parity in Dinfh and a reflection sign for Sigma. Compare EVERY
// sampled operation, not a visually inferred sigma/pi or one selected axis.
std::string match_linear(const Frame& f,const std::vector<double>& chars,std::size_t d){
    if(!f.valid||!f.linear||chars.size()!=f.ops.size()||(d!=1&&d!=2))return {};
    const char* symbols[]={"Sigma","Pi","Delta","Phi","Gamma"};
    std::string match;
    for(unsigned m=d==1?0:1;m<=(d==1?0:f.lmax);++m)for(int parity:{1,-1})for(int sign:{1,-1}){
        if(f.group=="Cinfv"&&parity<0)continue;if(m&&sign<0)continue;
        bool same=true;
        for(std::size_t i=0;i<f.ops.size();++i){const auto& op=f.ops[i].matrix;Vec transformed{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)transformed[a]+=op[3*a+b]*f.axis[b];
            const double direction=dot(f.axis,transformed);if(std::abs(std::abs(direction)-1)>metric_tolerance){same=false;break;}
            const bool reverse=direction<0;if(reverse&&f.group!="Dinfh"){same=false;break;}
            Mat r=op;if(reverse)for(auto& x:r)x=-x;
            const bool reflect=det(r)<0;
            double expected=m?(reflect?0:2*std::cos(m*std::acos(std::clamp((trace(r)-1)/2,-1.,1.)))):(reflect?sign:1);
            if(reverse)expected*=parity;
            if(!std::isfinite(chars[i])||std::abs(chars[i]-expected)>character_tolerance){same=false;break;}
        }
        if(same){std::string label=symbols[m];if(f.group=="Dinfh")label+=parity>0?"_g":"_u";if(!m)label+=sign>0?"+":"-";if(!match.empty()&&match!=label)return {};match=label;}
    }
    return match;
}
std::string match_simple(const Frame& f,const std::vector<double>& chars,std::size_t d){
    if(!f.valid||chars.size()!=f.ops.size()||d!=1)return {};
    const auto near=[&](std::size_t i,double target){return std::isfinite(chars[i])&&std::abs(chars[i]-target)<=character_tolerance;};
    if(!near(f.e,1))return {};
    if(f.group=="D2h")return orbital_characters::match_d2h(f.d2h_axes,f.ops,chars,d,character_tolerance);
    if(f.group=="C1")return "A";
    if(f.group=="Ci")return near(f.inv,1)?"Ag":near(f.inv,-1)?"Au":"";
    if(f.group=="Cs")return near(f.mirror0,1)?"A'":near(f.mirror0,-1)?"A''":"";
    if(f.group=="C2"||f.group=="C2h"){const auto base=near(f.rotation,1)?"A":near(f.rotation,-1)?"B":"";if(!*base)return {};if(f.group=="C2")return base;return near(f.inv,1)?std::string(base)+"g":near(f.inv,-1)?std::string(base)+"u":"";}
    if(f.group!="C2v")return {};
    if(near(f.rotation,1)&&near(f.mirror0,1)&&near(f.mirror1,1))return "A1";
    if(near(f.rotation,1)&&near(f.mirror0,-1)&&near(f.mirror1,-1))return "A2";
    if(!f.mirrors_named||!near(f.rotation,-1))return {};
    const auto other=f.yz==f.mirror0?f.mirror1:f.mirror0;
    if(near(f.yz,-1)&&near(other,1))return "B1";
    if(near(f.yz,1)&&near(other,-1))return "B2";
    return {};
}
std::string match_frame(const Frame& f,const std::vector<double>& chars,std::size_t d){
    auto label=match_simple(f,chars,d);if(!label.empty())return label;
    label=match_linear(f,chars,d);if(!label.empty())return label;
    // Axis-independent Cnv labels follow its real angular characters. Even-n
    // B1/B2 need a reflection-class convention and remain unresolved here.
    if(!f.valid||f.linear||f.rotation_order<3||f.group.front()!='C'||f.group.back()!='v'||chars.size()!=f.ops.size()||(d!=1&&d!=2))return {};
    for(unsigned m=d==1?0:1;m<=(d==1?0:(f.rotation_order-1)/2);++m)for(int sign:{1,-1}){if(m&&sign<0)continue;bool same=true;for(std::size_t i=0;i<chars.size();++i){const auto& op=f.ops[i].matrix;const bool reflect=det(op)<0;const double expected=m?(reflect?0:2*std::cos(m*std::acos(std::clamp((trace(op)-1)/2,-1.,1.)))):(reflect?sign:1);if(!std::isfinite(chars[i])||std::abs(chars[i]-expected)>character_tolerance)same=false;}if(same)return m?(f.rotation_order<=4?"E":"E"+std::to_string(m)):(sign>0?"A1":"A2");}
    return {};
}
double metric(const Wavefunction& w,const std::vector<double>& a,const std::vector<double>& b){
    const auto n=std::size_t(w.basis_count);double v=0;for(std::size_t i=0;i<n;++i){double row=0;for(std::size_t j=0;j<n;++j)row+=w.ao_overlap[i*n+j]*b[j];v+=a[i]*row;}return v;
}
const CanonicalActionCache& canonical_action(const Wavefunction& w,const Frame& f,std::size_t operation){
    const auto n=std::size_t(w.basis_count),m=w.orbitals.size();
    const auto metric_columns=[&](const std::vector<double>& packed){
        std::vector<double> result(n*m,0);
        // Contiguous column updates share each S entry and let the compiler
        // vectorise the independent columns of the complete immutable set.
        for(std::size_t row=0;row<n;++row)for(std::size_t j=0;j<n;++j){const auto s=w.ao_overlap[row*n+j];for(std::size_t col=0;col<m;++col)result[row*m+col]+=s*packed[j*m+col];}
        return result;
    };
    if(!f.canonical_action||f.canonical_action->wavefunction!=&w){
        f.canonical_action=std::make_shared<CanonicalActionCache>();auto& cache=*f.canonical_action;cache.wavefunction=&w;cache.packed.assign(n*m,0);
        for(std::size_t col=0;col<m;++col)if(w.orbitals[col].coefficients.size()==n)for(std::size_t row=0;row<n;++row)cache.packed[row*m+col]=w.orbitals[col].coefficients[row];
        cache.metric_packed=metric_columns(cache.packed);cache.transformed.resize(f.ops.size());cache.metric_transformed.resize(f.ops.size());
    }
    auto& cache=*f.canonical_action;
    if(cache.transformed[operation].empty()){
        cache.transformed[operation]=apply_orbital_symmetry_operation(w,f.ops[operation],cache.packed,m);
        if(cache.transformed[operation].size()==n*m)cache.metric_transformed[operation]=metric_columns(cache.transformed[operation]);
    }
    return cache;
}
// Complete S-orthonormal subspace character and leakage test. This examines
// the actual canonical coefficients, including core and all virtual rows.
std::vector<double> characters(const Wavefunction& w,const Frame& f,const std::vector<std::size_t>& indices){
    const auto n=std::size_t(w.basis_count),k=indices.size();if(!f.valid||!n||!k||w.ao_overlap.size()!=n*n)return {};
    std::vector<std::vector<double>> q;std::vector<double> packed(n*k);
    for(std::size_t col=0;col<k;++col){const auto index=indices[col];if(index>=w.orbitals.size()||w.orbitals[index].coefficients.size()!=n)return {};q.push_back(w.orbitals[index].coefficients);for(std::size_t row=0;row<n;++row)packed[row*k+col]=q.back()[row];}
    const auto m=w.orbitals.size();const auto& action=canonical_action(w,f,0);
    std::vector<std::vector<double>> sq(k,std::vector<double>(n,0));
    for(std::size_t col=0;col<k;++col)for(std::size_t row=0;row<n;++row)sq[col][row]=action.metric_packed[row*m+indices[col]];
    for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<k;++b){const double g=std::inner_product(q[a].begin(),q[a].end(),sq[b].begin(),0.);if(!std::isfinite(g)||std::abs(g-(a==b?1.:0.))>metric_tolerance)return {};}
    std::vector<double> result;
    for(std::size_t g=0;g<f.ops.size();++g){const auto& cache=canonical_action(w,f,g);const auto& tq=cache.transformed[g];const auto& stq=cache.metric_transformed[g];if(tq.size()!=n*m||stq.size()!=n*m)return {};double chi=0;
        for(std::size_t col=0;col<k;++col){std::vector<double> t(n),st(n);for(std::size_t row=0;row<n;++row){t[row]=tq[row*m+indices[col]];st[row]=stq[row*m+indices[col]];}const double norm=std::inner_product(t.begin(),t.end(),st.begin(),0.);if(!std::isfinite(norm)||std::abs(norm-1)>metric_tolerance)return {};auto residual=t,residual_metric=st;
            for(std::size_t a=0;a<k;++a){const double v=std::inner_product(sq[a].begin(),sq[a].end(),t.begin(),0.);if(a==col)chi+=v;for(std::size_t row=0;row<n;++row){residual[row]-=v*q[a][row];residual_metric[row]-=v*sq[a][row];}}
            const double leakage=std::inner_product(residual.begin(),residual.end(),residual_metric.begin(),0.);if(!std::isfinite(leakage)||leakage< -1e-10||leakage>character_tolerance*character_tolerance)return {};
        }result.push_back(chi);
    }return result;
}
// Energy supplies candidates only. Edges require measured operation coupling;
// the resulting blocks still have to pass the complete leakage test above.
std::vector<std::vector<std::size_t>> connected_blocks(const Wavefunction& w,const Frame& f,const std::vector<std::size_t>& rows){
    const auto n=std::size_t(w.basis_count),k=rows.size();std::vector<std::vector<std::size_t>> result;
    if(!k)return result;
    const auto singles=[&](){std::vector<std::vector<std::size_t>> v;for(auto i:rows)v.push_back({i});return v;};
    if(k==1)return singles();
    if(!f.valid||!n||w.ao_overlap.size()!=n*n)return singles();
    std::vector<double> packed(n*k),sq(n*k,0);
    for(std::size_t col=0;col<k;++col){if(rows[col]>=w.orbitals.size()||w.orbitals[rows[col]].coefficients.size()!=n)return singles();for(std::size_t a=0;a<n;++a)packed[a*k+col]=w.orbitals[rows[col]].coefficients[a];}
    const auto m=w.orbitals.size();const auto& action=canonical_action(w,f,0);
    for(std::size_t a=0;a<n;++a)for(std::size_t col=0;col<k;++col)sq[a*k+col]=action.metric_packed[a*m+rows[col]];
    for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<k;++b){double g=0;for(std::size_t row=0;row<n;++row)g+=packed[row*k+a]*sq[row*k+b];if(!std::isfinite(g)||std::abs(g-(a==b?1.:0.))>metric_tolerance)return singles();}
    std::vector<std::size_t> parent(k);std::iota(parent.begin(),parent.end(),0);auto root=[&](std::size_t a){while(parent[a]!=a)a=parent[a];return a;};
    for(std::size_t g=0;g<f.ops.size();++g){const auto& cache=canonical_action(w,f,g);const auto& tq=cache.transformed[g];if(tq.size()!=n*m)return singles();for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<a;++b){double ab=0,ba=0;for(std::size_t row=0;row<n;++row){ab+=sq[row*k+a]*tq[row*m+rows[b]];ba+=sq[row*k+b]*tq[row*m+rows[a]];}if(std::max(std::abs(ab),std::abs(ba))>character_tolerance)parent[root(b)]=root(a);}}
    std::map<std::size_t,std::vector<std::size_t>> groups;for(std::size_t a=0;a<k;++a)groups[root(a)].push_back(rows[a]);for(auto& [key,members]:groups)result.push_back(std::move(members));return result;
}
std::vector<std::vector<std::size_t>> canonical_blocks(const Wavefunction& w,const Frame& f){
    std::vector<std::vector<std::size_t>> result;std::vector<std::size_t> rows(w.orbitals.size());std::iota(rows.begin(),rows.end(),0);
    std::stable_sort(rows.begin(),rows.end(),[&](auto a,auto b){const auto& x=w.orbitals[a];const auto& y=w.orbitals[b];if(x.spin!=y.spin)return x.spin<y.spin;const bool xf=std::isfinite(x.energy_hartree),yf=std::isfinite(y.energy_hartree);if(xf!=yf)return xf;if(xf&&x.energy_hartree!=y.energy_hartree)return x.energy_hartree<y.energy_hartree;return a<b;});
    for(std::size_t lo=0;lo<rows.size();){std::size_t hi=lo+1;const auto& first=w.orbitals[rows[lo]];while(hi<rows.size()&&w.orbitals[rows[hi]].spin==first.spin&&std::isfinite(first.energy_hartree)&&std::abs(w.orbitals[rows[hi]].energy_hartree-first.energy_hartree)<=1e-5)++hi;
        auto blocks=connected_blocks(w,f,{rows.begin()+lo,rows.begin()+hi});result.insert(result.end(),blocks.begin(),blocks.end());lo=hi;
    }return result;
}
// Reuse the native finite/linear classifier on an independent grouping view.
// Synthetic energies are private block IDs, never displayed or used as data.
std::vector<std::string> native_labels(const Wavefunction& w,const Frame& f,const std::vector<std::vector<std::size_t>>& blocks){
    std::vector<std::string> labels(blocks.size());if(!f.valid)return labels;
    Wavefunction view;view.atoms=w.atoms;view.shells=w.shells;view.primitives=w.primitives;view.basis_count=w.basis_count;view.ao_overlap=w.ao_overlap;
    std::vector<std::size_t> offsets;std::vector<bool> valid;
    for(std::size_t b=0;b<blocks.size();++b){offsets.push_back(view.orbitals.size());bool ok=!blocks[b].empty();for(auto i:blocks[b])if(i>=w.orbitals.size()||w.orbitals[i].coefficients.size()!=w.basis_count)ok=false;valid.push_back(ok);if(!ok)continue;for(auto i:blocks[b]){MolecularOrbital mo;mo.coefficients=w.orbitals[i].coefficients;mo.energy_hartree=double(b);view.orbitals.push_back(std::move(mo));}}
    OrbitalSymmetryOptions options;options.character_tolerance=character_tolerance;options.minimum_subspace_retention=1-metric_tolerance;
    if(derive_orbital_symmetry(view,options).point_group!=f.group)return labels;
    for(std::size_t b=0;b<blocks.size();++b){if(!valid[b])continue;const auto& label=view.orbitals[offsets[b]].symmetry;bool same=!label.empty();for(std::size_t j=0;j<blocks[b].size();++j)if(view.orbitals[offsets[b]+j].symmetry!=label)same=false;if(same)labels[b]=label;}return labels;
}
struct Unit {
    std::vector<std::size_t> members;std::string irrep,scope,detail;
    double energy=0;bool energy_available=true;std::size_t copies=1;
    bool copies_resolved=true,counting_only=false;
    double lower=std::numeric_limits<double>::quiet_NaN(),upper=std::numeric_limits<double>::quiet_NaN();
};
void certify_salc_energy_bounds(Unit& unit,const Wavefunction& side,const Wavefunction& canonical,const NboSalcModel& model){
    if(unit.members.empty())return;const auto spin=model.orbitals[unit.members.front()].spin;
    if(spin==NboSpin::Total&&std::any_of(canonical.orbitals.begin(),canonical.orbitals.end(),[](const auto& mo){return mo.spin==Spin::Beta;}))return;
    const NboSalcEnergyEvidence* proof=nullptr;for(const auto& energy:model.energies)if(energy.spin==spin&&energy.available&&energy.electronic_symmetry_verified&&energy.status=="verified_same_operator")proof=&energy;
    if(!proof)return;
    const auto n=std::size_t(canonical.basis_count),k=unit.members.size();if(!n||canonical.ao_overlap.size()!=n*n)return;
    std::vector<std::size_t> columns;std::map<std::size_t,std::size_t> column_map;
    for(std::size_t i=0;i<canonical.orbitals.size();++i)if(canonical.orbitals[i].spin==(spin==NboSpin::Beta?Spin::Beta:Spin::Alpha)){if(canonical.orbitals[i].coefficients.size()!=n||!std::isfinite(canonical.orbitals[i].energy_hartree))return;column_map[i]=columns.size();columns.push_back(i);}
    if(columns.empty()||columns.size()!=proof->canonical_columns_checked)return;
    const auto m=columns.size();std::vector<double> projections(k*m,0);std::vector<bool> seen(k*m,false);
    for(const auto& link:model.links){auto member=std::find(unit.members.begin(),unit.members.end(),link.side_index);if(member==unit.members.end())continue;const auto col=column_map.find(link.canonical_index);if(col==column_map.end())continue;const auto at=std::size_t(member-unit.members.begin())*m+col->second;if(seen[at]||!std::isfinite(link.coefficient))return;seen[at]=true;projections[at]=link.coefficient;}
    if(std::find(seen.begin(),seen.end(),false)!=seen.end())return;
    // The producer-validated canonical eigenbasis must actually span each
    // retained SALC, in the same S metric and spin. Missing virtual columns
    // cannot silently turn a projected operator into a certified spectrum.
    for(std::size_t a=0;a<k;++a){const auto row=unit.members[a];if(row>=side.orbitals.size()||side.orbitals[row].coefficients.size()!=n)return;auto residual=side.orbitals[row].coefficients;for(std::size_t j=0;j<m;++j)for(std::size_t r=0;r<n;++r)residual[r]-=projections[a*m+j]*canonical.orbitals[columns[j]].coefficients[r];const double loss=metric(canonical,residual,residual);if(!std::isfinite(loss)||loss< -1e-10||loss>character_tolerance*character_tolerance)return;}
    std::vector<double> f(k*k,0);double gram_error=0;
    for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<k;++b){double gram=0;for(std::size_t j=0;j<m;++j){const double product=projections[a*m+j]*projections[b*m+j];gram+=product;f[a*k+b]+=product*canonical.orbitals[columns[j]].energy_hartree;}gram_error=std::max(gram_error,std::abs(gram-(a==b?1.:0.)));}
    if(!std::isfinite(gram_error)||gram_error>metric_tolerance||double(k)*gram_error>=.1)return;
    double lower=std::numeric_limits<double>::infinity(),upper=-lower;
    for(std::size_t a=0;a<k;++a){const auto energy=model.orbitals[unit.members[a]].energy_hartree;if(!energy||!std::isfinite(*energy)||std::abs(f[a*k+a]-*energy)>energy_order_tolerance)return;double radius=0;for(std::size_t b=0;b<k;++b)if(a!=b)radius+=std::abs(f[a*k+b]);lower=std::min(lower,f[a*k+a]-radius);upper=std::max(upper,f[a*k+a]+radius);}
    // Gershgorin bounds concern the complete projected Hermitian operator,
    // not the diagonal expectations or their mean. Enlarge for metric and
    // validated producer residuals; never rotate the displayed SALC basis.
    const double defect=double(k)*gram_error;const double padding=std::max(std::abs(lower),std::abs(upper))*defect/(1-defect)+energy_order_tolerance+std::abs(proof->eigenvalue_error_hartree);
    if(!std::isfinite(lower)||!std::isfinite(upper)||!std::isfinite(padding))return;unit.lower=lower-padding;unit.upper=upper+padding;
    std::ostringstream detail;detail<<"; same-spin complete-canonical projected Fock spectral bounds=["<<unit.lower<<','<<unit.upper<<"] Ha (Gershgorin with metric/error margin)";unit.detail+=detail.str();
}
std::string orbital_label(const std::string& irrep){
    auto label=format_symmetry_unicode(normalized(irrep));
    // These labels name one-electron orbitals. Preserve the machine irrep
    // identifier while using the spectroscopic lowercase orbital convention.
    const std::pair<const char*,const char*> symbols[]={{"Σ","σ"},{"Π","π"},{"Δ","δ"},{"Φ","φ"},{"Γ","γ"}};
    for(const auto& [upper,lower]:symbols)if(label.starts_with(upper)){label.replace(0,std::string(upper).size(),lower);break;}
    return label;
}
void assign_ordinals(std::vector<Unit> units,std::vector<NboAomoName>& names){
    for(auto& u:units)if(u.copies_resolved&&u.energy_available&&std::isfinite(u.energy)&&(!u.irrep.empty()||u.scope.starts_with("canonical ")))u.lower=u.upper=u.energy;
    for(std::size_t target=0;target<units.size();++target){const auto& u=units[target];if(u.irrep.empty()||u.counting_only)continue;
        std::size_t ordinal=1;bool resolved=u.copies_resolved&&u.energy_available&&std::isfinite(u.lower)&&std::isfinite(u.upper);
        for(std::size_t other=0;resolved&&other<units.size();++other){if(other==target)continue;const auto& v=units[other];if(v.scope!=u.scope||(!v.irrep.empty()&&v.irrep!=u.irrep))continue;
            if(std::isfinite(v.lower)&&v.lower>u.upper+energy_order_tolerance)continue;
            if(!v.irrep.empty()&&std::isfinite(v.upper)&&v.upper<u.lower-energy_order_tolerance){ordinal+=v.copies;continue;}
            if(!v.irrep.empty()&&v.copies_resolved&&v.energy_available&&std::isfinite(v.energy)){
                // Two independent occurrences at unresolved equal energy
                // have no physical copy order. Source row is not evidence.
                if(std::abs(v.energy-u.energy)<=energy_order_tolerance)resolved=false;
                else if(v.energy<u.energy)ordinal+=v.copies;
            }else resolved=false;
        }
        for(auto index:u.members){auto& name=names[index];name.irrep=u.irrep;name.verified=true;name.ordinal=resolved?ordinal:0;name.representation_multiplicity=u.copies;name.label=(name.ordinal?std::to_string(name.ordinal):std::string{})+orbital_label(u.irrep);name.detail=u.detail+"; ordinal scope="+u.scope+"; complete-set counts with certified subspace energy ordering"+(name.ordinal?"":"; irrep known; missing energy evidence or overlapping/unresolved earlier representation counts");}
    }
}
} // namespace

NboAomoNames build_nbo_aomo_names(const Wavefunction& w,const NboIntegration& data,const NboSalcModel* salc){
    NboAomoNames out;out.canonical.resize(w.orbitals.size());if(salc)out.salc.resize(salc->orbitals.size());
    const bool open=std::any_of(w.orbitals.begin(),w.orbitals.end(),[](const auto& mo){return mo.spin==Spin::Beta;});
    const bool associated=!salc || ((salc->dataset_id.empty()||data.id.empty()||salc->dataset_id==data.id)&&
        (salc->canonical_fingerprint.empty()||data.canonical_fingerprint.empty()||salc->canonical_fingerprint==data.canonical_fingerprint));
    const auto frame=salc?frame_for(w,associated?salc:nullptr):canonical_frame(w);std::vector<Unit> units;std::vector<bool> taken(w.orbitals.size(),false);
    const bool simple=frame.group=="C1"||frame.group=="Ci"||frame.group=="Cs"||frame.group=="C2"||frame.group=="C2h"||frame.group=="C2v"||frame.group=="D2h";
    struct NamedCharacters {std::string label;std::size_t dimension;std::vector<double> values;};std::vector<NamedCharacters> known;
    // Complete small Abelian tables are already defined by match_simple.
    // This also allows a mixed closed block to be counted when no separate
    // pure orbital happens to supply a named exemplar of one component.
    if(simple&&frame.valid&&frame.ops.size()<=8)for(std::size_t mask=0;mask<(std::size_t(1)<<frame.ops.size());++mask){std::vector<double> values(frame.ops.size());for(std::size_t g=0;g<values.size();++g)values[g]=(mask&(std::size_t(1)<<g))?-1:1;const auto label=match_simple(frame,values,1);bool representation=!label.empty();for(std::size_t a=0;representation&&a<values.size();++a)for(std::size_t b=0;b<values.size();++b){const auto product=multiply(frame.ops[a].matrix,frame.ops[b].matrix);for(std::size_t c=0;c<values.size();++c)if(matrix_error(product,frame.ops[c].matrix)<metric_tolerance&&values[c]!=values[a]*values[b])representation=false;}if(representation&&std::none_of(known.begin(),known.end(),[&](const auto& x){return x.label==label;}))known.push_back({label,1,std::move(values)});}
    const auto blocks=canonical_blocks(w,frame);
    const auto derived=(simple||frame.linear)?std::vector<std::string>(blocks.size()):native_labels(w,frame,blocks);
    for(std::size_t b=0;b<blocks.size();++b){const auto& members=blocks[b];const auto values=characters(w,frame,members);if(values.empty())continue;
        auto label=match_frame(frame,values,members.size());if(label.empty())label=derived[b];
        double norm=0;for(auto value:values)norm+=value*value/double(values.size());
        if(label.empty()||std::abs(norm-1)>1e-3)continue;
        Unit u;u.members=members;u.irrep=label;u.scope=w.orbitals[members.front()].spin==Spin::Beta?"canonical beta":"canonical alpha";u.detail="Calculated from complete canonical AO coefficients: S-orthonormal invariant irrep under every validated operation; "+frame.detail;
        for(auto i:members){taken[i]=true;u.energy+=w.orbitals[i].energy_hartree/double(members.size());u.energy_available=u.energy_available&&std::isfinite(w.orbitals[i].energy_hartree);}
        if(std::none_of(known.begin(),known.end(),[&](const auto& k){return k.label==label;}))known.push_back({label,members.size(),values});units.push_back(std::move(u));
    }
    for(std::size_t i=0;i<w.orbitals.size();++i){if(taken[i])continue;const auto& mo=w.orbitals[i];Unit u;u.members={i};u.scope=mo.spin==Spin::Beta?"canonical beta":"canonical alpha";u.energy=mo.energy_hartree;u.energy_available=std::isfinite(u.energy);
        auto evidence=molecular_orbital_symmetry(w,i);std::string group=evidence.point_group;
        if(evidence.origin==OrbitalSymmetryOrigin::Producer){group=evidence.producer_abelian_group.empty()?evidence.producer_detected_group:evidence.producer_abelian_group;}
        const auto existing=finite_irrep(group,evidence.label);const auto dim=dimension(group,existing);
        bool solid=evidence.origin==OrbitalSymmetryOrigin::Producer&&!evidence.source_path.empty()&&!existing.empty();
        if(evidence.origin==OrbitalSymmetryOrigin::MolecularOperations&&evidence.molecular_assignment){const auto& a=*evidence.molecular_assignment;solid=!existing.empty()&&std::isfinite(a.subspace_retention)&&a.subspace_retention>=.985&&a.orbital_indices.size()==dim;}
        // A producer's Abelian subgroup is not the full-group naming scope.
        // Keep its literal source record intact, but do not mix its labels
        // with calculated full-group labels and their occurrence counters.
        if(!frame.group.empty()&&group!=frame.group)solid=false;
        // Native derived assignments use their own retention threshold. They
        // must pass this naming layer's unchanged complete-span leakage gate
        // before being excluded from unknown-block representation counting.
        if(solid&&evidence.origin==OrbitalSymmetryOrigin::MolecularOperations&&dim==1&&characters(w,frame,{i}).empty())solid=false;
        if(solid&&dim>1){std::vector<std::size_t> members;
            if(evidence.molecular_assignment)members=evidence.molecular_assignment->orbital_indices;
            else for(std::size_t j=0;j<w.orbitals.size();++j)if(w.orbitals[j].spin==mo.spin&&std::abs(w.orbitals[j].energy_hartree-mo.energy_hartree)<=1e-5&&normalized(w.orbitals[j].symmetry)==normalized(mo.symmetry))members.push_back(j);
            if(members.size()!=dim)solid=false;
            for(auto j:members)if(j>=w.orbitals.size()||taken[j]||w.orbitals[j].spin!=mo.spin||normalized(w.orbitals[j].symmetry)!=normalized(mo.symmetry))solid=false;
            if(solid&&characters(w,frame,members).empty())solid=false;
            if(solid){u.members=members;u.energy=0;for(auto j:members)u.energy+=w.orbitals[j].energy_hartree/double(members.size());}
        }
        // Producer B1/B2 can use another axis convention. In the simple groups
        // verify agreement in this recorded frame before mixing its name with
        // labels derived for the other canonical rows or the SALCs.
        std::vector<double> simple_values;
        if(simple){simple_values=characters(w,frame,{i});if(solid&&dim==1&&match_simple(frame,simple_values,1)!=existing)solid=false;}
        if(solid){u.irrep=existing;u.detail="Existing "+std::string(orbital_symmetry_origin_name(evidence.origin))+" whole-orbital irrep; point group="+group;
            // A measured character signature transfers a known label to SALCs
            // in exactly the same input frame, never by dimension alone.
            if(!simple&&group==frame.group&&std::none_of(known.begin(),known.end(),[&](const auto& k){return k.label==existing;})){auto signature=characters(w,frame,u.members);if(!signature.empty())known.push_back({existing,u.members.size(),std::move(signature)});}
        }else{
            u.irrep=match_simple(frame,simple_values,1);
            u.detail=u.irrep.empty()?"No verified whole-orbital irrep (unsupported/ambiguous frame, mixed span, missing metric or operation action)":"Actual complete canonical AO coefficients pass every S-metric operation eigenfunction/leakage test; "+frame.detail;
        }
        for(auto j:u.members){taken[j]=true;out.canonical[j].label="?";out.canonical[j].detail=u.detail+"; canonical source row="+std::to_string(j+1);}
        units.push_back(std::move(u));
    }
    // A low-energy spectral cluster may have slightly mixed individual
    // eigenvectors while its COMPLETE span is accurately invariant. Count
    // its exact irrep content without labelling or rotating those members.
    std::vector<bool> replaced(units.size(),false);std::vector<Unit> count_units;
    for(const auto spin:{Spin::Alpha,Spin::Beta}){std::vector<std::size_t> unknown;for(const auto& u:units)if(u.irrep.empty())for(auto i:u.members)if(w.orbitals[i].spin==spin)unknown.push_back(i);
        for(const auto& members:connected_blocks(w,frame,unknown)){const auto values=characters(w,frame,members);if(values.empty())continue;std::vector<double> reconstructed(values.size(),0);std::vector<std::pair<const NamedCharacters*,std::size_t>> decomposition;std::size_t dimension_sum=0;bool valid=true;
            for(const auto& k:known){if(k.values.size()!=values.size())continue;double inner=0;for(std::size_t g=0;g<values.size();++g)inner+=values[g]*k.values[g]/double(values.size());if(!std::isfinite(inner)){valid=false;break;}const auto multiplicity=std::llround(inner);if(multiplicity<0||std::abs(inner-double(multiplicity))>character_tolerance){valid=false;break;}if(!multiplicity)continue;decomposition.push_back({&k,std::size_t(multiplicity)});dimension_sum+=k.dimension*std::size_t(multiplicity);for(std::size_t g=0;g<values.size();++g)reconstructed[g]+=double(multiplicity)*k.values[g];}
            if(!valid||dimension_sum!=members.size())continue;for(std::size_t g=0;g<values.size();++g)if(std::abs(values[g]-reconstructed[g])>character_tolerance)valid=false;if(!valid)continue;
            double low=std::numeric_limits<double>::infinity(),high=-low;for(auto i:members){const double energy=w.orbitals[i].energy_hartree;if(!std::isfinite(energy)){valid=false;break;}low=std::min(low,energy);high=std::max(high,energy);}if(!valid)continue;
            for(std::size_t i=0;i<units.size();++i)if(units[i].irrep.empty()&&std::all_of(units[i].members.begin(),units[i].members.end(),[&](auto member){return std::find(members.begin(),members.end(),member)!=members.end();}))replaced[i]=true;
            for(const auto& [k,copies]:decomposition){Unit count;count.members=members;count.irrep=k->label;count.scope=spin==Spin::Beta?"canonical beta":"canonical alpha";count.copies=copies;count.copies_resolved=false;count.counting_only=true;count.lower=low;count.upper=high;count_units.push_back(std::move(count));}
            for(auto i:members)out.canonical[i].detail+="; complete invariant spectral block has verified integer irrep counts and eigenvalue bounds; individual mixed member remains unclassified";
        }
    }
    std::vector<Unit> counted;for(std::size_t i=0;i<units.size();++i)if(!replaced[i])counted.push_back(std::move(units[i]));counted.insert(counted.end(),count_units.begin(),count_units.end());
    assign_ordinals(std::move(counted),out.canonical);
    for(auto& name:out.canonical) if(name.verified) name.point_group=frame.group;
    for(std::size_t i=0;i<out.canonical.size();++i){auto& name=out.canonical[i];if(name.verified&&!name.ordinal)name.label=canonical_mo_display_label(w,i,&name);else if(open)name.label+=w.orbitals[i].spin==Spin::Beta?" [beta]":" [alpha]";}
    if(!salc)return out;
    Wavefunction side;side.atoms=w.atoms;side.shells=w.shells;side.primitives=w.primitives;side.basis_count=w.basis_count;side.ao_overlap=w.ao_overlap;side.orbitals.resize(salc->orbitals.size());
    for(std::size_t i=0;i<salc->orbitals.size();++i){const auto& orbital=salc->orbitals[i];auto& mo=side.orbitals[i];mo.energy_hartree=orbital.energy_hartree.value_or(std::numeric_limits<double>::quiet_NaN());mo.spin=orbital.spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;
        bool valid=associated&&!orbital.terms.empty();mo.coefficients.assign(w.basis_count,0);
        for(const auto& term:orbital.terms){const auto* descriptor=nbo_orbital(data,term.orbital);if(!descriptor||descriptor->coefficients.size()!=w.basis_count||!std::isfinite(term.coefficient)){valid=false;break;}for(std::size_t a=0;a<w.basis_count;++a)mo.coefficients[a]+=term.coefficient*descriptor->coefficients[a];}
        if(!valid)mo.coefficients.clear();
    }
    std::vector<std::vector<std::size_t>> side_blocks;std::vector<std::size_t> side_subspaces;
    for(std::size_t s=0;s<salc->subspaces.size();++s){const auto& sub=salc->subspaces[s];bool have_coefficients=!sub.orbital_indices.empty();for(auto i:sub.orbital_indices)if(i>=side.orbitals.size()||side.orbitals[i].coefficients.size()!=w.basis_count)have_coefficients=false;
        const auto pieces=have_coefficients?connected_blocks(side,frame,sub.orbital_indices):std::vector<std::vector<std::size_t>>{sub.orbital_indices};for(const auto& piece:pieces){side_blocks.push_back(piece);side_subspaces.push_back(s);}}
    const auto side_derived=(simple||frame.linear)?std::vector<std::string>(side_blocks.size()):native_labels(side,frame,side_blocks);
    std::vector<bool> used(salc->orbitals.size(),false);units.clear();
    for(std::size_t b=0;b<side_blocks.size();++b){const auto& sub=salc->subspaces[side_subspaces[b]];Unit u;u.scope=sub.fragment_id+":"+nbo_spin_name(sub.spin);u.detail="Fixed fragment subspace; "+frame.detail;
        bool indices_ok=!side_blocks[b].empty();for(auto i:side_blocks[b])if(i>=out.salc.size()||used[i]||salc->orbitals[i].fragment_id!=sub.fragment_id||salc->orbitals[i].spin!=sub.spin)indices_ok=false;if(!indices_ok)continue;u.members=side_blocks[b];
        const bool verified=frame.valid&&sub.symmetry_verified&&sub.dimension==sub.orbital_indices.size()&&std::isfinite(sub.closure_error)&&sub.closure_error<=character_tolerance&&std::isfinite(sub.orthogonality_error)&&sub.orthogonality_error<=metric_tolerance;
        if(verified){auto values=characters(side,frame,u.members);bool actual=!values.empty();
            // Stored full-subspace characters remain useful even when the
            // attachment omitted AO columns. Never reuse them for a subset.
            if(!actual&&u.members==sub.orbital_indices)values=sub.characters;
            double norm=0;for(double chi:values)norm+=chi*chi/double(values.size());
            const auto copies=std::isfinite(norm)?std::size_t(std::llround(std::sqrt(std::max(0.,norm)))):0;
            if(copies&&std::abs(norm-double(copies*copies))<1e-3&&u.members.size()%copies==0){const auto dim=u.members.size()/copies;for(auto& value:values)value/=double(copies);
                u.irrep=match_frame(frame,values,dim);if(u.irrep.empty()&&copies==1&&actual)u.irrep=side_derived[b];
                if(u.irrep.empty())for(const auto& k:known){if(k.dimension!=dim||k.values.size()!=values.size())continue;bool same=true;for(std::size_t a=0;a<values.size();++a)if(!std::isfinite(values[a])||std::abs(k.values[a]-values[a])>character_tolerance)same=false;if(same){if(!u.irrep.empty()&&u.irrep!=k.label){u.irrep.clear();break;}u.irrep=k.label;}}
                u.copies=copies;u.copies_resolved=copies==1;
            }
        }
        if(u.irrep.empty())u.detail+="; unclassified: reducible/mixed span, unsupported action, ambiguous axes, or no measured named signature";
        else u.detail+="; full-subspace characters verified; irrep occurrence count="+std::to_string(u.copies)+(u.copies_resolved?"; one ordinal shared by true partners":"; repeated irrep known but individual copy numbering unresolved");
        for(auto i:u.members){used[i]=true;const auto& o=salc->orbitals[i];if(o.energy_hartree&&std::isfinite(*o.energy_hartree))u.energy+=*o.energy_hartree/double(u.members.size());else u.energy_available=false;out.salc[i].label="?";out.salc[i].detail=u.detail;
            if(verified && !u.irrep.empty() && u.copies_resolved) {
                out.salc[i].partner_block_id="salc-irrep-block:"+std::to_string(b);
                out.salc[i].partner_block_size=u.members.size();
            }
        }
        if(verified&&(u.irrep.empty()||!u.copies_resolved))certify_salc_energy_bounds(u,side,w,*salc);
        units.push_back(std::move(u));
    }
    for(std::size_t i=0;i<out.salc.size();++i)if(!used[i]){const auto& o=salc->orbitals[i];Unit u;u.members={i};u.scope=o.fragment_id+":"+nbo_spin_name(o.spin);u.energy=o.energy_hartree.value_or(0);u.energy_available=o.energy_hartree&&std::isfinite(*o.energy_hartree);units.push_back(u);out.salc[i].label="?";out.salc[i].detail="No validated containing SALC subspace";}
    assign_ordinals(std::move(units),out.salc);
    for(std::size_t i=0;i<out.salc.size();++i)if(out.salc[i].verified&&!out.salc[i].ordinal)out.salc[i].label+=" [SALC "+std::to_string(i+1)+"]";
    for(std::size_t i=0;i<out.salc.size();++i)if(open && !salc->orbitals[i].spatial_spin)out.salc[i].label+=salc->orbitals[i].spin==NboSpin::Beta?" [beta]":salc->orbitals[i].spin==NboSpin::Alpha?" [alpha]":" [total]";
    for(auto& name:out.salc) if(name.verified) name.point_group=frame.group;
    (void)data; // Identity is immutable and belongs to the caller's attachment.
    return out;
}

std::shared_ptr<const NboAomoNames> canonical_mo_names(const Wavefunction& w){
    struct Cache {
        const Wavefunction* wavefunction=nullptr;
        const MolecularOrbital* orbitals=nullptr;
        const Atom* atoms=nullptr;
        const DerivedOrbitalSymmetryAssignment* assignments=nullptr;
        std::size_t orbital_count=0,atom_count=0,assignment_count=0,revision=0;
        std::string group,enrichment;
        std::shared_ptr<const NboAomoNames> names;
    };
    static thread_local Cache cache;
    if(!cache.names||cache.revision!=canonical_name_revision||cache.wavefunction!=&w||cache.orbitals!=w.orbitals.data()||
       cache.atoms!=w.atoms.data()||cache.assignments!=w.derived_orbital_symmetry_assignments.data()||
       cache.orbital_count!=w.orbitals.size()||cache.atom_count!=w.atoms.size()||
       cache.assignment_count!=w.derived_orbital_symmetry_assignments.size()||
       cache.group!=w.point_group_detected||cache.enrichment!=w.enrichment_source){
        cache.wavefunction=&w;cache.orbitals=w.orbitals.data();cache.atoms=w.atoms.data();
        cache.assignments=w.derived_orbital_symmetry_assignments.data();
        cache.orbital_count=w.orbitals.size();cache.atom_count=w.atoms.size();
        cache.assignment_count=w.derived_orbital_symmetry_assignments.size();
        cache.group=w.point_group_detected;cache.enrichment=w.enrichment_source;
        cache.revision=canonical_name_revision;
        OpenProfile profile;
        cache.names=std::make_shared<const NboAomoNames>(build_nbo_aomo_names(w,NboIntegration{},nullptr));
        profile.stage("standalone-orbital-names");
    }
    return cache.names;
}

void invalidate_canonical_mo_names_cache(){++canonical_name_revision;}

std::string canonical_mo_source_label(const Wavefunction& w,std::size_t index){
    if(index>=w.orbitals.size())return "MO ?";
    const auto& mo=w.orbitals[index];
    const bool source=mo.source_orbital_index!=std::numeric_limits<std::size_t>::max();
    std::string label=std::string(source?"MO ":"MO [list] ")+std::to_string((source?mo.source_orbital_index:index)+1);
    const bool explicit_spin=w.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin||
        std::any_of(w.orbitals.begin(),w.orbitals.end(),[](const auto& o){return o.spin==Spin::Beta;});
    if(explicit_spin)label+=mo.spin==Spin::Beta?" [beta]":" [alpha]";
    return label;
}

std::string canonical_mo_display_label(const Wavefunction& w,std::size_t index,const NboAomoName* name){
    if(index>=w.orbitals.size())return "MO ?";
    const auto standalone=name?std::shared_ptr<const NboAomoNames>{}:canonical_mo_names(w);
    if(!name&&standalone&&index<standalone->canonical.size())name=&standalone->canonical[index];
    if(name&&name->verified&&name->ordinal&&!name->irrep.empty()){
        auto label=std::to_string(name->ordinal)+orbital_label(name->irrep);
        const bool explicit_spin=w.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin||
            std::any_of(w.orbitals.begin(),w.orbitals.end(),[](const auto& o){return o.spin==Spin::Beta;});
        if(explicit_spin)
            label+=w.orbitals[index].spin==Spin::Beta?" [beta]":" [alpha]";
        return label;
    }
    if(name&&name->verified&&!name->irrep.empty()){
        auto source=canonical_mo_source_label(w,index);std::string suffix;
        for(const std::string spin:{" [alpha]"," [beta]"})if(source.ends_with(spin)){source.resize(source.size()-spin.size());suffix=spin;break;}
        return orbital_label(name->irrep)+" ["+source+"]"+suffix;
    }
    return canonical_mo_source_label(w,index);
}

std::string canonical_mo_current_irrep(const Wavefunction& w,std::size_t index,const NboAomoName* name){
    if(index>=w.orbitals.size())return "?";
    const auto standalone=name?std::shared_ptr<const NboAomoNames>{}:canonical_mo_names(w);
    if(!name&&standalone&&index<standalone->canonical.size())name=&standalone->canonical[index];
    return name?orbital_irrep_display_label(*name):"?";
}

std::string orbital_irrep_display_label(const NboAomoName& name){
    return name.verified&&!name.irrep.empty()?orbital_label(name.irrep):"?";
}
} // namespace cov::ui
