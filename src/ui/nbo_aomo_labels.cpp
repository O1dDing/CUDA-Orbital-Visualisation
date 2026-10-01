#include "cov/nbo_aomo_labels.hpp"
#include "cov/d2h_orbital_characters.hpp"
#include "cov/open_profile.hpp"
#include "cov/orbital_symmetry_scope.hpp"
#include "cov/orbital_symmetry.hpp"
#include "cov/molecular_point_group_frame.hpp"
#include "cov/mo_diagram.hpp"
#include <Eigen/Core>

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
#include <iomanip>

namespace cov::ui {
namespace {
constexpr double character_tolerance=2e-4;
constexpr double metric_tolerance=2e-5;
constexpr double energy_order_tolerance=2e-5;
thread_local std::size_t canonical_name_revision=0;
using Vec=std::array<double,3>;
using Mat=std::array<double,9>;
using Dense=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
double dot(const Vec& a,const Vec& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double det(const Mat& a){return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6]);}
double trace(const Mat& a){return a[0]+a[4]+a[8];}
double matrix_error(const Mat& a,const Mat& b){double e=0;for(std::size_t i=0;i<9;++i)e=std::max(e,std::abs(a[i]-b[i]));return e;}
Mat multiply(const Mat& a,const Mat& b){Mat c{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)c[3*i+j]+=a[3*i+k]*b[3*k+j];return c;}
const Mat identity{1,0,0,0,1,0,0,0,1};
const Mat inversion{-1,0,0,0,-1,0,0,0,-1};
std::string normalized(std::string s){std::string o;for(unsigned char c:s)if(!std::isspace(c)&&c!='_')o+=char(std::tolower(c));return o;}
struct CanonicalActionCache {
    const Wavefunction* wavefunction=nullptr;
    std::vector<double> packed,metric_packed;
    std::vector<std::vector<double>> transformed,metric_transformed;
};
struct Frame {
    std::string group, detail;
    std::vector<SymmetryOperation> ops;
    bool valid=false, linear=false;
    Vec axis{};
    unsigned lmax=0;
    PointGroupIrrepTable table;
    mutable std::shared_ptr<CanonicalActionCache> canonical_action;
};
Frame frame_for(const Wavefunction& w,const NboSalcModel* model){
    Frame f;
    if(!model||!model->group_verified)return f;
    f.group=model->point_group;f.linear=f.group=="Dinfh"||f.group=="Cinfv";
    if(model->used_group!=f.group&&(!f.linear||model->used_group!="finite sampling subgroup of "+f.group))return f;
    f.ops=model->operations;const auto n=f.ops.size();
    const auto expected_order=finite_point_group_order(f.group);
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
    bool have_e=false;
    for(std::size_t i=0;i<n;++i){const auto& op=f.ops[i];
        if(op.atom_permutation.size()!=w.atoms.size())return f;
        std::set<std::size_t> seen;for(std::size_t a=0;a<w.atoms.size();++a){const auto b=op.atom_permutation[a];if(b>=w.atoms.size()||w.atoms[a].atomic_number!=w.atoms[b].atomic_number||!seen.insert(b).second)return f;}
        Mat trans{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)trans[3*a+b]=op.matrix[3*b+a];
        if(matrix_error(multiply(trans,op.matrix),identity)>metric_tolerance)return f;
        if(matrix_error(op.matrix,identity)<metric_tolerance)have_e=true;
    }
    if(!have_e)return f;
    // Recheck matrix closure, not just the advertised point-group string.
    for(const auto& a:f.ops)for(const auto& b:f.ops){const auto ab=multiply(a.matrix,b.matrix);if(std::none_of(f.ops.begin(),f.ops.end(),[&](const auto& c){return matrix_error(ab,c.matrix)<metric_tolerance;}))return f;}
    f.valid=true;
    if(!f.linear&&f.valid){
        auto group=analyse_molecular_symmetry(w);group.point_group=f.group;group.operations=f.ops;
        f.table=molecular_point_group_irreps(w,group);
        f.valid=f.table.valid;
        f.detail=f.table.axis_detail+(f.table.valid?"":"; "+f.table.reason);
    }
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
        const auto complete=complete_molecular_point_group(w,geometry);
        frame.operations=complete.operations;frame.group_verified=!complete.operations.empty();
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
double metric(const Wavefunction& w,const std::vector<double>& a,const std::vector<double>& b){
    const auto n=std::size_t(w.basis_count);double v=0;for(std::size_t i=0;i<n;++i){double row=0;for(std::size_t j=0;j<n;++j)row+=w.ao_overlap[i*n+j]*b[j];v+=a[i]*row;}return v;
}
const CanonicalActionCache& canonical_action(const Wavefunction& w,const Frame& f,std::size_t operation){
    const auto n=std::size_t(w.basis_count),m=w.orbitals.size();
    const auto metric_columns=[&](const std::vector<double>& packed){
        std::vector<double> result(n*m,0);
        Eigen::Map<Dense>(result.data(),n,m).noalias()=
            Eigen::Map<const Dense>(w.ao_overlap.data(),n,n)*Eigen::Map<const Dense>(packed.data(),n,m);
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
struct CharacterFailure {std::string status,detail;};
std::vector<double> characters(const Wavefunction& w,const Frame& f,const std::vector<std::size_t>& indices,CharacterFailure* failure=nullptr){
    const auto fail=[&](const char* status,const char* quantity,double value,double limit){
        if(failure){failure->status=status;std::ostringstream text;text<<std::setprecision(12)<<quantity<<'='<<value<<"; limit="<<limit;failure->detail=text.str();}
        return std::vector<double>{};
    };
    const auto n=std::size_t(w.basis_count),k=indices.size(),m=w.orbitals.size();
    if(!f.valid||!n||!k||w.ao_overlap.size()!=n*n)return fail("symmetry_evidence_unavailable","AO metric dimension",double(w.ao_overlap.size()),double(n*n));
    const auto& action=canonical_action(w,f,0);Dense q(n,k),sq(n,k);
    for(std::size_t col=0;col<k;++col){const auto i=indices[col];if(i>=m||w.orbitals[i].coefficients.size()!=n)return {};
        for(std::size_t r=0;r<n;++r){q(r,col)=w.orbitals[i].coefficients[r];sq(r,col)=action.metric_packed[r*m+i];}}
    const Dense gram=q.transpose()*sq;
    if(!gram.allFinite()||(gram-Dense::Identity(k,k)).cwiseAbs().maxCoeff()>metric_tolerance)
        return fail("source_metric_not_orthonormal","source Gram error",(gram-Dense::Identity(k,k)).cwiseAbs().maxCoeff(),metric_tolerance);
    std::vector<double> result;
    for(std::size_t g=0;g<f.ops.size();++g){const auto& cache=canonical_action(w,f,g);
        if(cache.transformed[g].size()!=n*m||cache.metric_transformed[g].size()!=n*m)return {};
        Dense tq(n,k),stq(n,k);for(std::size_t c=0;c<k;++c)for(std::size_t r=0;r<n;++r){tq(r,c)=cache.transformed[g][r*m+indices[c]];stq(r,c)=cache.metric_transformed[g][r*m+indices[c]];}
        const Dense d=sq.transpose()*tq;
        const Dense residual=tq-q*d,residual_metric=stq-sq*d;
        for(std::size_t c=0;c<k;++c){const double norm=tq.col(c).dot(stq.col(c)),loss=residual.col(c).dot(residual_metric.col(c));
            if(!std::isfinite(norm)||std::abs(norm-1)>metric_tolerance)
                return fail("symmetry_action_not_isometric","transformed norm error",std::abs(norm-1),metric_tolerance);
            if(!std::isfinite(loss)||loss< -1e-10||loss>character_tolerance*character_tolerance)
                return fail("subspace_not_closed","S-metric leakage squared",loss,character_tolerance*character_tolerance);}
        result.push_back(d.trace());
    }return result;
}
// Residual of the real central idempotent on each unchanged source column.
// Both canonical MOs and reconstructed source SALCs use this same S-metric
// test after their complete containing span passes the closure/Gram gates.
std::vector<double> projection_residuals(const Wavefunction& w,const Frame& f,
    const std::vector<std::size_t>& members,std::size_t dimension,double norm,
    const std::vector<double>& row_characters){
    const auto n=std::size_t(w.basis_count),m=w.orbitals.size(),k=members.size();
    if(!f.valid||!n||!k||w.ao_overlap.size()!=n*n||row_characters.size()!=f.ops.size()||norm<=0)return {};
    Dense projected=Dense::Zero(n,k);
    for(std::size_t g=0;g<f.ops.size();++g){const auto& cache=canonical_action(w,f,g);
        if(cache.transformed[g].size()!=n*m)return {};
        const double factor=double(dimension)*row_characters[g]/(double(f.ops.size())*norm);
        for(std::size_t c=0;c<k;++c){if(members[c]>=m||w.orbitals[members[c]].coefficients.size()!=n)return {};
            for(std::size_t a=0;a<n;++a)projected(a,c)+=factor*cache.transformed[g][a*m+members[c]];}}
    Dense residual(n,k);for(std::size_t c=0;c<k;++c)for(std::size_t a=0;a<n;++a)residual(a,c)=w.orbitals[members[c]].coefficients[a]-projected(a,c);
    const Dense sr=Eigen::Map<const Dense>(w.ao_overlap.data(),n,n)*residual;
    std::vector<double> losses;for(std::size_t c=0;c<k;++c)losses.push_back(residual.col(c).dot(sr.col(c)));
    return losses;
}
bool pure_projection(double loss){return std::isfinite(loss)&&loss>=-1e-10&&loss<=character_tolerance*character_tolerance;}
void retain_projection_residual(NboAomoName& name,double loss){
    if(std::isfinite(loss)&&loss>=-1e-10&&(!name.projection_residual||loss<*name.projection_residual))
        name.projection_residual=std::max(0.,loss);
}
// Every edge is measured from the actual group action in the AO metric.
// The scaled edge bound prevents many individually small couplings from
// evading the complete-column leakage test. Energy never excludes partners.
std::vector<std::vector<std::size_t>> connected_blocks(const Wavefunction& w,const Frame& f,const std::vector<std::size_t>& rows){
    const auto n=std::size_t(w.basis_count),k=rows.size(),m=w.orbitals.size();
    const auto singles=[&](){std::vector<std::vector<std::size_t>> v;for(auto i:rows)v.push_back({i});return v;};
    if(k<2||!f.valid||!n||w.ao_overlap.size()!=n*n)return singles();
    const auto& action=canonical_action(w,f,0);Dense q(n,k),sq(n,k);
    for(std::size_t c=0;c<k;++c){if(rows[c]>=m||w.orbitals[rows[c]].coefficients.size()!=n)return singles();
        for(std::size_t r=0;r<n;++r){q(r,c)=w.orbitals[rows[c]].coefficients[r];sq(r,c)=action.metric_packed[r*m+rows[c]];}}
    const Dense gram=q.transpose()*sq;
    if(!gram.allFinite()||(gram-Dense::Identity(k,k)).cwiseAbs().maxCoeff()>metric_tolerance)return singles();
    std::vector<std::size_t> parent(k);std::iota(parent.begin(),parent.end(),0);
    auto root=[&](std::size_t a){while(parent[a]!=a){parent[a]=parent[parent[a]];a=parent[a];}return a;};
    const double edge_bound=character_tolerance/std::sqrt(double(k));
    for(std::size_t g=0;g<f.ops.size();++g){const auto& cache=canonical_action(w,f,g);if(cache.transformed[g].size()!=n*m)return singles();
        Dense tq(n,k);for(std::size_t c=0;c<k;++c)for(std::size_t r=0;r<n;++r)tq(r,c)=cache.transformed[g][r*m+rows[c]];
        const Dense d=sq.transpose()*tq;if(!d.allFinite())return singles();
        for(std::size_t a=0;a<k;++a)for(std::size_t b=0;b<a;++b)
            if(std::max(std::abs(d(a,b)),std::abs(d(b,a)))>edge_bound)parent[root(b)]=root(a);
    }
    std::map<std::size_t,std::vector<std::size_t>> groups;
    for(std::size_t a=0;a<k;++a)groups[root(a)].push_back(rows[a]);
    std::vector<std::vector<std::size_t>> result;for(auto& [key,members]:groups)result.push_back(std::move(members));return result;
}
std::vector<std::vector<std::size_t>> canonical_blocks(const Wavefunction& w,const Frame& f){
    std::vector<std::vector<std::size_t>> result;
    for(const auto spin:{Spin::Alpha,Spin::Beta}){std::vector<std::size_t> rows;
        for(std::size_t i=0;i<w.orbitals.size();++i)if(w.orbitals[i].spin==spin)rows.push_back(i);
        auto blocks=connected_blocks(w,f,rows);result.insert(result.end(),blocks.begin(),blocks.end());
    }return result;
}
// Restrict the measured action to unchanged, individually pure source columns.
// A graph component is only a proposal: its dimension, complete character row,
// S Gram, transformed norms and leakage must all certify one actual copy.
// Energy never removes a coupling or supplies missing partners.
std::vector<std::vector<std::size_t>> verified_source_copies(const Wavefunction& w,const Frame& f,
    const std::vector<std::size_t>& pure,std::size_t dimension,double norm,
    const std::vector<double>& row_characters,std::size_t maximum_copies){
    std::vector<std::vector<std::size_t>> result;
    for(const auto& candidate:connected_blocks(w,f,pure)){
        if(candidate.size()!=dimension)continue;
        const auto values=characters(w,f,candidate);
        if(values.size()!=row_characters.size())continue;
        bool matches=true;for(std::size_t g=0;g<values.size();++g)
            if(std::abs(values[g]-row_characters[g])>character_tolerance)matches=false;
        const auto losses=projection_residuals(w,f,candidate,dimension,norm,row_characters);
        if(losses.size()!=candidate.size()||!std::all_of(losses.begin(),losses.end(),pure_projection))matches=false;
        if(matches)result.push_back(candidate);
    }
    // Connected components are disjoint. Reject inconsistent counting evidence
    // instead of choosing a subset of equally plausible source copies.
    if(result.size()>maximum_copies)result.clear();
    return result;
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
    struct NamedCharacters {std::string label;std::size_t dimension;double norm=1;std::vector<double> values;};
    std::vector<NamedCharacters> known;
    for(const auto& row:frame.table.rows)known.push_back({row.label,row.dimension,row.character_norm,row.characters});
    if(frame.valid&&frame.linear){
        const char* symbols[]={"Sigma","Pi","Delta","Phi","Gamma"};
        for(unsigned l=0;l<=frame.lmax;++l)for(int parity:{1,-1})for(int sign:{1,-1}){
            if(frame.group=="Cinfv"&&parity<0)continue;if(l&&sign<0)continue;
            NamedCharacters row;row.label=symbols[l];row.dimension=l?2:1;
            if(frame.group=="Dinfh")row.label+=parity>0?"_g":"_u";
            if(!l)row.label+=sign>0?"+":"-";
            for(const auto& op:frame.ops){Vec z{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)z[a]+=op.matrix[3*a+b]*frame.axis[b];
                const bool reverse=dot(z,frame.axis)<0;Mat r=op.matrix;if(reverse)for(auto& x:r)x=-x;
                const bool reflect=det(r)<0;
                double chi=l?(reflect?0:2*std::cos(l*std::acos(std::clamp((trace(r)-1)/2,-1.,1.)))):(reflect?sign:1);
                if(reverse)chi*=parity;row.values.push_back(chi);
            }known.push_back(std::move(row));
        }
    }
    const auto decompose=[&](const std::vector<double>& values,std::size_t dimension){
        std::vector<std::pair<std::size_t,std::size_t>> content;std::vector<double> reconstructed(values.size(),0);std::size_t dimensions=0;
        if(values.empty()||known.empty())return content;
        for(std::size_t r=0;r<known.size();++r){const auto& row=known[r];if(row.values.size()!=values.size()||row.norm<=0)return decltype(content){};
            double inner=std::inner_product(values.begin(),values.end(),row.values.begin(),0.)/(double(values.size())*row.norm);
            if(!std::isfinite(inner))return decltype(content){};const auto copies=std::llround(inner);
            if(copies<0||std::abs(inner-double(copies))>character_tolerance)return decltype(content){};
            if(!copies)continue;content.emplace_back(r,std::size_t(copies));dimensions+=row.dimension*std::size_t(copies);
            for(std::size_t g=0;g<values.size();++g)reconstructed[g]+=double(copies)*row.values[g];
        }
        if(dimensions!=dimension)return decltype(content){};
        for(std::size_t g=0;g<values.size();++g)if(std::abs(values[g]-reconstructed[g])>character_tolerance)return decltype(content){};
        return content;
    };
    const auto blocks=canonical_blocks(w,frame);
    for(std::size_t b=0;b<blocks.size();++b){const auto& members=blocks[b];if(members.empty())continue;
        CharacterFailure failure;const auto values=characters(w,frame,members,&failure);const auto content=decompose(values,members.size());
        const auto scope=w.orbitals[members.front()].spin==Spin::Beta?"canonical beta":"canonical alpha";
        double low=std::numeric_limits<double>::infinity(),high=-low;bool have_energy=true;
        for(auto i:members){const double e=w.orbitals[i].energy_hartree;have_energy=have_energy&&std::isfinite(e);low=std::min(low,e);high=std::max(high,e);
            auto& name=out.canonical[i];name.label="?";name.containing_members=members;
            name.status=!frame.valid?"operation_frame_unavailable":values.empty()?(failure.status.empty()?"symmetry_evidence_unavailable":failure.status):content.empty()?"character_decomposition_unavailable":"mixed_source_member";
            name.detail="Actual full AO coefficients; "+frame.detail+"; "+name.status+(failure.detail.empty()?"":"; "+failure.detail);
            for(const auto& [r,copies]:content)name.containing_irreps.push_back({known[r].label,known[r].dimension,copies});
        }
        if(content.empty()){
            for(auto i:members){Unit u;u.members={i};u.scope=scope;u.energy=w.orbitals[i].energy_hartree;u.energy_available=std::isfinite(u.energy);units.push_back(std::move(u));}
            continue;
        }
        std::vector<std::size_t> proved_counts(known.size(),0);std::set<std::size_t> proved_members;
        for(const auto& [r,copies]:content){const auto& row=known[r];std::vector<std::size_t> pure;
            // Central idempotent projector, with real conjugate-pair norm.
            // It tests each *unchanged source column*. No rotated orbital is
            // substituted for a mixed source MO merely to obtain a name.
            const auto losses=projection_residuals(w,frame,members,row.dimension,row.norm,row.values);
            for(std::size_t c=0;c<losses.size();++c){retain_projection_residual(out.canonical[members[c]],losses[c]);if(pure_projection(losses[c]))pure.push_back(members[c]);}
            for(auto i:pure){auto& name=out.canonical[i];name.verified=true;name.irrep=row.label;name.point_group=frame.group;name.representation_multiplicity=copies;
                name.status="verified_isotypic_member";name.detail="Source column lies in verified isotypic projection; independent occurrence membership unresolved; "+frame.detail;}
            const auto source_copies=verified_source_copies(w,frame,pure,row.dimension,row.norm,row.values,copies);
            for(std::size_t copy=0;copy<source_copies.size();++copy){const auto& partners=source_copies[copy];
                if(std::any_of(partners.begin(),partners.end(),[&](auto i){return proved_members.contains(i);}))continue;
                Unit u;u.members=partners;u.irrep=row.label;u.scope=scope;u.energy=0;u.energy_available=true;
                u.detail="Verified complete irreducible source span under every operation; "+frame.detail;
                for(auto i:partners){u.energy+=w.orbitals[i].energy_hartree/double(partners.size());u.energy_available=u.energy_available&&std::isfinite(w.orbitals[i].energy_hartree);proved_members.insert(i);
                    auto& name=out.canonical[i];name.status="verified_irrep";name.partner_block_id="canonical-irrep:"+std::to_string(b)+":"+row.label+":"+std::to_string(copy);name.partner_block_size=partners.size();}
                ++proved_counts[r];
                units.push_back(std::move(u));
            }
        }
        // Remove only certified complete source copies. A narrower energy
        // range for the remaining counts requires independent closure and an
        // exact character decomposition equal to the original minus these
        // copies. Failed residual evidence retains the original wide bounds.
        std::vector<std::size_t> remaining;for(auto i:members)if(!proved_members.contains(i))remaining.push_back(i);
        std::vector<std::pair<std::size_t,std::size_t>> expected;
        for(const auto& [r,copies]:content)if(copies>proved_counts[r])expected.emplace_back(r,copies-proved_counts[r]);
        const bool residual_verified=!remaining.empty()&&decompose(characters(w,frame,remaining),remaining.size())==expected;
        if(residual_verified){low=std::numeric_limits<double>::infinity();high=-low;have_energy=true;
            for(auto i:remaining){const double e=w.orbitals[i].energy_hartree;have_energy=have_energy&&std::isfinite(e);low=std::min(low,e);high=std::max(high,e);}}
        for(const auto& [r,copies]:expected){Unit count;count.members=residual_verified?remaining:members;count.irrep=known[r].label;count.scope=scope;count.copies=copies;count.copies_resolved=false;count.counting_only=true;
            if(have_energy){count.lower=low;count.upper=high;}units.push_back(std::move(count));}
    }
    assign_ordinals(std::move(units),out.canonical);
    for(auto& name:out.canonical){if(name.verified)name.point_group=frame.group;name.complete_set_ordinal=name.ordinal;name.ordinal_scope="complete canonical set";}
    for(std::size_t i=0;i<out.canonical.size();++i){auto& name=out.canonical[i];if(name.verified)name.label=canonical_mo_display_label(w,i,&name);else if(open)name.label+=w.orbitals[i].spin==Spin::Beta?" [beta]":" [alpha]";}
    if(!salc)return out;
    Wavefunction side;side.atoms=w.atoms;side.shells=w.shells;side.primitives=w.primitives;side.basis_count=w.basis_count;side.ao_overlap=w.ao_overlap;side.orbitals.resize(salc->orbitals.size());
    for(std::size_t i=0;i<salc->orbitals.size();++i){const auto& orbital=salc->orbitals[i];auto& mo=side.orbitals[i];mo.energy_hartree=orbital.energy_hartree.value_or(std::numeric_limits<double>::quiet_NaN());mo.spin=orbital.spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;
        bool valid=associated&&!orbital.terms.empty();mo.coefficients.assign(w.basis_count,0);
        for(const auto& term:orbital.terms){const auto* descriptor=nbo_orbital(data,term.orbital);if(!descriptor||descriptor->coefficients.size()!=w.basis_count||!std::isfinite(term.coefficient)){valid=false;break;}for(std::size_t a=0;a<w.basis_count;++a)mo.coefficients[a]+=term.coefficient*descriptor->coefficients[a];}
        if(!valid)mo.coefficients.clear();
    }
    std::vector<std::vector<std::size_t>> side_blocks;std::vector<std::size_t> side_subspaces;
    std::vector<bool> side_actual(salc->subspaces.size(),false),side_complete(salc->subspaces.size(),false);
    std::vector<CharacterFailure> side_failures(salc->subspaces.size());
    for(std::size_t s=0;s<salc->subspaces.size();++s){const auto& sub=salc->subspaces[s];bool have_coefficients=!sub.orbital_indices.empty();for(auto i:sub.orbital_indices)if(i>=side.orbitals.size()||side.orbitals[i].coefficients.size()!=w.basis_count)have_coefficients=false;
        side_actual[s]=have_coefficients;
        if(have_coefficients){const auto values=characters(side,frame,sub.orbital_indices,&side_failures[s]);
            side_complete[s]=!decompose(values,sub.orbital_indices.size()).empty();}
        const auto pieces=side_complete[s]?connected_blocks(side,frame,sub.orbital_indices):std::vector<std::vector<std::size_t>>{sub.orbital_indices};
        for(const auto& piece:pieces){side_blocks.push_back(piece);side_subspaces.push_back(s);}}

    std::vector<bool> used(salc->orbitals.size(),false);units.clear();
    for(std::size_t b=0;b<side_blocks.size();++b){const auto sub_index=side_subspaces[b];const auto& sub=salc->subspaces[sub_index];Unit u;u.scope=sub.fragment_id+":"+nbo_spin_name(sub.spin);u.detail="Fixed fragment subspace; "+frame.detail;
        bool indices_ok=!side_blocks[b].empty();for(auto i:side_blocks[b])if(i>=out.salc.size()||used[i]||salc->orbitals[i].fragment_id!=sub.fragment_id||salc->orbitals[i].spin!=sub.spin)indices_ok=false;if(!indices_ok)continue;u.members=side_blocks[b];
        const bool verified=frame.valid&&sub.symmetry_verified&&sub.dimension==sub.orbital_indices.size()&&std::isfinite(sub.closure_error)&&sub.closure_error<=character_tolerance&&std::isfinite(sub.orthogonality_error)&&sub.orthogonality_error<=metric_tolerance;
        const bool actual=side_actual[sub_index];
        std::vector<double> values;
        if(verified&&actual&&side_complete[sub_index])values=characters(side,frame,u.members);
        // A stored signature can stand in for omitted AO columns only for its
        // complete recorded span. Failed actual closure/metric evidence never
        // falls back to a stored signature that happens to look irreducible.
        if(verified&&!actual&&u.members==sub.orbital_indices)values=sub.characters;
        const auto content=decompose(values,u.members.size());
        for(auto i:u.members){used[i]=true;const auto& o=salc->orbitals[i];if(o.energy_hartree&&std::isfinite(*o.energy_hartree))u.energy+=*o.energy_hartree/double(u.members.size());else u.energy_available=false;out.salc[i].label="?";out.salc[i].detail=u.detail;
            auto& name=out.salc[i];name.containing_members=u.members;
            const auto& failure=side_failures[sub_index];
            name.status=!frame.valid?"operation_frame_unavailable":actual&&!failure.status.empty()?failure.status:!verified?"subspace_not_verified":values.empty()?(actual?"subspace_not_closed":"stored_characters_unavailable"):content.empty()?"character_decomposition_unavailable":actual?"mixed_source_member":"stored_span_unclassified";
            name.detail+="; "+name.status;
            if(!failure.detail.empty())name.detail+="; "+failure.detail;
            for(const auto& [r,copies]:content)name.containing_irreps.push_back({known[r].label,known[r].dimension,copies});
        }
        if(content.empty()){if(verified)certify_salc_energy_bounds(u,side,w,*salc);units.push_back(std::move(u));continue;}
        if(!actual){
            if(content.size()==1){const auto& [r,copies]=content.front();u.irrep=known[r].label;u.copies=copies;u.copies_resolved=copies==1;
                u.detail+="; complete stored characters verified; irrep occurrence count="+std::to_string(copies)+(u.copies_resolved?"; one ordinal shared by stored partners":"; repeated irrep span; individual copy numbering unresolved");
                for(auto i:u.members){auto& name=out.salc[i];name.status=u.copies_resolved?"verified_stored_irrep_span":"verified_stored_isotypic_span";
                    if(u.copies_resolved){name.partner_block_id="salc-irrep-block:"+std::to_string(b);name.partner_block_size=u.members.size();}}
                if(!u.copies_resolved)certify_salc_energy_bounds(u,side,w,*salc);
                units.push_back(std::move(u));
            }else{
                // Full stored reducible characters provide counts but cannot
                // identify any particular source column without AO evidence.
                for(const auto& [r,copies]:content){Unit count=u;count.irrep=known[r].label;count.copies=copies;count.copies_resolved=false;count.counting_only=true;certify_salc_energy_bounds(count,side,w,*salc);units.push_back(std::move(count));}
            }
            continue;
        }
        std::vector<std::size_t> proved_counts(known.size(),0);std::set<std::size_t> proved_members;
        for(const auto& [r,copies]:content){const auto& row=known[r];std::vector<std::size_t> pure;
            const auto losses=projection_residuals(side,frame,u.members,row.dimension,row.norm,row.values);
            for(std::size_t c=0;c<losses.size();++c){auto& name=out.salc[u.members[c]];retain_projection_residual(name,losses[c]);if(pure_projection(losses[c]))pure.push_back(u.members[c]);}
            Unit evidence=u;certify_salc_energy_bounds(evidence,side,w,*salc);
            for(auto i:pure){auto& name=out.salc[i];name.verified=true;name.irrep=row.label;name.label=orbital_label(row.label);name.representation_multiplicity=copies;
                name.status="verified_isotypic_member";name.detail="Unchanged source SALC column lies in verified S-metric isotypic projection; independent occurrence membership unresolved; "+evidence.detail;}
            const auto source_copies=verified_source_copies(side,frame,pure,row.dimension,row.norm,row.values,copies);
            for(std::size_t copy=0;copy<source_copies.size();++copy){const auto& partners=source_copies[copy];
                if(std::any_of(partners.begin(),partners.end(),[&](auto i){return proved_members.contains(i);}))continue;
                Unit named;named.members=partners;named.irrep=row.label;named.scope=u.scope;
                named.detail="Verified irreducible source SALC span under every operation; "+frame.detail;
                for(auto i:partners){const auto& energy=salc->orbitals[i].energy_hartree;if(energy&&std::isfinite(*energy))named.energy+=*energy/double(partners.size());else named.energy_available=false;proved_members.insert(i);
                    auto& name=out.salc[i];name.status="verified_irrep";name.partner_block_id="salc-irrep:"+std::to_string(b)+":"+row.label+":"+std::to_string(copy);name.partner_block_size=partners.size();}
                ++proved_counts[r];
                units.push_back(std::move(named));
            }
        }
        std::vector<std::size_t> remaining;for(auto i:u.members)if(!proved_members.contains(i))remaining.push_back(i);
        std::vector<std::pair<std::size_t,std::size_t>> expected;
        for(const auto& [r,copies]:content)if(copies>proved_counts[r])expected.emplace_back(r,copies-proved_counts[r]);
        const bool residual_verified=!remaining.empty()&&decompose(characters(side,frame,remaining),remaining.size())==expected;
        for(const auto& [r,copies]:expected){Unit count=u;count.members=residual_verified?remaining:u.members;count.irrep=known[r].label;count.copies=copies;count.copies_resolved=false;count.counting_only=true;
            certify_salc_energy_bounds(count,side,w,*salc);units.push_back(std::move(count));}
    }
    for(std::size_t i=0;i<out.salc.size();++i)if(!used[i]){const auto& o=salc->orbitals[i];Unit u;u.members={i};u.scope=o.fragment_id+":"+nbo_spin_name(o.spin);u.energy=o.energy_hartree.value_or(0);u.energy_available=o.energy_hartree&&std::isfinite(*o.energy_hartree);units.push_back(u);out.salc[i].label="?";out.salc[i].status="subspace_unavailable";out.salc[i].detail="No validated containing SALC subspace";}
    assign_ordinals(std::move(units),out.salc);
    for(std::size_t i=0;i<out.salc.size();++i)if(out.salc[i].verified&&!out.salc[i].ordinal)out.salc[i].label+=" [SALC "+std::to_string(i+1)+"]";
    for(std::size_t i=0;i<out.salc.size();++i)if(open && !salc->orbitals[i].spatial_spin)out.salc[i].label+=salc->orbitals[i].spin==NboSpin::Beta?" [beta]":salc->orbitals[i].spin==NboSpin::Alpha?" [alpha]":" [total]";
    for(auto& name:out.salc){if(name.verified)name.point_group=frame.group;name.complete_set_ordinal=name.ordinal;name.ordinal_scope="complete fragment spin set";}
    (void)data; // Identity is immutable and belongs to the caller's attachment.
    return out;
}

NboAomoNames nbo_aomo_names_for_view(const Wavefunction& w,const NboAomoNames& source,
    const std::vector<std::size_t>& canonical_indices,const std::vector<std::size_t>& salc_indices,
    const NboSalcModel* model,const std::string& scope){
    auto out=source;
    const auto number=[&](std::vector<NboAomoName>& names,const std::vector<std::size_t>& indices,bool side){
        struct Occurrence {std::string family,block;std::vector<std::size_t> members;double energy=0;bool available=true;};
        std::map<std::string,Occurrence> groups;std::set<std::size_t> selected(indices.begin(),indices.end());
        for(auto i:selected){if(i>=names.size()||(!side&&i>=w.orbitals.size())||(side&&(!model||i>=model->orbitals.size())))continue;
            auto& name=names[i];name.ordinal=0;name.ordinal_scope=scope;
            if(!name.verified||name.irrep.empty()||name.partner_block_id.empty()||name.representation_multiplicity!=1)continue;
            if(side&&model->orbitals[i].atoms.size()<=1){name.ordinal=source.salc[i].ordinal;name.ordinal_scope=source.salc[i].ordinal_scope;continue;}
            const std::string family=side?model->orbitals[i].fragment_id+":"+nbo_spin_name(model->orbitals[i].spin):
                w.orbitals[i].spin==Spin::Beta?"canonical beta":"canonical alpha";
            auto& group=groups[family+":"+name.irrep+":"+name.partner_block_id];group.family=family+":"+name.irrep;group.block=name.partner_block_id;group.members.push_back(i);
            const auto energy=side?model->orbitals[i].energy_hartree:std::optional<double>(w.orbitals[i].energy_hartree);
            group.available=group.available&&energy&&std::isfinite(*energy);if(energy&&std::isfinite(*energy))group.energy+=*energy;
        }
        std::vector<Occurrence> ordered;
        for(auto& [key,g]:groups){if(!g.available||g.members.empty()||g.members.size()!=names[g.members.front()].partner_block_size)continue;
            g.energy/=double(g.members.size());ordered.push_back(std::move(g));}
        std::stable_sort(ordered.begin(),ordered.end(),[](const auto& a,const auto& b){if(a.family!=b.family)return a.family<b.family;if(a.energy!=b.energy)return a.energy<b.energy;return a.members.front()<b.members.front();});
        std::map<std::string,std::size_t> counters;
        for(const auto& g:ordered){const auto ordinal=++counters[g.family];for(auto i:g.members){auto& name=names[i];name.ordinal=ordinal;
            name.detail+="; display ordinal in "+scope+": verified visible occurrences ordered by mean energy; exact ties use stable source identity (display convention)";}}
        for(auto i:selected){if(i>=names.size())continue;auto& name=names[i];if(!name.verified)continue;
            if(!side)name.label=canonical_mo_display_label(w,i,&name);
            else {name.label=(name.ordinal?std::to_string(name.ordinal):"")+orbital_label(name.irrep);if(!name.ordinal)name.label+=" [SALC "+std::to_string(i+1)+"]";
                if(model&&model->orbitals[i].spin!=NboSpin::Total)name.label+=" ["+std::string(nbo_spin_name(model->orbitals[i].spin))+"]";}
        }
    };
    number(out.canonical,canonical_indices,false);number(out.salc,salc_indices,true);return out;
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

std::string serialize_orbital_name_json(const NboAomoName& name){
    const auto quote=[](const std::string& text){std::ostringstream s;s<<'"';for(unsigned char c:text){
        if(c=='"'||c=='\\')s<<'\\'<<char(c);else if(c<32)s<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<unsigned(c);else s<<char(c);}s<<'"';return s.str();};
    std::ostringstream out;out<<std::setprecision(17)<<"{\"label\":"<<quote(name.label)
        <<",\"irrep\":"<<quote(name.irrep)<<",\"point_group\":"<<quote(name.point_group)
        <<",\"ordinal\":"<<name.ordinal<<",\"ordinal_scope\":"<<quote(name.ordinal_scope)
        <<",\"complete_set_ordinal\":"<<name.complete_set_ordinal
        <<",\"verified\":"<<(name.verified?"true":"false")<<",\"status\":"<<quote(name.status)
        <<",\"representation_multiplicity\":"<<name.representation_multiplicity
        <<",\"partner_block_id\":"<<quote(name.partner_block_id)<<",\"partner_block_size\":"<<name.partner_block_size
        <<",\"projection_residual_squared\":";
    if(name.projection_residual)out<<*name.projection_residual;else out<<"null";
    out<<",\"containing_members\":[";for(std::size_t i=0;i<name.containing_members.size();++i){if(i)out<<',';out<<name.containing_members[i];}
    out<<"],\"containing_irreps\":[";for(std::size_t i=0;i<name.containing_irreps.size();++i){if(i)out<<',';const auto& r=name.containing_irreps[i];
        out<<"{\"irrep\":"<<quote(r.irrep)<<",\"dimension\":"<<r.dimension<<",\"multiplicity\":"<<r.multiplicity<<'}';}
    out<<"],\"detail\":"<<quote(name.detail)<<'}';return out.str();
}
std::string serialize_orbital_names_json(const NboAomoNames& names){
    std::ostringstream out;out<<"{\"schema\":\"cov.orbital.display-names.v2\",\"canonical\":[";
    const auto rows=[&](const std::vector<NboAomoName>& entries){for(std::size_t i=0;i<entries.size();++i){if(i)out<<',';
        auto row=serialize_orbital_name_json(entries[i]);out<<"{\"index\":"<<i<<','<<row.substr(1);}};
    rows(names.canonical);out<<"],\"salc\":[";rows(names.salc);out<<"]}";return out.str();
}
} // namespace cov::ui
