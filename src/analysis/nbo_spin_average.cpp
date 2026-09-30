#include "cov/nbo_salc.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <set>
#include <sstream>

namespace cov {
namespace {
using M=Eigen::MatrixXd;
using V=Eigen::VectorXd;
using RM=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
constexpr double spatial_tolerance=2e-5,metric_tolerance=5e-5;
double error(const M& a){return a.size()?a.cwiseAbs().maxCoeff():0;}
std::vector<double> flat(const M& a){RM b=a;return {b.data(),b.data()+b.size()};}
M dense(const std::vector<double>& v,std::size_t n){
    if(!n||v.size()!=n*n||!std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);}))return {};
    return Eigen::Map<const RM>(v.data(),n,n);
}
const NboMatrix* matrix(const NboArchive& a,const char* kind,NboSpin spin){
    const NboMatrix* found=nullptr;
    for(const auto& m:a.matrices)if(m.kind==kind&&m.spin==spin){if(found)return nullptr;found=&m;}
    return found;
}
std::string quote(const std::string& s){std::ostringstream out;out<<'"';
    for(unsigned char c:s){if(c=='"'||c=='\\')out<<'\\'<<char(c);else if(c=='\n')out<<"\\n";
        else if(c=='\r')out<<"\\r";else if(c=='\t')out<<"\\t";else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else out<<char(c);}
    out<<'"';return out.str();}
void optional(std::ostream& out,const std::optional<double>& x){if(x&&std::isfinite(*x))out<<*x;else out<<"null";}
template<class T,class F>void array(std::ostream& out,const std::vector<T>& a,F f){out<<'[';bool comma=false;for(const auto& x:a){if(comma)out<<',';comma=true;f(x);}out<<']';}
M columns(const Wavefunction& w,const NboIntegration& d,const NboSalcModel& model,
          const std::vector<std::size_t>& members){
    M c=M::Zero(w.basis_count,members.size());
    for(std::size_t j=0;j<members.size();++j){if(members[j]>=model.orbitals.size())return {};
        for(const auto& term:model.orbitals[members[j]].terms){const auto* o=nbo_orbital(d,term.orbital);
            if(!o||o->coefficients.size()!=w.basis_count||!std::isfinite(term.coefficient))return {};
            c.col(j)+=term.coefficient*Eigen::Map<const V>(o->coefficients.data(),w.basis_count);}}
    return c.allFinite()?c:M{};
}
struct Correspondence {bool valid=false;M beta_from_alpha;double residual=0,metric_error=0;};
Correspondence compare(const M& a,const M& b,const M& s){
    Correspondence out;if(!a.size()||a.rows()!=b.rows()||a.cols()!=b.cols())return out;
    const M aa=a.transpose()*s*a,bb=b.transpose()*s*b,ba=b.transpose()*s*a;
    out.metric_error=std::max(error(aa-M::Identity(a.cols(),a.cols())),error(bb-M::Identity(b.cols(),b.cols())));
    if(out.metric_error>metric_tolerance)return out;
    const Eigen::LDLT<M> ad(aa),bd(bb);if(ad.info()!=Eigen::Success||bd.info()!=Eigen::Success||!ad.isPositive()||!bd.isPositive())return out;
    out.beta_from_alpha=bd.solve(ba);
    const M ra=a-b*out.beta_from_alpha,rb=b-a*ad.solve(ba.transpose());
    const M qa=ra.transpose()*s*ra,qb=rb.transpose()*s*rb;
    if(!qa.allFinite()||!qb.allFinite())return out;
    out.residual=std::sqrt(std::max({0.0,qa.diagonal().maxCoeff(),qb.diagonal().maxCoeff()}));
    out.valid=out.residual<=spatial_tolerance;return out;
}
struct OperatorBlocks {M fock,density;std::string energy_status,density_status;};
OperatorBlocks project_operators(const Wavefunction& w,const NboIntegration& d,
    const NboSalcModel& model,NboSpin spin,const M& c,const M& s){
    OperatorBlocks out;const NboSalcSpinOperator* op=nullptr;
    if(!c.size()){out.energy_status="No verified spatial columns";out.density_status=out.energy_status;return out;}
    for(const auto& item:model.spin_operators)if(item.spin==spin){if(op)return out;op=&item;}
    if(!op){out.energy_status="No verified spin Fock operator";out.density_status="No verified spin density";return out;}
    out.energy_status=op->energy_status;out.density_status=op->density_status;
    M b(w.basis_count,op->basis.size());
    for(std::size_t i=0;i<op->basis.size();++i){const auto* o=nbo_orbital(d,op->basis[i]);
        if(!o||o->coefficients.size()!=w.basis_count)return out;b.col(i)=Eigen::Map<const V>(o->coefficients.data(),w.basis_count);}
    const M gram=b.transpose()*s*b;const Eigen::LDLT<M> solver(gram);
    if(!b.size()||error(gram-M::Identity(b.cols(),b.cols()))>metric_tolerance||solver.info()!=Eigen::Success||!solver.isPositive())return out;
    const M x=solver.solve(b.transpose()*s*c),residual=c-b*x,q=residual.transpose()*s*residual;
    if(!q.allFinite()||std::sqrt(std::max(0.0,q.diagonal().maxCoeff()))>spatial_tolerance){
        out.energy_status="Common spatial direction is outside the verified spin operator basis";out.density_status=out.energy_status;return out;}
    const auto f=dense(op->fock,op->basis.size()),p=dense(op->density,op->basis.size());
    if(f.size())out.fock=x.transpose()*f*x;
    if(p.size())out.density=x.transpose()*p*x;
    return out;
}
std::string without_spin(std::string label){for(const auto* suffix:{" [alpha]"," [beta]"}){
    const std::string ending(suffix);if(label.ends_with(ending))label.resize(label.size()-ending.size());}return label;}
} // namespace

NboRestrictedOpenShellEvidence verify_nbo_restricted_open_shell(const Wavefunction& w,const NboIntegration& d){
    NboRestrictedOpenShellEvidence out;out.alpha_electrons=w.alpha_electrons;out.beta_electrons=w.beta_electrons;
    const auto reject=[&](const char* status,const char* detail){out.status=status;out.detail=detail;return out;};
    // Read the producer's method field, never the arbitrary title or NBO OPEN.
    std::istringstream header(w.source_route);std::string job;header>>job>>out.method;
    std::transform(out.method.begin(),out.method.end(),out.method.begin(),[](unsigned char c){return char(std::toupper(c));});
    if(w.source!=WavefunctionSource::Fchk||out.method.size()<3||out.method.rfind("RO",0)!=0)
        return reject("method_not_verified_ro","Producer method does not identify a restricted open-shell calculation");
    if(w.electron_counts_provenance!=DataProvenance::Producer||w.alpha_electrons<=w.beta_electrons||
       (w.multiplicity_provenance!=DataProvenance::Unavailable&&w.multiplicity!=w.alpha_electrons-w.beta_electrons+1))
        return reject("inconsistent_ro_electrons","Producer electron counts and open-shell multiplicity do not agree");
    const bool shared=w.orbital_occupation_model==OrbitalOccupationModel::CanonicalShared||
        w.orbital_occupation_model==OrbitalOccupationModel::SharedIntegerDeterminant;
    if(!shared&&w.orbital_occupation_model!=OrbitalOccupationModel::ExplicitSpin)
        return reject("occupation_model_unavailable","No verifiable integer shared or explicit-spin occupation model");
    if(!d.dataset.association.compatible||d.canonical_fingerprint!=nbo_canonical_fingerprint(w)||!d.dataset.archive)
        return reject("association_unavailable","Current canonical wavefunction and NBO archive are not associated");
    for(const auto spin:{NboSpin::Alpha,NboSpin::Beta}){
        const NboCanonicalEvidence* ev=nullptr;for(const auto& e:d.dataset.association.canonical_evidence)if(e.spin==spin){if(ev)return reject("ambiguous_association","Multiple same-spin association records");ev=&e;}
        if(!ev||!ev->density_verified||(spin==NboSpin::Alpha&&!ev->direct_fchk_coefficients))
            return reject("spatial_evidence_missing","Matched alpha canonical coefficients and both spin densities are required");
    }
    const auto n=std::size_t(w.basis_count);const auto& archive=*d.dataset.archive;
    const auto* am=matrix(archive,"LCAOMO",NboSpin::Alpha);const auto* bm=matrix(archive,"LCAOMO",NboSpin::Beta);
    const auto* sm=matrix(archive,"OVERLAP",NboSpin::Total);
    if(!am||!bm||!sm||am->rows!=n||bm->rows!=n||am->columns!=n||bm->columns!=n||sm->rows!=n||sm->columns!=n)
        return reject("shared_coefficients_unavailable","Complete same-source alpha/beta spatial columns and metric are required");
    const M a=dense(am->values,n),b=dense(bm->values,n),s=dense(sm->values,n);
    if(!a.size()||!b.size()||!s.size())return reject("invalid_coefficients","Nonfinite or incomplete spatial coefficient evidence");
    double na=0,nb=0;std::set<std::pair<Spin,std::size_t>> seen;
    for(const auto& mo:w.orbitals){
        if(!seen.insert({mo.spin,mo.source_orbital_index}).second||mo.source_orbital_index>=n||
           mo.occupation_provenance==DataProvenance::Unavailable||!std::isfinite(mo.occupation))
            return reject("occupation_identity_invalid","Missing or repeated source orbital occupation identity");
        const double expected=shared?(double(mo.source_orbital_index<w.alpha_electrons)+double(mo.source_orbital_index<w.beta_electrons)):
            double(mo.source_orbital_index<(mo.spin==Spin::Alpha?w.alpha_electrons:w.beta_electrons));
        out.occupation_error=std::max(out.occupation_error,std::abs(mo.occupation-expected));
        if(shared){if(mo.spin!=Spin::Alpha)return reject("shared_model_has_beta_rows","Shared model unexpectedly contains independent beta rows");
            na+=std::min(1.0,mo.occupation);nb+=std::max(0.0,mo.occupation-1);}
        else if(mo.spin==Spin::Alpha)na+=mo.occupation;else nb+=mo.occupation;
    }
    if(out.occupation_error>1e-8||std::abs(na-w.alpha_electrons)>1e-8||std::abs(nb-w.beta_electrons)>1e-8)
        return reject("occupation_mismatch","Canonical occupations do not reproduce both producer electron counts");
    for(std::size_t j=0;j<n;++j){
        const double an=(a.col(j).transpose()*s*a.col(j))(0,0),bn=(b.col(j).transpose()*s*b.col(j))(0,0);
        if(an<1e-18&&bn<1e-18){if(j<w.alpha_electrons)return reject("occupied_null_column","Occupied canonical source column is null");continue;}
        if(std::abs(an-1)>metric_tolerance||std::abs(bn-1)>metric_tolerance)
            return reject("spatial_norm_mismatch","Alpha/beta canonical columns have different active/null spatial dimensions");
        const double phase=(a.col(j).transpose()*s*b.col(j))(0,0)<0?-1:1;
        const V delta=a.col(j)-phase*b.col(j);const double residual=std::sqrt(std::max(0.0,(delta.transpose()*s*delta)(0,0)));
        out.coefficient_residual=std::max(out.coefficient_residual,residual);
        if(residual>out.tolerance||delta.cwiseAbs().maxCoeff()>3e-6)
            return reject("different_canonical_spaces","Corresponding source alpha/beta canonical spatial columns differ");
        ++out.shared_columns;
    }
    if(out.shared_columns<w.alpha_electrons)return reject("incomplete_occupied_space","Shared canonical space omits occupied directions");
    for(const auto spin:{NboSpin::Alpha,NboSpin::Beta}){
        const auto* pm=matrix(archive,"DENSITY",spin);
        if(!pm||pm->rows!=n||pm->columns!=n)return reject("spin_density_missing","Complete producer spin densities are required");
        const M p=dense(pm->values,n);if(!p.size())return reject("spin_density_invalid","Nonfinite or incomplete spin density");
        const auto count=spin==NboSpin::Alpha?w.alpha_electrons:w.beta_electrons;
        if(count>n)return reject("electron_count_out_of_range","More occupied spatial columns than basis directions");
        const M occupied=a.leftCols(count);M expected=occupied*occupied.transpose();
        if(!archive.density_is_bond_order)expected=(s*expected*s).eval();
        out.density_error=std::max(out.density_error,error(expected-p));
        if(out.density_error>2e-5)return reject("shared_spin_density_mismatch","Common canonical columns and their spin occupations do not reproduce both producer spin densities");
    }
    out.verified=true;out.status="verified_restricted_open_shell";
    out.detail="Producer RO method, integer occupations and electron counts agree; associated alpha/beta canonical spatial columns match up to phase, and both spin densities are verified";
    return out;
}

NboSalcModel build_nbo_spin_averaged_model(const Wavefunction& w,const NboIntegration& d,const NboSalcModel& raw){
    NboSalcModel out=raw;out.restricted_open_shell=verify_nbo_restricted_open_shell(w,d);
    out.cache_key+="|spatial-spin-average-v1";
    if(raw.dataset_id!=d.id||raw.canonical_fingerprint!=d.canonical_fingerprint){
        out.restricted_open_shell.verified=false;out.restricted_open_shell.status="stale_source_model";
        out.restricted_open_shell.detail="Source side model belongs to another integration or canonical wavefunction";return out;}
    if(!out.restricted_open_shell.verified||!raw.available)return out;
    const M s=dense(w.ao_overlap,w.basis_count);if(!s.size())return out;
    std::vector<std::size_t> aa,bb;std::set<std::size_t> used;
    for(std::size_t i=0;i<raw.orbitals.size();++i){if(raw.orbitals[i].spin==NboSpin::Alpha)aa.push_back(i);else if(raw.orbitals[i].spin==NboSpin::Beta)bb.push_back(i);}
    struct Block {std::vector<std::size_t> a,b;M mapping;double residual=0,metric_error=0;};
    std::vector<Block> blocks;
    const auto compatible=[&](std::size_t ai,std::size_t bi){const auto& a=raw.orbitals[ai];const auto& b=raw.orbitals[bi];
        return a.atoms==b.atoms&&a.type==b.type&&a.fragment_id==b.fragment_id;};
    // Mutual uniqueness avoids source-order or energy-order tie breaking.
    std::map<std::size_t,std::vector<std::size_t>> candidates,reverse;
    std::map<std::pair<std::size_t,std::size_t>,Correspondence> matches;
    for(auto ai:aa)for(auto bi:bb)if(compatible(ai,bi)){
        auto match=compare(columns(w,d,raw,{ai}),columns(w,d,raw,{bi}),s);
        if(match.valid){candidates[ai].push_back(bi);reverse[bi].push_back(ai);matches[{ai,bi}]=std::move(match);}}
    for(auto ai:aa)if(candidates[ai].size()==1){const auto bi=candidates[ai].front();if(reverse[bi].size()!=1)continue;
        const auto& m=matches.at({ai,bi});blocks.push_back({{ai},{bi},m.beta_from_alpha,m.residual,m.metric_error});used.insert(ai);used.insert(bi);}
    // Remaining complete invariant spaces, including rotated repeated irreps,
    // can share a common basis. Local single-atom families use the same proof.
    std::vector<std::vector<std::size_t>> ag,bg;
    for(const auto& sub:raw.subspaces){std::vector<std::size_t> remaining;
        for(auto i:sub.orbital_indices)if(!used.contains(i))remaining.push_back(i);
        if(remaining.size()<2)continue;
        if(!sub.symmetry_verified&&raw.orbitals[remaining.front()].atoms.size()!=1)continue;
        if(sub.spin==NboSpin::Alpha)ag.push_back(remaining);else if(sub.spin==NboSpin::Beta)bg.push_back(remaining);}
    candidates.clear();reverse.clear();matches.clear();
    for(std::size_t ai=0;ai<ag.size();++ai)for(std::size_t bi=0;bi<bg.size();++bi)
        if(ag[ai].size()==bg[bi].size()&&compatible(ag[ai].front(),bg[bi].front())){
            auto match=compare(columns(w,d,raw,ag[ai]),columns(w,d,raw,bg[bi]),s);
            if(match.valid){candidates[ai].push_back(bi);reverse[bi].push_back(ai);matches[{ai,bi}]=std::move(match);}}
    for(std::size_t ai=0;ai<ag.size();++ai)if(candidates[ai].size()==1){auto bi=candidates[ai].front();if(reverse[bi].size()!=1)continue;
        const auto& m=matches.at({ai,bi});blocks.push_back({ag[ai],bg[bi],m.beta_from_alpha,m.residual,m.metric_error});
        used.insert(ag[ai].begin(),ag[ai].end());used.insert(bg[bi].begin(),bg[bi].end());}
    std::map<std::size_t,NboSpatialSpinInfo> merged;std::set<std::size_t> removed;
    // Each old side becomes one or more common-basis components. Used only to
    // transform signed projection coefficients; channel weights are never added.
    std::map<std::size_t,std::vector<std::pair<std::size_t,double>>> remap;
    // Project each verified operator once, then slice source-member blocks.
    // Re-factoring the full NAO metric for every singleton is unnecessary.
    const auto alpha_ops=project_operators(w,d,raw,NboSpin::Alpha,columns(w,d,raw,aa),s);
    const auto beta_ops=project_operators(w,d,raw,NboSpin::Beta,columns(w,d,raw,bb),s);
    std::map<std::size_t,std::size_t> alpha_position,beta_position;
    for(std::size_t i=0;i<aa.size();++i)alpha_position[aa[i]]=i;
    for(std::size_t i=0;i<bb.size();++i)beta_position[bb[i]]=i;
    const auto slice=[](const OperatorBlocks& all,const std::map<std::size_t,std::size_t>& positions,const std::vector<std::size_t>& indices){
        OperatorBlocks result;result.energy_status=all.energy_status;result.density_status=all.density_status;
        const auto select=[&](const M& full){if(!full.size())return M{};M part(indices.size(),indices.size());
            for(std::size_t i=0;i<indices.size();++i)for(std::size_t j=0;j<indices.size();++j)part(i,j)=full(positions.at(indices[i]),positions.at(indices[j]));return part;};
        result.fock=select(all.fock);result.density=select(all.density);return result;};
    for(const auto& block:blocks){const auto k=block.a.size();
        const auto fa=slice(alpha_ops,alpha_position,block.a),fb=slice(beta_ops,beta_position,block.b);
        const std::string block_id="ro-spatial:"+raw.orbitals[block.a.front()].id;
        for(std::size_t j=0;j<k;++j){NboSpatialSpinInfo info;
            info.id=block_id+":member"+std::to_string(j);info.block_id=block_id;info.member_index=j;info.dimension=k;
            info.correspondence=k==1?"same_spatial_orbital_up_to_phase":"verified_equal_complete_spatial_subspace";
            info.detail="Common basis is the original alpha spatial basis; beta operators and density are expressed in that basis. Source coefficients are unchanged";
            info.ro_method=out.restricted_open_shell.method;info.tolerance=spatial_tolerance;
            info.metric_error=block.metric_error;info.spatial_residual=block.residual;
            for(int ch=0;ch<2;++ch){NboSpinChannelProjection channel;channel.spin=ch?NboSpin::Beta:NboSpin::Alpha;
                const auto& indices=ch?block.b:block.a;const auto& op=ch?fb:fa;V mapping=V::Zero(k);
                if(ch)mapping=block.mapping.col(j);else mapping[j]=1;
                channel.mapping.assign(mapping.data(),mapping.data()+mapping.size());
                for(auto index:indices){const auto& source=raw.orbitals[index];channel.members.push_back({source.id,source.label,source.terms,source.energy_hartree,source.occupation});}
                if(op.fock.size()){channel.fock_matrix=flat(op.fock);channel.energy_hartree=(mapping.transpose()*op.fock*mapping)(0,0);channel.energy_status="verified_common_basis_spin_fock_expectation";}
                else if(k==1&&channel.members[0].energy_hartree){channel.energy_hartree=channel.members[0].energy_hartree;channel.energy_status="verified_single_spatial_orbital_expectation";}
                else channel.energy_status=op.energy_status.empty()?"Spin Fock unavailable":op.energy_status;
                if(op.density.size()){channel.density_matrix=flat(op.density);channel.occupation=(mapping.transpose()*op.density*mapping)(0,0);channel.occupation_status="verified_common_basis_spin_density_expectation";}
                else if(k==1&&channel.members[0].occupation){channel.occupation=channel.members[0].occupation;channel.occupation_status="verified_single_spatial_orbital_occupation";}
                else channel.occupation_status=op.density_status.empty()?"Spin density unavailable":op.density_status;
                info.channels.push_back(std::move(channel));}
            const auto& ca=info.channels[0];const auto& cb=info.channels[1];
            if(ca.energy_hartree&&cb.energy_hartree){info.energy_hartree=(*ca.energy_hartree+*cb.energy_hartree)/2;info.energy_status="spin_average_of_verified_common_basis_expectations";}
            else info.energy_status="Mean unavailable: alpha: "+ca.energy_status+"; beta: "+cb.energy_status;
            if(ca.occupation&&cb.occupation){info.occupation=*ca.occupation+*cb.occupation;info.occupation_status="sum_of_two_common_basis_spin_density_expectations";}
            else info.occupation_status="Total unavailable: alpha: "+ca.occupation_status+"; beta: "+cb.occupation_status;
            merged[block.a[j]]=std::move(info);remap[block.a[j]]={{block.a[j],1}};
            for(std::size_t i=0;i<k;++i)remap[block.b[i]].push_back({block.a[j],block.mapping(i,j)});
        }removed.insert(block.b.begin(),block.b.end());}
    // Total is an alternative representation, never a third spin population.
    // Compare complete blocks first so a rotated Total basis is not duplicated.
    ag.clear();bg.clear();candidates.clear();reverse.clear();matches.clear();
    for(const auto& sub:raw.subspaces){
        if(sub.orbital_indices.size()<2)continue;
        if(sub.spin==NboSpin::Alpha&&std::all_of(sub.orbital_indices.begin(),sub.orbital_indices.end(),[&](auto i){return merged.contains(i);}))ag.push_back(sub.orbital_indices);
        else if(sub.spin==NboSpin::Total)bg.push_back(sub.orbital_indices);}
    for(std::size_t ai=0;ai<ag.size();++ai)for(std::size_t bi=0;bi<bg.size();++bi)
        if(ag[ai].size()==bg[bi].size()&&compatible(ag[ai].front(),bg[bi].front())){
            auto match=compare(columns(w,d,raw,ag[ai]),columns(w,d,raw,bg[bi]),s);
            if(match.valid){candidates[ai].push_back(bi);reverse[bi].push_back(ai);matches[{ai,bi}]=std::move(match);}}
    for(std::size_t ai=0;ai<ag.size();++ai)if(candidates[ai].size()==1){const auto bi=candidates[ai].front();if(reverse[bi].size()!=1)continue;
        const auto& m=matches.at({ai,bi});
        for(std::size_t j=0;j<ag[ai].size();++j){auto& info=merged.at(ag[ai][j]);
            info.total_alias_residual=m.residual;info.total_alias_metric_error=m.metric_error;
            for(std::size_t k=0;k<bg[bi].size();++k){const auto& o=raw.orbitals[bg[bi][k]];
                info.total_aliases.push_back({o.id,o.label,o.terms,o.energy_hartree,o.occupation});info.total_alias_mapping.push_back(m.beta_from_alpha(k,j));}}
        removed.insert(bg[bi].begin(),bg[bi].end());}
    for(std::size_t i=0;i<raw.orbitals.size();++i)if(raw.orbitals[i].spin==NboSpin::Total&&!removed.contains(i)){
        std::vector<std::size_t> aliases;
        for(const auto& [alpha,info]:merged)if(info.total_aliases.empty()&&compatible(alpha,i)&&
            compare(columns(w,d,raw,{alpha}),columns(w,d,raw,{i}),s).valid)aliases.push_back(alpha);
        if(aliases.size()==1){const auto& o=raw.orbitals[i];
            const auto match=compare(columns(w,d,raw,{aliases.front()}),columns(w,d,raw,{i}),s);
            merged.at(aliases.front()).total_aliases.push_back({o.id,o.label,o.terms,o.energy_hartree,o.occupation});
            merged.at(aliases.front()).total_alias_mapping.push_back(match.beta_from_alpha(0,0));
            merged.at(aliases.front()).total_alias_residual=match.residual;merged.at(aliases.front()).total_alias_metric_error=match.metric_error;
            removed.insert(i);}
    }
    out.orbitals.clear();out.subspaces.clear();out.links.clear();
    std::map<std::size_t,std::size_t> new_index;
    for(std::size_t i=0;i<raw.orbitals.size();++i){if(removed.contains(i))continue;auto orbital=raw.orbitals[i];
        if(const auto it=merged.find(i);it!=merged.end()){
            orbital.spatial_spin=it->second;orbital.id=it->second.id;orbital.subspace_id+="|spatial";
            orbital.spin=NboSpin::Total;orbital.label=without_spin(orbital.label);
            orbital.energy_hartree=it->second.energy_hartree;orbital.occupation=it->second.occupation;
            orbital.energy_semantics="spin-averaged side-orbital Fock expectation; (alpha+beta)/2 in one verified spatial basis";
            orbital.detail+="; "+it->second.energy_status;orbital.spin_correspondence_status=it->second.correspondence;
            ++out.merged_spatial_count;
        }else if(orbital.spin!=NboSpin::Total){orbital.spin_correspondence_status="kept_separate: no unique equal spatial orbital or complete subspace verified";orbital.detail+="; "+orbital.spin_correspondence_status;++out.separate_spin_count;}
        new_index[i]=out.orbitals.size();out.orbitals.push_back(std::move(orbital));}
    for(const auto& sub:raw.subspaces){for(bool spatial:{false,true}){NboSalcSubspace converted=sub;converted.orbital_indices.clear();
        for(auto i:sub.orbital_indices)if(new_index.contains(i)&&merged.contains(i)==spatial)converted.orbital_indices.push_back(new_index[i]);
        if(converted.orbital_indices.empty())continue;
        if(spatial){converted.id+="|spatial";converted.spin=NboSpin::Total;}
        if(converted.orbital_indices.size()!=sub.orbital_indices.size()){converted.symmetry_verified=false;converted.characters.clear();converted.multiplicity=0;converted.irrep_dimension=0;}
        converted.dimension=converted.orbital_indices.size();out.subspaces.push_back(std::move(converted));}}
    std::map<std::pair<std::size_t,std::size_t>,NboSalcLink> links;
    for(const auto& link:raw.links){
        if(link.side_index>=raw.orbitals.size()||link.canonical_index>=w.orbitals.size())continue;
        const auto spin=raw.orbitals[link.side_index].spin;
        const auto target_spin=w.orbitals[link.canonical_index].spin==Spin::Beta?NboSpin::Beta:NboSpin::Alpha;
        // Shared canonical MOs have one projection, not alpha+beta populations.
        if(spin!=NboSpin::Total&&spin!=target_spin)continue;
        const auto append=[&](std::size_t target,double mapping){auto& value=links[{new_index.at(target),link.canonical_index}];
            value.side_index=new_index.at(target);value.canonical_index=link.canonical_index;value.source_spin=spin;
            value.coefficient+=mapping*link.coefficient;value.source_side_indices.push_back(link.side_index);
            value.source_coefficients.push_back(link.coefficient);value.basis_mapping.push_back(mapping);};
        if(const auto it=remap.find(link.side_index);it!=remap.end())for(const auto& [target,coefficient]:it->second)append(target,coefficient);
        else if(new_index.contains(link.side_index))append(link.side_index,1);
    }
    for(auto& [identity,link]:links){link.weight=link.coefficient*link.coefficient;out.links.push_back(std::move(link));}
    out.spin_averaged=out.merged_spatial_count>0;
    return out;
}

std::string serialize_nbo_spatial_spin_json(const NboSpatialSpinInfo& info){
    std::ostringstream out;out<<std::setprecision(17)<<"{\"id\":"<<quote(info.id)<<",\"block_id\":"<<quote(info.block_id)
        <<",\"display_mode\":\"spin_averaged_spatial\",\"correspondence\":"<<quote(info.correspondence)<<",\"detail\":"<<quote(info.detail)
        <<",\"ro_method\":"<<quote(info.ro_method)<<",\"common_basis\":\"original alpha spatial members\",\"member_index\":"<<info.member_index
        <<",\"dimension\":"<<info.dimension<<",\"tolerance\":"<<info.tolerance<<",\"metric_error\":"<<info.metric_error<<",\"spatial_residual\":"<<info.spatial_residual
        <<",\"energy_definition\":\"chi^T (F_alpha+F_beta) chi / 2\",\"occupation_definition\":\"chi^T (P_alpha+P_beta) chi\",\"energy_hartree\":";
    optional(out,info.energy_hartree);out<<",\"occupation\":";optional(out,info.occupation);
    out<<",\"energy_status\":"<<quote(info.energy_status)<<",\"occupation_status\":"<<quote(info.occupation_status)<<",\"channels\":";
    array(out,info.channels,[&](const auto& channel){out<<"{\"spin\":"<<quote(nbo_spin_name(channel.spin))<<",\"energy_hartree\":";optional(out,channel.energy_hartree);
        out<<",\"occupation\":";optional(out,channel.occupation);out<<",\"energy_status\":"<<quote(channel.energy_status)<<",\"occupation_status\":"<<quote(channel.occupation_status)
            <<",\"mapping\":";array(out,channel.mapping,[&](double x){out<<x;});out<<",\"source_fock_matrix_hartree\":";array(out,channel.fock_matrix,[&](double x){out<<x;});
        out<<",\"source_density_matrix\":";array(out,channel.density_matrix,[&](double x){out<<x;});out<<",\"members\":";
        array(out,channel.members,[&](const auto& member){out<<"{\"id\":"<<quote(member.id)<<",\"label\":"<<quote(member.label)<<",\"energy_hartree\":";
            optional(out,member.energy_hartree);out<<",\"occupation\":";optional(out,member.occupation);out<<",\"terms\":";
            array(out,member.terms,[&](const auto& term){out<<"{\"kind\":"<<quote(nbo_orbital_kind_name(term.orbital.kind))<<",\"spin\":"<<quote(nbo_spin_name(term.orbital.spin))
                <<",\"index\":"<<term.orbital.index<<",\"coefficient\":"<<term.coefficient<<'}';});out<<'}';});out<<'}';});
    out<<",\"total_aliases_not_added\":";array(out,info.total_aliases,[&](const auto& member){
        out<<"{\"id\":"<<quote(member.id)<<",\"label\":"<<quote(member.label)<<",\"energy_hartree\":";optional(out,member.energy_hartree);
        out<<",\"occupation\":";optional(out,member.occupation);out<<",\"terms\":";
        array(out,member.terms,[&](const auto& term){out<<"{\"kind\":"<<quote(nbo_orbital_kind_name(term.orbital.kind))<<",\"spin\":"<<quote(nbo_spin_name(term.orbital.spin))
            <<",\"index\":"<<term.orbital.index<<",\"coefficient\":"<<term.coefficient<<'}';});out<<'}';});
    out<<",\"total_alias_mapping\":";array(out,info.total_alias_mapping,[&](double x){out<<x;});
    out<<",\"total_alias_residual\":"<<info.total_alias_residual<<",\"total_alias_metric_error\":"<<info.total_alias_metric_error;
    out<<'}';return out.str();
}

std::string serialize_nbo_orbital_selection_json(const NboOrbitalSelection& selection){
    std::ostringstream out;out<<std::setprecision(17)<<"{\"dataset_id\":"<<quote(selection.dataset_id)
        <<",\"label\":"<<quote(selection.label)<<",\"semantic_kind\":"<<quote(selection.semantic_kind)
        <<",\"group_id\":"<<quote(selection.group_id)<<",\"source_id\":"<<quote(selection.source_id)
        <<",\"mode\":"<<static_cast<int>(selection.mode)<<",\"normalize\":"<<(selection.normalize?"true":"false")
        <<",\"target_canonical_index\":";
    if(selection.target_canonical_index)out<<*selection.target_canonical_index;else out<<"null";
    out<<",\"terms\":";array(out,selection.terms,[&](const auto& term){out<<"{\"kind\":"<<quote(nbo_orbital_kind_name(term.orbital.kind))
        <<",\"spin\":"<<quote(nbo_spin_name(term.orbital.spin))<<",\"index\":"<<term.orbital.index<<",\"coefficient\":"<<term.coefficient<<'}';});
    out<<",\"spatial_spin\":"<<(selection.spatial_spin?serialize_nbo_spatial_spin_json(*selection.spatial_spin):"null")<<'}';return out.str();
}
} // namespace cov
