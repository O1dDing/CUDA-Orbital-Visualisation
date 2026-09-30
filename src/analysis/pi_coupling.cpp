#include "cov/pi_coupling.hpp"
#include "cov/nbo_salc.hpp"
#include "cov/local_angular_projection.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <sstream>
#include <span>

namespace cov { namespace {
using M=Eigen::MatrixXd;
using V=Eigen::VectorXd;
using RM=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using V3=Eigen::Vector3d;
M mat(const NboMatrix& x){return Eigen::Map<const RM>(x.values.data(),x.rows,x.columns);}
double maxabs(const M& x){return x.size()?x.cwiseAbs().maxCoeff():0;}
const NboMatrix* unique_matrix(const std::vector<NboMatrix>& rows,const char* kind,NboSpin spin){
    const NboMatrix* out=nullptr;
    for(const auto& row:rows)if(row.kind==kind&&row.spin==spin){if(out)return nullptr;out=&row;}
    return out;
}
bool full(const NboMatrix* x,std::size_t n){
    return x&&x->rows==n&&x->columns==n&&x->values.size()==n*n&&
        std::all_of(x->values.begin(),x->values.end(),[](double v){return std::isfinite(v);});
}
bool valence(const NboNao& x,char angular){
    return x.type.find("Val")!=std::string::npos &&
        x.type.find(angular)!=std::string::npos;
}
int axis_of(const std::string& angular){
    if(angular=="px")return 0;if(angular=="py")return 1;
    if(angular=="pz")return 2;return -1;
}
struct Shells {
    std::map<std::size_t,std::array<std::size_t,3>> p;
    std::map<std::size_t,std::vector<std::size_t>> d;
};
Shells shells(const NboIntegration& integration,NboSpin spin){
    Shells out;std::map<std::size_t,std::array<bool,3>> seen;
    for(const auto& row:integration.dataset.naos){
        if(row.spin!=spin||!row.atom||!row.id)continue;
        const auto atom=row.atom-1;
        if(valence(row,'p')){
            const int axis=axis_of(row.angular);
            if(axis>=0){
                if(seen[atom][axis])out.p.erase(atom);
                else {out.p[atom][axis]=row.id;seen[atom][axis]=true;}
            }
        }
        if(valence(row,'d') && row.angular.size() && row.angular[0]=='d')
            out.d[atom].push_back(row.id);
    }
    for(auto it=out.p.begin();it!=out.p.end();){
        const auto& flags=seen[it->first];
        if(!(flags[0]&&flags[1]&&flags[2]))it=out.p.erase(it);
        else ++it;
    }
    for(auto it=out.d.begin();it!=out.d.end();){
        if(it->second.size()!=5)it=out.d.erase(it);else ++it;
    }
    return out;
}
M d_columns(std::size_t n,const std::vector<std::size_t>& ids){
    M q=M::Zero(n,ids.size());
    for(std::size_t i=0;i<ids.size();++i)q(ids[i]-1,i)=1;
    return q;
}
V3 bond_axis(const Wavefunction& w,std::size_t a,std::size_t b){
    V3 v(w.atoms[b].x-w.atoms[a].x,w.atoms[b].y-w.atoms[a].y,
         w.atoms[b].z-w.atoms[a].z);
    return v.norm()>1e-9?v.normalized():V3::Zero();
}
std::array<double,9> axis_frame(const V3& z){
    V3 guide=std::abs(z.z())<0.85?V3::UnitZ():V3::UnitY();
    const V3 x=guide.cross(z).normalized(),y=z.cross(x);
    return {x.x(),y.x(),z.x(),x.y(),y.y(),z.y(),x.z(),y.z(),z.z()};
}
struct PiFrame{M q;double leakage=0,centre=0,minimum_p=0,residual=0;
    std::array<double,3> spectrum{};};
std::optional<PiFrame> derived_pi_frame(const Wavefunction& w,const NboIntegration& data,
                        std::size_t atom,const std::array<std::size_t,3>& ids,
                        NboSpin spin,const V3& axis,std::string& reason){
    const auto fail=[&](const char* why)->std::optional<PiFrame>{reason=why;return std::nullopt;};
    if(axis.norm()<0.5)return fail("Bond axis cannot be defined");
    const auto n=static_cast<std::size_t>(w.basis_count);
    std::array<const NboOrbitalDescriptor*,3> descriptors{};
    for(int component=0;component<3;++component){
        descriptors[component]=nbo_orbital(data,{NboOrbitalKind::NAO,spin,ids[component]-1});
        if(!descriptors[component]&&spin!=NboSpin::Total)
            descriptors[component]=nbo_orbital(data,{NboOrbitalKind::NAO,
                NboSpin::Total,ids[component]-1});
        if(!descriptors[component]||descriptors[component]->coefficients.size()!=n)
            return fail("NAO descriptor or complete coefficient column is missing");
    }
    Wavefunction view=w;view.orbitals.clear();
    const std::array<std::array<double,3>,6> samples{{
        {1,0,0},{0,1,0},{0,0,1},
        {std::numbers::sqrt2/2,std::numbers::sqrt2/2,0},
        {std::numbers::sqrt2/2,0,std::numbers::sqrt2/2},
        {0,std::numbers::sqrt2/2,std::numbers::sqrt2/2}}};
    for(const auto& sample:samples){
        MolecularOrbital mo;mo.spin=spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;
        mo.occupation=1;mo.coefficients.assign(n,0);
        for(int component=0;component<3;++component)
            for(std::size_t k=0;k<n;++k)
                mo.coefficients[k]+=sample[component]*
                    descriptors[component]->coefficients[k];
        view.orbitals.push_back(std::move(mo));
    }
    try{
        LocalAngularProjectionWorkspace workspace(view,atom,axis_frame(axis));
        std::array<double,6> sigma{},pi{};
        double centre=1,residual=0;
        const M metric=Eigen::Map<const RM>(w.ao_overlap.data(),n,n);
        for(std::size_t sample=0;sample<6;++sample){
            const auto projection=workspace.project(std::span<const std::size_t>(&sample,1));
            if(projection.status!=MetricSubspaceStatus::Available||
               projection.represented_spin_orbital_rank!=1||
               !std::isfinite(projection.angular_partition_residual)||
               std::abs(projection.angular_partition_residual)>2e-5)
                return fail("Local angular projection rank or partition residual failed");
            const V coeff=Eigen::Map<const V>(view.orbitals[sample].coefficients.data(),n);
            const double norm=coeff.dot(metric*coeff);
            if(!std::isfinite(norm)||norm<0.7||std::abs(norm-1)>5e-5)
                return fail("NAO metric normalization failed");
            sigma[sample]=norm*projection.component_projection_traces[1];
            pi[sample]=norm*(projection.component_projection_traces[2]+
                       projection.component_projection_traces[3]);
            centre=std::min(centre,projection.centre_mean_fraction);
            residual=std::max(residual,std::abs(projection.angular_partition_residual));
        }
        M gs=M::Zero(3,3),gp=M::Zero(3,3);
        for(int i=0;i<3;++i){gs(i,i)=sigma[i];gp(i,i)=sigma[i]+pi[i];}
        const std::array<std::pair<int,int>,3> pairs{{{0,1},{0,2},{1,2}}};
        for(int pair=0;pair<3;++pair){
            const auto [i,j]=pairs[pair];
            const double s=sigma[3+pair]-(sigma[i]+sigma[j])*0.5;
            const double p=sigma[3+pair]+pi[3+pair]-(gp(i,i)+gp(j,j))*0.5;
            gs(i,j)=gs(j,i)=s;gp(i,j)=gp(j,i)=p;
        }
        Eigen::SelfAdjointEigenSolver<M> ps(gp);
        if(ps.info()!=Eigen::Success||ps.eigenvalues()[0]<0.7||centre<0.7)
            return fail("Local p metric or centre projection is below the supported threshold");
        const M gpi=gp-gs;
        Eigen::GeneralizedSelfAdjointEigenSolver<M> angular(gpi,gp);
        if(angular.info()!=Eigen::Success||!angular.eigenvalues().allFinite()||
           angular.eigenvalues()[0]< -2e-4||angular.eigenvalues()[2]>1+2e-4||
           std::abs(angular.eigenvalues()[0])>2e-4||
           std::abs(1-angular.eigenvalues()[1])>2e-4||
           std::abs(1-angular.eigenvalues()[2])>2e-4)
            return fail("Generalized pi projector eigenvalues failed the sigma/pi separation gate");
        const M eig=angular.eigenvectors().rightCols(2);
        Eigen::HouseholderQR<M> qr(eig);
        const M local=qr.householderQ()*M::Identity(3,2);
        const M local_gp=local.transpose()*gp*local;
        const M local_gs=local.transpose()*gs*local;
        Eigen::GeneralizedSelfAdjointEigenSolver<M> verify(local_gs,local_gp);
        if(verify.info()!=Eigen::Success||!verify.eigenvalues().allFinite()||
           verify.eigenvalues().maxCoeff()>2e-4)return fail("Residual sigma leakage exceeds tolerance");
        PiFrame result;result.q=M::Zero(n,2);
        for(int row=0;row<3;++row)result.q.row(ids[row]-1)=local.row(row);
        result.leakage=std::max(std::abs(angular.eigenvalues()[0]),
            std::max(std::abs(1-angular.eigenvalues()[1]),
                     std::abs(1-angular.eigenvalues()[2])));
        result.centre=centre;result.minimum_p=ps.eigenvalues()[0];
        result.residual=residual;
        for(int i=0;i<3;++i)result.spectrum[i]=angular.eigenvalues()[i];
        return result;
    }catch(const std::exception& e){reason=std::string("Local angular projector unavailable: ")+e.what();return std::nullopt;}
}
std::string verified_symmetry(const Wavefunction& w,std::size_t index){
    for(const auto& assignment:w.derived_orbital_symmetry_assignments)
        if(std::find(assignment.orbital_indices.begin(),assignment.orbital_indices.end(),index)!=
           assignment.orbital_indices.end() &&
           std::isfinite(assignment.subspace_retention) && assignment.subspace_retention>0.999)
            return assignment.label;
    return w.orbitals[index].symmetry_provenance!=DataProvenance::Unavailable?
        w.orbitals[index].symmetry:"";
}
std::vector<std::vector<std::size_t>> canonical_groups(const Wavefunction& w,NboSpin spin){
    std::vector<std::size_t> indices;
    for(std::size_t i=0;i<w.orbitals.size();++i)
        if((spin==NboSpin::Beta)==(w.orbitals[i].spin==Spin::Beta))indices.push_back(i);
    std::sort(indices.begin(),indices.end(),[&](auto a,auto b){
        if(w.orbitals[a].energy_hartree!=w.orbitals[b].energy_hartree)
            return w.orbitals[a].energy_hartree<w.orbitals[b].energy_hartree;
        return a<b;
    });
    std::vector<std::vector<std::size_t>> groups;
    for(auto index:indices){
        if(groups.empty()||
           std::abs(w.orbitals[index].energy_hartree-
                    w.orbitals[groups.back().front()].energy_hartree)>1e-6 ||
           verified_symmetry(w,index)!=verified_symmetry(w,groups.back().front()))
            groups.push_back({});
        groups.back().push_back(index);
    }
    return groups;
}
struct E2Support{double donor=0,acceptor=0;std::vector<NboPiDirectionProjection> rows;};
E2Support directed_e2(const NboIntegration& data,const M& transform,
                 const M& qc,const M& ql,const std::vector<std::size_t>& centre,
                 const std::vector<std::size_t>& ligand,NboSpin spin,bool ligand_to_centre){
    E2Support support;std::set<std::size_t> donor_ids,acceptor_ids;
    for(const auto& row:data.structure){
        if(row.kind!="donor_acceptor"||row.spin!=spin||row.orbitals.size()<2)continue;
        const auto* donor=nbo_orbital(data,row.orbitals[0]);
        const auto* acceptor=nbo_orbital(data,row.orbitals[1]);
        if(!donor||!acceptor)continue;
        if(row.id.rfind("e2:",0)!=0)continue;
        const auto source_index=std::stoull(row.id.substr(3));
        if(source_index>=data.dataset.e2.size())continue;
        const auto& printed=data.dataset.e2[source_index];
        if(!(printed.value>0)||!(printed.energy_gap_hartree>0))continue;
        const auto& from=ligand_to_centre?ligand:centre;
        const auto& to=ligand_to_centre?centre:ligand;
        const bool source=std::any_of(donor->atoms.begin(),donor->atoms.end(),
            [&](auto a){return std::find(from.begin(),from.end(),a)!=from.end();});
        const bool target=std::any_of(acceptor->atoms.begin(),acceptor->atoms.end(),
            [&](auto a){return std::find(to.begin(),to.end(),a)!=to.end();});
        if(!source||!target||row.orbitals[0].index>=static_cast<std::size_t>(transform.cols())||
           row.orbitals[1].index>=static_cast<std::size_t>(transform.cols()))continue;
        const V bd=transform.col(row.orbitals[0].index),
                ba=transform.col(row.orbitals[1].index);
        if(std::abs(bd.squaredNorm()-1)>5e-5||std::abs(ba.squaredNorm()-1)>5e-5)
            continue;
        NboPiDirectionProjection evidence;
        evidence.e2_id=row.id;
        evidence.donor_nbo_id=row.orbitals[0].index+1;
        evidence.acceptor_nbo_id=row.orbitals[1].index+1;
        evidence.donor_ligand_weight=(ql.transpose()*bd).squaredNorm();
        evidence.donor_centre_weight=(qc.transpose()*bd).squaredNorm();
        evidence.acceptor_ligand_weight=(ql.transpose()*ba).squaredNorm();
        evidence.acceptor_centre_weight=(qc.transpose()*ba).squaredNorm();
        evidence.source=row.source;
        const double donor_target=ligand_to_centre?evidence.donor_ligand_weight:
            evidence.donor_centre_weight;
        const double acceptor_target=ligand_to_centre?evidence.acceptor_centre_weight:
            evidence.acceptor_ligand_weight;
        // A σ→π or π→σ E2 row cannot be combined with another row to
        // manufacture a π→π transition. The product is a gate, never an
        // allocated fraction of the printed E(2) energy.
        if(donor_target>=1e-4 && acceptor_target>=1e-4 &&
           donor_target*acceptor_target>1e-8){
            donor_ids.insert(row.orbitals[0].index);
            acceptor_ids.insert(row.orbitals[1].index);
        }
        support.rows.push_back(std::move(evidence));
    }
    for(auto id:donor_ids){const V column=transform.col(id);
        support.donor+=(ligand_to_centre?ql.transpose()*column:
                        qc.transpose()*column).squaredNorm();}
    for(auto id:acceptor_ids){const V column=transform.col(id);
        support.acceptor+=(ligand_to_centre?qc.transpose()*column:
                           ql.transpose()*column).squaredNorm();}
    support.donor/=std::max<Eigen::Index>(1,qc.cols());
    support.acceptor/=std::max<Eigen::Index>(1,qc.cols());
    return support;
}
void append_group(NboPiCoupling& out,const Wavefunction& w,const M& u,const M& f,
                  const M& qc,const M& ql,const M& pc,const M& pl,
                  const V* occupations,const std::vector<std::size_t>& members){
    M c(u.rows(),members.size());double occ=0;
    for(std::size_t j=0;j<members.size();++j){
        const auto& mo=w.orbitals[members[j]];
        if(mo.source_orbital_index>=static_cast<std::size_t>(u.cols()))return;
        c.col(j)=u.col(mo.source_orbital_index);
        if(occupations)occ+=(*occupations)[mo.source_orbital_index];
    }
    NboPiCanonicalGroup group;
    group.members=members;group.spin=out.spin;
    if(occupations)group.occupation_per_mo=occ/members.size();
    group.centre_weight=(qc.transpose()*c).squaredNorm()/members.size();
    group.ligand_weight=(ql.transpose()*c).squaredNorm()/members.size();
    if(group.centre_weight<0.01||group.ligand_weight<0.01)return;
    M cross=c.transpose()*(pc*f*pl+pl*f*pc)*c;
    cross=(cross+cross.transpose()).eval()*0.5;
    Eigen::SelfAdjointEigenSolver<M> eig(cross,Eigen::EigenvaluesOnly);
    if(eig.info()!=Eigen::Success)return;
    group.cross_fock_min_hartree=eig.eigenvalues()[0];
    group.cross_fock_max_hartree=eig.eigenvalues()[eig.eigenvalues().size()-1];
    group.cross_fock_mean_hartree=cross.trace()/members.size();
    const double gate=std::max(1e-5,10*u.rows()*out.operator_max_error_hartree);
    if(group.cross_fock_max_hartree< -gate)group.character="bonding_mixing";
    else if(group.cross_fock_min_hartree>gate)group.character="antibonding_mixing";
    out.groups.push_back(std::move(group));
}
} // namespace

NboPiCouplingAnalysis analyse_nbo_pi_couplings(const Wavefunction& w,
        const NboIntegration& data,
        const std::vector<std::pair<std::size_t,std::size_t>>& strong_connectivity){
    NboPiCouplingAnalysis out;
    std::vector<std::string> projection_failures;
    std::optional<NboSalcModel> physical_operator_model;
    if(!data.dataset.association.compatible||data.canonical_fingerprint!=nbo_canonical_fingerprint(w)){
        out.status="rejected";out.reason="NBO association does not match immutable canonical identity";
        return out;
    }
    const auto n=static_cast<std::size_t>(w.basis_count);
    if(!n||!data.dataset.archive){out.reason="Complete same-source archive is unavailable";return out;}
    const auto* s=unique_matrix(data.dataset.archive->matrices,"OVERLAP",NboSpin::Total);
    if(!full(s,n)){out.reason="Complete archive AO overlap is unavailable";return out;}
    const auto ro=verify_nbo_restricted_open_shell(w,data);
    const bool shared_ro=ro.verified&&
        (w.orbital_occupation_model==OrbitalOccupationModel::CanonicalShared||
         w.orbital_occupation_model==OrbitalOccupationModel::SharedIntegerDeterminant);
    for(auto spin:{NboSpin::Total,NboSpin::Alpha,NboSpin::Beta}){
        const auto* a=unique_matrix(data.dataset.matrices,"AONAO",spin);
        if(!a&&spin!=NboSpin::Total)
            a=unique_matrix(data.dataset.matrices,"AONAO",NboSpin::Total);
        const auto* u=unique_matrix(data.dataset.matrices,"NAOMO",spin);
        const auto* f=unique_matrix(data.dataset.archive->matrices,"FOCK",spin);
        if(!full(a,n)||!full(u,n)||!full(f,n))continue;
        const NboNaoValidation* validation=nullptr;
        for(const auto& candidate:data.dataset.nao_validation)
            if(candidate.spin==spin&&candidate.available&&
               (candidate.direct_fchk_coefficients||(shared_ro&&spin==NboSpin::Beta)))
                validation=&candidate;
        if(!validation||validation->effective_mo_columns!=n)continue;
        const M an=mat(*a),un=mat(*u),fn=an.transpose()*mat(*f)*an;
        const double orth=maxabs(an.transpose()*mat(*s)*an-M::Identity(n,n));
        if(orth>5e-5||maxabs(un.transpose()*un-M::Identity(n,n))>5e-5)continue;
        V energies(n),occupations(n);std::vector<bool> seen(n,false);
        bool occupation_provenance=true;
        for(const auto& mo:w.orbitals){
            if(!shared_ro&&(spin==NboSpin::Beta)!=(mo.spin==Spin::Beta))continue;
            if(shared_ro&&mo.spin!=Spin::Alpha)continue;
            const auto index=mo.source_orbital_index;
            if(index>=n||seen[index])continue;
            seen[index]=true;energies[index]=mo.energy_hartree;
            occupation_provenance=occupation_provenance&&
                (mo.occupation_provenance!=DataProvenance::Unavailable ||
                 w.electron_counts_provenance!=DataProvenance::Unavailable);
        }
        if(!std::all_of(seen.begin(),seen.end(),[](bool yes){return yes;}))continue;
        const double canonical_operator_residual=maxabs(fn-un*energies.asDiagonal()*un.transpose());
        double operator_error=canonical_operator_residual;
        std::string operator_kind="canonical-same-operator";
        double operator_validation_tolerance=2e-5;
        if(operator_error>2e-5) {
            if(!physical_operator_model)physical_operator_model=build_nbo_salc_model(w,data);
            const auto physical=std::find_if(physical_operator_model->energies.begin(),
                physical_operator_model->energies.end(),[&](const auto& e){
                    return e.spin==spin&&e.available&&e.printed_operator_verified;});
            const auto op=std::find_if(physical_operator_model->spin_operators.begin(),
                physical_operator_model->spin_operators.end(),[&](const auto& e){return e.spin==spin;});
            if(physical==physical_operator_model->energies.end()||
               op==physical_operator_model->spin_operators.end()||op->fock.size()!=n*n||
               op->basis.size()!=n) {
                projection_failures.push_back(std::string(nbo_spin_name(spin))+": canonical operator differs and independent printed physical-spin Fock validation is unavailable");
                continue;
            }
            bool same_basis=true;for(std::size_t i=0;i<n;++i)
                same_basis=same_basis&&op->basis[i].kind==NboOrbitalKind::NAO&&
                    op->basis[i].spin==spin&&op->basis[i].index==i;
            const M verified=Eigen::Map<const RM>(op->fock.data(),n,n);
            if(!same_basis||maxabs(fn-verified)>2e-5)continue;
            operator_error=std::max(maxabs(fn-verified),physical->hermiticity_error);
            operator_kind="independently-verified-physical-spin-fock";
            operator_validation_tolerance=5.1e-6;
        }
        bool density_validated=occupation_provenance &&
            validation->occupation_density_verified &&
            validation->canonical_occupations.size()==n &&
            std::all_of(validation->canonical_occupations.begin(),
                validation->canonical_occupations.end(),[](double x){return std::isfinite(x);});
        M density;
        if(density_validated){
            for(std::size_t col=0;col<n;++col)
                occupations[col]=validation->canonical_occupations[col];
            density=un*occupations.asDiagonal()*un.transpose();
        }
        if(!density_validated&&shared_ro&&physical_operator_model) {
            const auto op=std::find_if(physical_operator_model->spin_operators.begin(),
                physical_operator_model->spin_operators.end(),[&](const auto& e){return e.spin==spin;});
            if(op!=physical_operator_model->spin_operators.end()&&op->density.size()==n*n) {
                const std::size_t electrons=spin==NboSpin::Beta?w.beta_electrons:w.alpha_electrons;
                for(std::size_t i=0;i<n;++i)occupations[i]=i<electrons?1.0:0.0;
                const M reconstructed=un*occupations.asDiagonal()*un.transpose();
                const M independently_verified=Eigen::Map<const RM>(op->density.data(),n,n);
                if(maxabs(reconstructed-independently_verified)<=2e-5) {
                    density=independently_verified;density_validated=true;
                }
            }
        }
        std::optional<M> nbo_transform;
        const auto* naob=unique_matrix(data.dataset.matrices,"NAONBO",spin);
        const auto* aob=unique_matrix(data.dataset.matrices,"AONBO",spin);
        const auto* nbo_cap=nbo_capability(data,"nbo");
        if(nbo_cap&&nbo_cap->available()&&full(naob,n)&&full(aob,n)&&
           maxabs(an*mat(*naob)-mat(*aob))<=2e-5)
            nbo_transform=mat(*naob);
        const auto family=shells(data,spin);
        std::set<std::pair<std::size_t,std::size_t>> bonded;
        for(const auto& edge:strong_connectivity)bonded.emplace(std::minmax(edge.first,edge.second));
        struct Candidate{std::vector<std::size_t> centre,ligand,ci,li;M qc,ql;
            std::string ligand_family;std::vector<std::size_t> family_nbo_ids;
            std::vector<NboPiAngularEvidence> angular;};
        std::vector<Candidate> candidates;
        const auto add_angular=[&](Candidate& c,std::size_t atom,const PiFrame& frame){
            NboPiAngularEvidence x;x.atom=atom;
            x.p_metric_min_eigenvalue=frame.minimum_p;
            x.max_sigma_leakage=frame.leakage;
            x.centre_projection=frame.centre;
            x.partition_residual=frame.residual;
            x.pi_generalized_eigenvalues=frame.spectrum;
            c.angular.push_back(x);
        };
        for(const auto& [atom,d_ids]:family.d){
            Candidate candidate;candidate.centre={atom};candidate.ci=d_ids;
            candidate.qc=d_columns(n,d_ids);
            std::vector<M> parts;bool valid=true;
            for(const auto& [other,p_ids]:family.p){
                if(!bonded.count(std::minmax(atom,other))||atom==other)continue;
                const V3 axis=bond_axis(w,atom,other);if(axis.norm()<0.5)continue;
                std::string reason;
                auto frame=derived_pi_frame(w,data,other,p_ids,spin,axis,reason);
                if(!frame){projection_failures.push_back(std::string(nbo_spin_name(spin))+" atom "+std::to_string(other+1)+": "+reason);valid=false;break;}
                candidate.ligand.push_back(other);
                candidate.li.insert(candidate.li.end(),p_ids.begin(),p_ids.end());
                parts.push_back(frame->q);add_angular(candidate,other,*frame);
            }
            if(parts.empty()||!valid)continue;
            candidate.ql=M::Zero(n,2*parts.size());
            for(std::size_t k=0;k<parts.size();++k)
                candidate.ql.middleCols(2*k,2)=parts[k];
            candidates.push_back(std::move(candidate));
        }
        // Resolve the complete ligand's internal occupied pi and antibonding
        // pi families before reducing the d--ligand coupling. Atom-pi alone
        // combines these physically distinct spaces (e.g. both ends of any
        // heteronuclear multiple bond). This uses no ligand-name catalogue or
        // E2 printing threshold, and never assumes BD ordinal 1 means sigma.
        if(nbo_transform) {
            const std::size_t broad_count=candidates.size();
            for(std::size_t cindex=0;cindex<broad_count;++cindex) {
                const auto broad=candidates[cindex];
                std::set<std::size_t> fragment(broad.ligand.begin(),broad.ligand.end());
                bool changed=true;
                while(changed) {
                    changed=false;
                    for(const auto& [a,b]:bonded) {
                        if(family.d.count(a)||family.d.count(b))continue;
                        if(fragment.count(a)&&fragment.insert(b).second)changed=true;
                        if(fragment.count(b)&&fragment.insert(a).second)changed=true;
                    }
                }
                std::set<std::size_t> fragment_p_rows;
                for(const auto& row:data.dataset.naos)
                    if(row.spin==spin&&row.id&&row.atom&&fragment.count(row.atom-1)&&valence(row,'p'))
                        fragment_p_rows.insert(row.id-1);
                for(const bool antibonding:{false,true}) {
                    std::vector<V> columns;std::vector<std::size_t> ids;
                    for(const auto& local:data.dataset.orbitals) {
                        if(local.spin!=spin||!local.id||local.id>static_cast<std::size_t>(nbo_transform->cols())||
                           local.kind!=(antibonding?"BD*":"BD"))continue;
                        const auto* desc=nbo_orbital(data,{NboOrbitalKind::NBO,spin,local.id-1});
                        if(!desc||desc->atoms.size()<2||
                           !std::all_of(desc->atoms.begin(),desc->atoms.end(),[&](auto atom){return fragment.count(atom);})||
                           !std::any_of(desc->channels.begin(),desc->channels.end(),[](const auto& evidence){
                               return evidence.status=="available"&&evidence.channel=="pi";}))continue;
                        const V original=nbo_transform->col(local.id-1);
                        double direct_p=0;for(auto id:broad.li)direct_p+=original[id-1]*original[id-1];
                        const double transverse=(broad.ql.transpose()*original).squaredNorm();
                        // Establish that this internal pi family belongs to the
                        // metal-facing pi geometry, not a sigma-facing orbital.
                        if(direct_p<1e-4||transverse/direct_p<0.9)continue;
                        V projected=V::Zero(n);for(auto row:fragment_p_rows)projected[row]=original[row];
                        if(projected.squaredNorm()<0.7)continue;
                        columns.push_back(std::move(projected));ids.push_back(local.id);
                    }
                    if(columns.empty())continue;
                    M raw(n,columns.size());for(std::size_t i=0;i<columns.size();++i)raw.col(i)=columns[i];
                    Eigen::JacobiSVD<M> basis(raw,Eigen::ComputeThinU);
                    const double cutoff=std::max(1e-8,1e-8*basis.singularValues()[0]);
                    std::size_t rank=0;for(double value:basis.singularValues())if(value>cutoff)++rank;
                    if(!rank||rank!=columns.size())continue; // incomplete dependent family is not silently truncated
                    Candidate split=broad;split.ligand.assign(fragment.begin(),fragment.end());
                    split.li.clear();for(auto row:fragment_p_rows)split.li.push_back(row+1);
                    split.ql=basis.matrixU().leftCols(rank);
                    split.ligand_family=antibonding?"internal-pi-antibonding":"internal-pi-bonding";
                    split.family_nbo_ids=std::move(ids);
                    candidates.push_back(std::move(split));
                }
            }
        }
        for(const auto& [atom,other]:bonded){
            const auto ca=family.p.find(atom),lb=family.p.find(other);
            if(ca==family.p.end()||lb==family.p.end())continue;
            const V3 axis=bond_axis(w,atom,other);if(axis.norm()<0.5)continue;
            std::string centre_reason,ligand_reason;
            auto centre_frame=derived_pi_frame(w,data,atom,ca->second,spin,axis,centre_reason);
            auto ligand_frame=derived_pi_frame(w,data,other,lb->second,spin,axis,ligand_reason);
            if(!centre_frame||!ligand_frame){
                projection_failures.push_back(std::string(nbo_spin_name(spin))+" atoms "+
                    std::to_string(atom+1)+"-"+std::to_string(other+1)+": "+centre_reason+" "+ligand_reason);
                continue;
            }
            Candidate candidate;candidate.centre={atom};candidate.ligand={other};
            candidate.ci.assign(ca->second.begin(),ca->second.end());
            candidate.li.assign(lb->second.begin(),lb->second.end());
            candidate.qc=centre_frame->q;
            candidate.ql=ligand_frame->q;
            add_angular(candidate,atom,*centre_frame);
            add_angular(candidate,other,*ligand_frame);
            candidates.push_back(std::move(candidate));
        }
        for(const auto& candidate:candidates){
            if(maxabs(candidate.qc.transpose()*candidate.ql)>1e-8)continue;
            const M block=candidate.qc.transpose()*fn*candidate.ql;
            Eigen::JacobiSVD<M> svd(block,Eigen::ComputeThinU|Eigen::ComputeThinV);
            std::size_t rank=0;for(auto value:svd.singularValues())if(value>1e-6)++rank;
            if(!rank)continue;
            NboPiCoupling record;
            record.spin=spin;record.centre_atoms=candidate.centre;
            record.ligand_family=candidate.ligand_family;
            record.ligand_family_nbo_ids=candidate.family_nbo_ids;
            record.localized_family_verified=!candidate.ligand_family.empty();
            record.ligand_atoms=candidate.ligand;
            record.centre_nao_ids=candidate.ci;record.ligand_nao_ids=candidate.li;
            record.coupled_rank=rank;record.operator_max_error_hartree=operator_error;
            record.operator_kind=operator_kind;
            record.canonical_operator_residual_hartree=canonical_operator_residual;
            record.operator_validation_tolerance_hartree=operator_validation_tolerance;
            record.canonical_members_are_verified_shared_spatial=shared_ro;
            record.nao_orthogonality_error=orth;record.source=f->source;
            record.angular_projector_evidence=candidate.angular;
            record.minimum_centre_projection=1;
            for(const auto& evidence:candidate.angular){
                record.angular_leakage=std::max(record.angular_leakage,evidence.max_sigma_leakage);
                record.minimum_centre_projection=std::min(record.minimum_centre_projection,
                    evidence.centre_projection);
            }
            for(auto value:svd.singularValues())
                record.singular_values_hartree.push_back(value);
            std::ostringstream id;id<<"pi:"<<nbo_spin_name(spin)<<":";
            for(auto atom:candidate.centre)id<<atom<<',';id<<":";
            for(auto atom:candidate.ligand)id<<atom<<',';
            id<<":"<<candidate.ligand_family;record.id=id.str();
            const M qc=candidate.qc*svd.matrixU().leftCols(rank);
            const M ql=candidate.ql*svd.matrixV().leftCols(rank);
            const M pc=qc*qc.transpose(),pl=ql*ql.transpose();
            record.centre_onsite_hartree=(qc.transpose()*fn*qc).trace()/rank;
            record.ligand_onsite_hartree=(ql.transpose()*fn*ql).trace()/rank;
            const auto range=[&](const M& block){
                Eigen::SelfAdjointEigenSolver<M> solver((block+block.transpose())*0.5,
                    Eigen::EigenvaluesOnly);
                if(solver.info()!=Eigen::Success||!solver.eigenvalues().allFinite())
                    return std::optional<std::array<double,2>>{};
                return std::optional<std::array<double,2>>(
                    std::array<double,2>{solver.eigenvalues()[0],
                        solver.eigenvalues()[solver.eigenvalues().size()-1]});
            };
            const auto centre_fock_range=range(qc.transpose()*fn*qc),
                       ligand_fock_range=range(ql.transpose()*fn*ql);
            if(!centre_fock_range||!ligand_fock_range)continue;
            record.centre_onsite_range_hartree=*centre_fock_range;
            record.ligand_onsite_range_hartree=*ligand_fock_range;
            if(density_validated){
                const auto centre_density=qc.transpose()*density*qc,
                           ligand_density=ql.transpose()*density*ql;
                record.centre_occupation=centre_density.trace()/rank;
                record.ligand_occupation=ligand_density.trace()/rank;
                record.centre_occupation_range=range(centre_density);
                record.ligand_occupation_range=range(ligand_density);
                if(!record.centre_occupation_range||!record.ligand_occupation_range)
                    continue;
                record.occupation_status="available";
                record.occupation_reason="Validated same-spin canonical occupations reconstruct NAO density";
            }else record.occupation_reason=
                "Occupation provenance or independent canonical-to-NAO density validation unavailable; coupling retained without direction";
            for(const auto& group:canonical_groups(w,shared_ro?NboSpin::Alpha:spin))
                append_group(record,w,un,fn,qc,ql,pc,pl,
                    density_validated?&occupations:nullptr,group);
            if(record.groups.empty())continue;
            const double operator_gate=std::max(1e-5,10*n*operator_error);
            const bool ligand_lower=record.centre_onsite_range_hartree[0]>
                record.ligand_onsite_range_hartree[1]+operator_gate;
            const bool centre_lower=record.ligand_onsite_range_hartree[0]>
                record.centre_onsite_range_hartree[1]+operator_gate;
            const bool ligand_more=record.ligand_occupation_range&&
                record.centre_occupation_range&&
                (*record.ligand_occupation_range)[0]>
                (*record.centre_occupation_range)[1]+0.1;
            const bool centre_more=record.centre_occupation_range&&
                record.ligand_occupation_range&&
                (*record.centre_occupation_range)[0]>
                (*record.ligand_occupation_range)[1]+0.1;
            E2Support forward,reverse;
            if(nbo_transform){
                forward=directed_e2(data,*nbo_transform,qc,ql,candidate.centre,
                                    candidate.ligand,spin,true);
                reverse=directed_e2(data,*nbo_transform,qc,ql,candidate.centre,
                                    candidate.ligand,spin,false);
                record.direction_projection_evidence=forward.rows;
                record.direction_projection_evidence.insert(
                    record.direction_projection_evidence.end(),
                    reverse.rows.begin(),reverse.rows.end());
            }
            const bool symmetric_onsite=
                std::abs(record.centre_onsite_range_hartree[0]-
                         record.ligand_onsite_range_hartree[0])<1e-4 &&
                std::abs(record.centre_onsite_range_hartree[1]-
                         record.ligand_onsite_range_hartree[1])<1e-4;
            const bool symmetric_occupation=record.centre_occupation_range&&
                record.ligand_occupation_range&&
                std::abs((*record.centre_occupation_range)[0]-
                         (*record.ligand_occupation_range)[0])<0.05 &&
                std::abs((*record.centre_occupation_range)[1]-
                         (*record.ligand_occupation_range)[1])<0.05;
            if(candidate.ligand_family=="internal-pi-antibonding"&&centre_lower&&centre_more){
                record.direction="centre_to_ligand";
                record.direction_evidence="Verified full-ligand internal pi-antibonding family, same-operator Fock coupling and concordant occupied centre; independent of E2 printing threshold";
            }else if(candidate.ligand_family=="internal-pi-bonding"&&ligand_lower&&ligand_more){
                record.direction="ligand_to_centre";
                record.direction_evidence="Verified full-ligand internal pi-bonding family, same-operator Fock coupling and concordant occupied ligand; independent of E2 printing threshold";
            }else if(symmetric_onsite&&symmetric_occupation){
                record.direction="symmetric_coupled";
                record.direction_evidence="Coupled-projector on-site and occupation ranges are symmetric; no directional donor assignment";
            }else if(ligand_lower&&ligand_more&&
                     forward.donor>=0.7&&forward.acceptor>=0.7&&
                     (reverse.donor<0.3||reverse.acceptor<0.3)){
                record.direction="ligand_to_centre";
                record.direction_donor_weight=forward.donor;
                record.direction_acceptor_weight=forward.acceptor;
                record.direction_evidence="Lower occupied ligand coupled projector and direct same-channel localized E2 donor/acceptor spans";
            }else if(centre_lower&&centre_more&&
                     reverse.donor>=0.7&&reverse.acceptor>=0.7&&
                     (forward.donor<0.3||forward.acceptor<0.3)){
                record.direction="centre_to_ligand";
                record.direction_donor_weight=reverse.donor;
                record.direction_acceptor_weight=reverse.acceptor;
                record.direction_evidence="Lower occupied centre coupled projector and direct same-channel localized E2 donor/acceptor spans";
            }else record.direction_evidence=density_validated?
                "Reciprocal Fock coupling lacks concordant occupied/onsite ranges and mapped E2 direction evidence":
                "Reciprocal Fock coupling verified; occupation/density gate prevents directional inference";
            out.couplings.push_back(std::move(record));
        }
        out.evidence.push_back(f->source);
    }
    if(out.couplings.empty()){
        out.status=out.evidence.empty()||!projection_failures.empty()?"insufficient_evidence":"not_applicable";
        out.reason=out.evidence.empty()?"No spin has a complete matched NAO transform and independently consistent same-source Fock operator":
            !projection_failures.empty()?"Same-source Fock verified; local pi projection evidence is incomplete":
            "No supported bonded valence p-pi or d-to-ligand-p-pi coupling block";
    }else{
        out.status="available";
        out.reason="Independently verified local-operator NAO projectors; each channel declares canonical or physical-spin operator identity; canonical groups and eigenvalue intervals remain separate";
    }
    if(!projection_failures.empty()){
        out.reason+="; candidate pi subspaces could not be verified";
        for(const auto& failure:projection_failures)out.reason+="; "+failure;
    }
    return out;
}
} // namespace cov
