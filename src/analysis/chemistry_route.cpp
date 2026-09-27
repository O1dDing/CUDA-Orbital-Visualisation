#include "cov/chemistry_route.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <tuple>

namespace cov { namespace {

std::string quoted(const std::string& value) {
    std::ostringstream out;out<<'"';
    for(unsigned char c:value) {
        if(c=='"'||c=='\\')out<<'\\'<<c;
        else if(c=='\n')out<<"\\n";
        else if(c=='\r')out<<"\\r";
        else if(c=='\t')out<<"\\t";
        else if(c<0x20)out<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<int(c)<<std::dec;
        else out<<c;
    }
    return out.str()+'"';
}

RoutedStatus unavailable_state(const NboCapability* capability,bool attached) {
    if(!attached)return RoutedStatus::NotAnalysed;
    if(!capability)return RoutedStatus::Insufficient;
    switch(capability->state) {
        case NboCapabilityState::Rejected: return RoutedStatus::Rejected;
        case NboCapabilityState::Unsupported: return RoutedStatus::Unsupported;
        case NboCapabilityState::Missing:
        case NboCapabilityState::Ambiguous: return RoutedStatus::Insufficient;
        case NboCapabilityState::Available: return RoutedStatus::Insufficient;
    }
    return RoutedStatus::Insufficient;
}

template<class T> void unavailable(RoutedResult<T>& result,const NboCapability* cap,
                                    bool attached,const char* method) {
    result.status=unavailable_state(cap,attached);
    result.method=method;
    result.reason=cap?cap->detail:attached?"Required validated evidence is absent":
        "NBO analysis has not been attached";
    if(cap)result.evidence=cap->sources;
}

template<class T> void fallback(RoutedResult<T>& result,const std::string& reason) {
    result.fallback_reason=reason;
    result.status=RoutedStatus::Available;
    result.provider=RoutedProvider::Legacy;
}

std::optional<std::pair<std::vector<double>,std::vector<NboSource>>>
validated_atoms(const NboIntegration& integration,std::size_t atom_count,
                const char* kind,NboSpin spin) {
    std::vector<double> values(atom_count,0);
    std::vector<NboSource> sources;
    std::vector<bool> seen(atom_count,false);
    for(const auto& row:integration.structure) {
        if(row.kind!=kind || row.spin!=spin || row.atoms.size()!=1 ||
           !row.value || !std::isfinite(*row.value) || row.atoms.front()>=atom_count)
            continue;
        const auto atom=row.atoms.front();
        if(seen[atom])return std::nullopt;
        seen[atom]=true;values[atom]=*row.value;sources.push_back(row.source);
    }
    if(!std::all_of(seen.begin(),seen.end(),[](bool yes){return yes;}))
        return std::nullopt;
    return std::make_pair(std::move(values),std::move(sources));
}

RoutedMoComposition legacy_mo(const MolecularOrbital& mo,std::size_t index) {
    RoutedMoComposition out;out.canonical_index=index;
    out.semantics="legacy AO-reference projection; basis-dependent and not an NBO population";
    std::map<std::size_t,NboNaoGroupContribution> atoms;
    for(const auto& contribution:mo.chemistry.ao_contributions) {
        NboNaoGroupContribution shell;
        shell.atom=static_cast<std::size_t>(contribution.atom_index)+1;
        shell.label=contribution.label;
        shell.principal_n=contribution.principal_n;
        shell.angular_l=contribution.angular_momentum;
        shell.weight=contribution.weight;
        out.shells.push_back(shell);
        auto& atom=atoms[shell.atom];atom.atom=shell.atom;
        atom.label="atom "+std::to_string(shell.atom);
        atom.weight+=contribution.weight;
    }
    for(auto& [_,atom]:atoms)out.atoms.push_back(std::move(atom));
    out.legacy_coverage=std::clamp(1.0-mo.chemistry.unresolved_weight,0.0,1.0);
    out.legacy_unresolved_fraction=mo.chemistry.unresolved_weight;
    out.complete=false; // Legacy AO-reference coverage is not an S-metric reconstruction proof.
    return out;
}

std::optional<double> nao_reconstruction_error(const Wavefunction& canonical,
                                                const NboIntegration& integration,
                                                std::size_t mo_index) {
    if(mo_index>=canonical.orbitals.size())return std::nullopt;
    const auto n=static_cast<std::size_t>(canonical.basis_count);
    if(!n || canonical.ao_overlap.size()!=n*n ||
       canonical.orbitals[mo_index].coefficients.size()!=n)return std::nullopt;
    std::vector<double> residual=canonical.orbitals[mo_index].coefficients;
    std::size_t terms=0;
    for(const auto& link:integration.links) {
        if(link.canonical_index!=mo_index || link.orbital.kind!=NboOrbitalKind::NAO)
            continue;
        const auto* orbital=nbo_orbital(integration,link.orbital);
        if(!orbital || orbital->coefficients.size()!=n ||
           !std::isfinite(link.coefficient))return std::nullopt;
        for(std::size_t row=0;row<n;++row)
            residual[row]-=link.coefficient*orbital->coefficients[row];
        ++terms;
    }
    if(!terms)return std::nullopt;
    double squared=0;
    for(std::size_t row=0;row<n;++row) {
        double overlap_row=0;
        for(std::size_t col=0;col<n;++col)
            overlap_row+=canonical.ao_overlap[row*n+col]*residual[col];
        squared+=residual[row]*overlap_row;
    }
    if(!std::isfinite(squared) || squared< -1e-7)return std::nullopt;
    return std::sqrt(std::max(0.0,squared));
}

void result_header(std::ostringstream& out,const char* key,RoutedStatus status,
                   RoutedProvider provider,const std::string& method,
                   const std::string& reason,const std::string& fallback) {
    out<<quoted(key)<<":{\"status\":"<<quoted(routed_status_name(status))
       <<",\"provider\":"<<quoted(routed_provider_name(provider))
       <<",\"method\":"<<quoted(method)
       <<",\"reason\":"<<quoted(reason)
       <<",\"fallback_reason\":"<<quoted(fallback);
}

void source_json(std::ostringstream& out,const NboSource& source) {
    out<<"{\"path\":"<<quoted(source.path)<<",\"block\":"<<quoted(source.block)
       <<",\"line_begin\":"<<source.line_begin
       <<",\"line_end\":"<<source.line_end
       <<",\"analysis_segment\":"<<source.analysis_segment<<'}';
}

template<class T> void sources_json(std::ostringstream& out,const RoutedResult<T>& result) {
    out<<",\"evidence\":[";
    for(std::size_t i=0;i<result.evidence.size();++i) {
        if(i)out<<',';source_json(out,result.evidence[i]);
    }
    out<<']';
}

} // namespace

const char* routed_status_name(RoutedStatus s) noexcept {
    switch(s) {
        case RoutedStatus::Available:return "available";
        case RoutedStatus::NotAnalysed:return "not_analysed";
        case RoutedStatus::Insufficient:return "insufficient_evidence";
        case RoutedStatus::Rejected:return "rejected";
        case RoutedStatus::Unsupported:return "unsupported";
        case RoutedStatus::NotApplicable:return "not_applicable";
        case RoutedStatus::NotReportedAboveThreshold:return "not_reported_above_threshold";
    }
    return "insufficient_evidence";
}
const char* routed_provider_name(RoutedProvider p) noexcept {
    switch(p) {
        case RoutedProvider::None:return "none";
        case RoutedProvider::Legacy:return "legacy";
        case RoutedProvider::Nbo:return "nbo";
    }
    return "none";
}
const char* active_orbital_kind_name(ActiveOrbitalKind k) noexcept {
    switch(k) {
        case ActiveOrbitalKind::Canonical:return "canonical";
        case ActiveOrbitalKind::NboSet:return "nbo";
        case ActiveOrbitalKind::Inspection:return "inspection";
    }
    return "canonical";
}

RoutedAnalysis route_chemistry(const Wavefunction& canonical,
                                const NboIntegration* integration) {
    RoutedAnalysis out;out.canonical_fingerprint=nbo_canonical_fingerprint(canonical);
    if(integration)out.integration_id=integration->id;
    const bool associated=integration &&
        integration->canonical_fingerprint==out.canonical_fingerprint &&
        integration->dataset.association.compatible;
    const auto* charge_cap=associated?nbo_capability(*integration,"charges"):nullptr;
    const auto* spin_cap=associated?nbo_capability(*integration,"spin"):nullptr;
    const auto* aomo_cap=associated?nbo_capability(*integration,"aomo"):nullptr;
    unavailable(out.total_atomic_charge,charge_cap,integration!=nullptr,"NPA total charge");
    unavailable(out.atomic_spin,spin_cap,integration!=nullptr,"NPA alpha-minus-beta spin");
    if(integration && !associated) {
        const std::string reason="NBO association does not match this canonical wavefunction";
        out.total_atomic_charge.status=out.atomic_spin.status=RoutedStatus::Rejected;
        out.total_atomic_charge.reason=out.atomic_spin.reason=reason;
    }
    if(charge_cap && charge_cap->available()) {
        if(auto atoms=validated_atoms(*integration,canonical.atoms.size(),
                                      "atomic_charge",NboSpin::Total)) {
            out.total_atomic_charge.value=std::move(atoms->first);
            out.total_atomic_charge.evidence=std::move(atoms->second);
            out.total_atomic_charge.provider=RoutedProvider::Nbo;
            out.total_atomic_charge.status=RoutedStatus::Available;
            out.total_atomic_charge.reason="Complete independently validated total NPA charges";
        }else{
            out.total_atomic_charge.status=RoutedStatus::Insufficient;
            out.total_atomic_charge.reason="Validated total NPA charge is incomplete or ambiguous per atom";
        }
    }
    if(spin_cap && spin_cap->available()) {
        if(auto atoms=validated_atoms(*integration,canonical.atoms.size(),
                                      "atomic_spin",NboSpin::Total)) {
            out.atomic_spin.value=std::move(atoms->first);
            out.atomic_spin.evidence=std::move(atoms->second);
            out.atomic_spin.provider=RoutedProvider::Nbo;
            out.atomic_spin.status=RoutedStatus::Available;
            out.atomic_spin.reason="Complete independently validated NPA spin population";
        }else{
            out.atomic_spin.status=RoutedStatus::Insufficient;
            out.atomic_spin.reason="Validated spin capability lacks a complete per-atom value set";
        }
    }
    if(!out.total_atomic_charge.available() &&
       canonical.atomic_partial_charge_provenance!=DataProvenance::Unavailable &&
       canonical.atomic_partial_charges.size()==canonical.atoms.size() &&
       std::all_of(canonical.atomic_partial_charges.begin(),canonical.atomic_partial_charges.end(),
                   [](double x){return std::isfinite(x);})) {
        const auto why=out.total_atomic_charge.reason;
        out.total_atomic_charge.value=canonical.atomic_partial_charges;
        out.total_atomic_charge.method=canonical.atomic_partial_charge_scheme.empty()?
            "legacy atomic partial charge":canonical.atomic_partial_charge_scheme;
        fallback(out.total_atomic_charge,why);
        out.total_atomic_charge.reason="Producer atom-resolved partial charges; not NPA";
    }
    out.mo_composition.resize(canonical.orbitals.size());
    out.mo_relations.resize(canonical.orbitals.size());
    using LinkKey=std::tuple<std::size_t,int,int,std::size_t>;
    std::map<LinkKey,double> link_weights;
    if(associated){
        for(const auto& link:integration->links)if(link.weight)
            link_weights[{link.canonical_index,static_cast<int>(link.orbital.kind),
                static_cast<int>(link.orbital.spin),link.orbital.index}]=*link.weight;
        for(const auto& local:integration->structure){
            if(local.kind!="donor_acceptor" && local.kind!="localization_relation" &&
               local.kind!="multicentre")continue;
            RoutedLocalRelation relation;
            relation.id=local.id;relation.kind=local.kind;relation.label=local.label;
            relation.spin=local.spin;relation.atoms=local.atoms;
            relation.orbitals=local.orbitals;relation.printed_value=local.value;
            relation.units=local.units;relation.source=local.source;
            relation.semantics=local.kind=="donor_acceptor"?
                "Directed localized NBO E(2); independent donor and acceptor canonical projections, not allocated E(2)":
                local.kind=="localization_relation"?
                "Validated NLMO-parent relation and canonical projections; no inferred donor role":
                "Producer multicentre localized record related to canonical subspace; not a canonical relabel";
            out.local_relation_registry.push_back(std::move(relation));
        }
    }
    // Invert the actual source-orbital transforms once. Dense reports can
    // contain many E(2) rows; scanning every row for every canonical MO made
    // the route and its exported snapshot scale with all MO×relation pairs.
    using RefKey=std::tuple<int,int,std::size_t>;
    std::map<RefKey,std::vector<std::pair<std::size_t,std::size_t>>> relation_refs;
    for(std::size_t relation_index=0;relation_index<out.local_relation_registry.size();++relation_index)
        for(std::size_t component=0;
            component<out.local_relation_registry[relation_index].orbitals.size();++component){
            const auto& ref=out.local_relation_registry[relation_index].orbitals[component];
            relation_refs[{static_cast<int>(ref.kind),static_cast<int>(ref.spin),ref.index}]
                .push_back({relation_index,component});
        }
    std::vector<std::map<std::size_t,RoutedRelationProjection>> relation_hits(canonical.orbitals.size());
    for(const auto& [key,weight]:link_weights){
        const auto& [mo,kind,spin,index]=key;
        if(mo>=relation_hits.size() || weight<1e-4)continue;
        const auto found=relation_refs.find({kind,spin,index});
        if(found==relation_refs.end())continue;
        for(const auto& [relation_index,component]:found->second){
            auto& projection=relation_hits[mo][relation_index];
            if(projection.canonical_projection_weights.empty()){
                const auto& relation=out.local_relation_registry[relation_index];
                projection.relation_id=relation.id;
                projection.canonical_projection_weights.resize(relation.orbitals.size());
            }
            projection.canonical_projection_weights[component]=weight;
        }
    }
    // Once a relation is genuinely supported, retain its weaker measured
    // endpoint projections too; the display cutoff controls inclusion only.
    for(const auto& [key,weight]:link_weights){
        if(weight>=1e-4)continue;
        const auto& [mo,kind,spin,index]=key;
        if(mo>=relation_hits.size())continue;
        const auto found=relation_refs.find({kind,spin,index});
        if(found==relation_refs.end())continue;
        for(const auto& [relation_index,component]:found->second){
            const auto hit=relation_hits[mo].find(relation_index);
            if(hit!=relation_hits[mo].end())
                hit->second.canonical_projection_weights[component]=weight;
        }
    }
    for(std::size_t i=0;i<canonical.orbitals.size();++i) {
        auto& result=out.mo_composition[i];
        unavailable(result,aomo_cap,integration!=nullptr,"validated NAO projection");
        if(integration && !associated) {
            result.status=RoutedStatus::Rejected;
            result.reason="NBO association does not match this canonical wavefunction";
        }
        if(aomo_cap && aomo_cap->available()) {
            if(const auto* decomposition=nbo_mo_decomposition(integration->dataset,i);
               decomposition && decomposition->available) {
                RoutedMoComposition composition;
                composition.canonical_index=i;
                composition.rows=decomposition->rows;
                composition.atoms=decomposition->atoms;
                composition.shells=decomposition->shells;
                composition.retained_norm=decomposition->weight_sum;
                composition.residual_norm=decomposition->projection_residual_norm;
                if(!composition.residual_norm)
                    composition.residual_norm=nao_reconstruction_error(canonical,*integration,i);
                composition.complete=composition.residual_norm &&
                    *composition.residual_norm<=1e-5 &&
                    decomposition->normalization_error &&
                    *decomposition->normalization_error<=1e-5;
                composition.semantics="NAO squared coefficients in validated orthonormal active subspace; retained norm and residual are explicit";
                result.value=std::move(composition);
                result.status=RoutedStatus::Available;result.provider=RoutedProvider::Nbo;
                result.evidence={decomposition->matrix_source};
                result.reason=decomposition->detail;
            }
        }
        if(!result.available() && canonical.orbitals[i].chemistry.available) {
            const auto why=result.reason;
            result.value=legacy_mo(canonical.orbitals[i],i);
            result.method="legacy AO-reference projection";
            fallback(result,why);
            result.reason="Basis-dependent legacy AO reference; not NAO decomposition";
        }
        auto& relations=out.mo_relations[i];
        const auto* relation_cap=associated?nbo_capability(*integration,"structure"):nullptr;
        unavailable(relations,relation_cap,integration!=nullptr,
                    "verified NBO/NLMO-to-canonical projection relations");
        if(integration && !associated){relations.status=RoutedStatus::Rejected;
            relations.reason="NBO association does not match this canonical wavefunction";}
        if(associated){
            std::vector<RoutedRelationProjection> mapped;
            mapped.reserve(relation_hits[i].size());
            for(auto& [_,projection]:relation_hits[i])
                mapped.push_back(std::move(projection));
            if(!mapped.empty()){
                relations.status=RoutedStatus::Available;
                relations.provider=RoutedProvider::Nbo;
                relations.reason="Actual source-orbital transforms relate this canonical MO to localized records";
                relations.value=std::move(mapped);
            }else if(!out.local_relation_registry.empty()){
                relations.status=RoutedStatus::NotReportedAboveThreshold;
                relations.reason="No associated local relation passes the exported 1e-4 display-support cutoff for this canonical MO";
            }else if(relation_cap && relation_cap->available()){
                relations.status=RoutedStatus::NotReportedAboveThreshold;
                relations.reason="Validated local analysis reports no applicable relationship above its producer threshold";
            }
        }
    }
    InteractionEvidence graph_evidence;
    if(out.total_atomic_charge.available()) {
        graph_evidence.atomic_partial_charges=*out.total_atomic_charge.value;
        graph_evidence.atomic_charge_provenance=out.total_atomic_charge.provider==RoutedProvider::Nbo?
            DataProvenance::Derived:canonical.atomic_partial_charge_provenance;
    }
    out.interaction_graph.value=build_interaction_graph(canonical,graph_evidence);
    if(associated)if(const auto* cap=nbo_capability(*integration,"wiberg");
                   cap && cap->available()) {
        for(auto& edge:out.interaction_graph.value->edges)
            for(const auto& row:integration->structure)
                if(row.kind=="bond" && row.atoms.size()==2 && row.wiberg &&
                   std::minmax(row.atoms[0],row.atoms[1])==
                       std::minmax(static_cast<std::size_t>(edge.atom_a),
                                   static_cast<std::size_t>(edge.atom_b))) {
                    edge.wiberg_index=row.wiberg;
                    edge.wiberg_source_path=row.source.path;
                    break;
                }
    }
    out.interaction_graph.status=RoutedStatus::Available;
    out.interaction_graph.provider=out.total_atomic_charge.provider==RoutedProvider::Nbo?
        RoutedProvider::Nbo:RoutedProvider::Legacy;
    out.interaction_graph.method="legacy Mayer/geometry graph with routed atom-charge evidence";
    out.interaction_graph.reason="Mayer and Wiberg remain separate electronic observables";
    out.interaction_graph.fallback_reason=out.total_atomic_charge.fallback_reason;
    out.interaction_graph.evidence=out.total_atomic_charge.evidence;
    if(!integration){
        RoutedResult<NboPiCoupling> result;
        result.status=RoutedStatus::NotAnalysed;
        result.reason="NBO same-operator analysis has not been attached";
        result.method="verified NAO same-operator pi coupling";
        out.pi_couplings.push_back(std::move(result));
    }else if(!associated){
        RoutedResult<NboPiCoupling> result;
        result.status=RoutedStatus::Rejected;
        result.reason="NBO association does not match this canonical wavefunction";
        result.method="verified NAO same-operator pi coupling";
        out.pi_couplings.push_back(std::move(result));
    }else{
        std::vector<std::pair<std::size_t,std::size_t>> strong_connectivity;
        for(const auto& edge:out.interaction_graph.value->edges)
            if(edge.strength==InteractionStrength::StrongConnectivity)
                strong_connectivity.emplace_back(edge.atom_a,edge.atom_b);
        auto analysis=analyse_nbo_pi_couplings(canonical,*integration,strong_connectivity);
        for(auto& coupling:analysis.couplings){
            RoutedResult<NboPiCoupling> result;
            result.status=RoutedStatus::Available;result.provider=RoutedProvider::Nbo;
            result.method="same-operator NAO Fock block, SVD projector and canonical cross matrix";
            result.reason=analysis.reason;result.evidence={coupling.source};
            result.value=std::move(coupling);
            out.pi_couplings.push_back(std::move(result));
        }
        if(out.pi_couplings.empty()){
            RoutedResult<NboPiCoupling> result;
            result.status=analysis.status=="not_applicable"?RoutedStatus::NotApplicable:
                analysis.status=="rejected"?RoutedStatus::Rejected:RoutedStatus::Insufficient;
            result.method="verified NAO same-operator pi coupling";
            result.reason=analysis.reason;result.evidence=std::move(analysis.evidence);
            out.pi_couplings.push_back(std::move(result));
        }
    }
    return out;
}

std::string serialize_routed_analysis_json(const RoutedAnalysis& data) {
    std::ostringstream out;out<<std::setprecision(17)<<"{\"schema\":\"cov.chemistry.route.v1\""
       <<",\"canonical_fingerprint\":"<<quoted(data.canonical_fingerprint)
       <<",\"integration_id\":"<<quoted(data.integration_id)<<',';
    const auto atomic=[&](const char* key,const RoutedResult<std::vector<double>>& result) {
        result_header(out,key,result.status,result.provider,result.method,result.reason,
                      result.fallback_reason);
        sources_json(out,result);out<<",\"values\":";
        if(result.available()) {
            out<<'[';for(std::size_t i=0;i<result.value->size();++i) {
                if(i)out<<',';out<<(*result.value)[i];
            }out<<']';
        }else out<<"null";
        out<<'}';
    };
    atomic("total_atomic_charge",data.total_atomic_charge);out<<',';
    atomic("atomic_spin",data.atomic_spin);out<<",\"mo_composition\":[";
    for(std::size_t i=0;i<data.mo_composition.size();++i) {
        if(i)out<<',';const auto& result=data.mo_composition[i];
        out<<"{\"canonical_index\":"<<i<<",\"status\":"<<quoted(routed_status_name(result.status))
           <<",\"provider\":"<<quoted(routed_provider_name(result.provider))
           <<",\"method\":"<<quoted(result.method)
           <<",\"reason\":"<<quoted(result.reason)
           <<",\"fallback_reason\":"<<quoted(result.fallback_reason);
        sources_json(out,result);
        if(result.available()) {
            const auto& value=*result.value;
            out<<",\"retained_norm\":";
            if(value.retained_norm)out<<*value.retained_norm;else out<<"null";
            out<<",\"residual_norm\":";
            if(value.residual_norm)out<<*value.residual_norm;else out<<"null";
            out<<",\"legacy_coverage\":";
            if(value.legacy_coverage)out<<*value.legacy_coverage;else out<<"null";
            out<<",\"legacy_unresolved_fraction\":";
            if(value.legacy_unresolved_fraction)out<<*value.legacy_unresolved_fraction;else out<<"null";
            out<<",\"complete\":"<<(value.complete?"true":"false")
               <<",\"semantics\":"<<quoted(value.semantics);
            const auto groups=[&](const char* key,const auto& rows) {
                out<<",\""<<key<<"\":[";
                for(std::size_t j=0;j<rows.size();++j) {
                    if(j)out<<',';const auto& row=rows[j];
                    out<<"{\"atom\":"<<row.atom<<",\"label\":"<<quoted(row.label)
                       <<",\"weight\":"<<row.weight<<",\"principal_n\":";
                    if(row.principal_n)out<<*row.principal_n;else out<<"null";
                    out<<",\"angular_l\":";
                    if(row.angular_l)out<<*row.angular_l;else out<<"null";
                    out<<",\"nao_ids\":[";
                    for(std::size_t k=0;k<row.nao_ids.size();++k) {
                        if(k)out<<',';out<<row.nao_ids[k];
                    }
                    out<<"]}";
                }out<<']';
            };
            groups("atoms",value.atoms);groups("shells",value.shells);
            out<<",\"rows\":[";
            for(std::size_t j=0;j<value.rows.size();++j) {
                if(j)out<<',';const auto& row=value.rows[j];
                out<<"{\"nao_id\":"<<row.nao_id<<",\"atom\":"<<row.atom
                   <<",\"symbol\":"<<quoted(row.symbol)
                   <<",\"type\":"<<quoted(row.type)
                   <<",\"angular\":"<<quoted(row.angular)
                   <<",\"principal_n\":";
                if(row.principal_n)out<<*row.principal_n;else out<<"null";
                out<<",\"angular_l\":";
                if(row.angular_l)out<<*row.angular_l;else out<<"null";
                out<<",\"coefficient\":"<<row.coefficient
                   <<",\"weight\":"<<row.weight
                   <<",\"electron_contribution\":";
                if(row.electron_contribution)out<<*row.electron_contribution;else out<<"null";
                out<<'}';
            }
            out<<']';
        }else out<<",\"retained_norm\":null,\"residual_norm\":null,\"legacy_coverage\":null,\"legacy_unresolved_fraction\":null,\"complete\":null,\"atoms\":null,\"shells\":null,\"rows\":null";
        out<<'}';
    }
    out<<"],\"local_relation_registry\":[";
    for(std::size_t j=0;j<data.local_relation_registry.size();++j){
        if(j)out<<',';const auto& relation=data.local_relation_registry[j];
        out<<"{\"id\":"<<quoted(relation.id)<<",\"kind\":"<<quoted(relation.kind)
           <<",\"label\":"<<quoted(relation.label)<<",\"semantics\":"
           <<quoted(relation.semantics)<<",\"spin\":"<<quoted(nbo_spin_name(relation.spin))
           <<",\"atoms\":[";
        for(std::size_t k=0;k<relation.atoms.size();++k){if(k)out<<',';out<<relation.atoms[k];}
        out<<"],\"orbitals\":[";
        for(std::size_t k=0;k<relation.orbitals.size();++k){
            if(k)out<<',';const auto& ref=relation.orbitals[k];
            out<<"{\"kind\":"<<quoted(nbo_orbital_kind_name(ref.kind))
               <<",\"spin\":"<<quoted(nbo_spin_name(ref.spin))
               <<",\"index\":"<<ref.index<<'}';
        }
        out<<"],\"printed_value\":";
        if(relation.printed_value)out<<*relation.printed_value;else out<<"null";
        out<<",\"units\":"<<quoted(relation.units)<<",\"source\":";
        source_json(out,relation.source);out<<'}';
    }
    out<<"],\"mo_relations\":[";
    for(std::size_t i=0;i<data.mo_relations.size();++i){
        if(i)out<<',';
        const auto& result=data.mo_relations[i];
        out<<"{\"canonical_index\":"<<i<<",\"status\":"
           <<quoted(routed_status_name(result.status))<<",\"provider\":"
           <<quoted(routed_provider_name(result.provider))<<",\"method\":"
           <<quoted(result.method)<<",\"reason\":"<<quoted(result.reason)
           <<",\"support_cutoff\":0.0001";
        sources_json(out,result);out<<",\"relations\":";
        if(!result.available()){out<<"null}";continue;}
        out<<'[';
        for(std::size_t j=0;j<result.value->size();++j){
            if(j)out<<',';const auto& relation=(*result.value)[j];
            out<<"{\"relation_id\":"<<quoted(relation.relation_id)
               <<",\"canonical_projection_weights\":[";
            for(std::size_t k=0;k<relation.canonical_projection_weights.size();++k){
                if(k)out<<',';
                if(relation.canonical_projection_weights[k])
                    out<<*relation.canonical_projection_weights[k];else out<<"null";
            }
            out<<"]}";
        }
        out<<"]}";
    }
    out<<"],\"pi_couplings\":[";
    for(std::size_t i=0;i<data.pi_couplings.size();++i){
        if(i)out<<',';const auto& result=data.pi_couplings[i];
        out<<"{\"status\":"<<quoted(routed_status_name(result.status))
           <<",\"provider\":"<<quoted(routed_provider_name(result.provider))
           <<",\"method\":"<<quoted(result.method)
           <<",\"reason\":"<<quoted(result.reason);
        sources_json(out,result);out<<",\"value\":";
        if(!result.available()){out<<"null}";continue;}
        const auto& value=*result.value;
        out<<"{\"id\":"<<quoted(value.id)<<",\"spin\":"
           <<quoted(nbo_spin_name(value.spin))<<",\"channel\":"
           <<quoted(value.channel)<<",\"centre_atoms\":[";
        const auto indices=[&](const auto& rows){
            for(std::size_t j=0;j<rows.size();++j){if(j)out<<',';out<<rows[j];}
            out<<']';
        };
        indices(value.centre_atoms);out<<",\"ligand_atoms\":[";
        indices(value.ligand_atoms);out<<",\"centre_nao_ids\":[";
        indices(value.centre_nao_ids);out<<",\"ligand_nao_ids\":[";
        indices(value.ligand_nao_ids);out<<",\"singular_values_hartree\":[";
        indices(value.singular_values_hartree);
        out<<",\"coupled_rank\":"<<value.coupled_rank
           <<",\"centre_onsite_hartree\":"<<value.centre_onsite_hartree
           <<",\"ligand_onsite_hartree\":"<<value.ligand_onsite_hartree
           <<",\"centre_onsite_range_hartree\":["<<value.centre_onsite_range_hartree[0]
           <<','<<value.centre_onsite_range_hartree[1]<<']'
           <<",\"ligand_onsite_range_hartree\":["<<value.ligand_onsite_range_hartree[0]
           <<','<<value.ligand_onsite_range_hartree[1]<<']'
           <<",\"occupation_status\":"<<quoted(value.occupation_status)
           <<",\"occupation_reason\":"<<quoted(value.occupation_reason)
           <<",\"centre_occupation\":";
        if(value.centre_occupation)out<<*value.centre_occupation;else out<<"null";
        out<<",\"ligand_occupation\":";
        if(value.ligand_occupation)out<<*value.ligand_occupation;else out<<"null";
        const auto optional_range=[&](const std::optional<std::array<double,2>>& range){
            if(range)out<<'['<<(*range)[0]<<','<<(*range)[1]<<']';else out<<"null";
        };
        out<<",\"centre_occupation_range\":";optional_range(value.centre_occupation_range);
        out<<",\"ligand_occupation_range\":";optional_range(value.ligand_occupation_range);
        out<<",\"operator_max_error_hartree\":"<<value.operator_max_error_hartree
           <<",\"nao_orthogonality_error\":"<<value.nao_orthogonality_error
           <<",\"angular_leakage\":"<<value.angular_leakage
           <<",\"minimum_centre_projection\":"<<value.minimum_centre_projection
           <<",\"angular_projector_evidence\":[";
        for(std::size_t j=0;j<value.angular_projector_evidence.size();++j){
            if(j)out<<',';const auto& row=value.angular_projector_evidence[j];
            out<<"{\"atom\":"<<row.atom<<",\"method\":"<<quoted(row.method)
               <<",\"rank\":"<<row.rank
               <<",\"p_metric_min_eigenvalue\":"<<row.p_metric_min_eigenvalue
               <<",\"pi_generalized_eigenvalues\":[";
            for(std::size_t k=0;k<3;++k){if(k)out<<',';out<<row.pi_generalized_eigenvalues[k];}
            out<<"],\"max_sigma_leakage\":"<<row.max_sigma_leakage
               <<",\"centre_projection\":"<<row.centre_projection
               <<",\"partition_residual\":"<<row.partition_residual<<'}';
        }
        out<<"],\"groups\":[";
        for(std::size_t j=0;j<value.groups.size();++j){
            if(j)out<<',';const auto& group=value.groups[j];
            out<<"{\"members\":[";indices(group.members);
            out<<",\"spin\":"<<quoted(nbo_spin_name(group.spin))
               <<",\"centre_weight\":"<<group.centre_weight
               <<",\"ligand_weight\":"<<group.ligand_weight
               <<",\"occupation_per_mo\":";
            if(group.occupation_per_mo)out<<*group.occupation_per_mo;else out<<"null";
            out<<",\"cross_fock_min_hartree\":"<<group.cross_fock_min_hartree
               <<",\"cross_fock_max_hartree\":"<<group.cross_fock_max_hartree
               <<",\"cross_fock_mean_hartree\":"<<group.cross_fock_mean_hartree
               <<",\"character\":"<<quoted(group.character)<<'}';
        }
        out<<"],\"direction\":"<<quoted(value.direction)
           <<",\"direction_evidence\":"<<quoted(value.direction_evidence)
           <<",\"direction_donor_weight\":";
        if(value.direction_donor_weight)out<<*value.direction_donor_weight;else out<<"null";
        out<<",\"direction_acceptor_weight\":";
        if(value.direction_acceptor_weight)out<<*value.direction_acceptor_weight;else out<<"null";
        out<<",\"direction_projection_evidence\":[";
        for(std::size_t j=0;j<value.direction_projection_evidence.size();++j){
            if(j)out<<',';const auto& row=value.direction_projection_evidence[j];
            out<<"{\"e2_id\":"<<quoted(row.e2_id)
               <<",\"donor_nbo_id\":"<<row.donor_nbo_id
               <<",\"acceptor_nbo_id\":"<<row.acceptor_nbo_id
               <<",\"donor_ligand_weight\":"<<row.donor_ligand_weight
               <<",\"donor_centre_weight\":"<<row.donor_centre_weight
               <<",\"acceptor_ligand_weight\":"<<row.acceptor_ligand_weight
               <<",\"acceptor_centre_weight\":"<<row.acceptor_centre_weight
               <<",\"source\":";source_json(out,row.source);out<<'}';
        }
        out<<']'
           <<",\"source\":";source_json(out,value.source);out<<"}}";
    }
    out<<"],";
    const auto& graph=data.interaction_graph;
    result_header(out,"interaction_graph",graph.status,graph.provider,graph.method,
                  graph.reason,graph.fallback_reason);
    sources_json(out,graph);
    out<<",\"charge_scheme\":"<<quoted(data.total_atomic_charge.provider==RoutedProvider::Nbo?
        "validated NPA total":data.total_atomic_charge.provider==RoutedProvider::Legacy?
        data.total_atomic_charge.method:"none")
       <<",\"edges\":[";
    if(graph.available())for(std::size_t i=0;i<graph.value->edges.size();++i) {
        if(i)out<<',';const auto& edge=graph.value->edges[i];
        out<<"{\"atom_a\":"<<edge.atom_a<<",\"atom_b\":"<<edge.atom_b
           <<",\"kind\":"<<quoted(interaction_kind_name(edge.kind))
           <<",\"strength\":"<<quoted(edge.strength==InteractionStrength::StrongConnectivity?
                "strong_connectivity":"weak_contact")
           <<",\"mayer_order\":"<<edge.mayer_order
           <<",\"wiberg_index\":";
        if(edge.wiberg_index)out<<*edge.wiberg_index;else out<<"null";
        out<<",\"wiberg_source_path\":"<<quoted(edge.wiberg_source_path)
           <<",\"distance_bohr\":"<<edge.distance_bohr<<'}';
    }
    out<<"],\"fragments\":[";
    if(graph.available())for(std::size_t i=0;i<graph.value->fragment_analysis.fragments.size();++i) {
        if(i)out<<',';const auto& fragment=graph.value->fragment_analysis.fragments[i];
        out<<"{\"index\":"<<fragment.index<<",\"atoms\":[";
        for(std::size_t j=0;j<fragment.atoms.size();++j) {
            if(j)out<<',';out<<fragment.atoms[j];
        }
        out<<"],\"evidenced_charge\":";
        if(fragment.charge_provenance!=DataProvenance::Unavailable)
            out<<fragment.evidenced_charge;else out<<"null";
        out<<'}';
    }
    out<<"]}}";
    return out.str();
}

std::string serialize_active_orbital_view_json(const ActiveOrbitalView& view) {
    std::ostringstream out;out<<std::setprecision(17)
        <<"{\"schema\":\"cov.orbital.active-view.v1\",\"kind\":"
        <<quoted(active_orbital_kind_name(view.kind))
        <<",\"source_id\":"<<quoted(view.source_id)
        <<",\"label\":"<<quoted(view.label)
        <<",\"source_label\":"<<quoted(view.source_label)
        <<",\"display_name_evidence\":"<<quoted(view.display_name_evidence)
        <<",\"semantic_kind\":"<<quoted(view.semantic_kind)
        <<",\"group_id\":"<<quoted(view.group_id)
        <<",\"spin\":"<<quoted(nbo_spin_name(view.spin))
        <<",\"source_spin\":"<<quoted(nbo_spin_name(view.source_spin))
        <<",\"spin_semantics\":"<<quoted(view.spin_semantics)
        <<",\"canonical_index\":";
    if(view.canonical_index)out<<*view.canonical_index;else out<<"null";
    out<<",\"rendered_index\":";
    if(view.rendered_index)out<<*view.rendered_index;else out<<"null";
    out<<",\"selection\":";
    if(view.selection) {
        const auto& selection=*view.selection;
        out<<"{\"dataset_id\":"<<quoted(selection.dataset_id)
           <<",\"label\":"<<quoted(selection.label)
           <<",\"mode\":"<<static_cast<int>(selection.mode)
           <<",\"normalize\":"<<(selection.normalize?"true":"false")
           <<",\"target_canonical_index\":";
        if(selection.target_canonical_index)out<<*selection.target_canonical_index;else out<<"null";
        out<<",\"terms\":[";
        for(std::size_t i=0;i<selection.terms.size();++i) {
            if(i)out<<',';const auto& term=selection.terms[i];
            out<<"{\"kind\":"<<quoted(nbo_orbital_kind_name(term.orbital.kind))
               <<",\"spin\":"<<quoted(nbo_spin_name(term.orbital.spin))
               <<",\"index\":"<<term.orbital.index
               <<",\"coefficient\":"<<term.coefficient<<'}';
        }out<<"]}";
    }else out<<"null";
    out<<'}';return out.str();
}

} // namespace cov
