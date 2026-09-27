#pragma once
#include "cov/numerical_diagnostics.hpp"

namespace cov {
inline const char* pi_topology_name(DelocalisedPiTopology topology) noexcept {
    switch(topology) {
        case DelocalisedPiTopology::Path:return "path";
        case DelocalisedPiTopology::Cycle:return "cycle";
        case DelocalisedPiTopology::BranchedResonance:return "branched-resonance";
        case DelocalisedPiTopology::Spiro:return "spiro";
        case DelocalisedPiTopology::HapticMetal:return "haptic-metal";
        case DelocalisedPiTopology::SymmetryDirectSum:return "symmetry-direct-sum";
        case DelocalisedPiTopology::MultiChannel:return "multi-channel";
        default:return "unknown";
    }
}
inline const char* pi_topology_graph_source_name(PiTopologyGraphSource source) noexcept {
    switch(source) {
        case PiTopologyGraphSource::MayerDistanceModel:return "mayer-and-covalent-distance-model";
        case PiTopologyGraphSource::CovalentDistanceModel:return "covalent-distance-model";
        default:return "unavailable";
    }
}
namespace topology_json {
template<class Range> inline void indices(std::ostream& out,const Range& values) {
    out<<'[';bool first=true;for(const auto value:values) {
        if(!first)out<<',';first=false;out<<value;
    }out<<']';
}
inline void channels(std::ostream& out,const std::vector<PiOrientationChannel>& channels) {
    out<<'[';bool first=true;
    for(const auto& channel:channels) {
        if(!first)out<<',';first=false;
        out<<"{\"atoms\":";indices(out,channel.atoms);
        out<<",\"direction\":[";
        for(std::size_t i=0;i<3;++i) {if(i)out<<',';numerical_json::number(out,channel.direction[i]);}
        out<<"],\"coherence\":";numerical_json::number(out,channel.coherence);
        out<<",\"cyclic\":"<<(channel.cyclic?"true":"false")<<'}';
    }out<<']';
}
}

// This is a serializer of one production result, not an alternative detector.
inline void write_pi_topology_graph_json(std::ostream& out,const PiTopologyGraphEvidence& evidence) {
    using namespace numerical_json;
    if(evidence.source==PiTopologyGraphSource::Unavailable) {out<<"null";return;}
    out<<"{\"schema\":1,\"source\":";string(out,pi_topology_graph_source_name(evidence.source));
    out<<",\"scope\":\"full input atom graph; local ring paths with electronic-channel correspondence\","
          "\"atom_index_base\":0,\"channel_index_base\":0,\"closed_cycle_paths\":true,\"atom_count\":"<<evidence.atom_count;
    out<<",\"bond_order_provenance\":";string(out,data_provenance_name(evidence.bond_order_provenance));
    out<<",\"minimum_mayer_order\":";number(out,evidence.minimum_mayer_order);
    out<<",\"maximum_covalent_radius_factor\":";number(out,evidence.maximum_covalent_radius_factor);
    out<<",\"criteria_units\":\"dimensionless\",\"edges\":[";
    for(std::size_t i=0;i<evidence.edges.size();++i) {
        if(i)out<<',';out<<'['<<evidence.edges[i].first<<','<<evidence.edges[i].second<<']';
    }
    out<<"],\"association_policy\":\"two-exclusive-channel-atoms-per-ring; no other-channel atom except hub; intervening non-channel atoms allowed\","
          "\"absence_scope\":\"no matching witness under the declared graph and channel boundary; not a molecular acyclicity claim\","
          "\"channel_association_search_complete\":"<<(evidence.channel_association_search_complete?"true":"false");
    out<<",\"examined_hubs\":";topology_json::indices(out,evidence.examined_hubs);
    out<<",\"unscoped_cycle_union_hubs\":";topology_json::indices(out,evidence.unscoped_cycle_union_hubs);
    out<<",\"channel_ring_witnesses\":[";
    for(std::size_t i=0;i<evidence.channel_ring_witnesses.size();++i) {
        if(i)out<<',';const auto& record=evidence.channel_ring_witnesses[i];
        out<<"{\"channel_indices\":";topology_json::indices(out,record.channel_indices);
        out<<",\"hub\":"<<record.ring_union.hub<<",\"cycles\":[";
        topology_json::indices(out,record.ring_union.cycles[0]);out<<',';
        topology_json::indices(out,record.ring_union.cycles[1]);
        out<<"],\"additional_connection_without_hub\":"
            <<(record.ring_union.additional_connection_without_hub?"true":"false")<<'}';
    }
    out<<"]}";
}

inline void write_pi_topology_assignments_json(std::ostream& out,const Wavefunction& wf) {
    using namespace numerical_json;
    out<<'[';bool first=true;
    for(const auto& assignment:wf.delocalised_pi_assignments) {
        if(!first)out<<',';first=false;
        out<<"{\"family_id\":";string(out,assignment.family_id);
        out<<",\"atoms\":";topology_json::indices(out,assignment.atoms);
        out<<",\"orbitals\":";topology_json::indices(out,assignment.orbitals);
        out<<",\"electron_count\":";number(out,assignment.electron_count);
        out<<",\"topology\":";string(out,pi_topology_name(assignment.topology));
        out<<",\"cyclic_pi_topology\":"<<(assignment.cyclic_topology?"true":"false");
        out<<",\"orientation_channels\":";topology_json::channels(out,assignment.orientation_channels);
        out<<",\"topology_graph\":";write_pi_topology_graph_json(out,assignment.topology_graph);
        out<<",\"provenance\":";string(out,data_provenance_name(assignment.provenance));
        out<<'}';
    }
    out<<']';
}
}
