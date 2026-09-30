#include "cov/chemistry_route.hpp"
#include "cov/mo_diagram.hpp"
#include "cov/nbo_channels.hpp"
#include "cov/wavefunction_io.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>

// Small observation adapter for saved inputs; acceptance compares the actual
// channel/member records, not a process exit or a molecule-name lookup.
int main(int argc,char** argv) {try {
    if(argc!=4)throw std::runtime_error("Usage: cov_pi_diagram_probe canonical.fchk analysis_directory|- output.json");
    const auto w=cov::parse_wavefunction(argv[1]);
    const auto identity=cov::nbo_canonical_fingerprint(w);
    std::optional<cov::NboIntegration> integration;
    if(std::string(argv[2])!="-") {
        const auto found=cov::discover_nbo_inputs({std::filesystem::path(argv[2])});
        if(found.candidates.size()!=1)throw std::runtime_error("Exactly one NBO analysis is required");
        integration=cov::read_nbo_integration(w,found.candidates.front());
        cov::annotate_nbo_bond_channels(*integration,w);
    }
    const auto route=cov::route_chemistry(w,integration?&*integration:nullptr);
    std::ofstream out(argv[3]);
    if(!out)throw std::runtime_error("Cannot write output");
    out<<"{\"views\":[";
    for(int compact=0;compact<2;++compact) {
        cov::MODiagramOptions options;options.routed=&route;
        options.routed_identity=route.canonical_fingerprint;
        options.hide_ligand_centred_intermediates=compact!=0;
        const auto diagram=cov::build_mo_diagram_data(w,options);
        if(compact)out<<',';
        out<<"{\"compact\":"<<(compact?"true":"false")
           <<",\"pi_interactions\":"<<cov::orbital_energy_gap_array_json(diagram.pi_interactions)
           <<",\"candidates\":"<<cov::pi_partner_candidates_json(diagram.pi_partner_candidates)
           <<",\"visible_members\":[";
        bool comma=false;
        for(const auto& row:diagram.levels)for(const auto member:row.member_indices) {
            if(comma)out<<',';comma=true;out<<member;
        }
        out<<"]}";
    }
    out<<"],\"canonical_immutable\":"<<(identity==cov::nbo_canonical_fingerprint(w)?"true":"false")<<'}';
    return out?0:2;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
