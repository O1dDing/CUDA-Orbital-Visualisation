#include "cov/mo_diagram.hpp"
#include "cov/orbital_chemistry.hpp"
#include "cov/orbital_view.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

cov::Primitive primitive(const float exponent) {
    cov::Primitive value;
    value.exponent=exponent;
    value.coefficient=1.0f;
    return value;
}

cov::Wavefunction make_h3_reference_case() {
    cov::Wavefunction wf;
    wf.atoms={
        {"H",1,-1.0,0.0,0.0,1.0},
        {"H",1, 0.5,0.866025403784,0.0,1.0},
        {"H",1, 0.5,-0.866025403784,0.0,1.0},
    };

    for (std::uint32_t atom=0;atom<3u;++atom) {
        const std::uint32_t primitive_offset=
            static_cast<std::uint32_t>(wf.primitives.size());
        wf.primitives.push_back(primitive(1.2f));
        wf.shells.push_back({atom,primitive_offset,1u,wf.basis_count,0u,0u});
        wf.basis_count+=1u;

        const std::uint32_t diffuse_offset=
            static_cast<std::uint32_t>(wf.primitives.size());
        wf.primitives.push_back(primitive(0.010f));
        wf.shells.push_back({atom,diffuse_offset,1u,wf.basis_count,0u,0u});
        wf.basis_count+=1u;

        const std::uint32_t p_offset=
            static_cast<std::uint32_t>(wf.primitives.size());
        wf.primitives.push_back(primitive(0.050f));
        wf.shells.push_back({atom,p_offset,1u,wf.basis_count,1u,0u});
        wf.basis_count+=3u;
    }

    const std::size_t n=wf.basis_count;
    wf.ao_overlap.assign(n*n,0.0);
    for (std::size_t i=0;i<n;++i) wf.ao_overlap[i*n+i]=1.0;
    const std::array<std::size_t,3> main{0u,5u,10u};
    constexpr double s=0.20;
    for (std::size_t i=0;i<3u;++i) {
        for (std::size_t j=i+1u;j<3u;++j) {
            wf.ao_overlap[main[i]*n+main[j]]=s;
            wf.ao_overlap[main[j]*n+main[i]]=s;
        }
    }

    cov::MolecularOrbital a1;
    a1.energy_hartree=-1.022932;
    a1.occupation=2.0f;
    a1.symmetry="A1'";
    a1.coefficients.assign(n,0.0f);
    const double ca=1.0/std::sqrt(3.0+6.0*s);
    for (const auto index:main) a1.coefficients[index]=static_cast<float>(ca);

    cov::MolecularOrbital e1;
    e1.energy_hartree=-0.34;
    e1.occupation=0.0f;
    e1.symmetry="E'";
    e1.coefficients.assign(n,0.0f);
    const double ce1=1.0/std::sqrt(6.0*(1.0-s));
    e1.coefficients[main[0]]=static_cast<float>(2.0*ce1);
    e1.coefficients[main[1]]=static_cast<float>(-ce1);
    e1.coefficients[main[2]]=static_cast<float>(-ce1);

    cov::MolecularOrbital e2;
    e2.energy_hartree=-0.34;
    e2.occupation=0.0f;
    e2.symmetry="E'";
    e2.coefficients.assign(n,0.0f);
    const double ce2=1.0/std::sqrt(2.0*(1.0-s));
    e2.coefficients[main[1]]=static_cast<float>(ce2);
    e2.coefficients[main[2]]=static_cast<float>(-ce2);

    wf.orbitals={a1,e1,e2};
    for (std::size_t basis=0;basis<n;++basis) {
        if (basis==main[0] || basis==main[1] || basis==main[2]) continue;
        cov::MolecularOrbital extra;
        extra.energy_hartree=0.1+0.1*static_cast<double>(wf.orbitals.size());
        extra.occupation=0.0f;
        extra.coefficients.assign(n,0.0f);
        extra.coefficients[basis]=1.0f;
        wf.orbitals.push_back(std::move(extra));
    }

    wf.bond_orders={
        {0u,1u,0.5,cov::DataProvenance::Derived},
        {1u,2u,0.5,cov::DataProvenance::Derived},
        {0u,2u,0.5,cov::DataProvenance::Derived},
    };
    return wf;
}

} // namespace

int main() {
    auto wf=make_h3_reference_case();
    cov::derive_orbital_chemistry(wf);

    std::vector<std::size_t> selected;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        if (wf.orbitals[i].chemistry.valence_manifold) selected.push_back(i);
    }
    if (selected!=std::vector<std::size_t>{0u,1u,2u}) {
        std::cerr<<"H3 chemical-valence manifold must be MO1 + E' pair only\n";
        return EXIT_FAILURE;
    }
    if (wf.orbitals[0].chemistry.channel.dominant!=
            cov::OrbitalAngularFamily::Sigma ||
        wf.orbitals[0].chemistry.channel.status!=
            cov::ChemistryStatus::Determined) {
        std::cerr<<"H3 occupied valence MO was not classified as sigma\n";
        return EXIT_FAILURE;
    }
    if (wf.orbitals[0].chemistry.bonding.dominant!=
            cov::OrbitalBondingRole::Bonding) {
        std::cerr<<"H3 occupied valence MO was not classified as bonding\n";
        return EXIT_FAILURE;
    }
    if (wf.orbitals[0].chemistry.multicentre_label!="3c2e" ||
        wf.orbitals[0].chemistry.participating_atoms!=3u ||
        std::abs(wf.orbitals[0].chemistry.participating_electrons-2.0)>1.0e-8) {
        std::cerr<<"H3 3c2e family was not attached to the valence manifold\n";
        return EXIT_FAILURE;
    }

    cov::OrbitalFilterSettings filter;
    filter.mode=cov::OrbitalFilterMode::AutoReasonable;
    const auto frontier=cov::find_frontier_orbitals(wf.orbitals);
    const auto visible=cov::visible_orbital_indices(wf.orbitals,frontier,filter);
    if (visible!=selected) {
        std::cerr<<"Auto filter did not use the chemical-valence manifold\n";
        return EXIT_FAILURE;
    }

    cov::MODiagramOptions options;
    options.selected_index=0;
    const auto diagram=cov::build_mo_diagram_data(wf,options);
    if (diagram.levels.size()!=3u) {
        std::cerr<<"MO diagram did not reduce H3 to three chemical-valence MOs\n";
        return EXIT_FAILURE;
    }
    if (!diagram.levels[0].chemistry.available ||
        diagram.levels[0].chemistry.ao_contributions.size()!=3u) {
        std::cerr<<"MO diagram did not carry chemistry/AO composition\n";
        return EXIT_FAILURE;
    }

    std::cout<<"chemical valence H3 smoke test passed\n";
    return EXIT_SUCCESS;
}
