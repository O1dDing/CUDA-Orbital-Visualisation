#include "cov/mo_diagram.hpp"
#include "cov/orbital_chemistry.hpp"
#include "cov/wavefunction_io.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr std::array<std::size_t,5> kExpectedMembers{
    15u,16u,17u,41u,42u // raw MO16, MO17/MO18 and MO42/MO43
};

cov::Primitive primitive(const float exponent) {
    cov::Primitive value;
    value.exponent=exponent;
    value.coefficient=1.0f;
    return value;
}

void append_shell(cov::Wavefunction& wf,
                  const std::uint32_t atom,
                  const std::uint8_t angular_momentum,
                  const std::uint8_t pure,
                  const float exponent) {
    const auto primitive_offset=
        static_cast<std::uint32_t>(wf.primitives.size());
    wf.primitives.push_back(primitive(exponent));
    wf.shells.push_back({
        atom,primitive_offset,1u,wf.basis_count,angular_momentum,pure});
    wf.basis_count+=cov::shell_basis_count(wf.shells.back());
}

cov::Wavefunction make_planar_c5_case(const bool distort_out_of_plane) {
    cov::Wavefunction wf;
    constexpr double pi=3.1415926535897932384626433832795;
    constexpr double radius=2.25; // bohr; nearest-neighbour C--C = 2.65 bohr
    for (std::size_t atom=0;atom<5u;++atom) {
        const double angle=2.0*pi*static_cast<double>(atom)/5.0;
        wf.atoms.push_back({
            "C",6,radius*std::cos(angle),radius*std::sin(angle),0.0,6.0});
    }
    if (distort_out_of_plane) {
        // Deliberately far outside any chemically reasonable Cp ring plane.
        // Keeping all Mayer edges present makes this a direct planarity gate,
        // rather than accidentally passing because the graph was disconnected.
        wf.atoms[0].z=4.0;
    }

    // One minimal C 2s/2p shell set per centre. The p shell is Cartesian in
    // x,y,z order, so the five z functions are the ring p_perp reference.
    std::array<std::size_t,5> pz{};
    for (std::uint32_t atom=0;atom<5u;++atom) {
        append_shell(wf,atom,0u,0u,1.2f);
        append_shell(wf,atom,1u,0u,0.8f);
        pz[atom]=wf.shells.back().basis_offset+2u;
    }

    // Add 23 polarization/diffuse degrees of freedom so there are 43 complete
    // canonical MOs. These functions must not become extra carbon p-valence
    // references: three spherical f shells contribute 21, plus two diffuse s.
    for (std::uint32_t atom=0;atom<3u;++atom) {
        append_shell(wf,atom,3u,1u,0.35f+0.03f*static_cast<float>(atom));
    }
    append_shell(wf,0u,0u,0u,0.04f);
    append_shell(wf,1u,0u,0u,0.04f);

    if (wf.basis_count!=43u) {
        std::cerr<<"internal fixture error: expected 43 basis functions, got "
                 <<wf.basis_count<<'\n';
        std::abort();
    }
    const std::size_t n=wf.basis_count;
    wf.ao_overlap.assign(n*n,0.0);
    for (std::size_t i=0;i<n;++i) wf.ao_overlap[i*n+i]=1.0;
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;

    std::vector<std::size_t> complement;
    complement.reserve(n-5u);
    for (std::size_t basis=0;basis<n;++basis) {
        if (std::find(pz.begin(),pz.end(),basis)==pz.end()) {
            complement.push_back(basis);
        }
    }

    wf.orbitals.resize(n);
    std::size_t filler=0u;
    for (std::size_t i=0;i<n;++i) {
        auto& mo=wf.orbitals[i];
        mo.coefficients.assign(n,0.0f);
        mo.symmetry="A"+std::to_string(i+1u);
        if (i<15u) mo.energy_hartree=-1.50+0.05*static_cast<double>(i);
        else if (i==15u) mo.energy_hartree=-0.70;
        else if (i<=17u) mo.energy_hartree=-0.50;
        else if (i<41u) mo.energy_hartree=-0.30+
            0.018*static_cast<double>(i-18u);
        else mo.energy_hartree=0.30;
        mo.occupation=i<=17u?2.0f:0.0f;
        mo.occupation_provenance=cov::DataProvenance::Producer;

        if (std::find(kExpectedMembers.begin(),kExpectedMembers.end(),i)==
            kExpectedMembers.end()) {
            mo.coefficients[complement[filler++]]=1.0f;
        }
    }

    // Real orthonormal SALCs of a regular five-membered p_perp ring:
    // k=0; cos/sin(k=1); cos/sin(k=2). Their occupations are 2,2,2,0,0.
    for (std::size_t member=0;member<kExpectedMembers.size();++member) {
        auto& mo=wf.orbitals[kExpectedMembers[member]];
        for (std::size_t atom=0;atom<5u;++atom) {
            const double angle=2.0*pi*static_cast<double>(atom)/5.0;
            double coefficient=0.0;
            if (member==0u) coefficient=1.0/std::sqrt(5.0);
            else if (member==1u) coefficient=std::sqrt(2.0/5.0)*std::cos(angle);
            else if (member==2u) coefficient=std::sqrt(2.0/5.0)*std::sin(angle);
            else if (member==3u) coefficient=std::sqrt(2.0/5.0)*std::cos(2.0*angle);
            else coefficient=std::sqrt(2.0/5.0)*std::sin(2.0*angle);
            mo.coefficients[pz[atom]]=static_cast<float>(coefficient);
        }
    }
    wf.orbitals[15].symmetry="A2\"";
    wf.orbitals[16].symmetry="E1\"";
    wf.orbitals[17].symmetry="E1\"";
    wf.orbitals[41].symmetry="E2\"";
    wf.orbitals[42].symmetry="E2\"";

    // Explicit ring-only connectivity. Non-neighbouring carbon pairs remain
    // below the graph threshold even though the pentagon is compact.
    for (std::uint32_t atom=0;atom<5u;++atom) {
        wf.bond_orders.push_back({
            atom,static_cast<std::uint32_t>((atom+1u)%5u),
            1.20,cov::DataProvenance::Derived});
    }
    wf.bond_order_provenance=cov::DataProvenance::Derived;
    return wf;
}

std::vector<std::size_t> detected_members(const cov::Wavefunction& wf) {
    std::vector<std::size_t> result;
    for (std::size_t i=0;i<wf.orbitals.size();++i) {
        const auto& chemistry=wf.orbitals[i].chemistry;
        if (!chemistry.delocalised_family_id.empty()) result.push_back(i);
    }
    return result;
}

bool validate_cp_family(const cov::Wavefunction& wf,const char* context) {
    const std::vector<std::size_t> expected_indices(
        kExpectedMembers.begin(),kExpectedMembers.end());
    const std::vector<std::uint32_t> expected_orbitals{
        15u,16u,17u,41u,42u};
    const auto detected=detected_members(wf);
    if (detected!=expected_indices) {
        std::cerr<<context<<": expected delocalised-pi members "
                 <<"16,17,18,42,43; detected";
        for (const auto index:detected) std::cerr<<' '<<(index+1u);
        std::cerr<<'\n';
        return false;
    }

    if (wf.delocalised_pi_assignments.size()!=1u) {
        std::cerr<<context<<": expected one wavefunction-level delocalised family, got "
                 <<wf.delocalised_pi_assignments.size()<<'\n';
        return false;
    }
    const auto& assignment=wf.delocalised_pi_assignments.front();
    const std::vector<std::uint32_t> expected_atoms{0u,1u,2u,3u,4u};
    if (assignment.family_id.empty() || assignment.atoms!=expected_atoms ||
        assignment.orbitals!=expected_orbitals ||
        std::abs(assignment.electron_count-6.0)>1.0e-8 ||
        assignment.provenance!=cov::DataProvenance::Derived) {
        std::cerr<<context<<": wavefunction-level Π^5_6 assignment is incomplete\n";
        return false;
    }

    const auto& reference=wf.orbitals[kExpectedMembers.front()].chemistry;
    if (reference.delocalised_family_id.empty() ||
        reference.multicentre_label!="delocalised-pi" ||
        reference.family_symbol!="pi" ||
        reference.participating_atoms!=5u ||
        std::abs(reference.participating_electrons-6.0)>1.0e-8 ||
        reference.participating_atom_indices!=expected_atoms ||
        reference.delocalised_family_orbitals!=expected_orbitals) {
        std::cerr<<context<<": incomplete Π^5_6 family metadata\n";
        return false;
    }

    const std::string expected_label=
        "\xCE\xA0\xE2\x81\xB5\xE2\x82\x86"; // Π⁵₆
    for (const auto index:kExpectedMembers) {
        const auto& chemistry=wf.orbitals[index].chemistry;
        if (chemistry.delocalised_family_id!=reference.delocalised_family_id ||
            chemistry.delocalised_family_orbitals!=expected_orbitals ||
            chemistry.participating_atom_indices!=expected_atoms ||
            chemistry.channel.dominant!=cov::OrbitalAngularFamily::Pi) {
            std::cerr<<context<<": family membership differs at MO"
                     <<(index+1u)<<'\n';
            return false;
        }
        const auto annotation=cov::annotate_orbital(wf.orbitals[index]);
        if (!annotation.delocalised_pi.available ||
            annotation.delocalised_pi.label!=expected_label ||
            annotation.delocalised_pi.participating_atoms!=5u ||
            std::abs(annotation.delocalised_pi.participating_electrons-6.0)>1.0e-8 ||
            annotation.delocalised_pi.orbital_indices!=expected_indices ||
            annotation.delocalised_pi.family_id!=reference.delocalised_family_id) {
            std::cerr<<context<<": Π⁵₆ diagram descriptor missing at MO"
                     <<(index+1u)<<'\n';
            return false;
        }
    }
    return true;
}

bool validate_optional_real_input(const std::filesystem::path& path) {
    try {
        const auto wf=cov::parse_wavefunction(path);
        return validate_cp_family(wf,"real Cp input");
    } catch (const std::exception& error) {
        std::cerr<<"real Cp input failed: "<<error.what()<<'\n';
        return false;
    }
}

} // namespace

int main(const int argc,char** argv) {
    auto planar=make_planar_c5_case(false);
    cov::derive_orbital_chemistry(planar);
    if (!validate_cp_family(planar,"synthetic planar C5")) {
        return EXIT_FAILURE;
    }

    auto nonplanar=make_planar_c5_case(true);
    cov::derive_orbital_chemistry(nonplanar);
    if (!detected_members(nonplanar).empty()) {
        std::cerr<<"non-planar C5 negative control was classified as delocalised pi\n";
        return EXIT_FAILURE;
    }

    if (argc>1 && !validate_optional_real_input(argv[1])) {
        return EXIT_FAILURE;
    }

    std::cout<<"delocalised pi smoke test passed: Π⁵₆ = raw MO16/17/18/42/43\n";
    return EXIT_SUCCESS;
}
