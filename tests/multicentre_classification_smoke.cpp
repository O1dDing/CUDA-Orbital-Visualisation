#include "cov/bond_analysis.hpp"
#include "cov/density.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

cov::Wavefunction make_three_centre(const bool four_electron) {
    cov::Wavefunction wf;
    wf.atoms={
        {"X",1,-1.0,0.0,0.0},
        {"X",1, 0.0,0.0,0.0},
        {"X",1, 1.0,0.0,0.0},
    };
    wf.basis_count=3;
    for (std::uint32_t atom=0;atom<3;++atom) {
        const std::uint32_t primitive_offset=static_cast<std::uint32_t>(wf.primitives.size());
        wf.primitives.push_back({1.0f,1.0f});
        cov::Shell shell;
        shell.atom_index=atom;
        shell.primitive_offset=primitive_offset;
        shell.primitive_count=1;
        shell.basis_offset=atom;
        shell.angular_momentum=0;
        shell.pure=0;
        wf.shells.push_back(shell);
    }

    const double r3=std::sqrt(3.0);
    const double r2=std::sqrt(2.0);
    const double r6=std::sqrt(6.0);

    cov::MolecularOrbital bonding;
    bonding.energy_hartree=-1.0;
    bonding.occupation=2.0f;
    bonding.spin=cov::Spin::Alpha;
    bonding.coefficients={
        static_cast<float>(1.0/r3),
        static_cast<float>(1.0/r3),
        static_cast<float>(1.0/r3),
    };

    cov::MolecularOrbital nonbonding;
    nonbonding.energy_hartree=-0.2;
    nonbonding.occupation=four_electron?2.0f:0.0f;
    nonbonding.spin=cov::Spin::Alpha;
    nonbonding.coefficients={
        static_cast<float>(1.0/r2),0.0f,static_cast<float>(-1.0/r2)
    };

    cov::MolecularOrbital antibonding;
    antibonding.energy_hartree=0.4;
    antibonding.occupation=0.0f;
    antibonding.spin=cov::Spin::Alpha;
    antibonding.coefficients={
        static_cast<float>(1.0/r6),
        static_cast<float>(-2.0/r6),
        static_cast<float>(1.0/r6),
    };

    wf.orbitals={bonding,nonbonding,antibonding};
    wf.ao_overlap={1,0,0,0,1,0,0,0,1};
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;
    wf.total_density_packed=cov::reconstruct_total_density_packed(wf);
    wf.total_density_provenance=cov::DataProvenance::Derived;
    return wf;
}

bool check(const bool four_electron) {
    auto wf=make_three_centre(four_electron);
    cov::derive_bond_and_multicentre_analysis(wf);
    if (wf.multicentre_assignments.size()!=1u) {
        std::cerr << "expected exactly one strict multicentre assignment\n";
        return false;
    }
    const auto& assignment=wf.multicentre_assignments.front();
    const auto expected=four_electron
        ?cov::MulticentreKind::ThreeCentreFourElectron
        :cov::MulticentreKind::ThreeCentreTwoElectron;
    const double expected_e=four_electron?4.0:2.0;
    if (assignment.kind!=expected ||
        std::abs(assignment.electron_count-expected_e)>1.0e-6 ||
        assignment.atoms.size()!=3u || assignment.orbitals.size()!=3u ||
        assignment.provenance!=cov::DataProvenance::Derived ||
        assignment.confidence<0.70) {
        std::cerr << "strict multicentre assignment contents regression\n";
        return false;
    }
    return true;
}

} // namespace

int main() {
    if (!check(false) || !check(true)) return EXIT_FAILURE;

    // Separate spin blocks are deliberately not collapsed into an electron
    // count until a spin-resolved active-space implementation exists.
    auto open_shell=make_three_centre(false);
    auto beta=open_shell.orbitals.front();
    beta.spin=cov::Spin::Beta;
    beta.occupation=1.0f;
    open_shell.orbitals.push_back(beta);
    cov::derive_bond_and_multicentre_analysis(open_shell);
    if (!open_shell.multicentre_assignments.empty()) {
        std::cerr << "open-shell multicentre analysis should remain unclassified\n";
        return EXIT_FAILURE;
    }

    std::cout << "multicentre classification smoke test passed\n";
    return EXIT_SUCCESS;
}
