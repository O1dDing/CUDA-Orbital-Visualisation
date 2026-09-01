#include "cov/orbital_symmetry.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

cov::Atom atom(const char* symbol, int z, double x, double y, double zc) {
    return {symbol,z,x,y,zc};
}

void identity_overlap(cov::Wavefunction& wf) {
    wf.ao_overlap.assign(static_cast<std::size_t>(wf.basis_count)*wf.basis_count,0.0);
    for (std::size_t i=0;i<wf.basis_count;++i) wf.ao_overlap[i*wf.basis_count+i]=1.0;
    wf.ao_overlap_provenance=cov::DataProvenance::Derived;
}

void add_s_shell(cov::Wavefunction& wf, std::size_t atom_index) {
    const std::uint32_t primitive_offset=static_cast<std::uint32_t>(wf.primitives.size());
    wf.primitives.push_back({1.0f,1.0f});
    cov::Shell shell;
    shell.atom_index=static_cast<std::uint32_t>(atom_index);
    shell.primitive_offset=primitive_offset;
    shell.primitive_count=1;
    shell.basis_offset=wf.basis_count;
    shell.angular_momentum=0;
    shell.pure=0;
    wf.shells.push_back(shell);
    ++wf.basis_count;
}

void add_p_shell(cov::Wavefunction& wf, std::size_t atom_index) {
    const std::uint32_t primitive_offset=static_cast<std::uint32_t>(wf.primitives.size());
    wf.primitives.push_back({1.0f,1.0f});
    cov::Shell shell;
    shell.atom_index=static_cast<std::uint32_t>(atom_index);
    shell.primitive_offset=primitive_offset;
    shell.primitive_count=1;
    shell.basis_offset=wf.basis_count;
    shell.angular_momentum=1;
    shell.pure=0;
    wf.shells.push_back(shell);
    wf.basis_count+=3;
}

void add_mo(cov::Wavefunction& wf, double energy, std::vector<double> c) {
    cov::MolecularOrbital mo;
    mo.energy_hartree=energy;
    mo.spin=cov::Spin::Alpha;
    mo.occupation=0.0f;
    for (double value:c) mo.coefficients.push_back(static_cast<float>(value));
    wf.orbitals.push_back(std::move(mo));
}

bool all_label(const cov::Wavefunction& wf,
               std::size_t begin,
               std::size_t end,
               const std::string& expected) {
    for (std::size_t i=begin;i<end;++i) {
        if (wf.orbitals[i].symmetry!=expected ||
            wf.orbitals[i].symmetry_provenance!=cov::DataProvenance::Derived) {
            std::cerr << "expected " << expected << " at MO " << i
                      << ", got '" << wf.orbitals[i].symmetry << "'\n";
            return false;
        }
    }
    return true;
}

cov::Wavefunction d3h_ring() {
    cov::Wavefunction wf;
    const double r=1.6;
    for (int i=0;i<3;++i) {
        const double a=2.0*pi*static_cast<double>(i)/3.0;
        wf.atoms.push_back(atom("H",1,r*std::cos(a),r*std::sin(a),0.0));
        add_s_shell(wf,static_cast<std::size_t>(i));
    }
    const double s3=std::sqrt(3.0), s6=std::sqrt(6.0), s2=std::sqrt(2.0);
    add_mo(wf,-1.0,{1.0/s3,1.0/s3,1.0/s3});
    add_mo(wf,-0.4,{2.0/s6,-1.0/s6,-1.0/s6});
    add_mo(wf,-0.4,{0.0,1.0/s2,-1.0/s2});
    identity_overlap(wf);
    return wf;
}

cov::Wavefunction d5h_ring() {
    cov::Wavefunction wf;
    const double r=2.0;
    for (int i=0;i<5;++i) {
        const double a=2.0*pi*static_cast<double>(i)/5.0;
        wf.atoms.push_back(atom("C",6,r*std::cos(a),r*std::sin(a),0.0));
        add_s_shell(wf,static_cast<std::size_t>(i));
    }
    std::vector<double> a1(5,1.0/std::sqrt(5.0));
    std::vector<double> c1(5),s1(5),c2(5),s2(5);
    const double norm=std::sqrt(2.0/5.0);
    for (int i=0;i<5;++i) {
        const double a=2.0*pi*static_cast<double>(i)/5.0;
        c1[i]=norm*std::cos(a); s1[i]=norm*std::sin(a);
        c2[i]=norm*std::cos(2.0*a); s2[i]=norm*std::sin(2.0*a);
    }
    add_mo(wf,-1.0,a1);
    add_mo(wf,-0.5,c1); add_mo(wf,-0.5,s1);
    add_mo(wf, 0.1,c2); add_mo(wf, 0.1,s2);
    identity_overlap(wf);
    return wf;
}

cov::Wavefunction td_ligand_sigma() {
    cov::Wavefunction wf;
    wf.atoms={
        atom("Cl",17, 1, 1, 1), atom("Cl",17, 1,-1,-1),
        atom("Cl",17,-1, 1,-1), atom("Cl",17,-1,-1, 1),
    };
    for (std::size_t i=0;i<4;++i) add_s_shell(wf,i);
    add_mo(wf,-1.0,{0.5,0.5,0.5,0.5});
    add_mo(wf,-0.2,{1/std::sqrt(2.0),-1/std::sqrt(2.0),0,0});
    add_mo(wf,-0.2,{1/std::sqrt(6.0),1/std::sqrt(6.0),-2/std::sqrt(6.0),0});
    add_mo(wf,-0.2,{1/std::sqrt(12.0),1/std::sqrt(12.0),1/std::sqrt(12.0),-3/std::sqrt(12.0)});
    identity_overlap(wf);
    return wf;
}

cov::Wavefunction oh_central_p() {
    cov::Wavefunction wf;
    wf.atoms={
        atom("Ti",22,0,0,0),
        atom("F",9, 2,0,0),atom("F",9,-2,0,0),
        atom("F",9,0, 2,0),atom("F",9,0,-2,0),
        atom("F",9,0,0, 2),atom("F",9,0,0,-2),
    };
    add_p_shell(wf,0);
    add_mo(wf,-0.3,{1,0,0});
    add_mo(wf,-0.3,{0,1,0});
    add_mo(wf,-0.3,{0,0,1});
    identity_overlap(wf);
    return wf;
}

} // namespace

int main() {
    {
        auto wf=d3h_ring();
        const auto r=cov::derive_orbital_symmetry(wf);
        if (r.point_group!="D3h" || r.orbitals_labelled!=3u ||
            !all_label(wf,0,1,"A1'") || !all_label(wf,1,3,"E1'")) {
            std::cerr << "D3h derived-irrep regression\n";
            return EXIT_FAILURE;
        }
    }
    {
        auto wf=d5h_ring();
        const auto r=cov::derive_orbital_symmetry(wf);
        if (r.point_group!="D5h" || r.orbitals_labelled!=5u ||
            !all_label(wf,0,1,"A1'") || !all_label(wf,1,3,"E1'") ||
            !all_label(wf,3,5,"E2'")) {
            std::cerr << "D5h derived-irrep regression\n";
            return EXIT_FAILURE;
        }
    }
    {
        auto wf=td_ligand_sigma();
        const auto r=cov::derive_orbital_symmetry(wf);
        if (r.point_group!="Td" || r.orbitals_labelled!=4u ||
            !all_label(wf,0,1,"A1") || !all_label(wf,1,4,"T2")) {
            std::cerr << "Td derived-irrep regression\n";
            return EXIT_FAILURE;
        }
    }
    {
        auto wf=oh_central_p();
        const auto r=cov::derive_orbital_symmetry(wf);
        if (r.point_group!="Oh" || r.orbitals_labelled!=3u ||
            !all_label(wf,0,3,"T1u")) {
            std::cerr << "Oh AO-angular transform/derived-irrep regression"
                      << " point_group=" << r.point_group
                      << " groups_examined=" << r.groups_examined
                      << " groups_labelled=" << r.groups_labelled
                      << " orbitals_labelled=" << r.orbitals_labelled
                      << " worst_retention=" << r.worst_subspace_retention << '\n';
            for (std::size_t i=0;i<wf.orbitals.size();++i) {
                std::cerr << "  MO " << i << " label='" << wf.orbitals[i].symmetry
                          << "' provenance=" << static_cast<int>(wf.orbitals[i].symmetry_provenance)
                          << '\n';
            }
            return EXIT_FAILURE;
        }
    }
    {
        auto wf=d3h_ring();
        wf.orbitals[0].symmetry="producer";
        wf.orbitals[0].symmetry_provenance=cov::DataProvenance::Producer;
        (void)cov::derive_orbital_symmetry(wf);
        if (wf.orbitals[0].symmetry!="producer" ||
            wf.orbitals[0].symmetry_provenance!=cov::DataProvenance::Producer) {
            std::cerr << "producer MO symmetry was overwritten\n";
            return EXIT_FAILURE;
        }
    }
    std::cout << "orbital symmetry smoke test passed\n";
    return EXIT_SUCCESS;
}
