#include "cov/wavefunction_io.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

void prefix(std::ofstream& out,const std::string& label,const char type) {
    out << std::left << std::setw(40) << label << type;
}
void scalar_i(std::ofstream& out,const std::string& label,long long value) {
    prefix(out,label,'I'); out << std::right << std::setw(16) << value << '\n';
}
void array_i(std::ofstream& out,const std::string& label,const std::vector<long long>& v) {
    prefix(out,label,'I'); out << "   N=" << std::right << std::setw(12) << v.size() << '\n';
    for (std::size_t i=0;i<v.size();++i) {out << std::setw(12) << v[i]; if ((i+1u)%6u==0u||i+1u==v.size()) out << '\n';}
}
void array_r(std::ofstream& out,const std::string& label,const std::vector<double>& v) {
    prefix(out,label,'R'); out << "   N=" << std::right << std::setw(12) << v.size() << '\n';
    out << std::uppercase << std::scientific << std::setprecision(8);
    for (std::size_t i=0;i<v.size();++i) {out << std::setw(16) << v[i]; if ((i+1u)%5u==0u||i+1u==v.size()) out << '\n';}
}

} // namespace

int main() {
    const auto path=std::filesystem::temp_directory_path()/"cov_overlap_source.fchk";
    {
        std::ofstream out(path,std::ios::binary);
        if (!out) return EXIT_FAILURE;
        out << "Producer overlap regression\nSP        RHF STO-3G\n";
        scalar_i(out,"Number of atoms",2);
        scalar_i(out,"Number of alpha electrons",1);
        scalar_i(out,"Number of beta electrons",1);
        scalar_i(out,"Number of basis functions",2);
        scalar_i(out,"Number of independent functions",2);
        array_i(out,"Atomic numbers",{1,1});
        array_r(out,"Current cartesian coordinates",{0,0,-0.7,0,0,0.7});
        array_i(out,"Shell types",{0,0});
        array_i(out,"Number of primitives per shell",{1,1});
        array_i(out,"Shell to atom map",{1,2});
        array_r(out,"Primitive exponents",{1,1});
        array_r(out,"Contraction coefficients",{1,1});
        constexpr double s=0.2;
        const double cb=1.0/std::sqrt(2.0*(1.0+s));
        const double ca=1.0/std::sqrt(2.0*(1.0-s));
        array_r(out,"Alpha Orbital Energies",{-0.5,0.2});
        array_r(out,"Alpha MO coefficients",{cb,cb,ca,-ca});
        array_r(out,"Overlap Matrix",{1.0,s,1.0});
    }

    std::error_code ec;
    try {
        const auto wf=cov::parse_wavefunction(path);
        std::filesystem::remove(path,ec);
        if (wf.ao_overlap_provenance!=cov::DataProvenance::Producer ||
            wf.ao_overlap.size()!=4u ||
            std::abs(wf.ao_overlap[0]-1.0)>1.0e-10 ||
            std::abs(wf.ao_overlap[1]-0.2)>1.0e-10 ||
            std::abs(wf.ao_overlap[2]-0.2)>1.0e-10 ||
            std::abs(wf.ao_overlap[3]-1.0)>1.0e-10) {
            std::cerr << "producer FCHK overlap was not preserved\n";
            return EXIT_FAILURE;
        }
        std::cout << "FCHK producer overlap smoke test passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::filesystem::remove(path,ec);
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
