#include "cov/local_angular_generators.hpp"
#include <iomanip>
#include <iostream>
int main() {
    std::cout<<std::setprecision(17);
    int id=0,L=0,pure=0,l=0;
    while(std::cin>>id>>L>>pure>>l) {
        std::array<double,9> rotation{};
        for(auto& v:rotation)if(!(std::cin>>v))return 2;
        const auto x=cov::local_angular_generators(L,pure!=0,l,rotation);
        std::cout<<"{\"id\":"<<id<<",\"available\":"<<(x.available?"true":"false")
            <<",\"basis_count\":"<<x.basis_count<<",\"columns\":"<<x.columns
            <<",\"rotation_error\":"<<x.rotation_orthogonality_error
            <<",\"polynomial_error\":"<<x.polynomial_reconstruction_error<<",\"coefficients\":[";
        bool first=true;for(double v:x.coefficients){if(!first)std::cout<<',';first=false;std::cout<<v;}
        std::cout<<"]}\n";
    }
    return std::cin.eof()?0:2;
}
