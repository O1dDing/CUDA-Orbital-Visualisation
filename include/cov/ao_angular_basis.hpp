#pragma once
#include <array>
#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace cov::ao_angular {
// One internal AO order, phase and angular normalization for analytic
// integrals and local angular subspaces. This is not a sampled fit.
using Powers=std::array<int,3>;
struct Term { Powers powers{}; double coefficient=0.0; };
using Polynomial=std::vector<Term>;
inline double factorial(const int n) {
    double value=1.0;
    for (int i=2;i<=n;++i) value*=i;
    return value;
}
inline double odd_factorial(const int n) {
    double value=1.0;
    for (int i=n;i>1;i-=2) value*=i;
    return value;
}

inline std::vector<Powers> cartesian_components(const int l) {
    switch (l) {
        case 0:return {{0,0,0}};
        case 1:return {{1,0,0},{0,1,0},{0,0,1}};
        case 2:return {{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}};
        case 3:return {{3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},
                       {2,0,1},{1,0,2},{0,1,2},{0,2,1},{1,1,1}};
        case 4:return {{4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},
                       {1,3,0},{0,3,1},{1,0,3},{0,1,3},{2,2,0},
                       {2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}};
        default:throw std::invalid_argument("Analytic AO integral requires supported s through g shells");
    }
}

inline Polynomial solid_harmonic(const int l,const int component) {
    // Expand normalized real r^l Y_lm in the same Condon--Shortley convention
    // used by the production evaluators. Expand (x+i y)^m and r^(2k), without
    // fitting sampled values or using MO coefficients to define the basis.
    const int m=(component+1)/2;
    const bool imaginary=component>0 && component%2==0;
    const double normalization=std::sqrt((2*l+1)/(4*std::numbers::pi)*
        factorial(l-m)/factorial(l+m))*(m==0?1.0:std::sqrt(2.0));
    std::map<Powers,double> terms;
    for (int k=0;k<=(l-m)/2;++k) {
        const double legendre=((m+k)%2?-1.0:1.0)*factorial(2*l-2*k)/
            (std::pow(2.0,l)*factorial(k)*factorial(l-k)*factorial(l-m-2*k));
        for (int y=0;y<=m;++y) {
            if ((y%2!=0)!=imaginary) continue;
            const double complex_sign=(y/2)%2?-1.0:1.0;
            const double complex_coefficient=complex_sign*factorial(m)/(factorial(y)*factorial(m-y));
            for (int rx=0;rx<=k;++rx) {
                for (int ry=0;ry<=k-rx;++ry) {
                    const int rz=k-rx-ry;
                    const double radial=factorial(k)/(factorial(rx)*factorial(ry)*factorial(rz));
                    terms[{m-y+2*rx,y+2*ry,l-m-2*k+2*rz}]+=
                        normalization*legendre*complex_coefficient*radial;
                }
            }
        }
    }
    Polynomial result;
    for (const auto& [powers,coefficient]:terms) {
        if (coefficient!=0.0) result.push_back({powers,coefficient});
    }
    return result;
}

} // namespace cov::ao_angular
