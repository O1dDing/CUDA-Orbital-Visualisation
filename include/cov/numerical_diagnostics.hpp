#pragma once

#include "cov/model.hpp"

#include <cmath>
#include <iomanip>
#include <ostream>
#include <string_view>

namespace cov {
inline const char* numerical_status_name(const NumericalStatus status) noexcept {
    switch (status) {
        case NumericalStatus::Available:return "available";
        case NumericalStatus::MissingInput:return "missing-input";
        case NumericalStatus::InvalidInput:return "invalid-input";
        case NumericalStatus::Failed:return "failed";
        default:return "not-computed";
    }
}

namespace numerical_json {
inline void number(std::ostream& out,const double value) {
    if (std::isfinite(value)) out<<std::setprecision(17)<<value;
    else out<<"null";
}
inline void string(std::ostream& out,const std::string_view value) {
    out<<'"';
    constexpr char hex[]="0123456789abcdef";
    for (const unsigned char c:value) {
        if (c=='"' || c=='\\') out<<'\\'<<static_cast<char>(c);
        else if (c<32) out<<"\\u00"<<hex[c/16]<<hex[c%16];
        else out<<static_cast<char>(c);
    }
    out<<'"';
}
inline void metric(std::ostream& out,const AoMetricDiagnostics& data) {
    out<<"{\"status\":";string(out,numerical_status_name(data.status));
    out<<",\"dimension\":"<<data.dimension<<",\"finite\":"<<(data.finite?"true":"false")
       <<",\"positive_semidefinite\":"<<(data.positive_semidefinite?"true":"false")
       <<",\"symmetry_error\":";number(out,data.symmetry_error);
    out<<",\"minimum_eigenvalue\":";number(out,data.minimum_eigenvalue);
    out<<",\"maximum_eigenvalue\":";number(out,data.maximum_eigenvalue);
    out<<",\"condition_number\":";number(out,data.condition_number);
    out<<",\"eigen_residual\":";number(out,data.eigen_residual);
    out<<",\"numerical_rank\":"<<data.numerical_rank<<",\"numerical_rank_tolerance\":";
    number(out,data.numerical_rank_tolerance);
    out<<",\"detail\":";string(out,data.detail);
    out<<",\"spin_blocks\":[";
    bool first=true;
    for (const auto& block:data.orbital_blocks) {
        if (!first) out<<',';
        first=false;
        out<<"{\"spin\":";string(out,block.spin==Spin::Alpha?"alpha":"beta");
        out<<",\"orbital_count\":"<<block.orbital_count<<",\"status\":";
        string(out,numerical_status_name(block.status));
        out<<",\"maximum_diagonal_error\":";number(out,block.maximum_diagonal_error);
        out<<",\"maximum_off_diagonal_error\":";number(out,block.maximum_off_diagonal_error);
        out<<",\"orthonormality_tolerance\":";number(out,block.orthonormality_tolerance);
        out<<",\"detail\":";string(out,block.detail);out<<'}';
    }
    out<<"]}";
}
}

// Shared serialization of measured production data; it does not recompute or
// substitute independent science. Unknown/nonfinite numbers use JSON null.
inline void write_numerical_diagnostics_json(std::ostream& out,const Wavefunction& wf) {
    out<<"{\"schema\":1,\"scope\":\"full AO basis and every supplied spin MO block\","
         "\"units\":\"dimensionless inner products\",\"source\":\"analytic basis integrals\",\"basis\":";
    numerical_json::metric(out,wf.ao_metric_diagnostics);
    out<<",\"producer\":";numerical_json::metric(out,wf.producer_ao_metric_diagnostics);
    out<<",\"producer_basis_comparison\":{\"status\":";
    numerical_json::string(out,numerical_status_name(wf.producer_ao_overlap_basis_status));
    out<<",\"maximum_absolute_error\":";numerical_json::number(out,wf.producer_ao_overlap_basis_error);
    out<<",\"tolerance\":";numerical_json::number(out,wf.producer_ao_overlap_basis_tolerance);
    out<<"}}";
}
}
