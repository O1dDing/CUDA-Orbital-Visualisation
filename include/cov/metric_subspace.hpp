#pragma once
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace cov {
enum class MetricSubspaceStatus { Available, MissingInput, InvalidInput, Failed };
struct MetricSubspaceDiagnostics {
    std::size_t supplied_columns=0;
    std::size_t numerical_rank=0;
    std::size_t unresolved_norm_columns=0;
    double rank_cutoff=0;
    double gram_inverse_residual=0;
    double orthonormality_residual=0;
    double retained_condition_number=0;
};
struct MetricSubspaceOverlap {
    MetricSubspaceStatus status=MetricSubspaceStatus::MissingInput;
    std::size_t basis_dimension=0;
    std::size_t metric_numerical_rank=0;
    bool basis_norm_preconditioned=false;
    // Spectral diagnostics below refer to the congruently transformed metric
    // with each individual AO normalized. Original S and coefficient data are
    // not changed. X and C undergo the same inverse basis transformation.
    double metric_rank_cutoff=0;
    double metric_symmetry_error=0;
    double metric_eigen_residual=0;
    double metric_minimum_eigenvalue=std::numeric_limits<double>::quiet_NaN();
    double subspace_overlap_trace=std::numeric_limits<double>::quiet_NaN();
    double mean_orbital_fraction=std::numeric_limits<double>::quiet_NaN();
    MetricSubspaceDiagnostics reference;
    MetricSubspaceDiagnostics orbitals;
    std::string detail;
};

class MetricSubspaceContext;
class MetricPreparedSubspace {
public:
    MetricPreparedSubspace()=default;
    MetricSubspaceDiagnostics diagnostics() const;
    MetricSubspaceStatus status() const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
    friend class MetricSubspaceContext;
};

// The context owns its normalized metric and source norm factors. Each space
// owns its prepared columns. Mutating a caller's input vector cannot affect a
// prepared result, and spaces from different contexts cannot be silently mixed.
class MetricSubspaceContext {
public:
    MetricSubspaceContext(std::size_t dimension,const std::vector<double>& metric);
    MetricSubspaceOverlap diagnostics() const;
    MetricPreparedSubspace prepare_space(std::size_t columns,const std::vector<double>& coefficients) const;
    MetricSubspaceOverlap compare(const MetricPreparedSubspace& reference,const MetricPreparedSubspace& orbitals) const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};
// Row-major S, X and C. The reported trace is tr(P_X P_C), and its mean is
// divided by the represented rank of C. It is invariant under well-conditioned
// changes of basis inside either space. The unweighted amount is not an
// electron population, nor an additive partition over overlapping centres.
// All numerical rank decisions and residuals are returned. A zero-rank target
// has trace zero and an undefined mean, which is preserved as NaN.
MetricSubspaceOverlap metric_subspace_overlap(
    std::size_t dimension,const std::vector<double>& metric,
    std::size_t reference_columns,const std::vector<double>& reference,
    std::size_t orbital_columns,const std::vector<double>& orbitals);
}
