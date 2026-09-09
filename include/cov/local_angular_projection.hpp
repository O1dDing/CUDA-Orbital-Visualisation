#pragma once
#include "cov/metric_subspace.hpp"
#include "cov/model.hpp"
#include <array>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace cov {
struct LocalAngularProjectionComponent {
    int angular_degree = 0;
    std::size_t component = 0; // Real spherical order, including its phase convention.
    std::vector<std::size_t> source_shell_indices;
    MetricSubspaceOverlap metric;
    double fraction_of_local_projection = std::numeric_limits<double>::quiet_NaN();
};
struct LocalAngularSpinProjection {
    Spin spin = Spin::Alpha;
    std::vector<std::size_t> orbital_indices;
    MetricSubspaceOverlap centre;
    std::vector<LocalAngularProjectionComponent> components;
    double angular_partition_residual = std::numeric_limits<double>::quiet_NaN();
};
struct LocalAngularProjection {
    MetricSubspaceStatus status = MetricSubspaceStatus::MissingInput;
    std::size_t atom_index = 0;
    std::array<double,9> rotation_reference_to_input{};
    std::vector<std::size_t> source_shell_indices;
    std::vector<LocalAngularSpinProjection> spins;
    std::size_t represented_spin_orbital_rank = 0;
    double centre_projection_trace = std::numeric_limits<double>::quiet_NaN();
    double centre_mean_fraction = std::numeric_limits<double>::quiet_NaN();
    std::array<double,25> component_projection_traces = [] {
        std::array<double,25> values;
        values.fill(std::numeric_limits<double>::quiet_NaN());
        return values;
    }(); // l*l + component; unavailable values are not zero.
    double angular_partition_residual = std::numeric_limits<double>::quiet_NaN();
    std::string detail;
};

// Owns an immutable snapshot of the metric, AO definitions and MO coefficients.
// Prepare once outside group loops. Alpha and beta are separate orthogonal spin
// blocks; their traces and represented ranks are added after spatial projection.
class LocalAngularProjectionWorkspace {
public:
    LocalAngularProjectionWorkspace(const Wavefunction&,std::size_t atom_index,
        const std::array<double,9>& rotation_reference_to_input);
    LocalAngularProjection project(std::span<const std::size_t> orbital_indices) const;
private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
};
}
