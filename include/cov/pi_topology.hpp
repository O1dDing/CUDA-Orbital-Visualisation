#pragma once

#include "cov/model.hpp"

#include <array>
#include <cstdint>
#include <utility>
#include <vector>

namespace cov {

// Geometry only proposes compatible p directions. Electronic assignment is
// deliberately left to S-metric active-subspace analysis: a geometrical
// candidate must never be presented as proof of a pi bond by itself.
struct PiTopologyOptions {
    double linear_rank_tolerance = 0.08;
    double planar_rank_tolerance = 0.12;
    double common_direction_tolerance = 0.12;
    double adjacent_direction_cosine = 0.90;
    double cyclic_global_coherence = 0.86;
};

struct OrientedPiNetwork {
    std::vector<std::uint32_t> atoms;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    std::array<double, 3> representative_direction{0.0, 0.0, 1.0};
    double orientation_coherence = 0.0;
    double confidence = 0.0;
    bool cyclic = false;
    bool globally_oriented = false;
};

// bonded_pairs is the strong molecular graph, not a renderer distance graph.
// Hydrogens contribute to the local hybridisation frame but are never p
// centres. Separate orthogonal channels sharing an sp centre remain separate.
[[nodiscard]] std::vector<OrientedPiNetwork> infer_oriented_pi_networks(
    const Wavefunction& wavefunction,
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& bonded_pairs,
    const PiTopologyOptions& options = {});

} // namespace cov
