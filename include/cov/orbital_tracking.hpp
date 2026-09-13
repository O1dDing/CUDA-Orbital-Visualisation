#pragma once

#include "cov/model.hpp"

#include <cstddef>
#include <vector>

namespace cov {

enum class OrbitalTrackingSource {
    Unavailable = 0,
    SMetricChemicalDescriptor,
};

struct OrbitalTrackingOptions {
    double degeneracy_tolerance_hartree = 1.0e-5;
    double minimum_similarity = 0.55;
    // A union descriptor is deliberately broader than an individual orbital
    // descriptor because symmetry lowering can redistribute support among
    // equivalent atoms. Keep a separate, still conservative threshold.
    double minimum_composite_similarity = 0.20;
    double ambiguity_margin = 0.06;
    double occupation_weight = 0.20;
    double maximum_average_occupation_change = 0.75;
    double energy_weight = 0.05;
    double energy_scale_hartree = 0.25;
    // A temporary split/recombine union must remain in the same local energy
    // neighbourhood. This prevents chemically broad atom/angular descriptors
    // from joining an occupied core group to an unrelated valence degeneracy.
    double composite_energy_window_hartree = 1.0;
    // Split components may interleave with another crossing level while still
    // belonging to the same local neighbourhood; they need not be consecutive
    // in instantaneous energy order, but their internal span stays bounded.
    double composite_component_span_hartree = 0.25;
    // Per-orbital occupancy must remain locally compatible when a temporary
    // union is formed; otherwise a broad descriptor can exchange an occupied
    // subspace with an unrelated virtual one.
    double composite_occupation_window = 0.50;
    // Hard complexity bounds for dense near-degenerate spectra. Candidate
    // generation is exact inside the retained local pool, then keeps the best
    // chemically admissible unions for conflict-component optimisation.
    std::size_t maximum_split_components = 5u;
    std::size_t maximum_local_component_pool = 16u;
    std::size_t maximum_composite_candidates_per_anchor = 8u;
    std::size_t maximum_conflict_component_candidates = 128u;
    std::size_t maximum_optimizer_states = 250000u;
    double maximum_optimizer_milliseconds = 250.0;
};

// A match describes the identity of an entire orbital subspace between two
// geometries.  Members of a degenerate group are intentionally not paired one
// by one because an arbitrary orthogonal rotation inside the group has no
// chemical meaning.
struct OrbitalSubspaceMatch {
    std::vector<std::size_t> from_members;
    std::vector<std::size_t> to_members;
    double similarity = 0.0;
    double score = 0.0;
    bool ambiguous = false;
    OrbitalTrackingSource source = OrbitalTrackingSource::Unavailable;
};

struct OrbitalTrackingResult {
    std::vector<OrbitalSubspaceMatch> matches;
    std::vector<std::vector<std::size_t>> unmatched_from;
    std::vector<std::vector<std::size_t>> unmatched_to;
    bool atom_mapping_compatible = false;
    std::size_t composite_candidates_considered = 0u;
    std::size_t composite_matches_selected = 0u;
    std::size_t composite_optimizer_states = 0u;
    std::size_t composite_fallback_components = 0u;
    bool composite_optimisation_truncated = false;
};

// The present FCHK model does not contain cross-geometry AO integrals.  This
// fallback therefore matches phase-invariant S-metric chemical descriptors
// (atom/element/angular support, occupation and only weakly energy) rather
// than raw coefficients. It is suitable for preserving identities through
// ordinary scans and explicitly reports its source; it never claims to be an
// exact cross-AO maximum-overlap calculation.
[[nodiscard]] OrbitalTrackingResult track_orbital_subspaces(
    const Wavefunction& from,
    const Wavefunction& to,
    const OrbitalTrackingOptions& options = {});

[[nodiscard]] const char* orbital_tracking_source_name(
    OrbitalTrackingSource source) noexcept;

} // namespace cov
