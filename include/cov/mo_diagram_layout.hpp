#pragma once

#include <span>
#include <vector>

namespace cov {

// A row's actual stroke, arrow, label or hit-target footprint in display units.
// Layout never changes y or interprets pixel proximity as physical degeneracy.
struct DiagramRowFootprint {
    double y = 0.0;
    double above = 0.0;
    double below = 0.0;
    double width = 0.0;
};

struct DiagramLaneLayout {
    std::vector<double> centre_x;
    double width = 0.0;
};

// Overlapping vertical footprints receive disjoint horizontal lanes. Separate
// vertical clusters are centred independently; lanes can be reused after a
// preceding footprint ends. The canvas grows when its minimum width cannot
// contain every lane, rather than shrinking strokes or hiding members.
[[nodiscard]] DiagramLaneLayout layout_diagram_lanes(
    std::span<const DiagramRowFootprint> rows,
    double minimum_width,
    double horizontal_gap,
    double vertical_gap = 0.0);

} // namespace cov
