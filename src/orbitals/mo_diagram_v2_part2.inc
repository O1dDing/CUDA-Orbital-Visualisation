            return std::pair{centres, electrons};
        }
    }
    return std::nullopt;
}

std::string annotation_summary(const OrbitalAnnotation& annotation) {
    std::string out;
    if (annotation.family != "unavailable") out += annotation.family;
    if (annotation.bonding_class != BondingClass::Unclassified) {
        if (!out.empty()) out += " / ";
        out += bonding_class_name(annotation.bonding_class);
    }
    if (annotation.multicentre.available) {
        if (!out.empty()) out += " / ";
        out += annotation.multicentre.label;
    }
    if (annotation.delocalised_pi.available) {
        if (!out.empty()) out += " / ";
        out += annotation.delocalised_pi.label;
    }
    return out;
}

int group_member_index(const std::string& label) noexcept {
    const auto dash = label.find('-');
    if (dash == std::string::npos || dash + 1 >= label.size()) return 0;
    const char suffix = label[dash + 1];
    if (suffix < 'a' || suffix > 'z') return 0;
    return static_cast<int>(suffix - 'a');
}

std::size_t group_base_number(const OrbitalMetadata& metadata) noexcept {
    const int member = group_member_index(metadata.display_label);
    if (metadata.degeneracy_size <= 1 || member < 0) return metadata.raw_mo_number;
    const auto m = static_cast<std::size_t>(member);
    return metadata.raw_mo_number >= m ? metadata.raw_mo_number - m
                                       : metadata.raw_mo_number;
}

double group_layout_energy(const std::vector<OrbitalMetadata>& metadata,
                           const std::size_t index) {
    if (index >= metadata.size()) return 0.0;
    const auto& item = metadata[index];
    if (item.degeneracy_size <= 1) return item.energy_hartree;
    const std::size_t base = group_base_number(item);
    if (base == 0) return item.energy_hartree;
    const std::size_t begin = base - 1;
    const std::size_t end = std::min(metadata.size(), begin + item.degeneracy_size);
    if (begin >= end) return item.energy_hartree;
    double sum = 0.0;
    for (std::size_t i = begin; i < end; ++i) sum += metadata[i].energy_hartree;
    return sum / static_cast<double>(end - begin);
}

struct PositionedLevel {
    const MODiagramLevel* level = nullptr;
    double x = 0.0;
    double y = 0.0;
    double label_y = 0.0;
};

double coordinate_min(const EnergyTransform& transform) noexcept {
    if (transform.knots.empty()) return 0.0;
    return transform.knots.front().coordinate;
}

double coordinate_max(const EnergyTransform& transform) noexcept {
    if (transform.knots.empty()) return 1.0;
    return transform.knots.back().coordinate;
}

double map_energy_y(const double layout_energy_hartree,
                    const EnergyTransform& transform,
                    const int height,
                    const int top = 110,
                    const int bottom = 70) {
    const double c = energy_display_coordinate(layout_energy_hartree, transform);
    const double c0 = coordinate_min(transform);
    const double c1 = coordinate_max(transform);
    const double t = (c - c0) / std::max(kCoordinateEpsilon, c1 - c0);
    return static_cast<double>(height - bottom) -
           t * static_cast<double>(height - top - bottom);
}

double energy_at_axis_fraction(const double fraction,
                               const EnergyTransform& transform) noexcept {
    const double c0 = coordinate_min(transform);
    const double c1 = coordinate_max(transform);
    return energy_from_display_coordinate(c0 + fraction * (c1 - c0), transform);
}

std::vector<PositionedLevel> position_levels(const MODiagramData& data,
                                             const int width,
                                             const int height) {
    std::vector<PositionedLevel> positioned;
    positioned.reserve(data.levels.size());
    const double centre = static_cast<double>(width) * 0.50;

    // Levels with the same layout-energy coordinate (declared degeneracy or an
    // accidental exact coincidence) share y but are spread horizontally.
    std::map<long long, std::vector<std::size_t>> coincident;
    constexpr double kEnergyBucket = 1.0e12;
    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        const long long key = static_cast<long long>(
            std::llround(data.levels[i].layout_energy_hartree * kEnergyBucket));
        coincident[key].push_back(i);
    }

    std::vector<double> x(data.levels.size(), centre);
    for (const auto& [_, ids] : coincident) {
        const double spacing = std::clamp(static_cast<double>(width) * 0.105, 92.0, 138.0);
        for (std::size_t member = 0; member < ids.size(); ++member) {
            const double offset =
                (static_cast<double>(member) -
                 (static_cast<double>(ids.size()) - 1.0) * 0.5) * spacing;
            x[ids[member]] = centre + offset;
        }
    }

    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        const auto& level = data.levels[i];
        const double y = map_energy_y(level.layout_energy_hartree,
                                      data.energy_transform,
                                      height);
        positioned.push_back({&level, x[i], y, y - 7.0});
    }

    // Only text labels move for collision avoidance. The actual level line y
    // remains exactly the energy-axis coordinate.
    std::vector<std::size_t> ids(positioned.size());
    std::iota(ids.begin(), ids.end(), 0);
    std::sort(ids.begin(), ids.end(), [&](const std::size_t a, const std::size_t b) {
        if (positioned[a].label_y != positioned[b].label_y) {
            return positioned[a].label_y < positioned[b].label_y;
        }
        return positioned[a].x < positioned[b].x;
    });
    const double spacing = data.levels.size() > 28 ? 23.0 : 29.0;
    const double top_label = 108.0;
    const double bottom_label = static_cast<double>(height - 48);
    double previous = top_label - spacing;
    for (const std::size_t id : ids) {
        positioned[id].label_y = std::max(positioned[id].label_y, previous + spacing);
        previous = positioned[id].label_y;
    }
    const double overflow = previous - bottom_label;
    if (overflow > 0.0) {
        for (const std::size_t id : ids) positioned[id].label_y -= overflow;
    }
    double next = bottom_label + spacing;
    for (auto it = ids.rbegin(); it != ids.rend(); ++it) {
        positioned[*it].label_y = std::min(positioned[*it].label_y, next - spacing);
        next = positioned[*it].label_y;
    }
    return positioned;
}

struct Canvas {
    int width = 0;
    int height = 0;
    std::vector<std::uint8_t> rgba;

    Canvas(const int w, const int h)
        : width(std::max(1, w)),
          height(std::max(1, h)),
          rgba(static_cast<std::size_t>(std::max(1, w)) *
               static_cast<std::size_t>(std::max(1, h)) * 4u, 255u) {
        clear(250, 251, 253, 255);
    }

    void clear(const std::uint8_t r, const std::uint8_t g,
               const std::uint8_t b, const std::uint8_t a) {
        for (std::size_t i = 0; i < rgba.size(); i += 4) {
            rgba[i] = r; rgba[i + 1] = g; rgba[i + 2] = b; rgba[i + 3] = a;
        }
    }

    void pixel(const int x, const int y,
               const std::uint8_t r, const std::uint8_t g,
