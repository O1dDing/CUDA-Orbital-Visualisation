        zlib.push_back(static_cast<std::uint8_t>((nlen >> 8u) & 0xffu));
        zlib.insert(zlib.end(),
                    raw.begin() + static_cast<std::ptrdiff_t>(pos),
                    raw.begin() + static_cast<std::ptrdiff_t>(pos + block));
        pos += block;
    }
    push_u32_be(zlib, adler32(raw.data(), raw.size()));

    std::vector<std::uint8_t> png = {137,80,78,71,13,10,26,10};
    std::vector<std::uint8_t> ihdr;
    push_u32_be(ihdr, static_cast<std::uint32_t>(canvas.width));
    push_u32_be(ihdr, static_cast<std::uint32_t>(canvas.height));
    ihdr.push_back(8); ihdr.push_back(6); ihdr.push_back(0);
    ihdr.push_back(0); ihdr.push_back(0);
    append_png_chunk(png, "IHDR", ihdr);
    append_png_chunk(png, "IDAT", zlib);
    append_png_chunk(png, "IEND", {});

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        if (error) *error = "Unable to open PNG output path";
        return false;
    }
    file.write(reinterpret_cast<const char*>(png.data()),
               static_cast<std::streamsize>(png.size()));
    if (!file) {
        if (error) *error = "Failed while writing PNG";
        return false;
    }
    return true;
}

} // namespace

const char* energy_axis_mode_name(const EnergyAxisMode mode) noexcept {
    return mode == EnergyAxisMode::Linear ? "linear" : "nonlinear-focus";
}

const char* energy_transform_name(const EnergyAxisMode mode) noexcept {
    return mode == EnergyAxisMode::Linear ? "identity" : "adaptive-log-gap-v2";
}

EnergyTransform build_energy_transform(const std::vector<double>& energies_hartree,
                                       const EnergyAxisMode mode,
                                       const double minimum_gap_weight) {
    EnergyTransform transform;
    transform.mode = mode;
    transform.minimum_gap_weight = std::max(0.0, minimum_gap_weight);

    std::vector<double> energies;
    energies.reserve(energies_hartree.size());
    for (const double energy : energies_hartree) {
        if (std::isfinite(energy)) energies.push_back(energy);
    }
    if (energies.empty()) {
        transform.knots = {{0.0, 0.0}, {1.0, 1.0}};
        transform.focus_hartree = 0.5;
        transform.scale_hartree = 1.0;
        return transform;
    }

    std::sort(energies.begin(), energies.end());
    energies.erase(std::unique(energies.begin(), energies.end(),
                               [](const double a, const double b) {
                                   return std::abs(a - b) <= kCoordinateEpsilon;
                               }),
                   energies.end());

    transform.focus_hartree = energies[energies.size() / 2];
    transform.scale_hartree = median_positive_gap(energies);

    if (energies.size() == 1) {
        transform.knots.push_back({energies.front(), 0.0});
        return transform;
    }

    transform.knots.reserve(energies.size());
    double coordinate = mode == EnergyAxisMode::Linear ? energies.front() : 0.0;
    transform.knots.push_back({energies.front(), coordinate});

    for (std::size_t i = 1; i < energies.size(); ++i) {
        const double gap = energies[i] - energies[i - 1];
        if (mode == EnergyAxisMode::Linear) {
            coordinate = energies[i];
        } else {
            const double weight =
                std::log1p(gap / std::max(1.0e-12, transform.scale_hartree)) +
                transform.minimum_gap_weight;
            coordinate += std::max(kCoordinateEpsilon, weight);
        }
        transform.knots.push_back({energies[i], coordinate});
    }
    return transform;
}

double energy_display_coordinate(const double energy_hartree,
                                 const EnergyTransform& transform) noexcept {
    if (!std::isfinite(energy_hartree)) return 0.0;
    if (transform.knots.empty()) {
        if (transform.mode == EnergyAxisMode::Linear) return energy_hartree;
        const double scale = std::max(1.0e-12, std::abs(transform.scale_hartree));
        return std::asinh((energy_hartree - transform.focus_hartree) / scale);
    }
    if (transform.knots.size() == 1) return transform.knots.front().coordinate;
    if (energy_hartree <= transform.knots.front().energy_hartree) {
        const auto& a = transform.knots[0];
        const auto& b = transform.knots[1];
        return interpolate(energy_hartree, a.energy_hartree, b.energy_hartree,
                           a.coordinate, b.coordinate);
    }
    if (energy_hartree >= transform.knots.back().energy_hartree) {
        const auto& a = transform.knots[transform.knots.size() - 2];
        const auto& b = transform.knots.back();
        return interpolate(energy_hartree, a.energy_hartree, b.energy_hartree,
                           a.coordinate, b.coordinate);
    }
    const auto it = std::upper_bound(
        transform.knots.begin(), transform.knots.end(), energy_hartree,
        [](const double value, const EnergyAxisKnot& knot) {
            return value < knot.energy_hartree;
        });
    const std::size_t hi = static_cast<std::size_t>(
        std::distance(transform.knots.begin(), it));
    const std::size_t lo = hi - 1;
    const auto& a = transform.knots[lo];
    const auto& b = transform.knots[hi];
    return interpolate(energy_hartree, a.energy_hartree, b.energy_hartree,
                       a.coordinate, b.coordinate);
}

double energy_from_display_coordinate(const double coordinate,
                                      const EnergyTransform& transform) noexcept {
    if (!std::isfinite(coordinate)) return 0.0;
    if (transform.knots.empty()) {
        if (transform.mode == EnergyAxisMode::Linear) return coordinate;
        const double scale = std::max(1.0e-12, std::abs(transform.scale_hartree));
        return transform.focus_hartree + scale * std::sinh(coordinate);
    }
    if (transform.knots.size() == 1) return transform.knots.front().energy_hartree;
    if (coordinate <= transform.knots.front().coordinate) {
        const auto& a = transform.knots[0];
        const auto& b = transform.knots[1];
        return interpolate(coordinate, a.coordinate, b.coordinate,
                           a.energy_hartree, b.energy_hartree);
    }
    if (coordinate >= transform.knots.back().coordinate) {
        const auto& a = transform.knots[transform.knots.size() - 2];
        const auto& b = transform.knots.back();
        return interpolate(coordinate, a.coordinate, b.coordinate,
                           a.energy_hartree, b.energy_hartree);
    }
    const auto it = std::upper_bound(
        transform.knots.begin(), transform.knots.end(), coordinate,
        [](const double value, const EnergyAxisKnot& knot) {
            return value < knot.coordinate;
        });
    const std::size_t hi = static_cast<std::size_t>(
        std::distance(transform.knots.begin(), it));
    const std::size_t lo = hi - 1;
    const auto& a = transform.knots[lo];
    const auto& b = transform.knots[hi];
    return interpolate(coordinate, a.coordinate, b.coordinate,
                       a.energy_hartree, b.energy_hartree);
}

const char* annotation_source_name(const AnnotationSource source) noexcept {
    switch (source) {
        case AnnotationSource::Direct: return "direct";
        case AnnotationSource::ParsedLabel: return "parsed-label";
        case AnnotationSource::Heuristic: return "heuristic";
        default: return "unavailable";
    }
}

const char* bonding_class_name(const BondingClass value) noexcept {
    switch (value) {
        case BondingClass::Bonding: return "bonding";
        case BondingClass::Nonbonding: return "nonbonding";
        case BondingClass::Antibonding: return "antibonding";
        default: return "unclassified";
