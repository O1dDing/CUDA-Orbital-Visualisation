namespace {

int bitmap_text_width(const std::string& text, const int scale) {
    int width = 0;
    for (const char c : text) width += (c == ' ' ? 4 : 6) * scale;
    return width;
}

void draw_symmetry_png(Canvas& canvas,
                       const int centre_x,
                       const int baseline_y,
                       const std::string& raw) {
    const SymmetryNotation notation = parse_symmetry_notation(raw);
    if (notation.base.empty()) {
        const std::string na = "N/A";
        draw_text(canvas, centre_x - bitmap_text_width(na, 1) / 2,
                  baseline_y, na, 91, 105, 124, 1);
        return;
    }
    const int base_scale = 2;
    const int sub_scale = 1;
    const int base_w = bitmap_text_width(notation.base, base_scale);
    const int sub_w = bitmap_text_width(notation.subscript, sub_scale);
    const int sup_w = bitmap_text_width(notation.superscript, sub_scale);
    const int total = base_w + std::max(sub_w, sup_w);
    int x = centre_x - total / 2;
    draw_text(canvas, x, baseline_y, notation.base, 91, 105, 124, base_scale);
    x += base_w;
    if (!notation.subscript.empty()) {
        draw_text(canvas, x, baseline_y + 8, notation.subscript, 91, 105, 124, sub_scale);
    }
    if (!notation.superscript.empty()) {
        // Bitmap export uses ASCII prime/star fallback while preserving the
        // vertical superscript position. SVG carries the exact Unicode glyph.
        std::string sup = notation.superscript;
        for (char& c : sup) {
            if (static_cast<unsigned char>(c) >= 0x80) c = '\'';
        }
        draw_text(canvas, x, baseline_y - 4, sup, 91, 105, 124, sub_scale);
    }
}

} // namespace

bool write_mo_diagram_png(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    const int width = std::max(760, options.width);
    const int height = std::max(720, options.height);
    Canvas canvas(width, height);
    const auto positioned = position_levels(data, width, height);

    draw_text(canvas, 36, 28, "VALENCE MO DIAGRAM", 24, 34, 49, 2);
    draw_text(canvas, 36, 55,
              options.energy_axis_mode == EnergyAxisMode::Linear
                  ? "ENERGY AXIS: LINEAR"
                  : "ENERGY AXIS: ADAPTIVE NONLINEAR LOG-GAP V3",
              92, 106, 126, 1);

    const int axis_x = 70;
    canvas.line(axis_x, 110, axis_x, height - 62, 45, 55, 69, 2);
    canvas.line(axis_x, 110, axis_x - 5, 120, 45, 55, 69, 2);
    canvas.line(axis_x, 110, axis_x + 5, 120, 45, 55, 69, 2);

    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const double energy = energy_at_axis_fraction(fraction, data.energy_transform);
        const int y = static_cast<int>(std::lround(
            map_energy_y(energy, data.energy_transform, height)));
        canvas.line(axis_x - 5, y, axis_x + 5, y, 116, 130, 148, 1);
        const std::string tick = format_energy(
            energy, options.energy_unit,
            options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 3);
        draw_text(canvas, axis_x + 10, y - 4, tick, 105, 119, 138, 1);
    }

    for (const auto& item : positioned) {
        const auto& level = *item.level;
        std::uint8_t r = 41, g = 56, b = 76;
        if (level.metadata.region == OrbitalRegion::Virtual) { r = 119; g = 132; b = 150; }
        if (level.homo) { r = 19; g = 138; b = 98; }
        if (level.lumo) { r = 181; g = 106; b = 0; }
        if (level.metadata.selected) { r = 40; g = 109; b = 224; }

        const int x0 = static_cast<int>(std::lround(item.x - 45.0));
        const int x1 = static_cast<int>(std::lround(item.x + 45.0));
        const int y = static_cast<int>(std::lround(item.y));
        canvas.line(x0, y, x1, y, r, g, b, level.metadata.selected ? 4 : 2);
        if (level.electrons.alpha > 0) draw_electron(canvas, static_cast<int>(item.x) - 10, y - 2, true);
        if (level.electrons.beta > 0) draw_electron(canvas, static_cast<int>(item.x) + 10, y - 2, false);

        const std::string type = orbital_type_plain(level.annotation);
        const int type_x = static_cast<int>(item.x) - bitmap_text_width(type, 1) / 2;
        draw_text(canvas, type_x, y + 8, type, 67, 81, 102, 1);
        draw_symmetry_png(canvas, static_cast<int>(item.x), y + 20,
                          level.metadata.symmetry);

        if (level.annotation.delocalised_pi.available) {
            std::ostringstream pi;
            pi << "PI";
            if (level.annotation.delocalised_pi.participating_atoms > 0) {
                pi << "^" << level.annotation.delocalised_pi.participating_atoms;
            }
            if (level.annotation.delocalised_pi.participating_electrons > 0.0) {
                pi << "_" << static_cast<int>(std::lround(
                    level.annotation.delocalised_pi.participating_electrons));
            }
            const std::string value = pi.str();
            draw_text(canvas,
                      static_cast<int>(item.x) - bitmap_text_width(value, 1) / 2,
                      y + 31, value, 91, 101, 117, 1);
        }
    }

    draw_text(canvas, 36, height - 32,
              "FIGURE LABELS: ORBITAL TYPE / BONDING CLASS; MO NUMBERS ARE IN HOVER + JSON/CSV",
              113, 126, 144, 1);
    return write_png_rgba(canvas, path, error);
}

