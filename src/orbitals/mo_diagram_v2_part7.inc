namespace {

int bitmap_text_width(const std::string& text, const int scale) {
    int width = 0;
    for (const char c : text) width += (c == ' ' ? 4 : 6) * scale;
    return width;
}

void draw_delta_pi_png(Canvas& canvas,
                       const int right_x,
                       const int centre_y,
                       const ExportPiColour colour) {
    // Raster exports use vector strokes for the same compact Δπ symbol as the
    // interactive canvas.  This avoids depending on a platform bitmap font for
    // Greek glyphs and keeps the symbol legible at all supported export sizes.
    const int delta_x = right_x - 25;
    canvas.line(delta_x, centre_y + 5, delta_x + 5, centre_y - 5,
                colour.r, colour.g, colour.b, 2);
    canvas.line(delta_x + 5, centre_y - 5, delta_x + 10, centre_y + 5,
                colour.r, colour.g, colour.b, 2);
    canvas.line(delta_x, centre_y + 5, delta_x + 10, centre_y + 5,
                colour.r, colour.g, colour.b, 2);

    const int pi_x = delta_x + 14;
    canvas.line(pi_x, centre_y - 4, pi_x + 10, centre_y - 4,
                colour.r, colour.g, colour.b, 2);
    canvas.line(pi_x + 2, centre_y - 4, pi_x + 2, centre_y + 5,
                colour.r, colour.g, colour.b, 2);
    canvas.line(pi_x + 8, centre_y - 4, pi_x + 8, centre_y + 5,
                colour.r, colour.g, colour.b, 2);
}

void draw_pi_interactions_png(Canvas& canvas,
                              const MODiagramData& data,
                              const std::vector<PositionedLevel>& positioned) {
    std::size_t rendered = 0;
    for (const auto& interaction : data.pi_interactions) {
        const auto bracket = export_pi_bracket(interaction, positioned);
        if (!bracket) continue;
        const ExportPiColour colour = export_pi_colour(interaction.kind);
        const int x = canvas.width - 26 - static_cast<int>(rendered) * 48;
        const int top = static_cast<int>(std::lround(bracket->top));
        const int bottom = static_cast<int>(std::lround(bracket->bottom));
        constexpr int cap = 7;
        canvas.line(x, top, x, bottom, colour.r, colour.g, colour.b, 2);
        canvas.line(x - cap, top, x, top, colour.r, colour.g, colour.b, 2);
        canvas.line(x - cap, bottom, x, bottom, colour.r, colour.g, colour.b, 2);
        draw_delta_pi_png(canvas, x - cap - 5, (top + bottom) / 2, colour);
        ++rendered;
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

    const char* diagram_title="VALENCE MO DIAGRAM";
    if (data.mode==MODiagramMode::DelocalisedPiFamilyOnly) {
        diagram_title="DELOCALISED PI MO DIAGRAM";
    } else if (data.mode==MODiagramMode::MulticentreActiveSpaceOnly) {
        diagram_title="MULTICENTRE ACTIVE-SPACE MO DIAGRAM";
    }
    draw_text(canvas, 36, 28, diagram_title, 24, 34, 49, 2);
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

        // Human-facing export: orbital type directly below each level; MO
        // identifiers intentionally stay out of the figure. Collapse the
        // fully-unavailable pair to a single N/A to reduce visual noise.
        std::string type = orbital_type_plain(level.annotation);
        if (type == "N/A | N/A") type = "N/A";
        const int type_x = static_cast<int>(item.x) - bitmap_text_width(type, 1) / 2;
        draw_text(canvas, type_x, y + 13, type, 67, 81, 102, 1);
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
            const std::string topology_suffix=compact_pi_topology_suffix(
                level.annotation.delocalised_pi);
            if (!topology_suffix.empty()) {
                pi << " " << topology_suffix;
            }
            const std::string value = pi.str();
            draw_text(canvas,
                      static_cast<int>(item.x) - bitmap_text_width(value, 1) / 2,
                      y + 25, value, 91, 101, 117, 1);
        }
    }

    draw_pi_interactions_png(canvas, data, positioned);

    draw_text(canvas, 36, height - 32,
              "FIGURE LABELS: ORBITAL TYPE / BONDING CLASS; MO NUMBERS ARE IN HOVER + JSON/CSV",
              113, 126, 144, 1);
    return write_png_rgba(canvas, path, error);
}
