bool write_mo_diagram_png(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    const int width = std::max(720, options.width);
    const int height = std::max(520, options.height);
    Canvas canvas(width, height);
    canvas.line(62, height - 70, 62, 115, 20, 25, 32, 2);
    canvas.line(62, 115, 56, 127, 20, 25, 32, 2);
    canvas.line(62, 115, 68, 127, 20, 25, 32, 2);
    draw_text(canvas, 22, 92, "ENERGY", 70, 75, 85, 2);

    if (data.textbook_diatomic && data.levels.size() >= 2) {
        std::vector<MODiagramLevel> levels = data.levels;
        std::sort(levels.begin(), levels.end(), [](const auto& a, const auto& b) {
            return a.metadata.energy_hartree < b.metadata.energy_hartree;
        });
        const Layout layout = make_layout(data);
        const int y_bond = static_cast<int>(std::lround(
            energy_y(levels.front().metadata.energy_hartree, layout, height, 130, 80)));
        const int y_anti = static_cast<int>(std::lround(
            energy_y(levels.back().metadata.energy_hartree, layout, height, 130, 80)));
        const int y_ao = (y_bond + y_anti) / 2;
        const int x_left = 190;
        const int x_centre = width / 2;
        const int x_right = width - 190;
        const std::string species = data.species_label.empty() ? "H2" : data.species_label;

        draw_text(canvas, width / 2 - 120, 28, species + " MO DIAGRAM", 20, 25, 32, 3);
        draw_text(canvas, width / 2 - 250, 60,
                  "QUALITATIVE AO REFERENCE - MO ENERGIES FROM FILE", 105, 110, 120, 1);

        for (const int x : {x_left, x_right}) {
            canvas.line(x - 48, y_ao, x + 48, y_ao, 20, 25, 32, 3);
            draw_text(canvas, x - 30, y_ao + 22, "H 1S", 55, 60, 70, 2);
            if (data.textbook_h2) draw_arrow(canvas, x, y_ao - 4, true);
        }
        canvas.line(x_centre - 58, y_bond, x_centre + 58, y_bond, 20, 25, 32, 4);
        canvas.line(x_centre - 58, y_anti, x_centre + 58, y_anti, 20, 25, 32, 4);

        for (const int x : {x_left, x_right}) {
            const int source_x = x < x_centre ? x + 48 : x - 48;
            const int target_x = x < x_centre ? x_centre - 58 : x_centre + 58;
            for (int segment = 0; segment < 20; ++segment) {
                if ((segment % 2) != 0) continue;
                const double a = static_cast<double>(segment) / 20.0;
                const double b = static_cast<double>(segment + 1) / 20.0;
                canvas.line(static_cast<int>(source_x + (target_x - source_x) * a),
                            static_cast<int>(y_ao + (y_bond - y_ao) * a),
                            static_cast<int>(source_x + (target_x - source_x) * b),
                            static_cast<int>(y_ao + (y_bond - y_ao) * b),
                            145, 150, 160, 1);
                canvas.line(static_cast<int>(source_x + (target_x - source_x) * a),
                            static_cast<int>(y_ao + (y_anti - y_ao) * a),
                            static_cast<int>(source_x + (target_x - source_x) * b),
                            static_cast<int>(y_ao + (y_anti - y_ao) * b),
                            145, 150, 160, 1);
            }
        }

        if (levels.front().electrons.alpha > 0) draw_arrow(canvas, x_centre - 11, y_bond - 4, true);
        if (levels.front().electrons.beta > 0) draw_arrow(canvas, x_centre + 11, y_bond - 4, false);
        draw_text(canvas, x_centre - 45, y_bond + 26, "SIGMA 1S", 20, 25, 32, 2);
        draw_text(canvas, x_centre - 55, y_anti - 28, "SIGMA* 1S", 20, 25, 32, 2);
        draw_text(canvas, x_centre + 80, y_bond - 7,
                  format_energy(levels.front().metadata.energy_hartree, options.energy_unit, 4),
                  75, 80, 90, 1);
        draw_text(canvas, x_centre + 80, y_anti - 7,
                  format_energy(levels.back().metadata.energy_hartree, options.energy_unit, 4),
                  75, 80, 90, 1);
    } else {
        const Layout layout = make_layout(data);
        std::string title = data.plan.classification == DiagramClassification::SymmetryGrouped
            ? "SYMMETRY-GROUPED MO DIAGRAM"
            : data.plan.classification == DiagramClassification::SalcUnavailable
                ? "MO DIAGRAM - SALC UNAVAILABLE"
                : "MO ENERGY DIAGRAM";
        draw_text(canvas, width / 2 - 180, 28, title, 20, 25, 32, 2);

        const double column_span = static_cast<double>(width - 330) /
                                   static_cast<double>(std::max<std::size_t>(1, layout.columns.size()));
        for (std::size_t i = 0; i < layout.columns.size(); ++i) {
            const int x = static_cast<int>(150.0 + (static_cast<double>(i) + 0.5) * column_span) -
                          static_cast<int>(layout.columns[i].size() * 6);
            draw_text(canvas, x, 78, layout.columns[i], 80, 85, 95, 2);
        }

        for (const PositionedLevel& positioned : position_levels(data, layout, width, height)) {
            const MODiagramLevel& level = *positioned.level;
            std::array<std::uint8_t, 3> colour = level.metadata.region == OrbitalRegion::Virtual
                ? std::array<std::uint8_t, 3>{150, 155, 165}
                : std::array<std::uint8_t, 3>{35, 40, 48};
            if (level.homo) colour = {5, 150, 105};
            if (level.lumo) colour = {210, 120, 20};
            if (level.metadata.selected) colour = {37, 99, 235};
            const int x = static_cast<int>(std::lround(positioned.x));
            const int y = static_cast<int>(std::lround(positioned.y));
            canvas.line(x - 44, y, x + 44, y, colour[0], colour[1], colour[2],
                        level.metadata.selected ? 4 : 3);
            if (level.electrons.alpha > 0) draw_arrow(canvas, x - 11, y - 5, true);
            if (level.electrons.beta > 0) draw_arrow(canvas, x + 11, y - 5, false);
            const int label_y = static_cast<int>(std::lround(positioned.label_y));
            if (std::abs(label_y - y) > 3) {
                canvas.line(x + 44, y, x + 56, label_y - 5, 200, 205, 215, 1);
            }
            draw_text(canvas, x + 62, label_y - 7, level.metadata.display_label, 25, 30, 38, 2);
        }
    }

    return save_png(canvas, path, error);
}

