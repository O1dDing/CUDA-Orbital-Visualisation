MODiagramData build_mo_diagram_data(const Wavefunction& wavefunction,
                                    const MODiagramOptions& options) {
    MODiagramData data = cov_legacy_build_mo_diagram_data(wavefunction, options);
    data.textbook_diatomic = homonuclear_diatomic(wavefunction);
    data.textbook_h2 = data.textbook_diatomic && is_h2(wavefunction);
    data.species_label = species_label(wavefunction);
    if (data.textbook_diatomic) {
        data.plan.classification = DiagramClassification::Simple;
        data.plan.strict_salc_available = false;
        data.plan.machine_reason =
            "homonuclear-diatomic: qualitative AO-to-MO textbook presentation; "
            "molecular energies are quantitative";
    }
    return data;
}

bool write_mo_diagram_svg(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open SVG output path";
        return false;
    }

    const int width = std::max(720, options.width);
    const int height = std::max(520, options.height);
    out << "<?xml version='1.0' encoding='UTF-8'?>\n"
        << "<svg xmlns='http://www.w3.org/2000/svg' width='" << width
        << "' height='" << height << "' viewBox='0 0 " << width << ' ' << height << "'>\n"
        << "<rect width='100%' height='100%' fill='#ffffff'/>\n";
    svg_energy_axis(out, height);

    if (data.textbook_diatomic && data.levels.size() >= 2) {
        std::vector<MODiagramLevel> levels = data.levels;
        std::sort(levels.begin(), levels.end(), [](const auto& a, const auto& b) {
            return a.metadata.energy_hartree < b.metadata.energy_hartree;
        });
        const MODiagramLevel& bonding = levels.front();
        const MODiagramLevel& antibonding = levels.back();
        const Layout layout = make_layout(data);
        const double y_bond = energy_y(bonding.metadata.energy_hartree, layout, height, 130, 80);
        const double y_anti = energy_y(antibonding.metadata.energy_hartree, layout, height, 130, 80);
        const double y_ao = (y_bond + y_anti) * 0.5;
        const double x_left = 190.0;
        const double x_centre = static_cast<double>(width) * 0.5;
        const double x_right = static_cast<double>(width) - 190.0;
        const std::string species = data.species_label.empty() ? "H2" : data.species_label;

        out << "<text x='" << width / 2
            << "' y='42' text-anchor='middle' font-family='sans-serif' font-size='26' "
               "font-weight='700' fill='#111827'>"
            << xml_escape(species) << " molecular-orbital diagram</text>\n";
        out << "<text x='" << width / 2
            << "' y='68' text-anchor='middle' font-family='sans-serif' font-size='13' fill='#6b7280'>"
               "Qualitative AO reference; molecular-orbital energies are read from the wavefunction</text>\n";

        for (const double x : {x_left, x_right}) {
            out << "<line x1='" << x - 48 << "' y1='" << y_ao
                << "' x2='" << x + 48 << "' y2='" << y_ao
                << "' stroke='#111827' stroke-width='3'/>\n";
            out << "<text x='" << x << "' y='" << y_ao + 28
                << "' text-anchor='middle' font-family='sans-serif' font-size='16' fill='#374151'>H 1s</text>\n";
            if (data.textbook_h2) svg_arrow(out, x, y_ao - 4, true);
        }

        out << "<line x1='" << x_centre - 58 << "' y1='" << y_bond
            << "' x2='" << x_centre + 58 << "' y2='" << y_bond
            << "' stroke='#111827' stroke-width='4'/>\n";
        out << "<line x1='" << x_centre - 58 << "' y1='" << y_anti
            << "' x2='" << x_centre + 58 << "' y2='" << y_anti
            << "' stroke='#111827' stroke-width='4'/>\n";

        for (const double x : {x_left, x_right}) {
            const double source_x = x < x_centre ? x + 48 : x - 48;
            const double target_x = x < x_centre ? x_centre - 58 : x_centre + 58;
            out << "<line x1='" << source_x << "' y1='" << y_ao
                << "' x2='" << target_x << "' y2='" << y_bond
                << "' stroke='#9ca3af' stroke-width='2' stroke-dasharray='5 5'/>\n";
            out << "<line x1='" << source_x << "' y1='" << y_ao
                << "' x2='" << target_x << "' y2='" << y_anti
                << "' stroke='#9ca3af' stroke-width='2' stroke-dasharray='5 5'/>\n";
        }

        if (bonding.electrons.alpha > 0) svg_arrow(out, x_centre - 11, y_bond - 4, true);
        if (bonding.electrons.beta > 0) svg_arrow(out, x_centre + 11, y_bond - 4, false);
        out << "<text x='" << x_centre << "' y='" << y_bond + 31
            << "' text-anchor='middle' font-family='sans-serif' font-size='18' fill='#111827'>σ1s</text>\n";
        out << "<text x='" << x_centre << "' y='" << y_anti - 16
            << "' text-anchor='middle' font-family='sans-serif' font-size='18' fill='#111827'>σ*1s</text>\n";
        out << "<text x='" << x_centre + 78 << "' y='" << y_bond + 5
            << "' font-family='sans-serif' font-size='13' fill='#4b5563'>"
            << xml_escape(format_energy(bonding.metadata.energy_hartree, options.energy_unit, 5))
            << "</text>\n";
        out << "<text x='" << x_centre + 78 << "' y='" << y_anti + 5
            << "' font-family='sans-serif' font-size='13' fill='#4b5563'>"
            << xml_escape(format_energy(antibonding.metadata.energy_hartree, options.energy_unit, 5))
            << "</text>\n";
    } else {
        const Layout layout = make_layout(data);
        out << "<text x='" << width / 2
            << "' y='42' text-anchor='middle' font-family='sans-serif' font-size='25' "
               "font-weight='700' fill='#111827'>"
            << xml_escape(classification_title(data.plan.classification)) << "</text>\n";
        out << "<text x='" << width / 2
            << "' y='67' text-anchor='middle' font-family='sans-serif' font-size='13' fill='#6b7280'>"
            << xml_escape(data.plan.machine_reason) << "</text>\n";

        const double column_span = static_cast<double>(width - 330) /
                                   static_cast<double>(std::max<std::size_t>(1, layout.columns.size()));
        for (std::size_t i = 0; i < layout.columns.size(); ++i) {
            const double x = 150.0 + (static_cast<double>(i) + 0.5) * column_span;
            out << "<text x='" << x
                << "' y='96' text-anchor='middle' font-family='sans-serif' font-size='15' "
                   "font-weight='600' fill='#4b5563'>"
                << xml_escape(layout.columns[i]) << "</text>\n";
        }

        for (const PositionedLevel& positioned : position_levels(data, layout, width, height)) {
            const MODiagramLevel& level = *positioned.level;
            const char* colour = level.metadata.selected ? "#2563eb" :
                                 level.homo ? "#059669" :
                                 level.lumo ? "#d97706" :
                                 level.metadata.region == OrbitalRegion::Virtual ? "#9ca3af" : "#1f2937";
            out << "<line x1='" << positioned.x - 44 << "' y1='" << positioned.y
                << "' x2='" << positioned.x + 44 << "' y2='" << positioned.y
                << "' stroke='" << colour << "' stroke-width='"
                << (level.metadata.selected ? 4 : 2.5) << "'/>\n";
            if (level.electrons.alpha > 0) svg_arrow(out, positioned.x - 11, positioned.y - 5, true);
            if (level.electrons.beta > 0) svg_arrow(out, positioned.x + 11, positioned.y - 5, false);

            const double label_x = positioned.x + 62;
            if (std::abs(positioned.label_y - positioned.y) > 3.0) {
                out << "<line x1='" << positioned.x + 44 << "' y1='" << positioned.y
                    << "' x2='" << label_x - 8 << "' y2='" << positioned.label_y - 4
                    << "' stroke='#d1d5db' stroke-width='1'/>\n";
            }
            out << "<text x='" << label_x << "' y='" << positioned.label_y
                << "' font-family='sans-serif' font-size='13' fill='#111827'>"
                << xml_escape(level.metadata.display_label) << "</text>\n";
            out << "<text x='" << label_x << "' y='" << positioned.label_y + 16
                << "' font-family='sans-serif' font-size='11' fill='#6b7280'>"
                << xml_escape(format_energy(level.metadata.energy_hartree, options.energy_unit,
                                            options.energy_unit == EnergyUnit::JoulePerMol ? 1 : 5))
                << "</text>\n";
        }
    }

    out << "</svg>\n";
    return static_cast<bool>(out);
}

