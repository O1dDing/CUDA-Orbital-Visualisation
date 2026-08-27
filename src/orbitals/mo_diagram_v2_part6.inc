    data.plan.classification = DiagramClassification::Simple;
    data.plan.complex_system =
        wavefunction.atoms.size() > 12 || wavefunction.orbitals.size() > 80;
    data.plan.strict_salc_available = false;
    data.plan.machine_reason =
        "valence-central: direct Molden MO data; strict SALC reconstruction is not claimed";
    data.frontier =
        find_frontier_orbitals(wavefunction.orbitals,
                               options.filter.occupation_threshold);
    data.metadata =
        build_orbital_metadata(wavefunction,
                               options.selected_index,
                               options.degeneracy,
                               options.filter);

    data.annotations.reserve(wavefunction.orbitals.size());
    for (const auto& mo : wavefunction.orbitals) {
        data.annotations.push_back(annotate_orbital(mo));
    }

    data.selection =
        build_valence_selection_plan(wavefunction, options, data.metadata);
    data.levels.reserve(data.selection.included_indices.size());
    std::vector<double> axis_energies;
    axis_energies.reserve(data.selection.included_indices.size());

    for (const std::size_t index : data.selection.included_indices) {
        if (index >= data.metadata.size()) continue;
        MODiagramLevel level;
        level.metadata = data.metadata[index];
        level.annotation = data.annotations[index];
        level.electrons =
            electron_glyphs_for_orbital(
                wavefunction.orbitals[index],
                data.frontier.separate_spin_sets);
        level.homo = data.frontier.homo && *data.frontier.homo == index;
        level.lumo = data.frontier.lumo && *data.frontier.lumo == index;
        level.layout_energy_hartree =
            group_layout_energy(data.metadata, index);
        axis_energies.push_back(level.layout_energy_hartree);
        data.levels.push_back(std::move(level));
    }

    data.energy_transform =
        build_energy_transform(axis_energies,
                               options.energy_axis_mode,
                               options.nonlinear_minimum_gap_weight);
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

    const int width = std::max(760, options.width);
    const int height = std::max(620, options.height);
    const auto positioned = position_levels(data, width, height);

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 "
        << width << ' ' << height << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"#fbfcfe\"/>\n";
    out << "<text x=\"36\" y=\"40\" fill=\"#182231\" "
           "font-family=\"sans-serif\" font-size=\"24\" font-weight=\"600\">"
           "Valence MO diagram</text>\n";
    out << "<text x=\"36\" y=\"64\" fill=\"#68778b\" "
           "font-family=\"sans-serif\" font-size=\"12\">"
        << xml_escape(data.selection.summary) << "</text>\n";
    out << "<text x=\"36\" y=\"82\" fill=\"#7f8c9c\" "
           "font-family=\"sans-serif\" font-size=\"11\">"
           "Molden-derived MO energies/occupations/symmetry; classification "
           "is shown only when explicitly supported by producer text.</text>\n";

    const int axis_x = 70;
    out << "<line x1=\"" << axis_x << "\" y1=\"110\" x2=\"" << axis_x
        << "\" y2=\"" << height - 62
        << "\" stroke=\"#2d3745\" stroke-width=\"2\"/>\n";
    out << "<path d=\"M " << axis_x
        << " 98 l -5 10 h 10 z\" fill=\"#2d3745\"/>\n";
    out << "<text x=\"24\" y=\"124\" fill=\"#4d5a6c\" "
           "font-family=\"sans-serif\" font-size=\"12\" "
           "transform=\"rotate(-90 24 124)\">Energy · "
        << (options.energy_axis_mode == EnergyAxisMode::Linear
                ? "linear"
                : "adaptive nonlinear (log-gap)")
        << " (" << xml_escape(energy_unit_symbol(options.energy_unit))
        << ")</text>\n";

    if (options.energy_axis_mode == EnergyAxisMode::NonlinearFocus) {
        out << "<text x=\"82\" y=\"102\" fill=\"#8b6173\" "
               "font-family=\"sans-serif\" font-size=\"10\">"
               "nonlinear axis · ticks show exact physical energy</text>\n";
    }

    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const double energy =
            energy_at_axis_fraction(fraction, data.energy_transform);
        const double y =
            map_energy_y(energy, data.energy_transform, height);
        out << "<line x1=\"64\" y1=\"" << y << "\" x2=\"76\" y2=\"" << y
            << "\" stroke=\"#7c8795\" stroke-width=\"1\"/>\n";
        out << "<text x=\"80\" y=\"" << y + 4
            << "\" fill=\"#7c8795\" font-family=\"sans-serif\" "
               "font-size=\"10\">"
            << xml_escape(format_energy(
                   energy,
                   options.energy_unit,
                   options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 4))
            << "</text>\n";
    }

    for (const auto& item : positioned) {
        const auto& level = *item.level;
        const char* colour =
            level.metadata.selected ? "#286de0" :
            level.homo ? "#138a62" :
            level.lumo ? "#b56a00" :
            level.metadata.region == OrbitalRegion::Virtual
                ? "#778496" : "#29384c";
        const double x0 = item.x - 45.0;
        const double x1 = item.x + 45.0;
        out << "<line x1=\"" << x0 << "\" y1=\"" << item.y
            << "\" x2=\"" << x1 << "\" y2=\"" << item.y
            << "\" stroke=\"" << colour << "\" stroke-width=\""
            << (level.metadata.selected ? 4 : 2.2) << "\"/>\n";

        if (level.electrons.alpha > 0) {
            out << "<text x=\"" << item.x - 14 << "\" y=\"" << item.y - 5
                << "\" fill=\"#151d28\" font-family=\"sans-serif\" "
                   "font-size=\"19\">↑</text>\n";
        }
        if (level.electrons.beta > 0) {
            out << "<text x=\"" << item.x + 3 << "\" y=\"" << item.y - 5
                << "\" fill=\"#151d28\" font-family=\"sans-serif\" "
                   "font-size=\"19\">↓</text>\n";
        }

        const double label_x = static_cast<double>(width) * 0.72;
        if (std::abs(item.label_y - item.y) > 5.0) {
            out << "<line x1=\"" << x1 + 2 << "\" y1=\"" << item.y
                << "\" x2=\"" << label_x - 6 << "\" y2=\"" << item.label_y
                << "\" stroke=\"#b6c0cc\" stroke-width=\"1\"/>\n";
        }

        out << "<text x=\"" << label_x << "\" y=\"" << item.label_y
            << "\" fill=\"" << colour
            << "\" font-family=\"sans-serif\" font-size=\"12\" "
               "font-weight=\"600\">"
            << xml_escape(level.metadata.display_label)
            << "</text>\n";

        std::ostringstream secondary;
        if (!level.metadata.symmetry.empty()) secondary << level.metadata.symmetry;
        const std::string summary = annotation_summary(level.annotation);
        if (!summary.empty()) {
            if (secondary.tellp() > 0) secondary << " · ";
            secondary << summary;
        }
        if (secondary.tellp() > 0) {
            out << "<text x=\"" << label_x + 48 << "\" y=\"" << item.label_y
                << "\" fill=\"#667588\" font-family=\"sans-serif\" "
                   "font-size=\"10\">"
                << xml_escape(secondary.str()) << "</text>\n";
        }

        out << "<text x=\"" << label_x << "\" y=\"" << item.label_y + 14
            << "\" fill=\"#7c8795\" font-family=\"sans-serif\" "
               "font-size=\"10\">"
            << xml_escape(format_energy(
                   level.metadata.energy_hartree,
                   options.energy_unit,
                   options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 5))
            << "</text>\n";
