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

namespace {

std::string family_symbol_svg(const std::string& family) {
    if (family == "sigma") return "σ";
    if (family == "pi") return "π";
    if (family == "delta") return "δ";
    if (family == "phi") return "φ";
    return "N/A";
}

std::string family_plain(const std::string& family) {
    if (family == "sigma") return "SIGMA";
    if (family == "pi") return "PI";
    if (family == "delta") return "DELTA";
    if (family == "phi") return "PHI";
    return "N/A";
}

std::string bonding_plain(const BondingClass value) {
    switch (value) {
        case BondingClass::Bonding: return "BONDING";
        case BondingClass::Nonbonding: return "NONBONDING";
        case BondingClass::Antibonding: return "ANTIBONDING";
        default: return "N/A";
    }
}

std::string orbital_type_plain(const OrbitalAnnotation& annotation) {
    const std::string family = family_plain(annotation.family);
    const std::string bonding = bonding_plain(annotation.bonding_class);
    std::string out = family == "N/A" && bonding == "N/A"
                          ? "N/A"
                          : family + " | " + bonding;
    if (annotation.multicentre.available) out += " | " + annotation.multicentre.label;
    if (annotation.delocalised_pi.available) out += " | LARGE-PI";
    return out;
}

std::string orbital_type_svg(const OrbitalAnnotation& annotation) {
    const std::string family = family_symbol_svg(annotation.family);
    const std::string bonding = annotation.bonding_class == BondingClass::Unclassified
                                    ? "N/A"
                                    : bonding_class_name(annotation.bonding_class);
    std::string out = family == "N/A" && bonding == "N/A"
                          ? "N/A"
                          : family + " · " + bonding;
    if (annotation.multicentre.available) {
        out += " · " + annotation.multicentre.label;
    }
    return out;
}

void write_svg_symmetry(std::ostream& out,
                        const double x,
                        const double y,
                        const std::string& raw,
                        const char* colour,
                        const int font_size = 10) {
    const SymmetryNotation notation = parse_symmetry_notation(raw);
    if (notation.base.empty()) {
        out << "<text x=\"" << x << "\" y=\"" << y
            << "\" text-anchor=\"middle\" fill=\"" << colour
            << "\" font-family=\"sans-serif\" font-size=\"" << font_size
            << "\">N/A</text>\n";
        return;
    }
    out << "<text x=\"" << x << "\" y=\"" << y
        << "\" text-anchor=\"middle\" fill=\"" << colour
        << "\" font-family=\"sans-serif\" font-size=\"" << font_size << "\">"
        << xml_escape(notation.base);
    if (!notation.subscript.empty()) {
        out << "<tspan baseline-shift=\"sub\" font-size=\"70%\">"
            << xml_escape(notation.subscript) << "</tspan>";
    }
    if (!notation.superscript.empty()) {
        out << "<tspan baseline-shift=\"super\" font-size=\"70%\">"
            << xml_escape(notation.superscript) << "</tspan>";
    }
    out << "</text>\n";
}

void write_svg_electron(std::ostream& out,
                        const double x,
                        const double y,
                        const bool up) {
    constexpr double length = 11.0;
    constexpr double head_dx = 3.2;
    constexpr double head_dy = 3.5;
    const double tip = up ? y - length : y + length;
    const double wing_y = up ? tip + head_dy : tip - head_dy;
    out << "<line x1=\"" << x << "\" y1=\"" << y
        << "\" x2=\"" << x << "\" y2=\"" << tip
        << "\" stroke=\"#151d28\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>\n";
    out << "<line x1=\"" << x << "\" y1=\"" << tip
        << "\" x2=\"" << x - head_dx << "\" y2=\"" << wing_y
        << "\" stroke=\"#151d28\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>\n";
    out << "<line x1=\"" << x << "\" y1=\"" << tip
        << "\" x2=\"" << x + head_dx << "\" y2=\"" << wing_y
        << "\" stroke=\"#151d28\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>\n";
}

void write_svg_delocalised_pi(std::ostream& out,
                              const double x,
                              const double y,
                              const DelocalisedPiDescriptor& pi,
                              const char* colour) {
    if (!pi.available) return;
    out << "<text x=\"" << x << "\" y=\"" << y
        << "\" text-anchor=\"middle\" fill=\"" << colour
        << "\" font-family=\"sans-serif\" font-size=\"10\">Π";
    if (pi.participating_atoms > 0) {
        out << "<tspan baseline-shift=\"super\" font-size=\"70%\">"
            << pi.participating_atoms << "</tspan>";
    }
    if (pi.participating_electrons > 0.0) {
        out << "<tspan baseline-shift=\"sub\" font-size=\"70%\">"
            << static_cast<int>(std::lround(pi.participating_electrons)) << "</tspan>";
    }
    out << "</text>\n";
}

} // namespace

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
    const int height = std::max(720, options.height);
    const auto positioned = position_levels(data, width, height);

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 "
        << width << ' ' << height << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"#fbfcfe\"/>\n";
    out << "<text x=\"36\" y=\"40\" fill=\"#182231\" font-family=\"sans-serif\" "
           "font-size=\"24\" font-weight=\"600\">Valence MO diagram</text>\n";
    out << "<text x=\"36\" y=\"64\" fill=\"#68778b\" font-family=\"sans-serif\" "
           "font-size=\"12\">" << xml_escape(data.selection.summary) << "</text>\n";
    out << "<text x=\"36\" y=\"82\" fill=\"#7f8c9c\" font-family=\"sans-serif\" "
           "font-size=\"11\">Molden-derived MO energies/occupations/symmetry; type/class is N/A unless supported by source data.</text>\n";

    const int axis_x = 70;
    out << "<line x1=\"" << axis_x << "\" y1=\"110\" x2=\"" << axis_x
        << "\" y2=\"" << height - 62 << "\" stroke=\"#2d3745\" stroke-width=\"2\"/>\n";
    out << "<path d=\"M " << axis_x << " 98 l -5 10 h 10 z\" fill=\"#2d3745\"/>\n";
    out << "<text x=\"24\" y=\"124\" fill=\"#4d5a6c\" font-family=\"sans-serif\" "
           "font-size=\"12\" transform=\"rotate(-90 24 124)\">Energy · "
        << (options.energy_axis_mode == EnergyAxisMode::Linear ? "linear" : "adaptive nonlinear (log-gap v3)")
        << " (" << xml_escape(energy_unit_symbol(options.energy_unit)) << ")</text>\n";
    if (options.energy_axis_mode == EnergyAxisMode::NonlinearFocus) {
        out << "<text x=\"82\" y=\"102\" fill=\"#8b6173\" font-family=\"sans-serif\" "
               "font-size=\"10\">nonlinear axis · ticks show exact physical energy</text>\n";
    }

    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const double energy = energy_at_axis_fraction(fraction, data.energy_transform);
        const double y = map_energy_y(energy, data.energy_transform, height);
        out << "<line x1=\"64\" y1=\"" << y << "\" x2=\"76\" y2=\"" << y
            << "\" stroke=\"#7c8795\" stroke-width=\"1\"/>\n";
        out << "<text x=\"80\" y=\"" << y + 4
            << "\" fill=\"#7c8795\" font-family=\"sans-serif\" font-size=\"10\">"
            << xml_escape(format_energy(energy, options.energy_unit,
                   options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 4))
            << "</text>\n";
    }

    for (const auto& item : positioned) {
        const auto& level = *item.level;
        const char* colour = level.metadata.selected ? "#286de0" :
                             level.homo ? "#138a62" :
                             level.lumo ? "#b56a00" :
                             level.metadata.region == OrbitalRegion::Virtual ? "#778496" : "#29384c";
        const double x0 = item.x - 45.0;
        const double x1 = item.x + 45.0;
        out << "<line x1=\"" << x0 << "\" y1=\"" << item.y
            << "\" x2=\"" << x1 << "\" y2=\"" << item.y
            << "\" stroke=\"" << colour << "\" stroke-width=\""
            << (level.metadata.selected ? 4 : 2.2) << "\"/>\n";
        if (level.electrons.alpha > 0) write_svg_electron(out, item.x - 10.0, item.y - 2.0, true);
        if (level.electrons.beta > 0) write_svg_electron(out, item.x + 10.0, item.y - 2.0, false);

        // Export deliberately omits MO numbers. Human-readable orbital type is
        // attached directly below the physical level; raw numbering remains in
        // hover/JSON/CSV where it is useful without cluttering the figure.
        out << "<text x=\"" << item.x << "\" y=\"" << item.y + 16
            << "\" text-anchor=\"middle\" fill=\"#435166\" font-family=\"sans-serif\" "
               "font-size=\"10\">" << xml_escape(orbital_type_svg(level.annotation)) << "</text>\n";
        write_svg_symmetry(out, item.x, item.y + 29, level.metadata.symmetry, "#68778b", 9);
        if (level.annotation.delocalised_pi.available) {
            write_svg_delocalised_pi(out, item.x, item.y + 42,
                                     level.annotation.delocalised_pi, "#5b6575");
        }
    }

    out << "<text x=\"36\" y=\"" << height - 24
        << "\" fill=\"#8793a3\" font-family=\"sans-serif\" font-size=\"10\">"
           "MO numbering is intentionally omitted from the figure; exact/raw identifiers remain in interactive hover and machine metadata.</text>\n";
    out << "</svg>\n";
    if (!out) {
        if (error) *error = "Failed while writing SVG";
        return false;
    }
    return true;
}