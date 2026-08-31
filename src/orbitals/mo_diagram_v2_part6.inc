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

constexpr const char* kSvgFontFamily =
    "Segoe UI, Microsoft YaHei UI, Yu Gothic UI, Noto Sans, sans-serif";

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
            << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\"" << font_size
            << "\">N/A</text>\n";
        return;
    }

    const double suffix_size = std::max(6.0, static_cast<double>(font_size) * 0.70);
    const double base_width = std::max(1.0, static_cast<double>(notation.base.size())) *
                              static_cast<double>(font_size) * 0.58;
    const double suffix_width = std::max(notation.subscript.size(), notation.superscript.size()) *
                                suffix_size * 0.55;
    const double start_x = x - (base_width + suffix_width) * 0.5;
    const double suffix_x = start_x + base_width + 0.5;

    // Use separate text nodes instead of SVG baseline-shift. This renders
    // consistently in browsers, CairoSVG and print/PDF pipelines while still
    // giving textbook subscript/superscript placement.
    out << "<text x=\"" << start_x << "\" y=\"" << y
        << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\""
        << font_size << "\">" << xml_escape(notation.base) << "</text>\n";
    if (!notation.subscript.empty()) {
        out << "<text x=\"" << suffix_x << "\" y=\"" << y + 3.5
            << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\""
            << suffix_size << "\">" << xml_escape(notation.subscript) << "</text>\n";
    }
    if (!notation.superscript.empty()) {
        out << "<text x=\"" << suffix_x << "\" y=\"" << y - 5.0
            << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\""
            << suffix_size << "\">" << xml_escape(notation.superscript) << "</text>\n";
    }
}

void write_svg_electron(std::ostream& out,
                        const double x,
                        const double level_anchor_y,
                        const bool up) {
    // Textbook convention: both spin arrows occupy the space above the energy
    // line. Alpha points upward from the near-line anchor. Beta starts above
    // the line and points downward *towards* the same anchor, so neither arrow
    // can collide with the orbital-type/symmetry text below the level.
    constexpr double length = 11.0;
    constexpr double head_dx = 3.2;
    constexpr double head_dy = 3.5;
    const double tail = up ? level_anchor_y : level_anchor_y - length;
    const double tip = up ? level_anchor_y - length : level_anchor_y;
    const double wing_y = up ? tip + head_dy : tip - head_dy;
    out << "<line x1=\"" << x << "\" y1=\"" << tail
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
    const double base_x = x - 4.0;
    out << "<g data-pi-topology-available=\""
        << (pi.topology_available ? "true" : "false")
        << "\" data-pi-orientation-channels=\"";
    if (pi.topology_available) out << pi.orientation_channels;
    out << "\" data-pi-cyclic=\"";
    if (pi.topology_available) out << (pi.cyclic_topology ? "true" : "false");
    out << "\">\n<title>delocalised pi; ";
    if (pi.topology_available) {
        out << pi.orientation_channels << " orientation channel(s); "
            << (pi.cyclic_topology ? "cyclic" : "non-cyclic");
    } else {
        out << "topology unavailable";
    }
    out << "</title>\n";
    out << "<text x=\"" << base_x << "\" y=\"" << y
        << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\"10\">Π</text>\n";
    if (pi.participating_atoms > 0) {
        out << "<text x=\"" << base_x + 7.0 << "\" y=\"" << y - 5.0
            << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\"7\">"
            << pi.participating_atoms << "</text>\n";
    }
    if (pi.participating_electrons > 0.0) {
        out << "<text x=\"" << base_x + 7.0 << "\" y=\"" << y + 4.0
            << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily << "\" font-size=\"7\">"
            << static_cast<int>(std::lround(pi.participating_electrons)) << "</text>\n";
    }
    if (pi.topology_available) {
        out << "<text x=\"" << base_x + 16.0 << "\" y=\"" << y
            << "\" fill=\"" << colour << "\" font-family=\"" << kSvgFontFamily
            << "\" font-size=\"7\">" << pi.orientation_channels << "ch"
            << (pi.cyclic_topology ? " C" : " L") << "</text>\n";
    }
    out << "</g>\n";
}

struct ExportPiColour {
    const char* hex = "#e25ddc";
    std::uint8_t r = 226;
    std::uint8_t g = 93;
    std::uint8_t b = 220;
};

ExportPiColour export_pi_colour(const PiInteractionKind kind) noexcept {
    switch (kind) {
        case PiInteractionKind::Donor: return {"#f2a340", 242, 163, 64};
        case PiInteractionKind::Acceptor: return {"#46ceda", 70, 206, 218};
        case PiInteractionKind::WeakNearNonbonding:
            return {"#ae9abe", 174, 154, 190};
        default: return {};
    }
}

const char* export_pi_kind_id(const PiInteractionKind kind) noexcept {
    switch (kind) {
        case PiInteractionKind::Donor: return "donor";
        case PiInteractionKind::Acceptor: return "acceptor";
        case PiInteractionKind::WeakNearNonbonding: return "weak-near-nonbonding";
        default: return "coupled";
    }
}

struct ExportPiBracket {
    double top = 0.0;
    double bottom = 0.0;
};

std::optional<ExportPiBracket> export_pi_bracket(
    const PiInteractionDescriptor& interaction,
    const std::vector<PositionedLevel>& positioned) {
    if (interaction.retained_level >= positioned.size()) return std::nullopt;

    ExportPiBracket bracket;
    if (interaction.kind == PiInteractionKind::WeakNearNonbonding ||
        !interaction.lower_visible || !interaction.upper_visible ||
        interaction.lower_level >= positioned.size() ||
        interaction.upper_level >= positioned.size()) {
        const double y = positioned[interaction.retained_level].y;
        bracket.top = y - 7.0;
        bracket.bottom = y + 7.0;
        return bracket;
    }

    bracket.top = std::min(positioned[interaction.lower_level].y,
                           positioned[interaction.upper_level].y) - 5.0;
    bracket.bottom = std::max(positioned[interaction.lower_level].y,
                              positioned[interaction.upper_level].y) + 5.0;
    if (bracket.bottom - bracket.top < 14.0) {
        const double middle = 0.5 * (bracket.top + bracket.bottom);
        bracket.top = middle - 7.0;
        bracket.bottom = middle + 7.0;
    }
    return bracket;
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
         << width << ' ' << height << "\" data-diagram-row-count=\""
         << data.levels.size() << "\" data-pi-interaction-count=\""
         << data.pi_interactions.size() << "\">\n";
    out << "<rect width=\"100%\" height=\"100%\" fill=\"#fbfcfe\"/>\n";
    out << "<text x=\"36\" y=\"40\" fill=\"#182231\" font-family=\"" << kSvgFontFamily << "\" "
           "font-size=\"24\" font-weight=\"600\">Valence MO diagram</text>\n";
    out << "<text x=\"36\" y=\"64\" fill=\"#68778b\" font-family=\"" << kSvgFontFamily << "\" "
           "font-size=\"12\">" << xml_escape(data.selection.summary) << "</text>\n";
    out << "<text x=\"36\" y=\"82\" fill=\"#7f8c9c\" font-family=\"" << kSvgFontFamily << "\" "
           "font-size=\"11\">Molden-derived MO energies/occupations/symmetry; type/class is N/A unless supported by source data.</text>\n";

    const int axis_x = 70;
    out << "<line x1=\"" << axis_x << "\" y1=\"110\" x2=\"" << axis_x
        << "\" y2=\"" << height - 62 << "\" stroke=\"#2d3745\" stroke-width=\"2\"/>\n";
    out << "<path d=\"M " << axis_x << " 98 l -5 10 h 10 z\" fill=\"#2d3745\"/>\n";
    out << "<text x=\"24\" y=\"124\" fill=\"#4d5a6c\" font-family=\"" << kSvgFontFamily << "\" "
           "font-size=\"12\" transform=\"rotate(-90 24 124)\">Energy · "
        << (options.energy_axis_mode == EnergyAxisMode::Linear ? "linear" : "adaptive nonlinear (log-gap v3)")
        << " (" << xml_escape(energy_unit_symbol(options.energy_unit)) << ")</text>\n";
    if (options.energy_axis_mode == EnergyAxisMode::NonlinearFocus) {
        out << "<text x=\"82\" y=\"102\" fill=\"#8b6173\" font-family=\"" << kSvgFontFamily << "\" "
               "font-size=\"10\">nonlinear axis · ticks show exact physical energy</text>\n";
    }

    for (const double fraction : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        const double energy = energy_at_axis_fraction(fraction, data.energy_transform);
        const double y = map_energy_y(energy, data.energy_transform, height);
        out << "<line x1=\"64\" y1=\"" << y << "\" x2=\"76\" y2=\"" << y
            << "\" stroke=\"#7c8795\" stroke-width=\"1\"/>\n";
        out << "<text x=\"80\" y=\"" << y + 4
            << "\" fill=\"#7c8795\" font-family=\"" << kSvgFontFamily << "\" font-size=\"10\">"
            << xml_escape(format_energy(energy, options.energy_unit,
                   options.energy_unit == EnergyUnit::JoulePerMol ? 0 : 4))
            << "</text>\n";
    }

    for (std::size_t row_index = 0; row_index < positioned.size(); ++row_index) {
        const auto& item = positioned[row_index];
        const auto& level = *item.level;
        const char* colour = level.metadata.selected ? "#286de0" :
                             level.homo ? "#138a62" :
                             level.lumo ? "#b56a00" :
                             level.metadata.region == OrbitalRegion::Virtual ? "#778496" : "#29384c";
        const double x0 = item.x - 45.0;
        const double x1 = item.x + 45.0;
        out << "<line class=\"mo-level\" data-row-index=\"" << row_index
            << "\" x1=\"" << x0 << "\" y1=\"" << item.y
            << "\" x2=\"" << x1 << "\" y2=\"" << item.y
            << "\" stroke=\"" << colour << "\" stroke-width=\""
            << (level.metadata.selected ? 4 : 2.2) << "\"/>\n";
        if (level.electrons.alpha > 0) write_svg_electron(out, item.x - 10.0, item.y - 2.0, true);
        if (level.electrons.beta > 0) write_svg_electron(out, item.x + 10.0, item.y - 2.0, false);

        // Export deliberately omits MO numbers. Human-readable orbital type is
        // attached directly below the physical level; raw numbering remains in
        // hover/JSON/CSV where it is useful without cluttering the figure.
        out << "<text x=\"" << item.x << "\" y=\"" << item.y + 16
            << "\" text-anchor=\"middle\" fill=\"#435166\" font-family=\"" << kSvgFontFamily << "\" "
               "font-size=\"10\">" << xml_escape(orbital_type_svg(level.annotation)) << "</text>\n";
        write_svg_symmetry(out, item.x, item.y + 29, level.metadata.symmetry, "#68778b", 9);
        if (level.annotation.delocalised_pi.available) {
            write_svg_delocalised_pi(out, item.x, item.y + 42,
                                     level.annotation.delocalised_pi, "#5b6575");
        }
    }

    out << "<g id=\"pi-interactions\">\n";
    std::size_t rendered_interactions = 0;
    for (std::size_t interaction_index = 0;
         interaction_index < data.pi_interactions.size(); ++interaction_index) {
        const auto& interaction = data.pi_interactions[interaction_index];
        const auto bracket = export_pi_bracket(interaction, positioned);
        if (!bracket) continue;
        const ExportPiColour colour = export_pi_colour(interaction.kind);
        const double x = static_cast<double>(width) - 26.0 -
                         static_cast<double>(rendered_interactions) * 48.0;
        constexpr double cap = 7.0;
        const double middle = 0.5 * (bracket->top + bracket->bottom);
        out << "<g class=\"pi-interaction\" data-interaction-index=\""
            << interaction_index << "\" data-kind=\""
            << export_pi_kind_id(interaction.kind) << "\">\n";
        out << "<title>" << xml_escape(pi_interaction_kind_name(interaction.kind));
        if (!interaction.symmetry.empty()) {
            out << " · " << xml_escape(interaction.symmetry);
        }
        out << " · Δπ=" << std::setprecision(10)
            << interaction.splitting_hartree << " Ha · confidence="
            << std::setprecision(4) << interaction.confidence << "</title>\n";
        out << "<line x1=\"" << x << "\" y1=\"" << bracket->top
            << "\" x2=\"" << x << "\" y2=\"" << bracket->bottom
            << "\" stroke=\"" << colour.hex << "\" stroke-width=\"2\"/>\n";
        out << "<line x1=\"" << x - cap << "\" y1=\"" << bracket->top
            << "\" x2=\"" << x << "\" y2=\"" << bracket->top
            << "\" stroke=\"" << colour.hex << "\" stroke-width=\"2\"/>\n";
        out << "<line x1=\"" << x - cap << "\" y1=\"" << bracket->bottom
            << "\" x2=\"" << x << "\" y2=\"" << bracket->bottom
            << "\" stroke=\"" << colour.hex << "\" stroke-width=\"2\"/>\n";
        out << "<text x=\"" << x - cap - 5.0 << "\" y=\"" << middle + 4.0
            << "\" text-anchor=\"end\" fill=\"" << colour.hex
            << "\" font-family=\"" << kSvgFontFamily
            << "\" font-size=\"13\">Δπ</text>\n";
        out << "</g>\n";
        ++rendered_interactions;
    }
    out << "</g>\n";

    out << "<text x=\"36\" y=\"" << height - 24
        << "\" fill=\"#8793a3\" font-family=\"" << kSvgFontFamily << "\" font-size=\"10\">"
           "MO numbering is intentionally omitted from the figure; exact/raw identifiers remain in interactive hover and machine metadata.</text>\n";
    out << "</svg>\n";
    if (!out) {
        if (error) *error = "Failed while writing SVG";
        return false;
    }
    return true;
}
