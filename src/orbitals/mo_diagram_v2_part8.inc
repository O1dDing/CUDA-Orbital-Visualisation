namespace {

const char* machine_bonding_class(const BondingClass value) noexcept {
    switch (value) {
        case BondingClass::Bonding: return "bonding";
        case BondingClass::Nonbonding: return "nonbonding";
        case BondingClass::Antibonding: return "antibonding";
        default: return "unclassified";
    }
}

void write_json_index_array(std::ostream& out,
                            const std::vector<std::size_t>& values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        out << values[i];
    }
    out << ']';
}

std::string csv_index_list(const std::vector<std::size_t>& values) {
    std::ostringstream out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ';';
        out << values[i];
    }
    return out.str();
}

std::string csv_pi_interaction_summary(const MODiagramData& data) {
    std::ostringstream out;
    out << std::setprecision(14);
    for (std::size_t i = 0; i < data.pi_interactions.size(); ++i) {
        if (i) out << " || ";
        const auto& interaction = data.pi_interactions[i];
        out << "kind=" << export_pi_kind_id(interaction.kind)
            << ";symmetry=" << interaction.symmetry
            << ";lower_level=" << interaction.lower_level
            << ";upper_level=" << interaction.upper_level
            << ";retained_level=" << interaction.retained_level
            << ";lower_orbitals=" << csv_index_list(interaction.lower_orbitals)
            << ";upper_orbitals=" << csv_index_list(interaction.upper_orbitals)
            << ";splitting_hartree=" << interaction.splitting_hartree
            << ";confidence=" << interaction.confidence
            << ";lower_visible=" << (interaction.lower_visible ? 1 : 0)
            << ";upper_visible=" << (interaction.upper_visible ? 1 : 0);
    }
    return out.str();
}

} // namespace

bool write_mo_diagram_json(const MODiagramData& data,
                           const MODiagramOptions& options,
                           const std::filesystem::path& path,
                           std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open JSON output path";
        return false;
    }

    std::set<std::size_t> included(data.selection.included_indices.begin(),
                                   data.selection.included_indices.end());
    out << "{\n";
    out << "  \"schema\": \"cov.mo-diagram.v2\",\n";
    out << "  \"mode\": \"valence-central\",\n";
    out << "  \"strict_salc_claimed\": false,\n";
    out << "  \"energy_unit\": \"" << json_escape(energy_unit_symbol(options.energy_unit)) << "\",\n";
    out << "  \"energy_axis_mode\": \"" << energy_axis_mode_name(options.energy_axis_mode) << "\",\n";
    out << "  \"transform_name\": \"" << energy_transform_name(options.energy_axis_mode) << "\",\n";
    out << "  \"transform_parameters\": {\n";
    out << "    \"minimum_gap_fraction\": " << std::setprecision(12)
        << data.energy_transform.minimum_gap_weight << ",\n";
    out << "    \"scale_hartree\": " << data.energy_transform.scale_hartree << ",\n";
    out << "    \"focus_hartree\": " << data.energy_transform.focus_hartree << ",\n";
    out << "    \"knots\": [";
    for (std::size_t i = 0; i < data.energy_transform.knots.size(); ++i) {
        if (i) out << ',';
        const auto& knot = data.energy_transform.knots[i];
        out << "{\"energy_hartree\":" << knot.energy_hartree
            << ",\"coordinate\":" << knot.coordinate << '}';
    }
    out << "]\n  },\n";
    out << "  \"selection_summary\": \"" << json_escape(data.selection.summary) << "\",\n";
    out << "  \"shown_orbitals\": " << data.selection.included_indices.size() << ",\n";
    out << "  \"diagram_row_count\": " << data.levels.size() << ",\n";
    out << "  \"hidden_orbitals\": " << data.selection.hidden_count << ",\n";
    out << "  \"valence_occupied_count\": " << data.selection.valence_occupied_count << ",\n";
    out << "  \"frontier_virtual_count\": " << data.selection.frontier_virtual_count << ",\n";
    out << "  \"local_coordination\": {\n";
    out << "    \"available\": "
        << (!data.ligand_field_geometry_id.empty() &&
            data.ligand_field_coordination_number > 0 ? "true" : "false") << ",\n";
    out << "    \"point_group\": \""
        << json_escape(data.ligand_field_point_group) << "\",\n";
    out << "    \"geometry_id\": \""
        << json_escape(data.ligand_field_geometry_id) << "\",\n";
    out << "    \"geometry_name\": \""
        << json_escape(data.ligand_field_geometry_name) << "\",\n";
    out << "    \"coordination_number\": "
        << data.ligand_field_coordination_number << ",\n";
    out << "    \"metal_atom_index\": " << data.ligand_field_metal_atom << ",\n";
    out << "    \"ligand_atom_indices\": ";
    write_json_index_array(out, data.ligand_field_ligand_atoms);
    out << ",\n";
    out << "    \"confidence\": " << std::setprecision(14)
        << data.ligand_field_confidence << ",\n";
    out << "    \"angular_rms\": " << data.ligand_field_angular_rms << ",\n";
    out << "    \"shape_measure\": " << data.ligand_field_shape_measure << ",\n";
    out << "    \"radial_cv\": " << data.ligand_field_radial_cv << "\n";
    out << "  },\n";
    out << "  \"pi_interactions\": [\n";
    for (std::size_t i = 0; i < data.pi_interactions.size(); ++i) {
        const auto& interaction = data.pi_interactions[i];
        out << "    {\n";
        out << "      \"kind\": \""
            << export_pi_kind_id(interaction.kind) << "\",\n";
        out << "      \"description\": \""
            << json_escape(pi_interaction_kind_name(interaction.kind)) << "\",\n";
        out << "      \"symmetry\": \""
            << json_escape(interaction.symmetry) << "\",\n";
        out << "      \"lower_level\": " << interaction.lower_level << ",\n";
        out << "      \"upper_level\": " << interaction.upper_level << ",\n";
        out << "      \"retained_level\": " << interaction.retained_level << ",\n";
        out << "      \"lower_orbitals\": ";
        write_json_index_array(out, interaction.lower_orbitals);
        out << ",\n";
        out << "      \"upper_orbitals\": ";
        write_json_index_array(out, interaction.upper_orbitals);
        out << ",\n";
        out << "      \"splitting_hartree\": " << std::setprecision(14)
            << interaction.splitting_hartree << ",\n";
        out << "      \"splitting_display\": "
            << convert_hartree(interaction.splitting_hartree, options.energy_unit)
            << ",\n";
        out << "      \"confidence\": " << interaction.confidence << ",\n";
        out << "      \"lower_visible\": "
            << (interaction.lower_visible ? "true" : "false") << ",\n";
        out << "      \"upper_visible\": "
            << (interaction.upper_visible ? "true" : "false") << "\n";
        out << "    }" << (i + 1 < data.pi_interactions.size() ? "," : "") << "\n";
    }
    out << "  ],\n";
    out << "  \"diagram_rows\": [\n";
    for (std::size_t i = 0; i < data.levels.size(); ++i) {
        const auto& level = data.levels[i];
        out << "    {\"row_index\":" << i
            << ",\"layout_energy_hartree\":" << std::setprecision(14)
            << level.layout_energy_hartree
            << ",\"symmetry\":\"" << json_escape(level.metadata.symmetry)
            << "\",\"member_orbitals\":";
        if (level.member_indices.empty()) {
            write_json_index_array(out, std::vector<std::size_t>{level.metadata.orbital_index});
        } else {
            write_json_index_array(out, level.member_indices);
        }
        out << '}' << (i + 1 < data.levels.size() ? "," : "") << "\n";
    }
    out << "  ],\n";
    out << "  \"orbitals\": [\n";

    bool first = true;
    for (std::size_t i = 0; i < data.metadata.size(); ++i) {
        if (!options.include_hidden_in_metadata && !included.count(i)) continue;
        if (!first) out << ",\n";
        first = false;
        const auto& item = data.metadata[i];
        const auto& annotation = i < data.annotations.size() ? data.annotations[i] : OrbitalAnnotation{};
        const double coordinate = energy_display_coordinate(
            group_layout_energy(data.metadata, i), data.energy_transform);
        out << "    {\n";
        out << "      \"index\": " << item.orbital_index << ",\n";
        out << "      \"raw_mo\": " << item.raw_mo_number << ",\n";
        out << "      \"label\": \"" << json_escape(item.display_label) << "\",\n";
        out << "      \"energy_hartree\": " << std::setprecision(14) << item.energy_hartree << ",\n";
        out << "      \"source_energy_hartree\": " << item.energy_hartree << ",\n";
        out << "      \"energy_display\": " << convert_hartree(item.energy_hartree, options.energy_unit) << ",\n";
        out << "      \"display_coordinate\": " << coordinate << ",\n";
        out << "      \"occupation\": " << item.occupation << ",\n";
        out << "      \"spin\": \"" << spin_name(item.spin) << "\",\n";
        out << "      \"symmetry\": \"" << json_escape(item.symmetry) << "\",\n";
        out << "      \"symmetry_display\": \"" << json_escape(format_symmetry_unicode(item.symmetry)) << "\",\n";
        out << "      \"degeneracy_size\": " << item.degeneracy_size << ",\n";
        out << "      \"region\": \"" << region_name(item.region) << "\",\n";
        out << "      \"included_in_diagram\": " << (included.count(i) ? "true" : "false") << ",\n";
        out << "      \"visible\": " << (item.visible ? "true" : "false") << ",\n";
        out << "      \"selected\": " << (item.selected ? "true" : "false") << ",\n";
        out << "      \"orbital_family\": \"" << json_escape(annotation.family) << "\",\n";
        out << "      \"family_source\": \"" << annotation_source_name(annotation.family_source) << "\",\n";
        out << "      \"family_confidence\": " << annotation.family_confidence << ",\n";
        out << "      \"bonding_class\": \"" << machine_bonding_class(annotation.bonding_class) << "\",\n";
        out << "      \"bonding_class_source\": \"" << annotation_source_name(annotation.bonding_source) << "\",\n";
        out << "      \"bonding_confidence\": " << annotation.bonding_confidence << ",\n";
        out << "      \"multicentre_available\": " << (annotation.multicentre.available ? "true" : "false") << ",\n";
        out << "      \"multicentre_centres\": " << annotation.multicentre.centres << ",\n";
        out << "      \"multicentre_electrons\": " << annotation.multicentre.electrons << ",\n";
        out << "      \"multicentre_label\": \"" << json_escape(annotation.multicentre.label) << "\",\n";
        out << "      \"multicentre_source\": \"" << annotation_source_name(annotation.multicentre.source) << "\",\n";
        out << "      \"multicentre_confidence\": " << annotation.multicentre.confidence << ",\n";
        out << "      \"multicentre_heuristic\": " << (annotation.multicentre.heuristic ? "true" : "false") << ",\n";
        out << "      \"delocalised_pi_available\": " << (annotation.delocalised_pi.available ? "true" : "false") << ",\n";
        out << "      \"delocalised_pi_participating_atoms\": " << annotation.delocalised_pi.participating_atoms << ",\n";
        out << "      \"delocalised_pi_participating_electrons\": " << annotation.delocalised_pi.participating_electrons << ",\n";
        out << "      \"delocalised_pi_label\": \"" << json_escape(annotation.delocalised_pi.label) << "\",\n";
        out << "      \"delocalised_pi_source\": \"" << annotation_source_name(annotation.delocalised_pi.source) << "\",\n";
        out << "      \"delocalised_pi_confidence\": " << annotation.delocalised_pi.confidence << ",\n";
        out << "      \"delocalised_pi_heuristic\": " << (annotation.delocalised_pi.heuristic ? "true" : "false") << ",\n";
        out << "      \"heuristic\": " << (annotation.heuristic ? "true" : "false") << "\n";
        out << "    }";
    }
    out << "\n  ]\n}\n";
    if (!out) {
        if (error) *error = "Failed while writing JSON";
        return false;
    }
    return true;
}

bool write_mo_diagram_csv(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        if (error) *error = "Unable to open CSV output path";
        return false;
    }
    std::set<std::size_t> included(data.selection.included_indices.begin(),
                                   data.selection.included_indices.end());
    out << "index,raw_mo,label,energy_hartree,energy_display,energy_unit,display_coordinate,"
           "energy_axis_mode,transform_name,symmetry,symmetry_display,occupation,spin,"
           "degeneracy_size,region,included_in_diagram,visible,selected,orbital_family,"
           "family_source,family_confidence,bonding_class,bonding_source,bonding_confidence,"
           "multicentre_available,multicentre_centres,multicentre_electrons,multicentre_label,"
           "multicentre_source,multicentre_confidence,multicentre_heuristic,"
           "delocalised_pi_available,delocalised_pi_participating_atoms,"
           "delocalised_pi_participating_electrons,delocalised_pi_label,delocalised_pi_source,"
           "delocalised_pi_confidence,delocalised_pi_heuristic,heuristic,"
           "local_point_group,local_geometry_id,local_geometry_name,local_coordination_number,"
           "local_metal_atom_index,local_ligand_atom_indices,local_confidence,local_angular_rms,"
           "local_shape_measure,local_radial_cv,diagram_row_count,pi_interaction_count,"
           "pi_interactions\n";
    const std::string ligand_atoms = csv_index_list(data.ligand_field_ligand_atoms);
    const std::string pi_interactions = csv_pi_interaction_summary(data);
    for (std::size_t i = 0; i < data.metadata.size(); ++i) {
        if (!options.include_hidden_in_metadata && !included.count(i)) continue;
        const auto& item = data.metadata[i];
        const auto& annotation = i < data.annotations.size() ? data.annotations[i] : OrbitalAnnotation{};
        const double coordinate = energy_display_coordinate(
            group_layout_energy(data.metadata, i), data.energy_transform);
        out << item.orbital_index << ',' << item.raw_mo_number << ','
            << csv_escape(item.display_label) << ',' << std::setprecision(14) << item.energy_hartree << ','
            << convert_hartree(item.energy_hartree, options.energy_unit) << ','
            << csv_escape(energy_unit_symbol(options.energy_unit)) << ',' << coordinate << ','
            << energy_axis_mode_name(options.energy_axis_mode) << ','
            << energy_transform_name(options.energy_axis_mode) << ','
            << csv_escape(item.symmetry) << ',' << csv_escape(format_symmetry_unicode(item.symmetry)) << ','
            << item.occupation << ',' << spin_name(item.spin) << ',' << item.degeneracy_size << ','
            << region_name(item.region) << ',' << (included.count(i) ? 1 : 0) << ','
            << (item.visible ? 1 : 0) << ',' << (item.selected ? 1 : 0) << ','
            << csv_escape(annotation.family) << ',' << annotation_source_name(annotation.family_source) << ','
            << annotation.family_confidence << ',' << machine_bonding_class(annotation.bonding_class) << ','
            << annotation_source_name(annotation.bonding_source) << ',' << annotation.bonding_confidence << ','
            << (annotation.multicentre.available ? 1 : 0) << ',' << annotation.multicentre.centres << ','
            << annotation.multicentre.electrons << ',' << csv_escape(annotation.multicentre.label) << ','
            << annotation_source_name(annotation.multicentre.source) << ',' << annotation.multicentre.confidence << ','
            << (annotation.multicentre.heuristic ? 1 : 0) << ','
            << (annotation.delocalised_pi.available ? 1 : 0) << ','
            << annotation.delocalised_pi.participating_atoms << ','
            << annotation.delocalised_pi.participating_electrons << ','
            << csv_escape(annotation.delocalised_pi.label) << ','
            << annotation_source_name(annotation.delocalised_pi.source) << ','
            << annotation.delocalised_pi.confidence << ','
            << (annotation.delocalised_pi.heuristic ? 1 : 0) << ','
            << (annotation.heuristic ? 1 : 0) << ','
            << csv_escape(data.ligand_field_point_group) << ','
            << csv_escape(data.ligand_field_geometry_id) << ','
            << csv_escape(data.ligand_field_geometry_name) << ','
            << data.ligand_field_coordination_number << ','
            << data.ligand_field_metal_atom << ','
            << csv_escape(ligand_atoms) << ','
            << data.ligand_field_confidence << ','
            << data.ligand_field_angular_rms << ','
            << data.ligand_field_shape_measure << ','
            << data.ligand_field_radial_cv << ','
            << data.levels.size() << ','
            << data.pi_interactions.size() << ','
            << csv_escape(pi_interactions) << '\n';
    }
    if (!out) {
        if (error) *error = "Failed while writing CSV";
        return false;
    }
    return true;
}

MODiagramExportResult export_mo_diagram_bundle(const Wavefunction& wavefunction,
                                               const MODiagramOptions& options,
                                               const std::filesystem::path& base_path) {
    MODiagramExportResult result;
    const MODiagramData data = build_mo_diagram_data(wavefunction, options);
    const std::filesystem::path svg = base_path.string() + ".mo.svg";
    const std::filesystem::path png = base_path.string() + ".mo.png";
    const std::filesystem::path json = base_path.string() + ".mo.json";
    const std::filesystem::path csv = base_path.string() + ".mo.csv";
    std::string error;
    result.svg = write_mo_diagram_svg(data, options, svg, &error);
    if (!result.svg) { result.error = error; return result; }
    result.png = write_mo_diagram_png(data, options, png, &error);
    if (!result.png) { result.error = error; return result; }
    result.json = write_mo_diagram_json(data, options, json, &error);
    if (!result.json) { result.error = error; return result; }
    result.csv = write_mo_diagram_csv(data, options, csv, &error);
    if (!result.csv) result.error = error;
    return result;
}

} // namespace cov
