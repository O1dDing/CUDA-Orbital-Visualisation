bool write_mo_diagram_json(const MODiagramData& data,
                           const MODiagramOptions& options,
                           const std::filesystem::path& path,
                           std::string* error) {
    return cov_legacy_write_mo_diagram_json(data, options, path, error);
}

bool write_mo_diagram_csv(const MODiagramData& data,
                          const MODiagramOptions& options,
                          const std::filesystem::path& path,
                          std::string* error) {
    return cov_legacy_write_mo_diagram_csv(data, options, path, error);
}

MODiagramExportResult export_mo_diagram_bundle(const Wavefunction& wavefunction,
                                                const MODiagramOptions& options,
                                                const std::filesystem::path& base_path) {
    MODiagramExportResult result;
    const MODiagramData data = build_mo_diagram_data(wavefunction, options);
    std::filesystem::path base = base_path;
    if (base.has_extension()) base.replace_extension();
    std::string error;
    result.svg = write_mo_diagram_svg(data, options, base.string() + ".mo.svg", &error);
    if (!result.svg) result.error = error;
    error.clear();
    result.png = write_mo_diagram_png(data, options, base.string() + ".mo.png", &error);
    if (!result.png && result.error.empty()) result.error = error;
    error.clear();
    result.json = write_mo_diagram_json(data, options, base.string() + ".mo.json", &error);
    if (!result.json && result.error.empty()) result.error = error;
    error.clear();
    result.csv = write_mo_diagram_csv(data, options, base.string() + ".mo.csv", &error);
    if (!result.csv && result.error.empty()) result.error = error;
    return result;
}

} // namespace cov
