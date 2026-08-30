#include "cov/mo_diagram.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

namespace {

std::string read_all(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
}

std::size_t count_token(const std::string& text, const std::string& token) {
    std::size_t count = 0;
    std::size_t offset = 0;
    while ((offset = text.find(token, offset)) != std::string::npos) {
        ++count;
        offset += token.size();
    }
    return count;
}

cov::OrbitalMetadata metadata(const std::size_t index,
                              const double energy,
                              const char* symmetry) {
    cov::OrbitalMetadata result;
    result.orbital_index = index;
    result.raw_mo_number = index + 1;
    result.display_label = std::to_string(index + 1);
    result.energy_hartree = energy;
    result.occupation = index == 0 ? 2.0f : 0.0f;
    result.symmetry = symmetry;
    result.region = index == 0 ? cov::OrbitalRegion::Valence
                               : cov::OrbitalRegion::Virtual;
    return result;
}

} // namespace

int main() {
    cov::MODiagramData data;
    data.selection.included_indices = {0, 1};
    data.selection.summary = "compact validation rows";
    data.metadata = {metadata(0, -0.20, "T2g"),
                     metadata(1, 0.10, "T2g*")};
    data.annotations.resize(2);
    data.annotations[0].family = "pi";
    data.annotations[1].family = "pi";

    for (std::size_t i = 0; i < data.metadata.size(); ++i) {
        cov::MODiagramLevel level;
        level.metadata = data.metadata[i];
        level.annotation = data.annotations[i];
        level.layout_energy_hartree = data.metadata[i].energy_hartree;
        level.member_indices = {i};
        level.electrons = i == 0 ? cov::ElectronGlyphs{1, 1}
                                 : cov::ElectronGlyphs{};
        data.levels.push_back(level);
    }
    data.energy_transform = cov::build_energy_transform(
        {-0.20, 0.10}, cov::EnergyAxisMode::NonlinearFocus, 0.055);

    data.ligand_field_point_group = "Oh";
    data.ligand_field_geometry_id = "OC-6";
    data.ligand_field_geometry_name = "Octahedron";
    data.ligand_field_coordination_number = 6;
    data.ligand_field_metal_atom = 4;
    data.ligand_field_ligand_atoms = {0, 1, 2, 3, 5, 6};
    data.ligand_field_confidence = 0.93;
    data.ligand_field_angular_rms = 0.0125;
    data.ligand_field_shape_measure = 0.018;
    data.ligand_field_radial_cv = 0.021;

    cov::PiInteractionDescriptor pi;
    pi.lower_level = 0;
    pi.upper_level = 1;
    pi.retained_level = 0;
    pi.lower_orbitals = {0, 2, 4};
    pi.upper_orbitals = {1, 3, 5};
    pi.symmetry = "T2g";
    pi.kind = cov::PiInteractionKind::Acceptor;
    pi.splitting_hartree = 0.30;
    pi.confidence = 0.87;
    data.pi_interactions.push_back(pi);

    cov::MODiagramOptions options;
    options.energy_unit = cov::EnergyUnit::ElectronVolt;
    options.width = 760;
    options.height = 720;
    options.include_hidden_in_metadata = false;

    const auto base = std::filesystem::temp_directory_path() /
                      "cov_mo_diagram_export_semantics_smoke";
    const auto svg_path = std::filesystem::path(base.string() + ".svg");
    const auto png_path = std::filesystem::path(base.string() + ".png");
    const auto json_path = std::filesystem::path(base.string() + ".json");
    const auto csv_path = std::filesystem::path(base.string() + ".csv");
    std::string error;
    if (!cov::write_mo_diagram_svg(data, options, svg_path, &error) ||
        !cov::write_mo_diagram_png(data, options, png_path, &error) ||
        !cov::write_mo_diagram_json(data, options, json_path, &error) ||
        !cov::write_mo_diagram_csv(data, options, csv_path, &error)) {
        std::cerr << "export failed: " << error << '\n';
        return 1;
    }

    const std::string svg = read_all(svg_path);
    if (svg.find("data-diagram-row-count=\"2\"") == std::string::npos ||
        count_token(svg, "class=\"mo-level\"") != data.levels.size() ||
        count_token(svg, "class=\"pi-interaction\"") != 1 ||
        svg.find("data-kind=\"acceptor\"") == std::string::npos ||
        svg.find("splitting; approximately") != std::string::npos) {
        std::cerr << "SVG compact-row or pi-marker semantics mismatch\n";
        return 2;
    }

    const std::string png = read_all(png_path);
    const std::string acceptor_rgb{
        static_cast<char>(70), static_cast<char>(206),
        static_cast<char>(218), static_cast<char>(255)};
    if (png.find(acceptor_rgb) == std::string::npos) {
        std::cerr << "PNG acceptor side marker colour missing\n";
        return 3;
    }

    const std::string json = read_all(json_path);
    if (json.find("\"local_coordination\": {") == std::string::npos ||
        json.find("\"geometry_id\": \"OC-6\"") == std::string::npos ||
        json.find("\"geometry_name\": \"Octahedron\"") == std::string::npos ||
        json.find("\"coordination_number\": 6") == std::string::npos ||
        json.find("\"angular_rms\": 0.0125") == std::string::npos ||
        json.find("\"shape_measure\": 0.018") == std::string::npos ||
        json.find("\"radial_cv\": 0.021") == std::string::npos ||
        json.find("\"pi_interactions\": [") == std::string::npos ||
        json.find("\"kind\": \"acceptor\"") == std::string::npos ||
        json.find("\"splitting_hartree\": 0.3") == std::string::npos ||
        json.find("\"lower_orbitals\": [0,2,4]") == std::string::npos ||
        json.find("\"diagram_row_count\": 2") == std::string::npos) {
        std::cerr << "JSON local coordination or pi metadata missing\n";
        return 4;
    }

    const std::string csv = read_all(csv_path);
    if (csv.find("local_geometry_id,local_geometry_name,local_coordination_number") ==
            std::string::npos ||
        csv.find("local_angular_rms,local_shape_measure,local_radial_cv") ==
            std::string::npos ||
        csv.find("diagram_row_count,pi_interaction_count,pi_interactions") ==
            std::string::npos ||
        csv.find(",Oh,OC-6,Octahedron,6,") == std::string::npos ||
        csv.find("kind=acceptor;symmetry=T2g") == std::string::npos ||
        csv.find("splitting_hartree=0.3") == std::string::npos) {
        std::cerr << "CSV local coordination or pi metadata missing\n";
        return 5;
    }

    std::error_code ec;
    std::filesystem::remove(svg_path, ec);
    std::filesystem::remove(png_path, ec);
    std::filesystem::remove(json_path, ec);
    std::filesystem::remove(csv_path, ec);
    std::cout << "mo_diagram_export_semantics_smoke ok\n";
    return 0;
}
