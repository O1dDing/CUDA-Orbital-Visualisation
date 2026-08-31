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
    auto& derived_pi=data.annotations[0].delocalised_pi;
    derived_pi.available=true;
    derived_pi.participating_atoms=6u;
    derived_pi.participating_electrons=6.0;
    derived_pi.atom_indices={0,1,2,3,4,5};
    derived_pi.orbital_indices={0};
    derived_pi.family_id="pi-family-derived";
    derived_pi.topology_available=true;
    derived_pi.orientation_channels=2u;
    derived_pi.cyclic_topology=true;
    derived_pi.orientation_channel_details={
        {{0,1,2,3,4,5},{0.0,0.0,1.0},0.98,true},
        {{0,1},{1.0,0.0,0.0},0.94,false},
    };
    derived_pi.label="PI6_6";
    derived_pi.source=cov::AnnotationSource::Derived;
    auto& parsed_pi=data.annotations[1].delocalised_pi;
    parsed_pi.available=true;
    parsed_pi.family_id="pi-family-parsed-only";
    parsed_pi.label="delocalised-pi";
    parsed_pi.source=cov::AnnotationSource::ParsedLabel;
    auto& shared_multicentre=data.annotations[0].multicentre;
    shared_multicentre.available=true;
    shared_multicentre.label="2x3c2e (shared 4c/4e source)";
    shared_multicentre.centres=4u;
    shared_multicentre.electrons=4.0;
    shared_multicentre.atom_indices={0,1,2,3};
    shared_multicentre.channel_count=2u;
    shared_multicentre.source_subspace_id="mc-source-bridge-pair";
    shared_multicentre.source_subspace_electron_count=4.0;

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
    cov::LocalGeometryDiagramDescriptor local_geometry;
    local_geometry.centre_atom=4u;
    local_geometry.geometry_id="OC-6";
    local_geometry.geometry_name="Octahedron";
    local_geometry.point_group="Oh";
    local_geometry.neighbour_atoms={0,1,2,3,5,6};
    local_geometry.confidence=0.93;
    local_geometry.angular_rms=0.0125;
    local_geometry.shape_measure=0.018;
    local_geometry.radial_cv=0.021;
    data.local_geometries.push_back(local_geometry);
    data.electronic_state.source=cov::WavefunctionSource::Fchk;
    data.electronic_state.charge=-1;
    data.electronic_state.multiplicity=2u;
    data.electronic_state.alpha_electrons=5u;
    data.electronic_state.beta_electrons=4u;
    data.electronic_state.charge_provenance=cov::DataProvenance::Producer;
    data.electronic_state.multiplicity_provenance=cov::DataProvenance::Producer;
    data.electronic_state.electron_counts_provenance=cov::DataProvenance::Producer;
    data.electronic_state.scf_convergence=cov::ScfConvergenceStatus::Converged;
    data.electronic_state.scf_convergence_provenance=cov::DataProvenance::Producer;
    data.electronic_state.stability=cov::WavefunctionStabilityStatus::Stable;
    data.electronic_state.stability_provenance=cov::DataProvenance::Producer;
    data.electronic_state.stability_detail="stable internal test";
    data.electronic_state.spin_squared_before=0.76;
    data.electronic_state.spin_squared_after=0.75;
    data.electronic_state.spin_squared_provenance=cov::DataProvenance::Producer;
    data.electronic_state.atomic_partial_charge_scheme="Mulliken";
    data.electronic_state.atomic_partial_charges={-0.2,0.2};
    data.electronic_state.atomic_partial_charge_provenance=
        cov::DataProvenance::Producer;
    data.electronic_state.point_group_detected="Oh";
    data.electronic_state.point_group_used="Oh";
    data.electronic_state.point_group_provenance=cov::DataProvenance::Producer;
    data.electronic_state.source_title="export semantics fixture";
    data.electronic_state.source_route="#p test route";
    data.electronic_state.enrichment_source="fixture.log";

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
        svg.find("data-pi-topology-available=\"true\"") == std::string::npos ||
        svg.find("data-pi-orientation-channels=\"2\"") == std::string::npos ||
        svg.find("data-pi-cyclic=\"true\"") == std::string::npos ||
        svg.find("data-pi-topology-available=\"false\"") == std::string::npos ||
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
        json.find("\"local_molecular_geometries\": [") == std::string::npos ||
        json.find("\"centre_atom_index\":4") == std::string::npos ||
        json.find("\"neighbour_atom_indices\":[0,1,2,3,5,6]") ==
            std::string::npos ||
        json.find("\"electronic_state\": {") == std::string::npos ||
        json.find("\"charge\": -1") == std::string::npos ||
        json.find("\"multiplicity\": 2") == std::string::npos ||
        json.find("\"alpha_electrons\": 5") == std::string::npos ||
        json.find("\"beta_electrons\": 4") == std::string::npos ||
        json.find("\"scf_convergence\": \"converged\"") == std::string::npos ||
        json.find("\"wavefunction_stability\": \"stable\"") == std::string::npos ||
        json.find("\"atomic_partial_charges\": [-0.2,0.2]") == std::string::npos ||
        json.find("\"pi_interactions\": [") == std::string::npos ||
        json.find("\"kind\": \"acceptor\"") == std::string::npos ||
        json.find("\"splitting_hartree\": 0.3") == std::string::npos ||
        json.find("\"lower_orbitals\": [0,2,4]") == std::string::npos ||
        json.find("\"family_id\":\"pi-family-derived\"") == std::string::npos ||
        json.find("\"topology_available\":true") == std::string::npos ||
        json.find("\"orientation_channels\":2") == std::string::npos ||
        json.find("\"direction\":[0,0,1]") == std::string::npos ||
        json.find("\"delocalised_pi_topology_available\": false") == std::string::npos ||
        json.find("\"delocalised_pi_orientation_channels\": null") == std::string::npos ||
        json.find("\"delocalised_pi_cyclic_topology\": null") == std::string::npos ||
        json.find("\"multicentre_channel_count\": 2") == std::string::npos ||
        json.find("\"multicentre_source_subspace_id\": \"mc-source-bridge-pair\"") ==
            std::string::npos ||
        json.find("\"multicentre_source_subspace_electron_count\": 4") ==
            std::string::npos ||
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
        csv.find("delocalised_pi_participating_electrons,delocalised_pi_label") ==
            std::string::npos ||
        csv.find("pi_interactions,delocalised_pi_topology_available,") ==
            std::string::npos ||
        csv.find("delocalised_pi_orientation_channel_details") ==
            std::string::npos ||
        csv.find("local_geometry_count,local_geometries") == std::string::npos ||
        csv.find("centre=4;geometry_id=OC-6") == std::string::npos ||
        csv.find("wavefunction_source,charge,charge_provenance") ==
            std::string::npos ||
        csv.find("stability_detail,source_title,source_route,enrichment_source") ==
            std::string::npos ||
        csv.find(",FCHK,-1,producer,2,producer,5,4,producer,converged,producer,stable,") ==
            std::string::npos ||
        csv.find(",Oh,OC-6,Octahedron,6,") == std::string::npos ||
        csv.find("kind=acceptor;symmetry=T2g") == std::string::npos ||
        csv.find("splitting_hartree=0.3") == std::string::npos) {
        std::cerr << "CSV local coordination or pi metadata missing\n";
        return 5;
    }
    if (csv.find("multicentre_channel_count,multicentre_source_subspace_id,") ==
            std::string::npos ||
        csv.find(",2,mc-source-bridge-pair,4,stable internal test,"
                 "export semantics fixture,#p test route,fixture.log\n") ==
            std::string::npos) {
        std::cerr << "CSV shared multicentre source metadata missing\n";
        return 6;
    }

    std::error_code ec;
    std::filesystem::remove(svg_path, ec);
    std::filesystem::remove(png_path, ec);
    std::filesystem::remove(json_path, ec);
    std::filesystem::remove(csv_path, ec);
    std::cout << "mo_diagram_export_semantics_smoke ok\n";
    return 0;
}
