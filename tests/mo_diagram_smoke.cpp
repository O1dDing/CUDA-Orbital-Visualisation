#include "cov/mo_diagram.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

int main() {
    cov::Wavefunction wf;
    wf.atoms.resize(4);

    cov::MolecularOrbital core;
    core.energy_hartree = -3.0;
    core.occupation = 2.0f;
    core.symmetry = "core";
    wf.orbitals.push_back(core);

    cov::MolecularOrbital occupied;
    occupied.energy_hartree = -0.5;
    occupied.occupation = 2.0f;
    occupied.symmetry = "sigma_g bonding";
    wf.orbitals.push_back(occupied);

    cov::MolecularOrbital virtual_a;
    virtual_a.energy_hartree = 0.2;
    virtual_a.occupation = 0.0f;
    virtual_a.symmetry = "pi_u";
    wf.orbitals.push_back(virtual_a);

    cov::MolecularOrbital virtual_b = virtual_a;
    virtual_b.energy_hartree = 0.200000004;
    wf.orbitals.push_back(virtual_b);

    cov::MolecularOrbital producer_multicentre;
    producer_multicentre.energy_hartree = 0.45;
    producer_multicentre.occupation = 0.0f;
    producer_multicentre.symmetry = "sigma 3c2e antibonding";
    wf.orbitals.push_back(producer_multicentre);

    cov::MolecularOrbital very_high;
    very_high.energy_hartree = 7.0;
    very_high.occupation = 0.0f;
    very_high.symmetry = "?";
    wf.orbitals.push_back(very_high);

    cov::MODiagramOptions options;
    options.selected_index = 1;
    options.energy_unit = cov::EnergyUnit::ElectronVolt;
    options.width = 760;
    options.height = 720;
    options.filter.core_energy_cutoff_hartree = -1.5;
    options.filter.virtual_window_hartree = 1.5;
    options.degeneracy.tolerance_hartree = 1.0e-5;
    options.max_levels = 5;

    const auto data = cov::build_mo_diagram_data(wf, options);
    if (!data.frontier.homo || !data.frontier.lumo) return 1;
    if (data.mode != cov::MODiagramMode::ValenceCentral || data.levels.size() != 4) {
        std::cerr << "valence-central selection failed\n";
        return 2;
    }
    if (data.levels[1].metadata.display_label != "3-a" ||
        data.levels[2].metadata.display_label != "3-b") {
        std::cerr << "degenerate labels failed\n";
        return 3;
    }
    if (data.levels[0].electrons.alpha != 1 || data.levels[0].electrons.beta != 1) {
        std::cerr << "electron population failed\n";
        return 4;
    }
    if (data.levels[0].annotation.family != "sigma" ||
        data.levels[0].annotation.bonding_class != cov::BondingClass::Bonding ||
        data.levels[1].annotation.family != "pi") {
        std::cerr << "explicit family/bonding annotation failed\n";
        return 5;
    }
    if (!data.levels[3].annotation.multicentre.available ||
        data.levels[3].annotation.multicentre.label != "3c2e" ||
        data.levels[3].annotation.multicentre.source != cov::AnnotationSource::ParsedLabel) {
        std::cerr << "explicit 3c2e producer annotation failed\n";
        return 6;
    }

    const auto temp = std::filesystem::temp_directory_path() / "cov_mo_diagram_smoke";
    const auto result = cov::export_mo_diagram_bundle(wf, options, temp);
    if (!result.svg || !result.png || !result.json || !result.csv) {
        std::cerr << "diagram export failed: " << result.error << '\n';
        return 7;
    }

    const auto png_path = std::filesystem::path(temp.string() + ".mo.png");
    std::ifstream png(png_path, std::ios::binary);
    std::array<std::uint8_t, 8> signature{};
    png.read(reinterpret_cast<char*>(signature.data()), static_cast<std::streamsize>(signature.size()));
    const std::array<std::uint8_t, 8> expected{137,80,78,71,13,10,26,10};
    if (signature != expected) return 8;

    const auto svg_path = std::filesystem::path(temp.string() + ".mo.svg");
    std::ifstream svg_file(svg_path, std::ios::binary);
    std::ostringstream svg_buffer;
    svg_buffer << svg_file.rdbuf();
    const std::string svg = svg_buffer.str();
    if (svg.find("Valence MO diagram") == std::string::npos ||
        svg.find("MO numbering is intentionally omitted") == std::string::npos ||
        svg.find(">3-a<") != std::string::npos) {
        std::cerr << "human export policy failed\n";
        return 9;
    }

    const auto json_path = std::filesystem::path(temp.string() + ".mo.json");
    std::ifstream json_file(json_path, std::ios::binary);
    std::ostringstream json_buffer;
    json_buffer << json_file.rdbuf();
    const std::string json = json_buffer.str();
    if (json.find("\"mode\": \"valence-central\"") == std::string::npos ||
        json.find("\"label\": \"3-a\"") == std::string::npos ||
        json.find("\"multicentre_label\": \"3c2e\"") == std::string::npos ||
        json.find("\"strict_salc_claimed\": false") == std::string::npos) {
        std::cerr << "machine-readable markers missing\n";
        return 10;
    }

    std::error_code ec;
    std::filesystem::remove(svg_path, ec);
    std::filesystem::remove(png_path, ec);
    std::filesystem::remove(json_path, ec);
    std::filesystem::remove(std::filesystem::path(temp.string() + ".mo.csv"), ec);

    std::cout << "mo_diagram_smoke ok\n";
    return 0;
}
