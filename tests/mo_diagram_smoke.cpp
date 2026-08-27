#include "cov/mo_diagram.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

int main() {
    cov::Wavefunction wf;
    cov::Atom hydrogen;
    hydrogen.symbol = "H";
    hydrogen.atomic_number = 1;
    wf.atoms = {hydrogen, hydrogen};

    cov::MolecularOrbital occupied;
    occupied.energy_hartree = -0.5;
    occupied.occupation = 2.0f;
    occupied.symmetry = "sigma_g";
    wf.orbitals.push_back(occupied);

    cov::MolecularOrbital virtual_mo;
    virtual_mo.energy_hartree = 0.2;
    virtual_mo.occupation = 0.0f;
    virtual_mo.symmetry = "sigma_u";
    wf.orbitals.push_back(virtual_mo);

    cov::MODiagramOptions options;
    options.selected_index = 0;
    options.energy_unit = cov::EnergyUnit::ElectronVolt;
    options.width = 700;
    options.height = 520;

    const auto data = cov::build_mo_diagram_data(wf, options);
    if (data.levels.size() != 2 || !data.frontier.homo || !data.frontier.lumo) {
        std::cerr << "diagram model failed\n";
        return 1;
    }
    if (!data.textbook_diatomic || !data.textbook_h2 || data.species_label != "H2") {
        std::cerr << "H2 textbook diagram classification failed\n";
        return 2;
    }
    if (data.levels[0].electrons.alpha != 1 || data.levels[0].electrons.beta != 1) {
        std::cerr << "diagram electron population failed\n";
        return 3;
    }

    const auto temp = std::filesystem::temp_directory_path() / "cov_mo_diagram_smoke";
    const auto result = cov::export_mo_diagram_bundle(wf, options, temp);
    if (!result.svg || !result.png || !result.json || !result.csv) {
        std::cerr << "diagram export failed: " << result.error << '\n';
        return 4;
    }

    const auto png_path = std::filesystem::path(temp.string() + ".mo.png");
    std::ifstream png(png_path, std::ios::binary);
    std::array<std::uint8_t, 8> signature{};
    png.read(reinterpret_cast<char*>(signature.data()), static_cast<std::streamsize>(signature.size()));
    const std::array<std::uint8_t, 8> expected{137,80,78,71,13,10,26,10};
    if (signature != expected) {
        std::cerr << "PNG signature mismatch\n";
        return 5;
    }

    const auto svg_path = std::filesystem::path(temp.string() + ".mo.svg");
    std::ifstream svg_file(svg_path, std::ios::binary);
    std::ostringstream svg_buffer;
    svg_buffer << svg_file.rdbuf();
    const std::string svg = svg_buffer.str();
    if (svg.find("fill='#ffffff'") == std::string::npos ||
        svg.find("H2 molecular-orbital diagram") == std::string::npos ||
        svg.find("Qualitative AO reference") == std::string::npos) {
        std::cerr << "textbook SVG layout marker missing\n";
        return 6;
    }

    std::error_code ec;
    std::filesystem::remove(svg_path, ec);
    std::filesystem::remove(png_path, ec);
    std::filesystem::remove(std::filesystem::path(temp.string() + ".mo.json"), ec);
    std::filesystem::remove(std::filesystem::path(temp.string() + ".mo.csv"), ec);

    std::cout << "mo_diagram_smoke ok\n";
    return 0;
}
