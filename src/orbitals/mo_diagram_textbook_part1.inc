#include "cov/mo_diagram.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string_view>
#include <vector>

namespace cov {

// These functions are the original exporter/model implementation compiled
// under renamed symbols. Keeping them behind the new presentation layer lets
// the human-facing diagram evolve without perturbing the already-tested
// orbital metadata, filtering, JSON and CSV paths.
MODiagramData cov_legacy_build_mo_diagram_data(const Wavefunction&, const MODiagramOptions&);
bool cov_legacy_write_mo_diagram_json(const MODiagramData&, const MODiagramOptions&,
                                      const std::filesystem::path&, std::string*);
bool cov_legacy_write_mo_diagram_csv(const MODiagramData&, const MODiagramOptions&,
                                     const std::filesystem::path&, std::string*);

namespace {

std::string xml_escape(std::string_view value) {
    std::string out;
    for (const char c : value) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&apos;"; break;
            default: out += c; break;
        }
    }
    return out;
}

const char* classification_title(const DiagramClassification value) noexcept {
    switch (value) {
        case DiagramClassification::SymmetryGrouped:
            return "Symmetry-grouped MO diagram";
        case DiagramClassification::SalcUnavailable:
            return "MO diagram — SALC not confidently available";
        default:
            return "MO energy diagram";
    }
}

bool homonuclear_diatomic(const Wavefunction& wavefunction) {
    if (wavefunction.atoms.size() != 2 || wavefunction.orbitals.size() < 2 ||
        wavefunction.orbitals.size() > 8) {
        return false;
    }
    const Atom& a = wavefunction.atoms[0];
    const Atom& b = wavefunction.atoms[1];
    if (a.atomic_number > 0 && b.atomic_number > 0) {
        return a.atomic_number == b.atomic_number;
    }
    if (!a.symbol.empty() && !b.symbol.empty()) return a.symbol == b.symbol;
    return true;
}

bool is_h2(const Wavefunction& wavefunction) {
    if (wavefunction.atoms.size() != 2) return false;
    const Atom& a = wavefunction.atoms[0];
    const Atom& b = wavefunction.atoms[1];
    return (a.atomic_number == 1 && b.atomic_number == 1) ||
           (a.symbol == "H" && b.symbol == "H");
}

std::string species_label(const Wavefunction& wavefunction) {
    if (wavefunction.atoms.size() != 2) return {};
    std::string a = wavefunction.atoms[0].symbol.empty() ? "A" : wavefunction.atoms[0].symbol;
    std::string b = wavefunction.atoms[1].symbol.empty() ? a : wavefunction.atoms[1].symbol;
    return a == b ? a + "2" : a + b;
}

struct Layout {
    double e_min = -1.0;
    double e_max = 1.0;
    std::vector<std::string> columns;
    std::map<std::string, int> column_index;
};

Layout make_layout(const MODiagramData& data) {
    Layout layout;
    if (!data.levels.empty()) {
        layout.e_min = layout.e_max = data.levels.front().metadata.energy_hartree;
        for (const auto& level : data.levels) {
            layout.e_min = std::min(layout.e_min, level.metadata.energy_hartree);
            layout.e_max = std::max(layout.e_max, level.metadata.energy_hartree);
        }
        const double span = std::max(1.0e-6, layout.e_max - layout.e_min);
        layout.e_min -= 0.10 * span;
        layout.e_max += 0.10 * span;
    }
    if (data.plan.classification == DiagramClassification::SymmetryGrouped) {
        std::set<std::string> labels;
        for (const auto& level : data.levels) {
            labels.insert(level.metadata.symmetry.empty() ? "?" : level.metadata.symmetry);
        }
        layout.columns.assign(labels.begin(), labels.end());
    }
    if (layout.columns.empty()) layout.columns.push_back("MO");
    for (std::size_t i = 0; i < layout.columns.size(); ++i) {
        layout.column_index[layout.columns[i]] = static_cast<int>(i);
    }
    return layout;
}

double energy_y(const double energy, const Layout& layout, const int height,
                const int top = 118, const int bottom = 80) {
    const double t = (energy - layout.e_min) /
                     std::max(1.0e-12, layout.e_max - layout.e_min);
    return static_cast<double>(height - bottom) -
           t * static_cast<double>(height - top - bottom);
}

double column_x(const MODiagramLevel& level, const Layout& layout, const int width) {
    const std::string key = layout.columns.size() > 1
        ? (level.metadata.symmetry.empty() ? "?" : level.metadata.symmetry)
        : "MO";
    const auto it = layout.column_index.find(key);
    const int index = it == layout.column_index.end() ? 0 : it->second;
    const double span = static_cast<double>(width - 330);
    return 150.0 + (static_cast<double>(index) + 0.5) * span /
                   static_cast<double>(std::max<std::size_t>(1, layout.columns.size()));
}

int degeneracy_offset(const MODiagramLevel& level) {
    const auto dash = level.metadata.display_label.find('-');
    if (dash == std::string::npos || dash + 1 >= level.metadata.display_label.size()) return 0;
    const char suffix = level.metadata.display_label[dash + 1];
    if (suffix < 'a' || suffix > 'z') return 0;
    const double centred = static_cast<int>(suffix - 'a') -
        (static_cast<int>(level.metadata.degeneracy_size) - 1) * 0.5;
    return static_cast<int>(std::lround(centred * 26.0));
}

struct PositionedLevel {
    const MODiagramLevel* level = nullptr;
    double x = 0.0;
    double y = 0.0;
    double label_y = 0.0;
};

std::vector<PositionedLevel> position_levels(const MODiagramData& data,
                                             const Layout& layout,
                                             const int width,
                                             const int height) {
    std::vector<PositionedLevel> positioned;
    for (const auto& level : data.levels) {
        const double y = energy_y(level.metadata.energy_hartree, layout, height);
        positioned.push_back({&level, column_x(level, layout, width) + degeneracy_offset(level), y, y});
    }

    // Labels are allowed to move vertically, but the physical level line is not.
    // This preserves quantitative energy placement while avoiding unreadable
    // collisions for dense/near-degenerate regions.
    std::map<int, std::vector<std::size_t>> lanes;
    for (std::size_t i = 0; i < positioned.size(); ++i) {
        lanes[static_cast<int>(std::lround(positioned[i].x / 80.0))].push_back(i);
    }
    for (auto& [_, ids] : lanes) {
        std::sort(ids.begin(), ids.end(), [&](const std::size_t a, const std::size_t b) {
            return positioned[a].label_y < positioned[b].label_y;
        });
        double previous = 98.0;
        for (const std::size_t id : ids) {
            positioned[id].label_y = std::max(positioned[id].label_y, previous + 24.0);
            previous = positioned[id].label_y;
        }
        const double overflow = previous - static_cast<double>(height - 64);
        if (overflow > 0.0) {
            for (const std::size_t id : ids) positioned[id].label_y -= overflow;
        }
    }
    return positioned;
