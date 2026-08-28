#include "cov/gaussian_log.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace cov {
namespace {

std::string trim(std::string value) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(), value.end());
    return value;
}

std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool starts_with(const std::string& value, const char* prefix) {
    const std::string p(prefix);
    return value.size() >= p.size() && value.compare(0, p.size(), p) == 0;
}

std::vector<std::string> parenthesized_labels(const std::string& line) {
    std::vector<std::string> labels;
    std::size_t cursor = 0;
    while (true) {
        const auto open = line.find('(', cursor);
        if (open == std::string::npos) break;
        const auto close = line.find(')', open + 1u);
        if (close == std::string::npos) break;
        std::string label = trim(line.substr(open + 1u, close - open - 1u));
        if (!label.empty()) labels.push_back(std::move(label));
        cursor = close + 1u;
    }
    return labels;
}

std::string token_after(const std::string& line, const char* prefix) {
    const auto pos = line.find(prefix);
    if (pos == std::string::npos) return {};
    std::istringstream stream(line.substr(pos + std::string(prefix).size()));
    std::string token;
    stream >> token;
    return token;
}

struct SymmetryBlock {
    std::vector<std::string> alpha;
    std::vector<std::string> beta;
    bool explicit_spin = false;
};

} // namespace

GaussianLogEnrichmentResult enrich_from_gaussian_log(
    Wavefunction& wavefunction,
    const std::filesystem::path& path) {
    GaussianLogEnrichmentResult result;
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        result.warning = "Unable to open Gaussian enrichment file: " + path.string();
        return result;
    }
    result.opened = true;

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) lines.push_back(line);

    std::string point_group_detected;
    std::string point_group_used;
    SymmetryBlock last_block;

    for (std::size_t i = 0; i < lines.size(); ++i) {
        const std::string t = trim(lines[i]);
        if (t.find("Full point group") != std::string::npos) {
            const auto value = token_after(t, "Full point group");
            if (!value.empty()) point_group_detected = value;
        }
        if (t.find("Largest Abelian subgroup") != std::string::npos) {
            const auto value = token_after(t, "Largest Abelian subgroup");
            if (!value.empty()) point_group_used = value;
        }

        // Deliberately exclude "Initial guess orbital symmetries:". In an
        // optimization or converged SCF output there may be several final
        // blocks; the last complete producer block wins.
        if (t != "Orbital symmetries:") continue;

        SymmetryBlock block;
        Spin current_spin = Spin::Alpha;
        for (std::size_t j = i + 1u; j < lines.size() && j < i + 80u; ++j) {
            const std::string s = trim(lines[j]);
            if (s == "Alpha Orbitals:" || s == "Alpha orbitals:") {
                current_spin = Spin::Alpha;
                block.explicit_spin = true;
                continue;
            }
            if (s == "Beta Orbitals:" || s == "Beta orbitals:") {
                current_spin = Spin::Beta;
                block.explicit_spin = true;
                continue;
            }
            if (starts_with(s, "The electronic state") ||
                starts_with(s, "Alpha  occ. eigenvalues") ||
                starts_with(s, "Beta  occ. eigenvalues") ||
                starts_with(s, "SCF Done:") ||
                s.find("Leave Link") != std::string::npos) {
                break;
            }

            const auto labels = parenthesized_labels(s);
            if (labels.empty()) continue;
            auto& target = current_spin == Spin::Alpha ? block.alpha : block.beta;
            target.insert(target.end(), labels.begin(), labels.end());
        }

        if (!block.alpha.empty() || !block.beta.empty()) {
            last_block = std::move(block);
        }
    }

    if (!point_group_detected.empty() || !point_group_used.empty()) {
        wavefunction.point_group_detected = point_group_detected;
        wavefunction.point_group_used = point_group_used.empty()
                                            ? point_group_detected
                                            : point_group_used;
        wavefunction.point_group_provenance = DataProvenance::Producer;
        result.point_group_applied = true;
    }

    std::vector<std::size_t> alpha_indices;
    std::vector<std::size_t> beta_indices;
    for (std::size_t i = 0; i < wavefunction.orbitals.size(); ++i) {
        if (wavefunction.orbitals[i].spin == Spin::Beta) beta_indices.push_back(i);
        else alpha_indices.push_back(i);
    }

    bool symmetry_dimensions_ok = true;
    if (!last_block.beta.empty() || last_block.explicit_spin) {
        symmetry_dimensions_ok =
            last_block.alpha.size() == alpha_indices.size() &&
            last_block.beta.size() == beta_indices.size();
    } else if (!beta_indices.empty()) {
        // Never mirror one producer-reported spin block onto another.
        symmetry_dimensions_ok = false;
    } else {
        symmetry_dimensions_ok = last_block.alpha.size() == alpha_indices.size();
    }

    if (symmetry_dimensions_ok && !last_block.alpha.empty()) {
        for (std::size_t i = 0; i < alpha_indices.size(); ++i) {
            auto& mo = wavefunction.orbitals[alpha_indices[i]];
            mo.symmetry = last_block.alpha[i];
            mo.symmetry_provenance = DataProvenance::Producer;
            ++result.symmetry_labels_applied;
        }
        for (std::size_t i = 0; i < beta_indices.size(); ++i) {
            auto& mo = wavefunction.orbitals[beta_indices[i]];
            mo.symmetry = last_block.beta[i];
            mo.symmetry_provenance = DataProvenance::Producer;
            ++result.symmetry_labels_applied;
        }
    } else if (!last_block.alpha.empty() || !last_block.beta.empty()) {
        result.warning =
            "Gaussian orbital-symmetry count does not match the loaded MO blocks; "
            "symmetry enrichment was left unapplied";
    }

    if (result.applied()) wavefunction.enrichment_source = path.string();
    return result;
}

std::filesystem::path find_sibling_gaussian_log(
    const std::filesystem::path& wavefunction_path) {
    std::string filename = wavefunction_path.filename().string();
    const std::string folded = lower(filename);
    std::string root;
    const auto molden_pos = folded.find(".molden");
    if (molden_pos != std::string::npos) {
        root = filename.substr(0, molden_pos);
    } else {
        root = wavefunction_path.stem().string();
    }
    if (root.empty()) return {};

    const auto parent = wavefunction_path.parent_path();
    std::error_code ec;
    for (const char* extension : {".log", ".out"}) {
        const auto candidate = parent / (root + extension);
        if (std::filesystem::exists(candidate, ec) && !ec) return candidate;
        ec.clear();
    }
    return {};
}

} // namespace cov
