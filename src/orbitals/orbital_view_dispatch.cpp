#include "cov/orbital_view.hpp"

#include <algorithm>

namespace cov {

OrbitalRegion classify_orbital_region_legacy(
    const MolecularOrbital& orbital,
    const OrbitalFilterSettings& settings) noexcept;
bool orbital_visible_legacy(
    const std::vector<MolecularOrbital>& orbitals,
    std::size_t index,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings) noexcept;
std::vector<std::size_t> visible_orbital_indices_legacy(
    const std::vector<MolecularOrbital>& orbitals,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings);
std::vector<OrbitalMetadata> build_orbital_metadata_legacy(
    const Wavefunction& wavefunction,
    std::size_t selected_index,
    const DegeneracySettings& degeneracy,
    const OrbitalFilterSettings& filter);

namespace {

bool chemistry_available(const std::vector<MolecularOrbital>& orbitals) {
    return std::any_of(orbitals.begin(),orbitals.end(),
        [](const MolecularOrbital& orbital) {
            return orbital.chemistry.available;
        });
}

bool occupied(const MolecularOrbital& orbital,const double threshold) {
    return static_cast<double>(orbital.occupation)>threshold;
}

} // namespace

OrbitalRegion classify_orbital_region(
    const MolecularOrbital& orbital,
    const OrbitalFilterSettings& settings) noexcept {
    if (!orbital.chemistry.available) {
        return classify_orbital_region_legacy(orbital,settings);
    }
    if (orbital.chemistry.valence_manifold) return OrbitalRegion::Valence;
    if (!occupied(orbital,settings.occupation_threshold)) {
        return OrbitalRegion::Virtual;
    }
    return OrbitalRegion::Core;
}

bool orbital_visible(
    const std::vector<MolecularOrbital>& orbitals,
    const std::size_t index,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings) noexcept {
    if (index>=orbitals.size()) return false;
    if (!chemistry_available(orbitals)) {
        return orbital_visible_legacy(orbitals,index,frontier,settings);
    }

    const auto& orbital=orbitals[index];
    const auto region=classify_orbital_region(orbital,settings);
    switch (settings.mode) {
        case OrbitalFilterMode::All:
            return true;
        case OrbitalFilterMode::Occupied:
            return occupied(orbital,settings.occupation_threshold);
        case OrbitalFilterMode::Virtual:
            return !occupied(orbital,settings.occupation_threshold);
        case OrbitalFilterMode::Core:
            return region==OrbitalRegion::Core;
        case OrbitalFilterMode::Valence:
        case OrbitalFilterMode::AutoReasonable:
        default:
            return orbital.chemistry.available &&
                   orbital.chemistry.valence_manifold;
    }
}

std::vector<std::size_t> visible_orbital_indices(
    const std::vector<MolecularOrbital>& orbitals,
    const FrontierOrbitals& frontier,
    const OrbitalFilterSettings& settings) {
    std::vector<std::size_t> result;
    result.reserve(orbitals.size());
    for (std::size_t i=0;i<orbitals.size();++i) {
        if (orbital_visible(orbitals,i,frontier,settings)) {
            result.push_back(i);
        }
    }
    return result;
}

std::vector<OrbitalMetadata> build_orbital_metadata(
    const Wavefunction& wavefunction,
    const std::size_t selected_index,
    const DegeneracySettings& degeneracy,
    const OrbitalFilterSettings& filter) {
    auto result=build_orbital_metadata_legacy(
        wavefunction,selected_index,degeneracy,filter);
    const auto frontier=find_frontier_orbitals(
        wavefunction.orbitals,filter.occupation_threshold);
    for (std::size_t i=0;i<result.size();++i) {
        result[i].region=classify_orbital_region(
            wavefunction.orbitals[i],filter);
        result[i].visible=orbital_visible(
            wavefunction.orbitals,i,frontier,filter);
    }
    return result;
}

} // namespace cov
