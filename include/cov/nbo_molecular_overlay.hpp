#pragma once
#include "cov/molecular_overlay.hpp"
#include "cov/nbo_integration.hpp"
#include "cov/interaction_graph.hpp"
namespace cov {
MoleculeOverlay make_nbo_molecule_overlay(const NboIntegration&,const InteractionGraph&,
    std::size_t atom_count,const std::vector<std::size_t>& selected_atoms,
    std::optional<std::size_t> selected_structure,AtomScalarMode,bool show_indices,bool show_e2);
}
