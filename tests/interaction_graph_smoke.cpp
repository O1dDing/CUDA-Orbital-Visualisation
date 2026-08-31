#include "cov/interaction_graph.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

void add_record(cov::Wavefunction& wf,
                const std::uint32_t a,
                const std::uint32_t b,
                const double mayer) {
    wf.bond_order_provenance = cov::DataProvenance::Derived;
    wf.bond_orders.push_back({a, b, mayer, cov::DataProvenance::Derived});
}

std::size_t count_kind(const cov::InteractionGraph& graph,
                       const cov::InteractionKind kind) {
    return static_cast<std::size_t>(std::count_if(
        graph.edges.begin(), graph.edges.end(),
        [&](const auto& edge) { return edge.kind == kind; }));
}

bool has_pair(const cov::InteractionGraph& graph,
              const std::uint32_t a,
              const std::uint32_t b,
              const cov::InteractionKind kind) {
    return std::any_of(graph.edges.begin(), graph.edges.end(), [&](const auto& edge) {
        return edge.kind == kind &&
               ((edge.atom_a == a && edge.atom_b == b) ||
                (edge.atom_a == b && edge.atom_b == a));
    });
}

} // namespace

int main() {
    // Zero- and one-atom inputs are ordinary graph boundaries, not errors.
    if (!cov::build_interaction_graph({}).fragment_analysis.fragments.empty()) {
        std::cerr << "empty wavefunction produced a fragment\n";
        return EXIT_FAILURE;
    }
    cov::Wavefunction atom;
    atom.atoms.push_back({"H", 1, 0.0, 0.0, 0.0});
    const auto atomic_graph = cov::build_interaction_graph(atom);
    if (atomic_graph.fragment_analysis.fragments.size() != 1u ||
        atomic_graph.fragment_analysis.fragments[0].atoms !=
            std::vector<std::uint32_t>{0u} ||
        !atomic_graph.edges.empty()) {
        std::cerr << "single-atom fragment boundary regressed\n";
        return EXIT_FAILURE;
    }

    // Linear F-H-F with a trusted 3c4e active subspace: only the two adjacent
    // F-H pairs support the hyperedge. The terminal F...F pair is not a bond.
    cov::Wavefunction hf2;
    hf2.atoms = {
        {"F", 9, -1.14 * cov::kAngstromToBohr, 0.0, 0.0},
        {"H", 1, 0.0, 0.0, 0.0},
        {"F", 9, 1.14 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    add_record(hf2, 0, 1, 0.32);
    add_record(hf2, 1, 2, 0.32);
    hf2.multicentre_assignments.push_back({
        cov::MulticentreKind::ThreeCentreFourElectron,
        {0, 1, 2}, {0, 1, 2}, 4.0, 0.91,
        "synthetic active-subspace evidence", cov::DataProvenance::Derived});
    const auto hf2_graph = cov::build_interaction_graph(hf2);
    if (hf2_graph.fragment_analysis.fragments.size() != 1u ||
        hf2_graph.multicentre_groups.size() != 1u ||
        count_kind(hf2_graph, cov::InteractionKind::MulticentreSupport) != 2u ||
        has_pair(hf2_graph, 0, 2, cov::InteractionKind::MulticentreSupport)) {
        std::cerr << "HF2-type multicentre support graph regressed\n";
        return EXIT_FAILURE;
    }
    for (const auto& edge : hf2_graph.edges) {
        if (edge.kind == cov::InteractionKind::MulticentreSupport &&
            cov::interaction_is_ordinary_bond(edge.kind)) {
            std::cerr << "multicentre support was reduced to an ordinary bond\n";
            return EXIT_FAILURE;
        }
    }

    // Minimum-spanning support must not break true symmetry or add a remote
    // terminal--terminal bond. An equilateral 3c2e H3 group keeps all three
    // equivalent supports, whereas a bridged B-H-B group keeps only B-H.
    cov::Wavefunction h3;
    constexpr double root3_over_2 = 0.86602540378443864676;
    h3.atoms = {
        {"H", 1, 0.0, 0.0, 0.0},
        {"H", 1, 0.90 * cov::kAngstromToBohr, 0.0, 0.0},
        {"H", 1, 0.45 * cov::kAngstromToBohr,
         0.90 * root3_over_2 * cov::kAngstromToBohr, 0.0},
    };
    h3.multicentre_assignments.push_back({
        cov::MulticentreKind::ThreeCentreTwoElectron,
        {0, 1, 2}, {0, 1, 2}, 2.0, 0.90,
        "synthetic symmetric active subspace", cov::DataProvenance::Derived});
    const auto h3_graph = cov::build_interaction_graph(h3);
    if (count_kind(h3_graph, cov::InteractionKind::MulticentreSupport) != 3u) {
        std::cerr << "symmetric H3 multicentre supports were arbitrarily broken\n";
        return EXIT_FAILURE;
    }

    cov::Wavefunction bhb;
    const double bridge_height = std::sqrt(1.33 * 1.33 - 0.885 * 0.885);
    bhb.atoms = {
        {"B", 5, -0.885 * cov::kAngstromToBohr, 0.0, 0.0},
        {"H", 1, 0.0, bridge_height * cov::kAngstromToBohr, 0.0},
        {"B", 5, 0.885 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    bhb.multicentre_assignments.push_back({
        cov::MulticentreKind::ThreeCentreTwoElectron,
        {0, 1, 2}, {0, 1, 2}, 2.0, 0.90,
        "synthetic bridge active subspace", cov::DataProvenance::Derived});
    const auto bhb_graph = cov::build_interaction_graph(bhb);
    if (count_kind(bhb_graph, cov::InteractionKind::MulticentreSupport) != 2u ||
        has_pair(bhb_graph, 0, 2, cov::InteractionKind::MulticentreSupport) ||
        has_pair(bhb_graph, 0, 2, cov::InteractionKind::CovalentConnectivity)) {
        std::cerr << "B-H-B support added a terminal B-B ordinary bond\n";
        return EXIT_FAILURE;
    }

    // Octahedral SbF6-like connectivity is a p-block molecular fragment, not
    // an automatically fabricated transition-metal ligand-field graph.
    cov::Wavefunction sbf6;
    sbf6.atoms.push_back({"Sb", 51, 0.0, 0.0, 0.0});
    constexpr double r = 1.90 * cov::kAngstromToBohr;
    sbf6.atoms.insert(sbf6.atoms.end(), {
        {"F", 9, r, 0.0, 0.0}, {"F", 9, -r, 0.0, 0.0},
        {"F", 9, 0.0, r, 0.0}, {"F", 9, 0.0, -r, 0.0},
        {"F", 9, 0.0, 0.0, r}, {"F", 9, 0.0, 0.0, -r},
    });
    for (std::uint32_t f = 1; f <= 6; ++f) add_record(sbf6, 0, f, 0.45);
    const auto sbf6_graph = cov::build_interaction_graph(sbf6);
    if (sbf6_graph.fragment_analysis.fragments.size() != 1u ||
        count_kind(sbf6_graph, cov::InteractionKind::CovalentConnectivity) != 6u ||
        count_kind(sbf6_graph, cov::InteractionKind::CoordinationContact) != 0u ||
        sbf6_graph.edges.size() != 6u) {
        std::cerr << "SbF6-type internal connectivity regressed\n";
        return EXIT_FAILURE;
    }

    // A separated positively charged atom plus the intact anion exercises the
    // generic ion-pair path. The weak edge is explicit and does not merge the
    // two strong fragments.
    cov::Wavefunction ion_pair = sbf6;
    ion_pair.atoms.push_back({"N", 7, 4.90 * cov::kAngstromToBohr, 0.0, 0.0});
    cov::InteractionEvidence charges;
    charges.atomic_charge_provenance = cov::DataProvenance::Derived;
    charges.atomic_partial_charges = {0.0, -0.30, -0.14, -0.14, -0.14, -0.14, -0.14, 1.0};
    const auto ion_graph = cov::build_interaction_graph(ion_pair, charges);
    if (ion_graph.fragment_analysis.fragments.size() != 2u ||
        count_kind(ion_graph, cov::InteractionKind::IonicContact) == 0u ||
        ion_graph.fragment_analysis.weak_interfragment_edge_indices.empty()) {
        std::cerr << "evidence-backed ion-pair contact regressed\n";
        return EXIT_FAILURE;
    }
    for (const auto& edge : ion_graph.edges) {
        if (edge.kind == cov::InteractionKind::IonicContact &&
            (edge.strength != cov::InteractionStrength::WeakContact ||
             cov::interaction_merges_fragments(edge.kind) ||
             cov::interaction_is_ordinary_bond(edge.kind))) {
            std::cerr << "ionic contact leaked into ordinary connectivity\n";
            return EXIT_FAILURE;
        }
    }

    // Two water fragments with a nearly linear O-H...O contact.
    cov::Wavefunction water_dimer;
    water_dimer.atoms = {
        {"O", 8, 0.0, 0.0, 0.0},
        {"H", 1, 0.96 * cov::kAngstromToBohr, 0.0, 0.0},
        {"H", 1, -0.30 * cov::kAngstromToBohr, 0.91 * cov::kAngstromToBohr, 0.0},
        {"O", 8, 2.76 * cov::kAngstromToBohr, 0.0, 0.0},
        {"H", 1, 3.72 * cov::kAngstromToBohr, 0.0, 0.0},
        {"H", 1, 2.46 * cov::kAngstromToBohr, 0.91 * cov::kAngstromToBohr, 0.0},
    };
    add_record(water_dimer, 0, 1, 0.75);
    add_record(water_dimer, 0, 2, 0.75);
    add_record(water_dimer, 3, 4, 0.75);
    add_record(water_dimer, 3, 5, 0.75);
    const auto water_graph = cov::build_interaction_graph(water_dimer);
    if (water_graph.fragment_analysis.fragments.size() != 2u ||
        !has_pair(water_graph, 1, 3, cov::InteractionKind::HydrogenBond)) {
        std::cerr << "water-dimer hydrogen-bond layer regressed\n";
        return EXIT_FAILURE;
    }

    // Low-Mayer octahedral metal--ligand contacts remain strong coordination
    // connectivity even below the ordinary renderer floor.
    cov::Wavefunction tif6;
    tif6.atoms.push_back({"Ti", 22, 0.0, 0.0, 0.0});
    constexpr double tr = 1.90 * cov::kAngstromToBohr;
    tif6.atoms.insert(tif6.atoms.end(), {
        {"F", 9, tr, 0.0, 0.0}, {"F", 9, -tr, 0.0, 0.0},
        {"F", 9, 0.0, tr, 0.0}, {"F", 9, 0.0, -tr, 0.0},
        {"F", 9, 0.0, 0.0, tr}, {"F", 9, 0.0, 0.0, -tr},
    });
    for (std::uint32_t f = 1; f <= 6; ++f) add_record(tif6, 0, f, 0.03);
    cov::InteractionEvidence ti_charges;
    ti_charges.atomic_charge_provenance = cov::DataProvenance::Derived;
    ti_charges.atomic_partial_charges = {2.4, -0.4, -0.4, -0.4, -0.4, -0.4, -0.4};
    const auto ti_graph = cov::build_interaction_graph(tif6, ti_charges);
    if (ti_graph.fragment_analysis.fragments.size() != 1u ||
        count_kind(ti_graph, cov::InteractionKind::CoordinationContact) != 6u ||
        cov::interaction_is_ordinary_bond(cov::InteractionKind::CoordinationContact)) {
        std::cerr << "low-Mayer coordination connectivity regressed\n";
        return EXIT_FAILURE;
    }

    // The same charge evidence has a different outcome for an s-block salt:
    // low Mayer support does not turn an ion-pair contact into a covalent bond.
    cov::Wavefunction nacl;
    nacl.atoms = {
        {"Na", 11, 0.0, 0.0, 0.0},
        {"Cl", 17, 2.36 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    add_record(nacl, 0, 1, 0.03);
    cov::InteractionEvidence nacl_charges;
    nacl_charges.atomic_charge_provenance = cov::DataProvenance::Derived;
    nacl_charges.atomic_partial_charges = {0.9, -0.9};
    const auto nacl_graph = cov::build_interaction_graph(nacl, nacl_charges);
    if (nacl_graph.fragment_analysis.fragments.size() != 2u ||
        count_kind(nacl_graph, cov::InteractionKind::IonicContact) != 1u) {
        std::cerr << "low-order s-block ion pair became strong connectivity\n";
        return EXIT_FAILURE;
    }

    // Producer atom charges are used automatically. A borderline H...F bridge
    // between two internally strong, oppositely charged fragments remains an
    // ionic contact; a genuine isolated polar H-F bond is not split merely
    // because its two Mulliken charges have opposite signs.
    cov::Wavefunction polyatomic_salt;
    polyatomic_salt.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"H", 1, 1.09 * cov::kAngstromToBohr, 0.0, 0.0},
        {"F", 9, 2.37 * cov::kAngstromToBohr, 0.0, 0.0},
        {"Sb", 51, 4.25 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    add_record(polyatomic_salt, 0, 1, 0.82);
    add_record(polyatomic_salt, 1, 2, 0.17);
    add_record(polyatomic_salt, 2, 3, 0.64);
    polyatomic_salt.atomic_partial_charges = {0.40, 0.40, -0.40, -0.40};
    polyatomic_salt.atomic_partial_charge_scheme = "Mulliken";
    polyatomic_salt.atomic_partial_charge_provenance = cov::DataProvenance::Producer;
    const auto salt_graph = cov::build_interaction_graph(polyatomic_salt);
    if (salt_graph.fragment_analysis.fragments.size() != 2u ||
        !has_pair(salt_graph, 1, 2, cov::InteractionKind::IonicContact) ||
        has_pair(salt_graph, 1, 2, cov::InteractionKind::CovalentConnectivity)) {
        std::cerr << "polyatomic ionic H...F bridge became covalent connectivity\n";
        return EXIT_FAILURE;
    }

    cov::Wavefunction polar_hf;
    polar_hf.atoms = {
        {"H", 1, 0.0, 0.0, 0.0},
        {"F", 9, 0.92 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    add_record(polar_hf, 0, 1, 0.17);
    polar_hf.atomic_partial_charges = {0.70, -0.70};
    polar_hf.atomic_partial_charge_scheme = "Mulliken";
    polar_hf.atomic_partial_charge_provenance = cov::DataProvenance::Producer;
    const auto polar_hf_graph = cov::build_interaction_graph(polar_hf);
    if (polar_hf_graph.fragment_analysis.fragments.size() != 1u ||
        !has_pair(polar_hf_graph, 0, 1,
                  cov::InteractionKind::CovalentConnectivity)) {
        std::cerr << "isolated polar H-F bond was split into ions\n";
        return EXIT_FAILURE;
    }

    // Through-bond Mayer coupling to the second atom of a ligand must not make
    // it an additional first-shell donor.
    cov::Wavefunction carbonyl;
    carbonyl.atoms = {
        {"Cr", 24, 0.0, 0.0, 0.0},
        {"C", 6, 1.80 * cov::kAngstromToBohr, 0.0, 0.0},
        {"O", 8, 2.93 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    add_record(carbonyl, 0, 1, 0.30);
    add_record(carbonyl, 1, 2, 1.30);
    add_record(carbonyl, 0, 2, 0.04);
    const auto carbonyl_graph = cov::build_interaction_graph(carbonyl);
    if (!has_pair(carbonyl_graph, 0, 1, cov::InteractionKind::CoordinationContact) ||
        !has_pair(carbonyl_graph, 1, 2, cov::InteractionKind::CovalentConnectivity) ||
        has_pair(carbonyl_graph, 0, 2, cov::InteractionKind::CoordinationContact)) {
        std::cerr << "shielded through-bond donor became first-shell coordination\n";
        return EXIT_FAILURE;
    }

    // Weak electronic evidence remains explicit and uncertain; it must not
    // silently merge two close fragments into a covalent molecule.
    cov::Wavefunction ambiguous;
    ambiguous.atoms = {
        {"C", 6, 0.0, 0.0, 0.0},
        {"C", 6, 2.20 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    add_record(ambiguous, 0, 1, 0.03);
    const auto ambiguous_graph = cov::build_interaction_graph(ambiguous);
    if (ambiguous_graph.fragment_analysis.fragments.size() != 2u ||
        count_kind(ambiguous_graph, cov::InteractionKind::AmbiguousContact) != 1u) {
        std::cerr << "weak electronic contact was promoted to connectivity\n";
        return EXIT_FAILURE;
    }

    // A close noble-gas pair has no electronic or charge basis for a stronger
    // label and therefore remains a generic noncovalent contact.
    cov::Wavefunction neon_pair;
    neon_pair.atoms = {
        {"Ne", 10, 0.0, 0.0, 0.0},
        {"Ne", 10, 3.00 * cov::kAngstromToBohr, 0.0, 0.0},
    };
    // Electronic availability prevents the conservative geometry fallback
    // from manufacturing a structural bond.
    neon_pair.bond_order_provenance = cov::DataProvenance::Derived;
    const auto neon_graph = cov::build_interaction_graph(neon_pair);
    if (neon_graph.fragment_analysis.fragments.size() != 2u ||
        count_kind(neon_graph, cov::InteractionKind::NoncovalentContact) != 1u) {
        std::cerr << "generic noncovalent contact layer regressed\n";
        return EXIT_FAILURE;
    }

    // A 12-vertex electron-deficient framework remains one molecular fragment
    // without pretending its low-order cage edges are ordinary two-centre
    // covalent bonds.
    cov::Wavefunction cage;
    constexpr double phi = 1.6180339887498948482;
    constexpr double scale = 0.875 * cov::kAngstromToBohr;
    const std::array<std::array<double, 3>, 12> vertices{{
        {{0, 1, phi}}, {{0, -1, phi}}, {{0, 1, -phi}}, {{0, -1, -phi}},
        {{1, phi, 0}}, {{-1, phi, 0}}, {{1, -phi, 0}}, {{-1, -phi, 0}},
        {{phi, 0, 1}}, {{phi, 0, -1}}, {{-phi, 0, 1}}, {{-phi, 0, -1}},
    }};
    for (const auto& vertex : vertices) {
        cage.atoms.push_back({"B", 5, scale * vertex[0], scale * vertex[1],
                              scale * vertex[2], 5.0});
    }
    cage.bond_order_provenance = cov::DataProvenance::Derived;
    for (std::uint32_t a = 0u; a < cage.atoms.size(); ++a) {
        for (std::uint32_t b = a + 1u; b < cage.atoms.size(); ++b) {
            const double dx = cage.atoms[a].x - cage.atoms[b].x;
            const double dy = cage.atoms[a].y - cage.atoms[b].y;
            const double dz = cage.atoms[a].z - cage.atoms[b].z;
            if (std::sqrt(dx * dx + dy * dy + dz * dz) < 2.05 * scale) {
                // One real-looking cage edge is above the ordinary renderer
                // floor; the higher-level cage interpretation must replace it.
                add_record(cage, a, b, a == 0u && b == 1u ? 0.30 : 0.03);
            }
        }
    }
    // One triangular face also carries a local 3c2e assignment. Its three
    // pairs must retain both independent support semantics: the local
    // multicentre hyperedge and the global cage skeleton. Neither layer is an
    // ordinary two-centre bond, and either can be hidden independently by the
    // renderer.
    cage.multicentre_assignments.push_back({
        cov::MulticentreKind::ThreeCentreTwoElectron,
        {0, 1, 8}, {0}, 2.0, 0.90,
        "synthetic cage-face active subspace", cov::DataProvenance::Derived});
    const auto cage_graph = cov::build_interaction_graph(cage);
    if (cage_graph.polyhedral_cages.size() != 1u ||
        cage_graph.polyhedral_cages[0].support_edges.size() != 30u ||
        cage_graph.fragment_analysis.fragments.size() != 1u ||
        count_kind(cage_graph, cov::InteractionKind::PolyhedralCageSupport) != 30u ||
        count_kind(cage_graph, cov::InteractionKind::MulticentreSupport) != 3u ||
        !has_pair(cage_graph, 0, 1, cov::InteractionKind::MulticentreSupport) ||
        !has_pair(cage_graph, 0, 1, cov::InteractionKind::PolyhedralCageSupport) ||
        !has_pair(cage_graph, 0, 8, cov::InteractionKind::MulticentreSupport) ||
        !has_pair(cage_graph, 0, 8, cov::InteractionKind::PolyhedralCageSupport) ||
        !has_pair(cage_graph, 1, 8, cov::InteractionKind::MulticentreSupport) ||
        !has_pair(cage_graph, 1, 8, cov::InteractionKind::PolyhedralCageSupport) ||
        count_kind(cage_graph, cov::InteractionKind::CovalentConnectivity) != 0u) {
        std::cerr << "overlapping multicentre/polyhedral support regression\n";
        return EXIT_FAILURE;
    }

    std::cout << "interaction graph smoke test passed\n";
    return EXIT_SUCCESS;
}
